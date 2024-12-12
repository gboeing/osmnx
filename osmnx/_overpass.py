"""Tools to work with the Overpass API."""

from __future__ import annotations

import datetime as dt
import logging as lg
import time
from collections import OrderedDict
from typing import TYPE_CHECKING
from typing import Any

import numpy as np
import requests
from requests.exceptions import ConnectionError as RequestsConnectionError

from . import _http
from . import projection
from . import settings
from . import utils
from . import utils_geo
from ._errors import InsufficientResponseError

if TYPE_CHECKING:
    from collections.abc import Iterator

    from shapely import MultiPolygon
    from shapely import Polygon


def _get_network_filter(network_type: str) -> str:
    """
    Create a filter to query Overpass for the specified network type.

    The filter queries Overpass for every OSM way with a "highway" tag but
    excludes ways that are incompatible with the requested network type. You
    can choose from the following types:

    "all" retrieves all public and private-access ways currently in use.

    "all_public" retrieves all public ways currently in use.

    "bike" retrieves public bikeable ways and excludes foot ways, motor ways,
    and anything tagged biking=no.

    "drive" retrieves public drivable streets and excludes service roads,
    anything tagged motor=no, and certain non-service roads tagged as
    providing certain services (such as alleys or driveways).

    "drive_service" retrieves public drivable streets including service roads
    but excludes certain services (such as parking or emergency access).

    "walk" retrieves public walkable ways and excludes cycle ways, motor ways,
    and anything tagged foot=no. It includes service roads like parking lot
    aisles and alleys that you can walk on even if they are unpleasant walks.

    Parameters
    ----------
    network_type
        {"all", "all_public", "bike", "drive", "drive_service", "walk"}
        What type of street network to retrieve.

    Returns
    -------
    way_filter
        The Overpass query filter.
    """
    # define built-in queries to send to the API. specifying way["highway"]
    # means that all ways returned must have a highway tag. the filters then
    # remove ways by tag/value.
    filters = {}

    # driving: filter out un-drivable roads, service roads, private ways, and
    # anything tagged motor=no. also filter out any non-service roads that are
    # tagged as providing certain services
    filters["drive"] = (
        f'["highway"]["area"!~"yes"]{settings.default_access}'
        f'["highway"!~"abandoned|bridleway|bus_guideway|construction|corridor|cycleway|elevator|'
        f"escalator|footway|no|path|pedestrian|planned|platform|proposed|raceway|razed|service|"
        f'steps|track"]'
        f'["motor_vehicle"!~"no"]["motorcar"!~"no"]'
        f'["service"!~"alley|driveway|emergency_access|parking|parking_aisle|private"]'
    )

    # drive+service: allow ways tagged 'service' but filter out certain types
    filters["drive_service"] = (
        f'["highway"]["area"!~"yes"]{settings.default_access}'
        f'["highway"!~"abandoned|bridleway|bus_guideway|construction|corridor|cycleway|elevator|'
        f"escalator|footway|no|path|pedestrian|planned|platform|proposed|raceway|razed|steps|"
        f'track"]'
        f'["motor_vehicle"!~"no"]["motorcar"!~"no"]'
        f'["service"!~"emergency_access|parking|parking_aisle|private"]'
    )

    # walking: filter out cycle ways, motor ways, private ways, and anything
    # tagged foot=no. allow service roads, permitting things like parking lot
    # aisles, alleys, etc that you *can* walk on even if they're not exactly
    # pleasant walks. some cycleways may allow pedestrians, but this filter
    # ignores such cycleways.
    filters["walk"] = (
        f'["highway"]["area"!~"yes"]{settings.default_access}'
        f'["highway"!~"abandoned|bus_guideway|construction|cycleway|motor|no|planned|platform|'
        f'proposed|raceway|razed"]'
        f'["foot"!~"no"]["service"!~"private"]'
        f'["sidewalk"!~"separate"]["sidewalk:both"!~"separate"]'
        f'["sidewalk:left"!~"separate"]["sidewalk:right"!~"separate"]'
    )

    # biking: filter out foot ways, motor ways, private ways, and anything
    # tagged biking=no
    filters["bike"] = (
        f'["highway"]["area"!~"yes"]{settings.default_access}'
        f'["highway"!~"abandoned|bus_guideway|construction|corridor|elevator|escalator|footway|'
        f'motor|no|planned|platform|proposed|raceway|razed|steps"]'
        f'["bicycle"!~"no"]["service"!~"private"]'
    )

    # to download all public ways, just filter out everything not currently in
    # use or that is private-access only
    filters["all_public"] = (
        f'["highway"]["area"!~"yes"]{settings.default_access}'
        f'["highway"!~"abandoned|construction|no|planned|platform|proposed|raceway|razed"]'
        f'["service"!~"private"]'
    )

    # to download all ways, including private-access ones, just filter out
    # everything not currently in use
    filters["all"] = (
        '["highway"]["area"!~"yes"]["highway"!~"abandoned|construction|no|planned|platform|'
        'proposed|raceway|razed"]'
    )

    if network_type in filters:
        way_filter = filters[network_type]
    else:  # pragma: no cover
        msg = f"Unrecognized network_type {network_type!r}."
        raise ValueError(msg)

    return way_filter


def _get_overpass_pause(
    base_endpoint: str,
    *,
    recursive_delay: float = 5,
    default_duration: float = 60,
) -> float:
    """
    Retrieve a pause duration from the Overpass API status endpoint.

    Check the Overpass API status endpoint to determine how long to wait until
    the next slot is available. You can disable this via the `settings`
    module's `overpass_rate_limit` setting.

    Parameters
    ----------
    base_endpoint
        Base Overpass API URL (without "/status" at the end).
    recursive_delay
        How long to wait between recursive calls if the server is currently
        running a query.
    default_duration
        If a fatal error occurs, fall back on returning this value.

    Returns
    -------
    pause
        The current pause duration specified by the Overpass status endpoint.
    """
    if not settings.overpass_rate_limit:
        # if overpass rate limiting is False, then there is zero pause
        return 0

    url = base_endpoint.rstrip("/") + "/status"

    # try to retrieve the URL
    try:
        response = requests.get(
            url,
            headers=_http._get_http_headers(),
            timeout=settings.requests_timeout,
            **settings.requests_kwargs,
        )
        response_text = response.text
    except RequestsConnectionError as e:  # pragma: no cover
        # cannot reach status endpoint: log error and return default duration
        msg = f"Unable to query {url}, {e}"
        utils.log(msg, level=lg.ERROR)
        return default_duration

    # try to parse the output
    try:
        status = response_text.split("\n")[4]
        status_first_part = status.split(" ")[0]
    except (AttributeError, IndexError, ValueError):  # pragma: no cover
        # cannot parse output: log error and return default duration
        msg = f"Unable to parse {url} response: {response_text}"
        utils.log(msg, level=lg.ERROR)
        return default_duration

    # determine the current status of the server
    try:
        # if first token is numeric, it's how many slots you have available,
        # no wait required
        _ = int(status_first_part)  # number of available slots
        pause: float = 0

    except ValueError:  # pragma: no cover
        # if first token is 'Slot', it tells you when your slot will be free
        if status_first_part == "Slot":
            utc_time_str = status.split(" ")[3]
            pattern = "%Y-%m-%dT%H:%M:%SZ,"
            utc_time = dt.datetime.strptime(utc_time_str, pattern).astimezone(dt.timezone.utc)
            utc_now = dt.datetime.now(tz=dt.timezone.utc)
            seconds = int(np.ceil((utc_time - utc_now).total_seconds()))
            pause = max(seconds, 1)

        # if first token is 'Currently', it is currently running a query so
        # check back in recursive_delay seconds
        elif status_first_part == "Currently":
            time.sleep(recursive_delay)
            pause = _get_overpass_pause(base_endpoint)

        # any other status is unrecognized: log error, return default duration
        else:
            msg = f"Unrecognized server status: {status!r}"
            utils.log(msg, level=lg.ERROR)
            return default_duration

    return pause


def _make_overpass_settings() -> str:
    """
    Make settings string to send in Overpass query.

    Returns
    -------
    overpass_settings
        The `settings.overpass_settings` string formatted with "timeout" and
        "maxsize" values.
    """
    maxsize = "" if settings.overpass_memory is None else f"[maxsize:{settings.overpass_memory}]"
    return settings.overpass_settings.format(timeout=settings.requests_timeout, maxsize=maxsize)


def _make_overpass_polygon_coord_strs(polygon: Polygon | MultiPolygon) -> list[str]:
    """
    Subdivide query polygon and return list of coordinate strings.

    Project to UTM, divide `polygon` up into sub-polygons if area exceeds a
    max size (in meters), project back to lat-lon, then get a list of
    polygon(s) exterior coordinates. Ignore interior ("holes") coordinates.

    Parameters
    ----------
    polygon
        The (Multi)Polygon to convert to exterior coordinate strings.

    Returns
    -------
    coord_strs
        Exterior coordinates of polygon(s).
    """
    # first subdivide the polygon if its area exceeds max size
    # this results in a multipolygon of 1+ constituent polygons
    poly_proj, crs_proj = projection.project_geometry(polygon)
    multi_poly_proj = utils_geo._consolidate_subdivide_geometry(poly_proj)
    multi_poly, _ = projection.project_geometry(multi_poly_proj, crs=crs_proj, to_latlong=True)

    # then extract each's exterior coords to the string format Overpass
    # expects, rounding lats and lons to 6 decimals (approx 5 to 10 cm
    # resolution) so we can hash and cache URL strings consistently
    coord_strs = []
    for geom in multi_poly.geoms:
        x, y = geom.exterior.xy
        coord_list = [f"{xy[1]:.6f}{' '}{xy[0]:.6f}" for xy in zip(x, y)]
        coord_strs.append(" ".join(coord_list))

    return coord_strs


def _create_overpass_features_query(  # noqa: PLR0912
    polygon_coord_str: str,
    tags: dict[str, bool | str | list[str]],
) -> str:
    """
    Create an Overpass features query string based on tags.

    Parameters
    ----------
    polygon_coord_str
        The lat lon coordinates.
    tags
        Tags used for finding elements in the search area.

    Returns
    -------
    query
    """
    # create overpass settings string
    overpass_settings = _make_overpass_settings()

    # make sure every value in dict is bool, str, or list of str
    err_msg = "`tags` must be a dict with values of bool, str, or list of str."
    if not isinstance(tags, dict):  # pragma: no cover
        raise TypeError(err_msg)

    tags_dict: dict[str, bool | str | list[str]] = {}
    for key, value in tags.items():
        if isinstance(value, bool):
            tags_dict[key] = value

        elif isinstance(value, str):
            tags_dict[key] = [value]

        elif isinstance(value, list):
            if not all(isinstance(s, str) for s in value):  # pragma: no cover
                raise TypeError(err_msg)
            tags_dict[key] = value

        else:  # pragma: no cover
            raise TypeError(err_msg)

    # convert the tags dict into a list of {tag:value} dicts
    tags_list: list[dict[str, bool | str | list[str]]] = []
    for key, value in tags_dict.items():
        if isinstance(value, bool):
            tags_list.append({key: value})
        else:
            for value_item in value:
                tags_list.append({key: value_item})  # noqa: PERF401

    # add node/way/relation query components one at a time
    components = []
    for d in tags_list:
        for key, value in d.items():
            if isinstance(value, bool):
                # if bool (ie, True) just pass the key, no value
                tag_str = f"[{key!r}](poly:{polygon_coord_str!r});(._;>;);"
            else:
                # otherwise, pass "key"="value"
                tag_str = f"[{key!r}={value!r}](poly:{polygon_coord_str!r});(._;>;);"

            for kind in ("node", "way", "relation"):
                components.append(f"({kind}{tag_str});")  # noqa: PERF401

    # finalize query and return
    components_str = "".join(components)
    return f"{overpass_settings};({components_str});out;"


def _download_overpass_network(
    polygon: Polygon | MultiPolygon,
    network_type: str,
    custom_filter: str | list[str] | None,
) -> Iterator[dict[str, Any]]:
    """
    Retrieve networked ways and nodes within boundary from the Overpass API.

    Parameters
    ----------
    polygon
        The boundary to fetch the network ways/nodes within.
    network_type
        What type of street network to get if `custom_filter` is None.
    custom_filter
        A custom "ways" filter to be used instead of `network_type` presets.

    Yields
    ------
    response_json
        JSON response from the Overpass server.
    """
    # create filter(s) to exclude certain kinds of ways based on the requested
    # network_type, if provided, otherwise use custom_filter
    way_filters = []
    if isinstance(custom_filter, list):
        way_filters = custom_filter
    elif isinstance(custom_filter, str):
        way_filters = [custom_filter]
    else:
        way_filters = [_get_network_filter(network_type)]

    # create overpass settings string
    overpass_settings = _make_overpass_settings()

    # subdivide query polygon to get list of sub-divided polygon coord strings
    polygon_coord_strs = _make_overpass_polygon_coord_strs(polygon)
    msg = f"Requesting data from API in {len(polygon_coord_strs)} request(s)"
    utils.log(msg, level=lg.INFO)

    # pass exterior coordinates of each polygon in list to API, one at a time
    # the '>' makes it recurse so we get ways and the ways' nodes.
    for polygon_coord_str in polygon_coord_strs:
        for way_filter in way_filters:
            query_str = f"{overpass_settings};(way{way_filter}(poly:{polygon_coord_str!r});>;);out;"
            yield _overpass_request(OrderedDict(data=query_str))


def _download_overpass_features(
    polygon: Polygon,
    tags: dict[str, bool | str | list[str]],
) -> Iterator[dict[str, Any]]:
    """
    Retrieve OSM features within some boundary polygon from the Overpass API.

    Parameters
    ----------
    polygon
        Boundary to retrieve elements within.
    tags
        Tags used for finding elements in the selected area.

    Yields
    ------
    response_json
        JSON response from the Overpass server.
    """
    # subdivide query polygon to get list of sub-divided polygon coord strings
    polygon_coord_strs = _make_overpass_polygon_coord_strs(polygon)
    msg = f"Requesting data from API in {len(polygon_coord_strs)} request(s)"
    utils.log(msg, level=lg.INFO)

    # pass exterior coordinates of each polygon in list to API, one at a time
    for polygon_coord_str in polygon_coord_strs:
        query_str = _create_overpass_features_query(polygon_coord_str, tags)
        yield _overpass_request(OrderedDict(data=query_str))


def _overpass_request(
    data: OrderedDict[str, Any],
    *,
    pause: float | None = None,
    error_pause: float = 60,
) -> dict[str, Any]:
    """
    Send a HTTP POST request to the Overpass API and return response.

    Parameters
    ----------
    data
        Key-value pairs of parameters.
    pause
        How long to pause in seconds before request. If None, will query API
        status endpoint to find when next slot is available.
    error_pause
        How long to pause in seconds (in addition to `pause`) before re-trying
        request if error.

    Returns
    -------
    response_json
    """
    # resolve url to same IP even if there is server round-robin redirecting
    _http._config_dns(settings.overpass_url)

    # prepare the Overpass API URL and see if request already exists in cache
    url = settings.overpass_url.rstrip("/") + "/interpreter"
    prepared_url = str(requests.Request("GET", url, params=data).prepare().url)
    cached_response_json = _http._retrieve_from_cache(prepared_url)
    if isinstance(cached_response_json, dict):
        return cached_response_json

    # pause then request this URL
    if pause is None:
        this_pause = _get_overpass_pause(settings.overpass_url)
    hostname = _http._hostname_from_url(url)
    msg = f"Pausing {this_pause} second(s) before making HTTP POST request to {hostname!r}"
    utils.log(msg, level=lg.INFO)
    time.sleep(this_pause)

    # transmit the HTTP POST request
    msg = f"Post {prepared_url} with timeout={settings.requests_timeout}"
    utils.log(msg, level=lg.INFO)
    response = requests.post(
        url,
        data=data,
        timeout=settings.requests_timeout,
        headers=_http._get_http_headers(),
        **settings.requests_kwargs,
    )

    # handle 429 and 504 errors by pausing then recursively re-trying request
    if response.status_code in {429, 504}:  # pragma: no cover
        this_pause = error_pause + _get_overpass_pause(settings.overpass_url)
        msg = (
            f"{hostname!r} responded {response.status_code} {response.reason}: "
            f"we'll retry in {this_pause} secs"
        )
        utils.log(msg, level=lg.WARNING)
        time.sleep(this_pause)
        return _overpass_request(data, pause=pause, error_pause=error_pause)

    response_json = _http._parse_response(response)
    if not isinstance(response_json, dict):  # pragma: no cover
        msg = "Overpass API did not return a dict of results."
        raise InsufficientResponseError(msg)
    _http._save_to_cache(prepared_url, response_json, response.ok)
    return response_json
