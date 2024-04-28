"""Tools to work with the Overpass API."""

import datetime as dt
import logging as lg
import time
from warnings import warn

import numpy as np
import requests
from requests.exceptions import ConnectionError

from . import _downloader
from . import projection
from . import settings
from . import utils
from . import utils_geo


def _get_osm_filter(network_type):
    """
    Create a filter to query OSM for the specified network type.

    Parameters
    ----------
    network_type : string {"all", "all_public", "bike", "drive", "drive_service", "walk"}
        what type of street network to get

    Returns
    -------
    string
    """
    if network_type == "all_private":
        network_type = "all"
        msg = (
            "The 'all_private' network type has been renamed 'all'. The old "
            "'all_private' naming is deprecated and will be removed in the v2.0.0 release. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
        )
        warn(msg, FutureWarning, stacklevel=2)

    # define built-in queries to send to the API. specifying way["highway"]
    # means that all ways returned must have a highway tag. the filters then
    # remove ways by tag/value.
    filters = {}

    # driving: filter out un-drivable roads, service roads, private ways, and
    # anything specifying motor=no. also filter out any non-service roads that
    # are tagged as providing certain services
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
    # specifying foot=no. allow service roads, permitting things like parking
    # lot lanes, alleys, etc that you *can* walk on even if they're not
    # exactly pleasant walks. some cycleways may allow pedestrians, but this
    # filter ignores such cycleways.
    filters["walk"] = (
        f'["highway"]["area"!~"yes"]{settings.default_access}'
        f'["highway"!~"abandoned|bus_guideway|construction|cycleway|motor|no|planned|platform|'
        f'proposed|raceway|razed"]'
        f'["foot"!~"no"]["service"!~"private"]'
    )

    # biking: filter out foot ways, motor ways, private ways, and anything
    # specifying biking=no
    filters["bike"] = (
        f'["highway"]["area"!~"yes"]{settings.default_access}'
        f'["highway"!~"abandoned|bus_guideway|construction|corridor|elevator|escalator|footway|'
        f'motor|no|planned|platform|proposed|raceway|razed|steps"]'
        f'["bicycle"!~"no"]["service"!~"private"]'
    )

    # to download all ways, just filter out everything not currently in use or
    # that is private-access only
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
        osm_filter = filters[network_type]
    else:  # pragma: no cover
        msg = f"Unrecognized network_type {network_type!r}"
        raise ValueError(msg)

    return osm_filter


def _get_overpass_pause(base_endpoint, recursive_delay=5, default_duration=60):
    """
    Retrieve a pause duration from the Overpass API status endpoint.

    Check the Overpass API status endpoint to determine how long to wait until
    the next slot is available. You can disable this via the `settings`
    module's `overpass_rate_limit` setting.

    Parameters
    ----------
    base_endpoint : string
        base Overpass API url (without "/status" at the end)
    recursive_delay : int
        how long to wait between recursive calls if the server is currently
        running a query
    default_duration : int
        if fatal error, fall back on returning this value

    Returns
    -------
    pause : int
    """
    if settings.timeout is None:
        timeout = settings.requests_timeout
    else:
        timeout = settings.timeout
        msg = (
            "`settings.timeout` is deprecated and will be removed in the "
            "v2.0.0 release: use `settings.requests_timeout` instead. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
        )
        warn(msg, FutureWarning, stacklevel=2)

    if not settings.overpass_rate_limit:
        # if overpass rate limiting is False, then there is zero pause
        return 0

    try:
        url = base_endpoint.rstrip("/") + "/status"
        response = requests.get(
            url,
            headers=_downloader._get_http_headers(),
            timeout=timeout,
            **settings.requests_kwargs,
        )
        status = response.text.split("\n")[4]
        status_first_token = status.split(" ")[0]
    except ConnectionError:  # pragma: no cover
        # cannot reach status endpoint, log error and return default duration
        utils.log(f"Unable to query {url}, got status {response.status_code}", level=lg.ERROR)
        return default_duration
    except (AttributeError, IndexError, ValueError):  # pragma: no cover
        # cannot parse output, log error and return default duration
        utils.log(f"Unable to parse {url} response: {response.text}", level=lg.ERROR)
        return default_duration

    try:
        # if first token is numeric, it's how many slots you have available,
        # no wait required
        _ = int(status_first_token)  # number of available slots
        pause = 0

    except ValueError:  # pragma: no cover
        # if first token is 'Slot', it tells you when your slot will be free
        if status_first_token == "Slot":
            utc_time_str = status.split(" ")[3]
            pattern = "%Y-%m-%dT%H:%M:%SZ,"
            utc_time = dt.datetime.strptime(utc_time_str, pattern).astimezone(dt.timezone.utc)
            utc_now = dt.datetime.now(tz=dt.timezone.utc)
            seconds = int(np.ceil((utc_time - utc_now).total_seconds()))
            pause = max(seconds, 1)

        # if first token is 'Currently', it is currently running a query so
        # check back in recursive_delay seconds
        elif status_first_token == "Currently":
            time.sleep(recursive_delay)
            pause = _get_overpass_pause(base_endpoint)

        # any other status is unrecognized: log error, return default duration
        else:
            utils.log(f"Unrecognized server status: {status!r}", level=lg.ERROR)
            return default_duration

    return pause


def _make_overpass_settings():
    """
    Make settings string to send in Overpass query.

    Returns
    -------
    string
    """
    if settings.timeout is None:
        timeout = settings.requests_timeout
    else:
        timeout = settings.timeout
        msg = (
            "`settings.timeout` is deprecated and will be removed in the "
            "v2.0.0 release: use `settings.requests_timeout` instead. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
        )
        warn(msg, FutureWarning, stacklevel=2)

    if settings.memory is None:
        memory = settings.overpass_memory
    else:
        memory = settings.memory
        msg = (
            "`settings.memory` is deprecated and will be removed in the "
            " v2.0.0 release: use `settings.overpass_memory` instead. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
        )
        warn(msg, FutureWarning, stacklevel=2)

    maxsize = "" if memory is None else f"[maxsize:{memory}]"
    return settings.overpass_settings.format(timeout=timeout, maxsize=maxsize)


def _make_overpass_polygon_coord_strs(polygon):
    """
    Subdivide query polygon and return list of coordinate strings.

    Project to utm, divide polygon up into sub-polygons if area exceeds a
    max size (in meters), project back to lat-lon, then get a list of
    polygon(s) exterior coordinates. Ignore interior ("holes") coordinates.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        polygon to convert to exterior coordinate strings

    Returns
    -------
    coord_strs : list
        list of strings of exterior coordinates of polygon(s)
    """
    # first subdivide the polygon if its area exceeds max size
    # this results in a multipolygon of 1+ constituent polygons
    poly_proj, crs_proj = projection.project_geometry(polygon)
    multi_poly_proj = utils_geo._consolidate_subdivide_geometry(poly_proj)
    multi_poly, _ = projection.project_geometry(multi_poly_proj, crs=crs_proj, to_latlong=True)

    # then extract each's exterior coords to the string format Overpass
    # expects, rounding lats and lons to 6 decimals (ie, ~100 mm) so we
    # can hash and cache URL strings consistently
    coord_strs = []
    for geom in multi_poly.geoms:
        x, y = geom.exterior.xy
        coord_list = [f'{xy[1]:.6f}{" "}{xy[0]:.6f}' for xy in zip(x, y)]
        coord_strs.append(" ".join(coord_list))

    return coord_strs


def _create_overpass_query(polygon_coord_str, tags):
    """
    Create an Overpass features query string based on passed tags.

    Parameters
    ----------
    polygon_coord_str : list
        list of lat lon coordinates
    tags : dict
        dict of tags used for finding elements in the search area

    Returns
    -------
    query : string
    """
    # create overpass settings string
    overpass_settings = _make_overpass_settings()

    # make sure every value in dict is bool, str, or list of str
    err_msg = "tags must be a dict with values of bool, str, or list of str"
    if not isinstance(tags, dict):  # pragma: no cover
        raise TypeError(err_msg)

    tags_dict = {}
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
    tags_list = []
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
    components = "".join(components)
    return f"{overpass_settings};({components});out;"


def _download_overpass_network(polygon, network_type, custom_filter):
    """
    Retrieve networked ways and nodes within boundary from the Overpass API.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        boundary to fetch the network ways/nodes within
    network_type : string
        what type of street network to get if custom_filter is None
    custom_filter : string
        a custom "ways" filter to be used instead of the network_type presets

    Yields
    ------
    response_json : dict
        a generator of JSON responses from the Overpass server
    """
    # create a filter to exclude certain kinds of ways based on the requested
    # network_type, if provided, otherwise use custom_filter
    osm_filter = custom_filter if custom_filter is not None else _get_osm_filter(network_type)

    # create overpass settings string
    overpass_settings = _make_overpass_settings()

    # subdivide query polygon to get list of sub-divided polygon coord strings
    polygon_coord_strs = _make_overpass_polygon_coord_strs(polygon)
    utils.log(f"Requesting data from API in {len(polygon_coord_strs)} request(s)")

    # pass exterior coordinates of each polygon in list to API, one at a time
    # the '>' makes it recurse so we get ways and the ways' nodes.
    for polygon_coord_str in polygon_coord_strs:
        query_str = f"{overpass_settings};(way{osm_filter}(poly:{polygon_coord_str!r});>;);out;"
        yield _overpass_request(data={"data": query_str})


def _download_overpass_features(polygon, tags):
    """
    Retrieve OSM features within boundary from the Overpass API.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon
        boundaries to fetch elements within
    tags : dict
        dict of tags used for finding elements in the selected area

    Yields
    ------
    response_json : dict
        a generator of JSON responses from the Overpass server
    """
    # subdivide query polygon to get list of sub-divided polygon coord strings
    polygon_coord_strs = _make_overpass_polygon_coord_strs(polygon)
    utils.log(f"Requesting data from API in {len(polygon_coord_strs)} request(s)")

    # pass exterior coordinates of each polygon in list to API, one at a time
    for polygon_coord_str in polygon_coord_strs:
        query_str = _create_overpass_query(polygon_coord_str, tags)
        yield _overpass_request(data={"data": query_str})


def _overpass_request(data, pause=None, error_pause=60):
    """
    Send a HTTP POST request to the Overpass API and return response.

    Parameters
    ----------
    data : OrderedDict
        key-value pairs of parameters
    pause : float
        how long to pause in seconds before request, if None, will query API
        status endpoint to find when next slot is available
    error_pause : float
        how long to pause in seconds (in addition to `pause`) before re-trying
        request if error

    Returns
    -------
    response_json : dict
    """
    if settings.timeout is None:
        timeout = settings.requests_timeout
    else:
        timeout = settings.timeout
        msg = (
            "`settings.timeout` is deprecated and will be removed in the "
            "v2.0.0 release: use `settings.requests_timeout` instead. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
        )
        warn(msg, FutureWarning, stacklevel=2)

    if settings.overpass_endpoint is None:
        overpass_endpoint = settings.overpass_url
    else:
        overpass_endpoint = settings.overpass_endpoint
        msg = (
            "`settings.overpass_endpoint` is deprecated and will be removed in the "
            "v2.0.0 release: use `settings.overpass_url` instead. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
        )
        warn(msg, FutureWarning, stacklevel=2)

    # resolve url to same IP even if there is server round-robin redirecting
    _downloader._config_dns(overpass_endpoint)

    # prepare the Overpass API URL and see if request already exists in cache
    url = overpass_endpoint.rstrip("/") + "/interpreter"
    prepared_url = requests.Request("GET", url, params=data).prepare().url
    cached_response_json = _downloader._retrieve_from_cache(prepared_url)
    if cached_response_json is not None:
        return cached_response_json

    # pause then request this URL
    if pause is None:
        this_pause = _get_overpass_pause(overpass_endpoint)
    domain = _downloader._hostname_from_url(url)
    utils.log(f"Pausing {this_pause} second(s) before making HTTP POST request to {domain!r}")
    time.sleep(this_pause)

    # transmit the HTTP POST request
    utils.log(f"Post {prepared_url} with timeout={timeout}")
    response = requests.post(
        url,
        data=data,
        timeout=timeout,
        headers=_downloader._get_http_headers(),
        **settings.requests_kwargs,
    )

    # handle 429 and 504 errors by pausing then recursively re-trying request
    if response.status_code in {429, 504}:  # pragma: no cover
        this_pause = error_pause + _get_overpass_pause(overpass_endpoint)
        msg = (
            f"{domain!r} responded {response.status_code} {response.reason}: "
            f"we'll retry in {this_pause} secs"
        )
        utils.log(msg, level=lg.WARNING)
        time.sleep(this_pause)
        return _overpass_request(data, pause, error_pause)

    response_json = _downloader._parse_response(response)
    _downloader._save_to_cache(prepared_url, response_json, response.status_code)
    return response_json
