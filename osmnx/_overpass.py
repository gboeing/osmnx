"""Tools to work with the Overpass API."""

import datetime as dt
import logging as lg
import socket
import time

import numpy as np
import requests
from requests.exceptions import ConnectionError

from . import _downloader
from . import projection
from . import settings
from . import utils
from . import utils_geo

# capture getaddrinfo function to use original later after mutating it
_original_getaddrinfo = socket.getaddrinfo


def _get_osm_filter(network_type):
    """
    Create a filter to query OSM for the specified network type.

    Parameters
    ----------
    network_type : string {"all_private", "all", "bike", "drive", "drive_service", "walk"}
        what type of street network to get

    Returns
    -------
    string
    """
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
    filters["all"] = (
        f'["highway"]["area"!~"yes"]{settings.default_access}'
        f'["highway"!~"abandoned|construction|no|planned|platform|proposed|raceway|razed"]'
        f'["service"!~"private"]'
    )

    # to download all ways, including private-access ones, just filter out
    # everything not currently in use
    filters["all_private"] = (
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
    if not settings.overpass_rate_limit:
        # if overpass rate limiting is False, then there is zero pause
        return 0

    sc = None

    try:
        url = base_endpoint.rstrip("/") + "/status"
        response = requests.get(
            url,
            headers=_downloader._get_http_headers(),
            timeout=settings.timeout,
            **settings.requests_kwargs,
        )
        sc = response.status_code
        status = response.text.split("\n")[4]
        status_first_token = status.split(" ")[0]
    except ConnectionError:  # pragma: no cover
        # cannot reach status endpoint, log error and return default duration
        utils.log(f"Unable to query {url}, got status {sc}", level=lg.ERROR)
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
    if settings.memory is None:
        maxsize = ""
    else:
        maxsize = f"[maxsize:{settings.memory}]"
    return settings.overpass_settings.format(timeout=settings.timeout, maxsize=maxsize)


def _make_overpass_polygon_coord_strs(polygon):
    """
    Subdivide query polygon and return list of coordinate strings.

    Project to utm, divide polygon up into sub-polygons if area exceeds a
    max size (in meters), project back to lat-lng, then get a list of
    polygon(s) exterior coordinates

    Parameters
    ----------
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon

    Returns
    -------
    polygon_coord_strs : list
        list of exterior coordinate strings for smaller sub-divided polygons
    """
    geometry_proj, crs_proj = projection.project_geometry(polygon)
    gpcs = utils_geo._consolidate_subdivide_geometry(geometry_proj)
    geometry, _ = projection.project_geometry(gpcs, crs=crs_proj, to_latlong=True)
    return utils_geo._get_polygons_coordinates(geometry)


def _create_overpass_query(polygon_coord_str, tags):
    """
    Create an overpass query string based on passed tags.

    Parameters
    ----------
    polygon_coord_str : list
        list of lat lng coordinates
    tags : dict
        dict of tags used for finding elements in the selected area

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
                tags_list.append({key: value_item})

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
                components.append(f"({kind}{tag_str});")

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
        a custom ways filter to be used instead of the network_type presets

    Yields
    ------
    response_json : dict
        a generator of JSON responses from the Overpass server
    """
    # create a filter to exclude certain kinds of ways based on the requested
    # network_type, if provided, otherwise use custom_filter
    if custom_filter is not None:
        osm_filter = custom_filter
    else:
        osm_filter = _get_osm_filter(network_type)

    # create overpass settings string
    overpass_settings = _make_overpass_settings()

    # subdivide query polygon to get list of sub-divided polygon coord strings
    polygon_coord_strs = _make_overpass_polygon_coord_strs(polygon)
    utils.log(f"Requesting data from API in {len(polygon_coord_strs)} request(s)")

    # pass each polygon exterior coordinates in the list to the API, one at a
    # time. The '>' makes it recurse so we get ways and the ways' nodes.
    for polygon_coord_str in polygon_coord_strs:
        query_str = f"{overpass_settings};(way{osm_filter}(poly:{polygon_coord_str!r});>;);out;"
        yield _downloader._overpass_request(data={"data": query_str})


def _download_overpass_features(polygon, tags):
    """
    Retrieve OSM features within boundary from the Overpass API.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon
        boundaries to fetch elements within
    tags : dict
        dict of tags used for finding elements in the selected area

    Returns
    -------
    response_jsons : list
        list of JSON responses from the Overpass server
    """
    response_jsons = []

    # subdivide query polygon to get list of sub-divided polygon coord strings
    polygon_coord_strs = _make_overpass_polygon_coord_strs(polygon)
    utils.log(f"Requesting data from API in {len(polygon_coord_strs)} request(s)")

    # pass exterior coordinates of each polygon in list to API, one at a time
    for polygon_coord_str in polygon_coord_strs:
        query_str = _create_overpass_query(polygon_coord_str, tags)
        response_json = _downloader._overpass_request(data={"data": query_str})
        response_jsons.append(response_json)

    utils.log(
        f"Got all features data within polygon from API in {len(polygon_coord_strs)} request(s)"
    )

    return response_jsons
