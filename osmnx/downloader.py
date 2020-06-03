"""Interact with the OSM APIs."""

import datetime as dt
import hashlib
import json
import logging as lg
import math
import os
import re
import time
from collections import OrderedDict

import requests
from dateutil import parser as date_parser
from shapely.geometry import Polygon

from . import projection
from . import settings
from . import utils
from . import utils_geo
from ._errors import InsufficientNetworkQueryArguments
from ._errors import UnknownNetworkType


def _get_osm_filter(network_type):
    """
    Create a filter to query OSM for the specified network type.

    Parameters
    ----------
    network_type : string
        {'walk', 'bike', 'drive', 'drive_service', 'all', 'all_private'}
        what type of street or other network to get

    Returns
    -------
    string
    """
    # define preset queries to send to the API. specifying way["highway"]
    # means that all ways returned must have a highway tag. the filters then
    # remove ways by tag/value. the '>' makes it recurse so we get ways and
    # the ways' nodes.
    filters = {}

    # driving: filter out un-drivable roads, service roads, private ways, and
    # anything specifying motor=no. also filter out any non-service roads that
    # are tagged as providing parking, driveway, private, or emergency-access
    # services
    filters["drive"] = (
        f'["highway"]["area"!~"yes"]["highway"!~"cycleway|footway|path|pedestrian|steps|track|corridor|'
        f'elevator|escalator|proposed|construction|bridleway|abandoned|platform|raceway|service"]'
        f'["motor_vehicle"!~"no"]["motorcar"!~"no"]{settings.default_access}'
        f'["service"!~"parking|parking_aisle|driveway|private|emergency_access"]'
    )

    # drive+service: allow ways tagged 'service' but filter out certain types of
    # service ways
    filters["drive_service"] = (
        f'["highway"]["area"!~"yes"]["highway"!~"cycleway|footway|path|pedestrian|steps|track|'
        f'corridor|elevator|escalator|proposed|construction|bridleway|abandoned|platform|raceway"]'
        f'["motor_vehicle"!~"no"]["motorcar"!~"no"]{settings.default_access}'
        f'["service"!~"parking|parking_aisle|private|emergency_access"]'
    )

    # walking: filter out cycle ways, motor ways, private ways, and anything
    # specifying foot=no. allow service roads, permitting things like parking
    # lot lanes, alleys, etc that you *can* walk on even if they're not exactly
    # pleasant walks. some cycleways may allow pedestrians, but this filter ignores
    # such cycleways.
    filters["walk"] = (
        f'["highway"]["area"!~"yes"]["highway"!~"cycleway|motor|proposed|construction|abandoned|'
        f'platform|raceway"]["foot"!~"no"]["service"!~"private"]{settings.default_access}'
    )

    # biking: filter out foot ways, motor ways, private ways, and anything
    # specifying biking=no
    filters["bike"] = (
        f'["highway"]["area"!~"yes"]["highway"!~"footway|steps|corridor|elevator|'
        f'escalator|motor|proposed|construction|abandoned|platform|raceway"]'
        f'["bicycle"!~"no"]["service"!~"private"]{settings.default_access}'
    )

    # to download all ways, just filter out everything not currently in use or
    # that is private-access only
    filters["all"] = (
        f'["highway"]["area"!~"yes"]["highway"!~"proposed|construction|abandoned|platform|raceway"]'
        f'["service"!~"private"]{settings.default_access}'
    )

    # to download all ways, including private-access ones, just filter out
    # everything not currently in use
    filters[
        "all_private"
    ] = '["highway"]["area"!~"yes"]["highway"!~"proposed|construction|abandoned|platform|raceway"]'

    if network_type in filters:
        osm_filter = filters[network_type]
    else:
        raise UnknownNetworkType(f'unknown network_type "{network_type}"')

    return osm_filter


def _save_to_cache(url, response_json):
    """
    Save an HTTP response json object to the cache.

    If the request was sent to server via POST instead of GET, then URL should
    be a GET-style representation of request. Users should always pass
    OrderedDicts instead of dicts of parameters into request functions, so that
    the parameters stay in the same order each time, producing the same URL
    string, and thus the same hash. Otherwise the cache will eventually contain
    multiple saved responses for the same request because the URL's parameters
    appeared in a different order each time.

    Parameters
    ----------
    url : string
        the url of the request
    response_json : dict
        the json response

    Returns
    -------
    None
    """
    if settings.use_cache:
        if response_json is None:
            utils.log("Did not save to cache because response_json is None")
        else:
            # create the folder on the disk if it doesn't already exist
            if not os.path.exists(settings.cache_folder):
                os.makedirs(settings.cache_folder)

            # hash the url to make the filename succinct but unique
            filename = hashlib.md5(url.encode("utf-8")).hexdigest()
            cache_filepath = os.path.join(settings.cache_folder, os.extsep.join([filename, "json"]))

            # dump to json, and save to file
            json_str = str(json.dumps(response_json))
            with open(cache_filepath, "w", encoding="utf-8") as cache_file:
                cache_file.write(json_str)

            utils.log(f'Saved response to cache file "{cache_filepath}"')


def _url_in_cache(url):
    """
    Determine if a URL's response exists in the cache.

    Parameters
    ----------
    url : string
        the url to look for in the cache

    Returns
    -------
    filepath : string
        path to cached response for url if it exists, otherwise None
    """
    # hash the url to generate the cache filename
    filename = hashlib.md5(url.encode("utf-8")).hexdigest()
    filepath = os.path.join(settings.cache_folder, os.extsep.join([filename, "json"]))

    # if this file exists in the cache, return its full path
    if os.path.isfile(filepath):
        return filepath
    else:
        return None


def _get_from_cache(url, check_remark=False):
    """
    Retrieve a HTTP response json object from the cache.

    Parameters
    ----------
    url : string
        the url of the request
    check_remark : string
        if True, only return filepath if cached response does not have a
        remark key indicating a server warning

    Returns
    -------
    response_json : dict
        cached response for url if it exists in the cache, otherwise None
    """
    # if the tool is configured to use the cache
    if settings.use_cache:

        # return cached response for this url if it exists, otherwise return None
        cache_filepath = _url_in_cache(url)
        if cache_filepath is not None:
            with open(cache_filepath, encoding="utf-8") as cache_file:
                response_json = json.load(cache_file)

                # return None if check_remark is True and there is a server
                # remark in the cached response
                if check_remark and "remark" in response_json:
                    utils.log(f'Found remark, so ignoring cache file "{cache_filepath}"')
                    return None

            utils.log(
                f'Retrieved response from cache file "{cache_filepath}" for URL "{url[:200]}..."'
            )
            return response_json


def _get_http_headers(user_agent=None, referer=None, accept_language=None):
    """
    Update the default requests HTTP headers with OSMnx info.

    Parameters
    ----------
    user_agent : string
        the user agent string, if None will set with OSMnx default
    referer : string
        the referer string, if None will set with OSMnx default
    accept_language : string
        make accept-language explicit e.g. for consistent nominatim result sorting

    Returns
    -------
    headers : dict
    """
    if user_agent is None:
        user_agent = settings.default_user_agent
    if referer is None:
        referer = settings.default_referer
    if accept_language is None:
        accept_language = settings.default_accept_language

    headers = requests.utils.default_headers()
    headers.update(
        {"User-Agent": user_agent, "referer": referer, "Accept-Language": accept_language}
    )
    return headers


def _get_pause(recursive_delay=5, default_duration=10):
    """
    Get a pause duration from the Overpass API status endpoint.

    Check the Overpass API status endpoint to determine how long to wait until
    next slot is available.

    Parameters
    ----------
    recursive_delay : int
        how long to wait between recursive calls if server is currently running
        a query
    default_duration : int
        if fatal error, function falls back on returning this value

    Returns
    -------
    int
    """
    try:
        url = settings.overpass_endpoint.rstrip("/") + "/status"
        response = requests.get(url, headers=_get_http_headers())
        status = response.text.split("\n")[3]
        status_first_token = status.split(" ")[0]
    # if we cannot reach the status endpoint or parse its output, log an
    # error and return default duration
    except Exception:
        utils.log(f"Unable to query {url}", level=lg.ERROR)
        return default_duration

    try:
        # if first token is numeric, it's how many slots you have available - no
        # wait required
        _ = int(status_first_token)  # number of available slots
        pause = 0
    except Exception:  # pragma: no cover
        # if first token is 'Slot', it tells you when your slot will be free
        if status_first_token == "Slot":
            utc_time_str = status.split(" ")[3]
            utc_time = date_parser.parse(utc_time_str).replace(tzinfo=None)
            pause = math.ceil((utc_time - dt.datetime.utcnow()).total_seconds())
            pause = max(pause, 1)

        # if first token is 'Currently', it is currently running a query so
        # check back in recursive_delay seconds. any other status is unrecognized,
        # log an error and return default duration
        elif status_first_token == "Currently":
            time.sleep(recursive_delay)
            pause = _get_pause()
        else:
            utils.log(f'Unrecognized server status: "{status}"', level=lg.ERROR)
            return default_duration

    return pause


def _make_overpass_settings(custom_settings, timeout, memory):
    """
    Make settings string to send in Overpass query.

    Parameters
    ----------
    custom_settings : string
        custom settings to be used in the overpass query instead of defaults
    timeout : int
        the timeout interval for requests and to pass to API
    memory : int
        server memory allocation size for the query, in bytes. If none, server
        will use its default allocation size

    Returns
    -------
    overpass_settings : string
    """
    # pass server memory allocation in bytes for the query to the API
    # if None, pass nothing so the server will use its default allocation size
    # otherwise, define the query's maxsize parameter value as whatever the
    # caller passed in
    if memory is None:
        maxsize = ""
    else:
        maxsize = f"[maxsize:{memory}]"

    # use custom settings if delivered, otherwise just the default ones
    if custom_settings:
        overpass_settings = custom_settings
    else:
        overpass_settings = settings.default_overpass_query_settings.format(
            timeout=timeout, maxsize=maxsize
        )

    return overpass_settings


def _osm_net_download(
    polygon=None,
    north=None,
    south=None,
    east=None,
    west=None,
    network_type="all_private",
    timeout=180,
    memory=None,
    max_query_area_size=50 * 1000 * 50 * 1000,
    infrastructure=None,
    custom_filter=None,
    custom_settings=None,
):
    """
    Download OSM ways and nodes within some bounding box from the Overpass API.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        geographic boundaries to fetch the street network within
    north : float
        northern latitude of bounding box
    south : float
        southern latitude of bounding box
    east : float
        eastern longitude of bounding box
    west : float
        western longitude of bounding box
    network_type : string
        {'walk', 'bike', 'drive', 'drive_service', 'all', 'all_private'} what
        type of street network to get
    timeout : int
        the timeout interval for requests and to pass to API
    memory : int
        server memory allocation size for the query, in bytes. If none, server
        will use its default allocation size
    max_query_area_size : float
        max area for any part of the geometry in meters: any polygon bigger
        will get divided up for multiple queries to API (default 50km x 50km)
    infrastructure : string
        deprecated, use custom_filter instead
    custom_filter : string
        a custom network filter to be used instead of the network_type presets,
        e.g., '["power"~"line"]' or '["highway"~"motorway|trunk"]'
    custom_settings : string
        custom settings to be used in the overpass query instead of defaults

    Returns
    -------
    response_jsons : list
    """
    if infrastructure is not None:
        from warnings import warn

        msg = (
            "The `infrastructure` parameter has been deprecated and will be "
            "removed in the next release. Use the `custom_filter` parameter instead."
        )
        warn(msg)

    # check if we're querying by polygon or by bounding box based on which
    # argument(s) where passed into this function
    by_poly = polygon is not None
    by_bbox = not (north is None or south is None or east is None or west is None)
    if not (by_poly or by_bbox):
        raise InsufficientNetworkQueryArguments(
            "You must pass a polygon or north, south, east, and west"
        )

    # create a filter to exclude certain kinds of ways based on the requested
    # network_type
    if custom_filter is not None:
        osm_filter = custom_filter
    else:
        osm_filter = _get_osm_filter(network_type)
    response_jsons = []

    overpass_settings = _make_overpass_settings(custom_settings, timeout, memory)

    if by_bbox:
        # turn bbox into a polygon and project to local UTM
        polygon = Polygon([(west, south), (east, south), (east, north), (west, north)])
        geometry_proj, crs_proj = projection.project_geometry(polygon)

        # subdivide it if it exceeds the max area size (in meters), then project
        # back to lat-lng
        gpcs = utils_geo._consolidate_subdivide_geometry(
            geometry_proj, max_query_area_size=max_query_area_size
        )
        geometry, _ = projection.project_geometry(gpcs, crs=crs_proj, to_latlong=True)
        utils.log(
            f"Requesting network data within bounding box from API in {len(geometry)} request(s)"
        )

        # loop through each polygon rectangle in the geometry (there will only
        # be one if original bbox didn't exceed max area size)
        for poly in geometry:
            # represent bbox as south,west,north,east and round lat-lngs to 6
            # decimal places (ie, ~100 mm) so URL strings aren't different
            # due to float rounding issues (for consistent caching)
            west, south, east, north = poly.bounds
            query_str = (
                f"{overpass_settings};"
                f"(way{osm_filter}({south:.6f},{west:.6f},{north:.6f},{east:.6f});>;);out;"
            )
            response_json = overpass_request(data={"data": query_str}, timeout=timeout)
            response_jsons.append(response_json)
        utils.log(
            f"Got all network data within bounding box from API in {len(geometry)} request(s)"
        )

    elif by_poly:
        # project to utm, divide polygon up into sub-polygons if area exceeds a
        # max size (in meters), project back to lat-lng, then get a list of
        # polygon(s) exterior coordinates
        geometry_proj, crs_proj = projection.project_geometry(polygon)
        gpcs = utils_geo._consolidate_subdivide_geometry(
            geometry_proj, max_query_area_size=max_query_area_size
        )
        geometry, _ = projection.project_geometry(gpcs, crs=crs_proj, to_latlong=True)
        polygon_coord_strs = utils_geo._get_polygons_coordinates(geometry)
        utils.log(
            f"Requesting network data within polygon from API in {len(polygon_coord_strs)} request(s)"
        )

        # pass each polygon exterior coordinates in the list to the API, one at
        # a time
        for polygon_coord_str in polygon_coord_strs:
            query_str = f"{overpass_settings};(way{osm_filter}(poly:'{polygon_coord_str}');>;);out;"
            response_json = overpass_request(data={"data": query_str}, timeout=timeout)
            response_jsons.append(response_json)
        utils.log(
            f"Got all network data within polygon from API in {len(polygon_coord_strs)} request(s)"
        )

    return response_jsons


def _osm_polygon_download(query, limit=1, polygon_geojson=1):
    """
    Geocode a place and download its boundary geometry from OSM's Nominatim API.

    Parameters
    ----------
    query : string or dict
        query string or structured query dict to geocode/download
    limit : int
        max number of results to return
    polygon_geojson : int
        request the boundary geometry polygon from the API, 0=no, 1=yes

    Returns
    -------
    response_json : dict
    """
    # define the parameters
    params = OrderedDict()
    params["format"] = "json"
    params["limit"] = limit
    params[
        "dedupe"
    ] = 0  # prevent OSM from deduping results so we get precisely 'limit' # of results
    params["polygon_geojson"] = polygon_geojson

    # add the structured query dict (if provided) to params, otherwise query
    # with place name string
    if isinstance(query, str):
        params["q"] = query
    elif isinstance(query, dict):
        # add the query keys in alphabetical order so the URL is the same string
        # each time, for caching purposes
        for key in sorted(list(query.keys())):
            params[key] = query[key]
    else:
        raise TypeError("query must be a dict or a string")

    # request the URL, return the JSON
    response_json = nominatim_request(params=params, timeout=30)
    return response_json


def nominatim_request(params, request_type="search", pause=1, timeout=30, error_pause=180):
    """
    Send a request to the Nominatim API via HTTP GET and return JSON response.

    Parameters
    ----------
    params : dict or OrderedDict
        key-value pairs of parameters
    request_type : string
        Type of Nominatim query. One of: search, reverse, or lookup
    pause : int
        how long to pause before requests, in seconds
    timeout : int
        the timeout interval for the requests library
    error_pause : int
        how long to pause in seconds before re-trying requests if error

    Returns
    -------
    response_json : dict
    """
    known_requests = {"search", "reverse", "lookup"}
    if request_type not in known_requests:
        raise ValueError('Nominatim request_type must be "search", "reverse", or "lookup"')

    # prepare the Nominatim API URL and see if request already exists in the
    # cache
    url = settings.nominatim_endpoint.rstrip("/") + "/" + request_type
    prepared_url = requests.Request("GET", url, params=params).prepare().url
    cached_response_json = _get_from_cache(prepared_url)

    if settings.nominatim_key:
        params["key"] = settings.nominatim_key

    if cached_response_json is not None:
        # found this request in the cache, just return it instead of making a
        # new HTTP call
        return cached_response_json

    else:
        # if this URL is not already in the cache, pause, then request it
        utils.log(f"Pausing {pause} seconds before making HTTP GET request")
        time.sleep(pause)
        utils.log(f"Get {prepared_url} with timeout={timeout}")
        response = requests.get(url, params=params, timeout=timeout, headers=_get_http_headers())

        # get the response size and the domain, log result
        size_kb = len(response.content) / 1000.0
        domain = re.findall(r"(?s)//(.*?)/", url)[0]
        utils.log(f"Downloaded {size_kb:,.1f}KB from {domain}")

        try:
            response_json = response.json()
            _save_to_cache(prepared_url, response_json)
        except Exception:  # pragma: no cover
            # 429 is 'too many requests' and 504 is 'gateway timeout' from server
            # overload - handle these errors by recursively calling
            # nominatim_request until we get a valid response
            sc = response.status_code
            if sc in [429, 504]:
                # pause for error_pause seconds before re-trying request
                utils.log(
                    f"{domain} returned {sc} and no data: retrying in {error_pause:.2f} secs",
                    level=lg.WARNING,
                )
                time.sleep(error_pause)
                response_json = nominatim_request(params=params, pause=pause, timeout=timeout)

            # else, this was an unhandled status_code, throw an exception
            else:
                utils.log(f"{domain} returned {sc} and no data", level=lg.ERROR)
                raise Exception(
                    f"Server returned no JSON data\n{response} {response.reason}\n{response.text}"
                )

        return response_json


def overpass_request(data, pause=None, timeout=180, error_pause=None):
    """
    Send a request to the Overpass API via HTTP POST and return JSON response.

    Parameters
    ----------
    data : dict or OrderedDict
        key-value pairs of parameters to post to the API
    pause : int
        how long to pause in seconds before requests, if None, will query API
        status endpoint to find when next slot is available
    timeout : int
        the timeout interval for the requests library
    error_pause : int
        how long to pause in seconds before re-trying requests if error

    Returns
    -------
    dict
    """
    # define the Overpass API URL, then construct a GET-style URL as a string to
    # hash to look up/save to cache
    url = settings.overpass_endpoint.rstrip("/") + "/interpreter"
    prepared_url = requests.Request("GET", url, params=data).prepare().url
    cached_response_json = _get_from_cache(prepared_url, check_remark=True)

    if cached_response_json is not None:
        # found this request in the cache, just return it instead of making a
        # new HTTP call
        return cached_response_json

    else:
        # if this URL is not already in the cache, pause, then request it
        if pause is None:
            this_pause = _get_pause()
        utils.log(f"Pausing {this_pause} seconds before making HTTP POST request")
        time.sleep(this_pause)
        utils.log(f"Post {prepared_url} with timeout={timeout}")
        response = requests.post(url, data=data, timeout=timeout, headers=_get_http_headers())

        # get the response size and the domain, log result
        size_kb = len(response.content) / 1000.0
        domain = re.findall(r"(?s)//(.*?)/", url)[0]
        utils.log(f"Downloaded {size_kb:,.1f}KB from {domain}")

        try:
            response_json = response.json()
            if "remark" in response_json:
                utils.log(f'Server remark: "{response_json["remark"]}"', level=lg.WARNING)
            _save_to_cache(prepared_url, response_json)

        # 429 is 'too many requests' and 504 is 'gateway timeout' from server
        # overload - handle these errors by recursively calling overpass_request
        # after a pause until we get a valid response
        except Exception:  # pragma: no cover
            sc = response.status_code
            if sc in [429, 504]:
                if error_pause is None:
                    error_pause = _get_pause()
                utils.log(
                    f"{domain} returned {sc} and no data: retry in {error_pause} secs.",
                    level=lg.WARNING,
                )
                time.sleep(error_pause)
                response_json = overpass_request(data=data, pause=pause, timeout=timeout)
            # else, this was an unhandled status_code, throw an exception
            else:
                utils.log(f"{domain} returned {sc} and no data", level=lg.ERROR)
                raise Exception(
                    f"Server returned no JSON data\n{response} {response.reason}\n{response.text}"
                )

        return response_json
