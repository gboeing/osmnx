"""Interact with the OSM APIs."""

import bz2
import datetime as dt
import json
import logging as lg
import math
import os
import re
import time
import xml.sax
from collections import OrderedDict
from hashlib import md5

import requests
from dateutil import parser as date_parser

from . import projection
from . import settings
from . import utils
from . import utils_geo


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
    # remove ways by tag/value.
    filters = dict()

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
        raise ValueError(f'Unrecognized network_type "{network_type}"')

    return osm_filter


def _save_to_cache(url, response_json, sc):
    """
    Save a HTTP response JSON object to a file in the cache folder.

    Function calculates the checksum of url to generate the cache file's name.
    If the request was sent to server via POST instead of GET, then URL should
    be a GET-style representation of request. Response is only saved to a
    cache file if settings.use_cache is True, response_json is not None, and
    sc = 200.

    Users should always pass OrderedDicts instead of dicts of parameters into
    request functions, so the parameters remain in the same order each time,
    producing the same URL string, and thus the same hash. Otherwise the cache
    will eventually contain multiple saved responses for the same request
    because the URL's parameters appeared in a different order each time.

    Parameters
    ----------
    url : string
        the URL of the request
    response_json : dict
        the JSON response
    sc : int
        the response's HTTP status code

    Returns
    -------
    None
    """
    if settings.use_cache:

        if sc != 200:
            utils.log(f"Did not save to cache because status code is {sc}")

        elif response_json is None:
            utils.log("Did not save to cache because response_json is None")

        else:
            # create the folder on the disk if it doesn't already exist
            if not os.path.exists(settings.cache_folder):
                os.makedirs(settings.cache_folder)

            # hash the url to make the filename succinct but unique
            # md5 digest is 128 bits = 16 bytes = 32 hexadecimal characters
            filename = md5(url.encode("utf-8")).hexdigest()
            cache_filepath = os.path.join(settings.cache_folder, os.extsep.join([filename, "json"]))

            # dump to json, and save to file
            json_str = str(json.dumps(response_json))
            with open(cache_filepath, "w", encoding="utf-8") as cache_file:
                cache_file.write(json_str)

            utils.log(f'Saved response to cache file "{cache_filepath}"')


def _url_in_cache(url):
    """
    Determine if a URL's response exists in the cache.

    Calculates the checksum of url to determine the cache file's name.

    Parameters
    ----------
    url : string
        the URL to look for in the cache

    Returns
    -------
    filepath : string
        path to cached response for url if it exists, otherwise None
    """
    # hash the url to generate the cache filename
    filename = md5(url.encode("utf-8")).hexdigest()
    filepath = os.path.join(settings.cache_folder, os.extsep.join([filename, "json"]))

    # if this file exists in the cache, return its full path
    if os.path.isfile(filepath):
        return filepath
    else:
        return None


def _retrieve_from_cache(url, check_remark=False):
    """
    Retrieve a HTTP response JSON object from the cache, if it exists.

    Parameters
    ----------
    url : string
        the URL of the request
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

        # return cached response for this url if exists, otherwise return None
        cache_filepath = _url_in_cache(url)
        if cache_filepath is not None:
            with open(cache_filepath, encoding="utf-8") as cache_file:
                response_json = json.load(cache_file)

                # return None if check_remark is True and there is a server
                # remark in the cached response
                if check_remark and "remark" in response_json:
                    utils.log(f'Found remark, so ignoring cache file "{cache_filepath}"')
                    return None

            utils.log(f'Retrieved response from cache file "{cache_filepath}"')
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
        make accept-language explicit e.g. for consistent nominatim result
        sorting

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


def _get_pause(recursive_delay=5, default_duration=60):
    """
    Get a pause duration from the Overpass API status endpoint.

    Check the Overpass API status endpoint to determine how long to wait until
    next slot is available.

    Parameters
    ----------
    recursive_delay : int
        how long to wait between recursive calls if the server is currently
        running a query
    default_duration : int
        if fatal error, fall back on returning this value

    Returns
    -------
    pause : int
    """
    try:
        url = settings.overpass_endpoint.rstrip("/") + "/status"
        response = requests.get(url, headers=_get_http_headers())
        status = response.text.split("\n")[3]
        status_first_token = status.split(" ")[0]

    except Exception:
        # if we cannot reach the status endpoint or parse its output, log an
        # error and return default duration
        sc = response.status_code
        utils.log(f"Unable to query {url} got status {sc}", level=lg.ERROR)
        return default_duration

    try:
        # if first token is numeric, it's how many slots you have available,
        # no wait required
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
        # check back in recursive_delay seconds
        elif status_first_token == "Currently":
            time.sleep(recursive_delay)
            pause = _get_pause()

        # any other status is unrecognized: log error, return default duration
        else:
            utils.log(f'Unrecognized server status: "{status}"', level=lg.ERROR)
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
        geographic boundaries to fetch the OSM geometries within

    Returns
    -------
    polygon_coord_strs : list
        list of exterior coordinate strings for smaller sub-divided polygons
    """
    geometry_proj, crs_proj = projection.project_geometry(polygon)
    gpcs = utils_geo._consolidate_subdivide_geometry(geometry_proj)
    geometry, _ = projection.project_geometry(gpcs, crs=crs_proj, to_latlong=True)
    polygon_coord_strs = utils_geo._get_polygons_coordinates(geometry)
    utils.log(f"Requesting data within polygon from API in {len(polygon_coord_strs)} request(s)")
    return polygon_coord_strs


def _create_overpass_query(polygon_coord_str, tags):
    """
    Create an overpass query string based on passed tags.

    Parameters
    ----------
    polygon_coord_str : list
        list of lat lng coordinates
    tags : dict
        dict of tags used for finding geometry in the selected area

    Returns
    -------
    query : string
    """
    # create overpass settings string
    overpass_settings = _make_overpass_settings()

    # make sure every value in dict is bool, str, or list of str
    error_msg = "tags must be a dict with values of bool, str, or list of str"
    if not isinstance(tags, dict):
        raise TypeError(error_msg)

    tags_dict = dict()
    for key, value in tags.items():

        if isinstance(value, bool):
            tags_dict[key] = value

        elif isinstance(value, str):
            tags_dict[key] = [value]

        elif isinstance(value, list):
            if not all(isinstance(s, str) for s in value):
                raise TypeError(error_msg)
            tags_dict[key] = value

        else:
            raise TypeError(error_msg)

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
                tag_str = f"['{key}'](poly:'{polygon_coord_str}');(._;>;);"
            else:
                # otherwise, pass "key"="value"
                tag_str = f"['{key}'='{value}'](poly:'{polygon_coord_str}');(._;>;);"

            for kind in ("node", "way", "relation"):
                components.append(f"({kind}{tag_str});")

    # finalize query and return
    components = "".join(components)
    query = f"{overpass_settings};({components});out;"

    return query


def _osm_net_download(polygon, network_type, custom_filter):
    """
    Download OSM ways and nodes within some polygon from the Overpass API.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        geographic boundaries to fetch the street network within
    network_type : string
        what type of street network to get if custom_filter is not None
    custom_filter : string
        a custom network filter to be used instead of the network_type presets

    Returns
    -------
    response_jsons : list
    """
    # create a filter to exclude certain kinds of ways based on the requested
    # network_type, if provided, otherwise use custom_filter
    if custom_filter is not None:
        osm_filter = custom_filter
    else:
        osm_filter = _get_osm_filter(network_type)

    response_jsons = []

    # create overpass settings string
    overpass_settings = _make_overpass_settings()

    # subdivide query polygon to get list of sub-divided polygon coord strings
    polygon_coord_strs = _make_overpass_polygon_coord_strs(polygon)

    # pass each polygon exterior coordinates in the list to the API, one at a
    # time. The '>' makes it recurse so we get ways and the ways' nodes.
    for polygon_coord_str in polygon_coord_strs:
        query_str = f"{overpass_settings};(way{osm_filter}(poly:'{polygon_coord_str}');>;);out;"
        response_json = overpass_request(data={"data": query_str})
        response_jsons.append(response_json)
    utils.log(
        f"Got all network data within polygon from API in {len(polygon_coord_strs)} request(s)"
    )

    return response_jsons


def _osm_geometry_download(polygon, tags):
    """
    Download all OSM geometry within some polygon from the Overpass API.

    Note that if a polygon is passed-in, the query will be limited to the
    exterior ring only.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon
        geographic boundaries to fetch geometry within
    tags : dict
        dict of tags used for finding geometry in the selected area

    Returns
    -------
    response_jsons : list
        list of JSON responses from the Overpass server
    """
    response_jsons = []

    # subdivide query polygon to get list of sub-divided polygon coord strings
    polygon_coord_strs = _make_overpass_polygon_coord_strs(polygon)

    # pass exterior coordinates of each polygon in list to API, one at a time
    for polygon_coord_str in polygon_coord_strs:
        query_str = _create_overpass_query(polygon_coord_str, tags)
        response_json = overpass_request(data={"data": query_str})
        response_jsons.append(response_json)

    utils.log(
        f"Got all geometry data within polygon from API in {len(polygon_coord_strs)} request(s)"
    )

    return response_jsons


def _osm_polygon_download(query, limit=1, polygon_geojson=1):
    """
    Geocode a place and download its boundary geometry from Nominatim API.

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
    # prevent OSM from deduping results so we get precisely 'limit' # of results
    params["dedupe"] = 0
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
    response_json = nominatim_request(params=params)
    return response_json


def nominatim_request(params, request_type="search", pause=1, error_pause=60):
    """
    Send a HTTP GET request to the Nominatim API and return JSON response.

    Parameters
    ----------
    params : OrderedDict
        key-value pairs of parameters
    request_type : string
        Type of Nominatim query. One of: search, reverse, or lookup
    pause : int
        how long to pause before request, in seconds. per the nominatim usage
        policy: "an absolute maximum of 1 request per second" is allowed
    error_pause : int
        how long to pause in seconds before re-trying request if error

    Returns
    -------
    response_json : dict
    """
    if request_type not in {"search", "reverse", "lookup"}:
        raise ValueError('Nominatim request_type must be "search", "reverse", or "lookup"')

    # prepare Nominatim API URL and see if request already exists in cache
    url = settings.nominatim_endpoint.rstrip("/") + "/" + request_type
    prepared_url = requests.Request("GET", url, params=params).prepare().url
    cached_response_json = _retrieve_from_cache(prepared_url)

    if settings.nominatim_key:
        params["key"] = settings.nominatim_key

    if cached_response_json is not None:
        # found response in the cache, return it instead of calling server
        return cached_response_json

    else:
        # if this URL is not already in the cache, pause, then request it
        utils.log(f"Pausing {pause} seconds before making HTTP GET request")
        time.sleep(pause)

        # transmit the HTTP GET request
        utils.log(f"Get {prepared_url} with timeout={settings.timeout}")
        headers = _get_http_headers()
        response = requests.get(url, params=params, timeout=settings.timeout, headers=headers)
        sc = response.status_code

        # log the response size and domain
        size_kb = len(response.content) / 1000.0
        domain = re.findall(r"(?s)//(.*?)/", url)[0]
        utils.log(f"Downloaded {size_kb:,.1f}KB from {domain}")

        try:
            response_json = response.json()

        except Exception:  # pragma: no cover
            if sc in {429, 504}:
                # 429 is 'too many requests' and 504 is 'gateway timeout' from
                # server overload: handle these by pausing then recursively
                # re-trying until we get a valid response from the server
                utils.log(f"{domain} returned {sc}: retry in {error_pause} secs", level=lg.WARNING)
                time.sleep(error_pause)
                response_json = nominatim_request(params, request_type, pause, error_pause)

            else:
                # else, this was an unhandled status code, throw an exception
                utils.log(f"{domain} returned {sc}", level=lg.ERROR)
                raise Exception(f"Server returned:\n{response} {response.reason}\n{response.text}")

        _save_to_cache(prepared_url, response_json, sc)
        return response_json


def overpass_request(data, pause=None, error_pause=60):
    """
    Send a HTTP POST request to the Overpass API and return JSON response.

    Parameters
    ----------
    data : OrderedDict
        key-value pairs of parameters
    pause : int
        how long to pause in seconds before request, if None, will query API
        status endpoint to find when next slot is available
    error_pause : int
        how long to pause in seconds (in addition to `pause`) before re-trying
        request if error

    Returns
    -------
    response_json : dict
    """
    # define the Overpass API URL, then construct a GET-style URL as a string to
    # hash to look up/save to cache
    url = settings.overpass_endpoint.rstrip("/") + "/interpreter"
    prepared_url = requests.Request("GET", url, params=data).prepare().url
    cached_response_json = _retrieve_from_cache(prepared_url, check_remark=True)

    if cached_response_json is not None:
        # found response in the cache, return it instead of calling server
        return cached_response_json

    else:
        # if this URL is not already in the cache, pause, then request it
        if pause is None:
            this_pause = _get_pause()
        utils.log(f"Pausing {this_pause} seconds before making HTTP POST request")
        time.sleep(this_pause)

        # transmit the HTTP POST request
        utils.log(f"Post {prepared_url} with timeout={settings.timeout}")
        headers = _get_http_headers()
        response = requests.post(url, data=data, timeout=settings.timeout, headers=headers)
        sc = response.status_code

        # log the response size and domain
        size_kb = len(response.content) / 1000.0
        domain = re.findall(r"(?s)//(.*?)/", url)[0]
        utils.log(f"Downloaded {size_kb:,.1f}KB from {domain}")

        try:
            response_json = response.json()
            if "remark" in response_json:
                utils.log(f'Server remark: "{response_json["remark"]}"', level=lg.WARNING)

        except Exception:  # pragma: no cover
            if sc in {429, 504}:
                # 429 is 'too many requests' and 504 is 'gateway timeout' from
                # server overload: handle these by pausing then recursively
                # re-trying until we get a valid response from the server
                this_pause = error_pause + _get_pause()
                utils.log(f"{domain} returned {sc}: retry in {this_pause} secs", level=lg.WARNING)
                time.sleep(this_pause)
                response_json = overpass_request(data, pause, error_pause)

            else:
                # else, this was an unhandled status code, throw an exception
                utils.log(f"{domain} returned {sc}", level=lg.ERROR)
                raise Exception(f"Server returned\n{response} {response.reason}\n{response.text}")

        _save_to_cache(prepared_url, response_json, sc)
        return response_json


def _overpass_json_from_file(filepath):
    """
    Read OSM XML from file and return Overpass-like JSON.

    Parameters
    ----------
    filepath : string
        path to file containing OSM XML data

    Returns
    -------
    OSMContentHandler object
    """

    def _opener(filepath):
        _, ext = os.path.splitext(filepath)
        if ext == ".bz2":
            return bz2.BZ2File(filepath)
        else:
            # assume an unrecognized file extension is just XML
            return open(filepath, mode="rb")

    with _opener(filepath) as file:
        handler = _OSMContentHandler()
        xml.sax.parse(file, handler)
        return handler.object


class _OSMContentHandler(xml.sax.handler.ContentHandler):
    """
    SAX content handler for OSM XML.

    Used to build an Overpass-like response JSON object in self.object. For format
    notes, see http://wiki.openstreetmap.org/wiki/OSM_XML#OSM_XML_file_format_notes
    and http://overpass-api.de/output_formats.html#json
    """

    def __init__(self):
        self._element = None
        self.object = {"elements": []}

    def startElement(self, name, attrs):
        if name == "osm":
            self.object.update({k: attrs[k] for k in attrs.keys() if k in {"version", "generator"}})

        elif name in {"node", "way"}:
            self._element = dict(type=name, tags={}, nodes=[], **attrs)
            self._element.update({k: float(attrs[k]) for k in attrs.keys() if k in {"lat", "lon"}})
            self._element.update(
                {
                    k: int(attrs[k])
                    for k in attrs.keys()
                    if k in {"id", "uid", "version", "changeset"}
                }
            )

        elif name == "relation":
            self._element = dict(type=name, tags={}, members=[], **attrs)
            self._element.update(
                {
                    k: int(attrs[k])
                    for k in attrs.keys()
                    if k in {"id", "uid", "version", "changeset"}
                }
            )

        elif name == "tag":
            self._element["tags"].update({attrs["k"]: attrs["v"]})

        elif name == "nd":
            self._element["nodes"].append(int(attrs["ref"]))

        elif name == "member":
            self._element["members"].append(
                {k: (int(attrs[k]) if k == "ref" else attrs[k]) for k in attrs.keys()}
            )

    def endElement(self, name):
        if name in {"node", "way", "relation"}:
            self.object["elements"].append(self._element)
