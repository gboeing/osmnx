################################################################################
# Module: core.py
# Description: Retrieve and construct spatial geometries and street networks
#              from OpenStreetMap
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

import datetime as dt
import geopandas as gpd
import hashlib
import io
import json
import logging as lg
import math
import networkx as nx
import numpy as np
import os
import pandas as pd
import re
import requests
import time

from collections import OrderedDict
from dateutil import parser as date_parser
from itertools import groupby
from shapely.geometry import LineString
from shapely.geometry import MultiPolygon
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.ops import unary_union

from . import settings
from .projection import project_geometry
from .projection import project_gdf
from .simplify import simplify_graph
from .utils import log
from .utils import make_str
from .utils import get_largest_component
from .utils import great_circle_vec
from .utils import get_nearest_node
from .utils import geocode
from .utils import count_streets_per_node
from .utils import overpass_json_from_file


class EmptyOverpassResponse(ValueError):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


class InvalidDistanceType(ValueError):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


class UnknownNetworkType(ValueError):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


class InsufficientNetworkQueryArguments(ValueError):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


def save_to_cache(url, response_json):
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
            log('Saved nothing to cache because response_json is None')
        else:
            # create the folder on the disk if it doesn't already exist
            if not os.path.exists(settings.cache_folder):
                os.makedirs(settings.cache_folder)

            # hash the url (to make filename shorter than the often extremely
            # long url)
            filename = hashlib.md5(url.encode('utf-8')).hexdigest()
            cache_path_filename = os.path.join(settings.cache_folder, os.extsep.join([filename, 'json']))

            # dump to json, and save to file
            json_str = make_str(json.dumps(response_json))
            with io.open(cache_path_filename, 'w', encoding='utf-8') as cache_file:
                cache_file.write(json_str)

            log('Saved response to cache file "{}"'.format(cache_path_filename))


def get_from_cache(url):
    """
    Retrieve a HTTP response json object from the cache.

    Parameters
    ----------
    url : string
        the url of the request

    Returns
    -------
    response_json : dict
    """
    # if the tool is configured to use the cache
    if settings.use_cache:
        # determine the filename by hashing the url
        filename = hashlib.md5(url.encode('utf-8')).hexdigest()

        cache_path_filename = os.path.join(settings.cache_folder, os.extsep.join([filename, 'json']))
        # open the cache file for this url hash if it already exists, otherwise
        # return None
        if os.path.isfile(cache_path_filename):
            with io.open(cache_path_filename, encoding='utf-8') as cache_file:
                response_json = json.load(cache_file)
            log('Retrieved response from cache file "{}" for URL "{}"'.format(cache_path_filename, url))
            return response_json


def get_http_headers(user_agent=None, referer=None, accept_language=None):
    """
    Update the default requests HTTP headers with OSMnx info.

    Parameters
    ----------
    user_agent : str
        the user agent string, if None will set with OSMnx default
    referer : str
        the referer string, if None will set with OSMnx default
    accept_language : str
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
    headers.update({'User-Agent': user_agent, 'referer': referer, 'Accept-Language': accept_language})
    return headers


def get_pause_duration(recursive_delay=5, default_duration=10):
    """
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
        response = requests.get('http://overpass-api.de/api/status', headers=get_http_headers())
        status = response.text.split('\n')[3]
        status_first_token = status.split(' ')[0]
    except Exception:
        # if we cannot reach the status endpoint or parse its output, log an
        # error and return default duration
        log('Unable to query http://overpass-api.de/api/status', level=lg.ERROR)
        return default_duration

    try:
        # if first token is numeric, it's how many slots you have available - no
        # wait required
        available_slots = int(status_first_token)
        pause_duration = 0
    except Exception:
        # if first token is 'Slot', it tells you when your slot will be free
        if status_first_token == 'Slot':
            utc_time_str = status.split(' ')[3]
            utc_time = date_parser.parse(utc_time_str).replace(tzinfo=None)
            pause_duration = math.ceil((utc_time - dt.datetime.utcnow()).total_seconds())
            pause_duration = max(pause_duration, 1)

        # if first token is 'Currently', it is currently running a query so
        # check back in recursive_delay seconds
        elif status_first_token == 'Currently':
            time.sleep(recursive_delay)
            pause_duration = get_pause_duration()

        else:
            # any other status is unrecognized - log an error and return default
            # duration
            log('Unrecognized server status: "{}"'.format(status), level=lg.ERROR)
            return default_duration

    return pause_duration


def nominatim_request(params, type = "search", pause_duration=1, timeout=30, error_pause_duration=180):
    """
    Send a request to the Nominatim API via HTTP GET and return the JSON
    response.

    Parameters
    ----------
    params : dict or OrderedDict
        key-value pairs of parameters
    type : string
        Type of Nominatim query. One of the following: search, reverse or lookup
    pause_duration : int
        how long to pause before requests, in seconds
    timeout : int
        the timeout interval for the requests library
    error_pause_duration : int
        how long to pause in seconds before re-trying requests if error

    Returns
    -------
    response_json : dict
    """

    known_requests = {"search", "reverse", "lookup"}
    if type not in known_requests:
        raise ValueError("The type of Nominatim request is invalid. Please choose one of {{'search', 'reverse', 'lookup'}}")

    # prepare the Nominatim API URL and see if request already exists in the
    # cache
    url = 'https://nominatim.openstreetmap.org/{}'.format(type)
    prepared_url = requests.Request('GET', url, params=params).prepare().url
    cached_response_json = get_from_cache(prepared_url)

    if cached_response_json is not None:
        # found this request in the cache, just return it instead of making a
        # new HTTP call
        return cached_response_json

    else:
        # if this URL is not already in the cache, pause, then request it
        log('Pausing {:,.2f} seconds before making API GET request'.format(pause_duration))
        time.sleep(pause_duration)
        start_time = time.time()
        log('Requesting {} with timeout={}'.format(prepared_url, timeout))
        response = requests.get(url, params=params, timeout=timeout, headers=get_http_headers())

        # get the response size and the domain, log result
        size_kb = len(response.content) / 1000.
        domain = re.findall(r'(?s)//(.*?)/', url)[0]
        log('Downloaded {:,.1f}KB from {} in {:,.2f} seconds'.format(size_kb, domain, time.time()-start_time))

        try:
            response_json = response.json()
            save_to_cache(prepared_url, response_json)
        except Exception:
            #429 is 'too many requests' and 504 is 'gateway timeout' from server
            # overload - handle these errors by recursively calling
            # nominatim_request until we get a valid response
            if response.status_code in [429, 504]:
                # pause for error_pause_duration seconds before re-trying request
                log('Server at {} returned status code {} and no JSON data. Re-trying request in {:.2f} seconds.'.format(domain,
                                                                                                                         response.status_code,
                                                                                                                         error_pause_duration),
                                                                                                                         level=lg.WARNING)
                time.sleep(error_pause_duration)
                response_json = nominatim_request(params=params, pause_duration=pause_duration, timeout=timeout)

            # else, this was an unhandled status_code, throw an exception
            else:
                log('Server at {} returned status code {} and no JSON data'.format(domain, response.status_code), level=lg.ERROR)
                raise Exception('Server returned no JSON data.\n{} {}\n{}'.format(response, response.reason, response.text))

        return response_json


def overpass_request(data, pause_duration=None, timeout=180, error_pause_duration=None):
    """
    Send a request to the Overpass API via HTTP POST and return the JSON
    response.

    Parameters
    ----------
    data : dict or OrderedDict
        key-value pairs of parameters to post to the API
    pause_duration : int
        how long to pause in seconds before requests, if None, will query API
        status endpoint to find when next slot is available
    timeout : int
        the timeout interval for the requests library
    error_pause_duration : int
        how long to pause in seconds before re-trying requests if error

    Returns
    -------
    dict
    """

    # define the Overpass API URL, then construct a GET-style URL as a string to
    # hash to look up/save to cache
    url = 'http://overpass-api.de/api/interpreter'
    prepared_url = requests.Request('GET', url, params=data).prepare().url
    cached_response_json = get_from_cache(prepared_url)

    if cached_response_json is not None:
        # found this request in the cache, just return it instead of making a
        # new HTTP call
        return cached_response_json

    else:
        # if this URL is not already in the cache, pause, then request it
        if pause_duration is None:
            this_pause_duration = get_pause_duration()
        log('Pausing {:,.2f} seconds before making API POST request'.format(this_pause_duration))
        time.sleep(this_pause_duration)
        start_time = time.time()
        log('Posting to {} with timeout={}, "{}"'.format(url, timeout, data))
        response = requests.post(url, data=data, timeout=timeout, headers=get_http_headers())

        # get the response size and the domain, log result
        size_kb = len(response.content) / 1000.
        domain = re.findall(r'(?s)//(.*?)/', url)[0]
        log('Downloaded {:,.1f}KB from {} in {:,.2f} seconds'.format(size_kb, domain, time.time()-start_time))

        try:
            response_json = response.json()
            if 'remark' in response_json:
                log('Server remark: "{}"'.format(response_json['remark'], level=lg.WARNING))
            save_to_cache(prepared_url, response_json)
        except Exception:
            #429 is 'too many requests' and 504 is 'gateway timeout' from server
            # overload - handle these errors by recursively calling
            # overpass_request until we get a valid response
            if response.status_code in [429, 504]:
                # pause for error_pause_duration seconds before re-trying request
                if error_pause_duration is None:
                    error_pause_duration = get_pause_duration()
                log('Server at {} returned status code {} and no JSON data. Re-trying request in {:.2f} seconds.'.format(domain,
                                                                                                                         response.status_code,
                                                                                                                         error_pause_duration),
                                                                                                                         level=lg.WARNING)
                time.sleep(error_pause_duration)
                response_json = overpass_request(data=data, pause_duration=pause_duration, timeout=timeout)

            # else, this was an unhandled status_code, throw an exception
            else:
                log('Server at {} returned status code {} and no JSON data'.format(domain, response.status_code), level=lg.ERROR)
                raise Exception('Server returned no JSON data.\n{} {}\n{}'.format(response, response.reason, response.text))

        return response_json


def osm_polygon_download(query, limit=1, polygon_geojson=1):
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
    dict
    """
    # define the parameters
    params = OrderedDict()
    params['format'] = 'json'
    params['limit'] = limit
    params['dedupe'] = 0 #this prevents OSM from de-duping results so we're guaranteed to get precisely 'limit' number of results
    params['polygon_geojson'] = polygon_geojson

    # add the structured query dict (if provided) to params, otherwise query
    # with place name string
    if isinstance(query, str):
        params['q'] = query
    elif isinstance(query, dict):
        # add the query keys in alphabetical order so the URL is the same string
        # each time, for caching purposes
        for key in sorted(list(query.keys())):
            params[key] = query[key]
    else:
        raise TypeError('query must be a dict or a string')

    # request the URL, return the JSON
    response_json = nominatim_request(params=params, timeout=30)
    return response_json


def gdf_from_place(query, gdf_name=None, which_result=1, buffer_dist=None):
    """
    Create a GeoDataFrame from a single place name query.

    Parameters
    ----------
    query : string or dict
        query string or structured query dict to geocode/download
    gdf_name : string
        name attribute metadata for GeoDataFrame (this is used to save shapefile
        later)
    which_result : int
        max number of results to return and which to process upon receipt
    buffer_dist : float
        distance to buffer around the place geometry, in meters

    Returns
    -------
    GeoDataFrame
    """
    # if no gdf_name is passed, just use the query
    assert (isinstance(query, dict) or isinstance(query, str)), 'query must be a dict or a string'
    if (gdf_name is None) and isinstance(query, dict):
        gdf_name = ', '.join(list(query.values()))
    elif (gdf_name is None) and isinstance(query, str):
        gdf_name = query

    # get the data from OSM
    data = osm_polygon_download(query, limit=which_result)
    if len(data) >= which_result:

        # extract data elements from the JSON response
        result = data[which_result - 1]
        bbox_south, bbox_north, bbox_west, bbox_east = [float(x) for x in result['boundingbox']]
        geometry = result['geojson']
        place = result['display_name']
        features = [{'type': 'Feature',
                     'geometry': geometry,
                     'properties': {'place_name': place,
                                    'bbox_north': bbox_north,
                                    'bbox_south': bbox_south,
                                    'bbox_east': bbox_east,
                                    'bbox_west': bbox_west}}]

        # if we got an unexpected geometry type (like a point), log a warning
        if geometry['type'] not in ['Polygon', 'MultiPolygon']:
            log('OSM returned a {} as the geometry.'.format(geometry['type']), level=lg.WARNING)

        # create the GeoDataFrame, name it, and set its original CRS to default_crs
        gdf = gpd.GeoDataFrame.from_features(features)
        gdf.gdf_name = gdf_name
        gdf.crs = settings.default_crs

        # if buffer_dist was passed in, project the geometry to UTM, buffer it
        # in meters, then project it back to lat-long
        if buffer_dist is not None:
            gdf_utm = project_gdf(gdf)
            gdf_utm['geometry'] = gdf_utm['geometry'].buffer(buffer_dist)
            gdf = project_gdf(gdf_utm, to_latlong=True)
            log('Buffered the GeoDataFrame "{}" to {} meters'.format(gdf.gdf_name, buffer_dist))

        # return the gdf
        log('Created GeoDataFrame with {} row for query "{}"'.format(len(gdf), query))
        return gdf
    else:
        # if there was no data returned (or fewer results than which_result
        # specified)
        log('OSM returned no results (or fewer than which_result) for query "{}"'.format(query), level=lg.WARNING)
        gdf = gpd.GeoDataFrame()
        gdf.gdf_name = gdf_name
        return gdf


def gdf_from_places(queries, gdf_name='unnamed', buffer_dist=None):
    """
    Create a GeoDataFrame from a list of place names to query.

    Parameters
    ----------
    queries : list
        list of query strings or structured query dicts to geocode/download, one
        at a time
    gdf_name : string
        name attribute metadata for GeoDataFrame (this is used to save shapefile
        later)
    buffer_dist : float
        distance to buffer around the place geometry, in meters

    Returns
    -------
    GeoDataFrame
    """
    # create an empty GeoDataFrame then append each result as a new row
    gdf = gpd.GeoDataFrame()
    for query in queries:
        gdf = gdf.append(gdf_from_place(query, buffer_dist=buffer_dist))

    # reset the index, name the GeoDataFrame
    gdf = gdf.reset_index().drop(labels='index', axis=1)
    gdf.gdf_name = gdf_name

    # set the original CRS of the GeoDataFrame to default_crs, and return it
    gdf.crs = settings.default_crs
    log('Finished creating GeoDataFrame with {} rows from {} queries'.format(len(gdf), len(queries)))
    return gdf


def get_osm_filter(network_type):
    """
    Create a filter to query OSM for the specified network type.

    Parameters
    ----------
    network_type : string
        {'walk', 'bike', 'drive', 'drive_service', 'all', 'all_private', 'none'}
        what type of street or other network to get

    Returns
    -------
    string
    """
    filters = {}

    # driving: filter out un-drivable roads, service roads, private ways, and
    # anything specifying motor=no. also filter out any non-service roads that
    # are tagged as providing parking, driveway, private, or emergency-access
    # services
    filters['drive'] = ('["area"!~"yes"]["highway"!~"cycleway|footway|path|pedestrian|steps|track|corridor|'
                        'proposed|construction|bridleway|abandoned|platform|raceway|service"]'
                        '["motor_vehicle"!~"no"]["motorcar"!~"no"]{}'
                        '["service"!~"parking|parking_aisle|driveway|private|emergency_access"]').format(settings.default_access)

    # drive+service: allow ways tagged 'service' but filter out certain types of
    # service ways
    filters['drive_service'] = ('["area"!~"yes"]["highway"!~"cycleway|footway|path|pedestrian|steps|track|corridor|'
                                'proposed|construction|bridleway|abandoned|platform|raceway"]'
                                '["motor_vehicle"!~"no"]["motorcar"!~"no"]{}'
                                '["service"!~"parking|parking_aisle|private|emergency_access"]').format(settings.default_access)

    # walking: filter out cycle ways, motor ways, private ways, and anything
    # specifying foot=no. allow service roads, permitting things like parking
    # lot lanes, alleys, etc that you *can* walk on even if they're not exactly
    # pleasant walks. some cycleways may allow pedestrians, but this filter ignores
    # such cycleways.
    filters['walk'] = ('["area"!~"yes"]["highway"!~"cycleway|motor|proposed|construction|abandoned|platform|raceway"]'
                       '["foot"!~"no"]["service"!~"private"]{}').format(settings.default_access)

    # biking: filter out foot ways, motor ways, private ways, and anything
    # specifying biking=no
    filters['bike'] = ('["area"!~"yes"]["highway"!~"footway|steps|corridor|motor|proposed|construction|abandoned|platform|raceway"]'
                       '["bicycle"!~"no"]["service"!~"private"]{}').format(settings.default_access)

    # to download all ways, just filter out everything not currently in use or
    # that is private-access only
    filters['all'] = ('["area"!~"yes"]["highway"!~"proposed|construction|abandoned|platform|raceway"]'
                      '["service"!~"private"]{}').format(settings.default_access)

    # to download all ways, including private-access ones, just filter out
    # everything not currently in use
    filters['all_private'] = '["area"!~"yes"]["highway"!~"proposed|construction|abandoned|platform|raceway"]'

    # no filter, needed for infrastructures other than "highway"
    filters['none'] = ''

    if network_type in filters:
        osm_filter = filters[network_type]
    else:
        raise UnknownNetworkType('unknown network_type "{}"'.format(network_type))

    return osm_filter


def osm_net_download(polygon=None, north=None, south=None, east=None, west=None,
                     network_type='all_private', timeout=180, memory=None,
                     max_query_area_size=50*1000*50*1000, infrastructure='way["highway"]',
                     custom_filter=None):
    """
    Download OSM ways and nodes within some bounding box from the Overpass API.

    Parameters
    ----------
    polygon : shapely Polygon or MultiPolygon
        geographic shape to fetch the street network within
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
        max area for any part of the geometry, in the units the geometry is in:
        any polygon bigger will get divided up for multiple queries to API
        (default is 50,000 * 50,000 units [ie, 50km x 50km in area, if units are
        meters])
    infrastructure : string
        download infrastructure of given type. default is streets, ie,
        'way["highway"]') but other infrastructures may be selected like power
        grids, ie, 'way["power"~"line"]'
    custom_filter : string
        a custom network filter to be used instead of the network_type presets

    Returns
    -------
    response_jsons : list
    """

    # check if we're querying by polygon or by bounding box based on which
    # argument(s) where passed into this function
    by_poly = polygon is not None
    by_bbox = not (north is None or south is None or east is None or west is None)
    if not (by_poly or by_bbox):
        raise InsufficientNetworkQueryArguments(
            'You must pass a polygon or north, south, east, and west')

    # create a filter to exclude certain kinds of ways based on the requested
    # network_type
    if custom_filter:
        osm_filter = custom_filter
    else:
        osm_filter = get_osm_filter(network_type)
    response_jsons = []

    # pass server memory allocation in bytes for the query to the API
    # if None, pass nothing so the server will use its default allocation size
    # otherwise, define the query's maxsize parameter value as whatever the
    # caller passed in
    if memory is None:
        maxsize = ''
    else:
        maxsize = '[maxsize:{}]'.format(memory)

    # define the query to send the API
    # specifying way["highway"] means that all ways returned must have a highway
    # key. the {filters} then remove ways by key/value. the '>' makes it recurse
    # so we get ways and way nodes. maxsize is in bytes.
    if by_bbox:
        # turn bbox into a polygon and project to local UTM
        polygon = Polygon([(west, south), (east, south), (east, north), (west, north)])
        geometry_proj, crs_proj = project_geometry(polygon)

        # subdivide it if it exceeds the max area size (in meters), then project
        # back to lat-long
        geometry_proj_consolidated_subdivided = consolidate_subdivide_geometry(geometry_proj, max_query_area_size=max_query_area_size)
        geometry, _ = project_geometry(geometry_proj_consolidated_subdivided, crs=crs_proj, to_latlong=True)
        log('Requesting network data within bounding box from API in {:,} request(s)'.format(len(geometry)))
        start_time = time.time()

        # loop through each polygon rectangle in the geometry (there will only
        # be one if original bbox didn't exceed max area size)
        for poly in geometry:
            # represent bbox as south,west,north,east and round lat-longs to 6
            # decimal places (ie, ~100 mm) so URL strings aren't different
            # due to float rounding issues (for consistent caching)
            west, south, east, north = poly.bounds
            query_template = '[out:json][timeout:{timeout}]{maxsize};({infrastructure}{filters}({south:.6f},{west:.6f},{north:.6f},{east:.6f});>;);out;'
            query_str = query_template.format(north=north, south=south,
                                              east=east, west=west,
                                              infrastructure=infrastructure,
                                              filters=osm_filter,
                                              timeout=timeout, maxsize=maxsize)
            response_json = overpass_request(data={'data':query_str}, timeout=timeout)
            response_jsons.append(response_json)
        log('Got all network data within bounding box from API in {:,} request(s) and {:,.2f} seconds'.format(len(geometry), time.time()-start_time))

    elif by_poly:
        # project to utm, divide polygon up into sub-polygons if area exceeds a
        # max size (in meters), project back to lat-long, then get a list of
        # polygon(s) exterior coordinates
        geometry_proj, crs_proj = project_geometry(polygon)
        geometry_proj_consolidated_subdivided = consolidate_subdivide_geometry(geometry_proj, max_query_area_size=max_query_area_size)
        geometry, _ = project_geometry(geometry_proj_consolidated_subdivided, crs=crs_proj, to_latlong=True)
        polygon_coord_strs = get_polygons_coordinates(geometry)
        log('Requesting network data within polygon from API in {:,} request(s)'.format(len(polygon_coord_strs)))
        start_time = time.time()

        # pass each polygon exterior coordinates in the list to the API, one at
        # a time
        for polygon_coord_str in polygon_coord_strs:
            query_template = '[out:json][timeout:{timeout}]{maxsize};({infrastructure}{filters}(poly:"{polygon}");>;);out;'
            query_str = query_template.format(polygon=polygon_coord_str, infrastructure=infrastructure, filters=osm_filter, timeout=timeout, maxsize=maxsize)
            response_json = overpass_request(data={'data':query_str}, timeout=timeout)
            response_jsons.append(response_json)
        log('Got all network data within polygon from API in {:,} request(s) and {:,.2f} seconds'.format(len(polygon_coord_strs), time.time()-start_time))

    return response_jsons


def consolidate_subdivide_geometry(geometry, max_query_area_size):
    """
    Consolidate a geometry into a convex hull, then subdivide it into smaller sub-polygons if its area exceeds max size (in geometry's units).

    Parameters
    ----------
    geometry : shapely Polygon or MultiPolygon
        the geometry to consolidate and subdivide
    max_query_area_size : float
        max area for any part of the geometry, in the units the geometry is in.
        any polygon bigger will get divided up for multiple queries to API (
        default is 50,000 * 50,000 units (ie, 50km x 50km in area, if units are meters))

    Returns
    -------
    geometry : Polygon or MultiPolygon
    """

    # let the linear length of the quadrats (with which to subdivide the
    # geometry) be the square root of max area size
    quadrat_width = math.sqrt(max_query_area_size)

    if not isinstance(geometry, (Polygon, MultiPolygon)):
        raise TypeError('Geometry must be a shapely Polygon or MultiPolygon')

    # if geometry is a MultiPolygon OR a single Polygon whose area exceeds the
    # max size, get the convex hull around the geometry
    if isinstance(geometry, MultiPolygon) or (isinstance(geometry, Polygon) and geometry.area > max_query_area_size):
        geometry = geometry.convex_hull

    # if geometry area exceeds max size, subdivide it into smaller sub-polygons
    if geometry.area > max_query_area_size:
        geometry = quadrat_cut_geometry(geometry, quadrat_width=quadrat_width)

    if isinstance(geometry, Polygon):
        geometry = MultiPolygon([geometry])

    return geometry


def get_polygons_coordinates(geometry):
    """
    Extract exterior coordinates from polygon(s) to pass to OSM in a query by
    polygon. Ignore the interior ("holes") coordinates.

    Parameters
    ----------
    geometry : shapely Polygon or MultiPolygon
        the geometry to extract exterior coordinates from

    Returns
    -------
    polygon_coord_strs : list
    """

    # extract the exterior coordinates of the geometry to pass to the API later
    polygons_coords = []
    if isinstance(geometry, Polygon):
        x, y = geometry.exterior.xy
        polygons_coords.append(list(zip(x, y)))
    elif isinstance(geometry, MultiPolygon):
        for polygon in geometry:
            x, y = polygon.exterior.xy
            polygons_coords.append(list(zip(x, y)))
    else:
        raise TypeError('Geometry must be a shapely Polygon or MultiPolygon')

    # convert the exterior coordinates of the polygon(s) to the string format
    # the API expects
    polygon_coord_strs = []
    for coords in polygons_coords:
        s = ''
        separator = ' '
        for coord in list(coords):
            # round floating point lats and longs to 6 decimal places (ie, ~100 mm),
            # so we can hash and cache strings consistently
            s = '{}{}{:.6f}{}{:.6f}'.format(s, separator, coord[1], separator, coord[0])
        polygon_coord_strs.append(s.strip(separator))

    return polygon_coord_strs


def get_node(element):
    """
    Convert an OSM node element into the format for a networkx node.

    Parameters
    ----------
    element : dict
        an OSM node element

    Returns
    -------
    dict
    """

    node = {}
    node['y'] = element['lat']
    node['x'] = element['lon']
    node['osmid'] = element['id']
    if 'tags' in element:
        for useful_tag in settings.useful_tags_node:
            if useful_tag in element['tags']:
                node[useful_tag] = element['tags'][useful_tag]
    return node


def get_path(element):
    """
    Convert an OSM way element into the format for a networkx graph path.

    Parameters
    ----------
    element : dict
        an OSM way element

    Returns
    -------
    dict
    """

    path = {}
    path['osmid'] = element['id']

    # remove any consecutive duplicate elements in the list of nodes
    grouped_list = groupby(element['nodes'])
    path['nodes'] = [group[0] for group in grouped_list]

    if 'tags' in element:
        for useful_tag in settings.useful_tags_path:
            if useful_tag in element['tags']:
                path[useful_tag] = element['tags'][useful_tag]
    return path


def parse_osm_nodes_paths(osm_data):
    """
    Construct dicts of nodes and paths with key=osmid and value=dict of
    attributes.

    Parameters
    ----------
    osm_data : dict
        JSON response from from the Overpass API

    Returns
    -------
    nodes, paths : tuple
    """

    nodes = {}
    paths = {}
    for element in osm_data['elements']:
        if element['type'] == 'node':
            key = element['id']
            nodes[key] = get_node(element)
        elif element['type'] == 'way': #osm calls network paths 'ways'
            key = element['id']
            paths[key] = get_path(element)

    return nodes, paths


def remove_isolated_nodes(G):
    """
    Remove from a graph all the nodes that have no incident edges (ie, node
    degree = 0).

    Parameters
    ----------
    G : networkx multidigraph
        the graph from which to remove nodes

    Returns
    -------
    networkx multidigraph
    """

    isolated_nodes = [node for node, degree in dict(G.degree()).items() if degree < 1]
    G.remove_nodes_from(isolated_nodes)
    log('Removed {:,} isolated nodes'.format(len(isolated_nodes)))
    return G


def truncate_graph_dist(G, source_node, max_distance=1000, weight='length', retain_all=False):
    """
    Remove everything further than some network distance from a specified node
    in graph.

    Parameters
    ----------
    G : networkx multidigraph
    source_node : int
        the node in the graph from which to measure network distances to other
        nodes
    max_distance : int
        remove every node in the graph greater than this distance from the
        source_node
    weight : string
        how to weight the graph when measuring distance (default 'length' is
        how many meters long the edge is)
    retain_all : bool
        if True, return the entire graph even if it is not connected

    Returns
    -------
    networkx multidigraph
    """

    # get the shortest distance between the node and every other node, then
    # remove every node further than max_distance away
    start_time = time.time()
    G = G.copy()
    distances = nx.shortest_path_length(G, source=source_node, weight=weight)
    distant_nodes = {key:value for key, value in dict(distances).items() if value > max_distance}
    G.remove_nodes_from(distant_nodes.keys())
    log('Truncated graph by weighted network distance in {:,.2f} seconds'.format(time.time()-start_time))

    # remove any isolated nodes and retain only the largest component (if
    # retain_all is True)
    if not retain_all:
        G = remove_isolated_nodes(G)
        G = get_largest_component(G)

    return G


def truncate_graph_bbox(G, north, south, east, west, truncate_by_edge=False, retain_all=False):
    """
    Remove every node in graph that falls outside a bounding box.

    Needed because overpass returns entire ways that also include nodes outside
    the bbox if the way (that is, a way with a single OSM ID) has a node inside
    the bbox at some point.

    Parameters
    ----------
    G : networkx multidigraph
    north : float
        northern latitude of bounding box
    south : float
        southern latitude of bounding box
    east : float
        eastern longitude of bounding box
    west : float
        western longitude of bounding box
    truncate_by_edge : bool
        if True retain node if it's outside bbox but at least one of node's
        neighbors are within bbox
    retain_all : bool
        if True, return the entire graph even if it is not connected

    Returns
    -------
    networkx multidigraph
    """

    start_time = time.time()
    G = G.copy()
    nodes_outside_bbox = []

    for node, data in G.nodes(data=True):
        if data['y'] > north or data['y'] < south or data['x'] > east or data['x'] < west:
            # this node is outside the bounding box
            if not truncate_by_edge:
                # if we're not truncating by edge, add node to list of nodes
                # outside the bounding box
                nodes_outside_bbox.append(node)
            else:
                # if we're truncating by edge, see if any of node's neighbors
                # are within bounding box
                any_neighbors_in_bbox = False
                neighbors = list(G.successors(node)) + list(G.predecessors(node))
                for neighbor in neighbors:
                    x = G.nodes[neighbor]['x']
                    y = G.nodes[neighbor]['y']
                    if y < north and y > south and x < east and x > west:
                        any_neighbors_in_bbox = True
                        break

                # if none of its neighbors are within the bounding box, add node
                # to list of nodes outside the bounding box
                if not any_neighbors_in_bbox:
                    nodes_outside_bbox.append(node)

    G.remove_nodes_from(nodes_outside_bbox)
    log('Truncated graph by bounding box in {:,.2f} seconds'.format(time.time()-start_time))

    # remove any isolated nodes and retain only the largest component (if
    # retain_all is True)
    if not retain_all:
        G = remove_isolated_nodes(G)
        G = get_largest_component(G)

    return G


def quadrat_cut_geometry(geometry, quadrat_width, min_num=3, buffer_amount=1e-9):
    """
    Split a Polygon or MultiPolygon up into sub-polygons of a specified size,
    using quadrats.

    Parameters
    ----------
    geometry : shapely Polygon or MultiPolygon
        the geometry to split up into smaller sub-polygons
    quadrat_width : numeric
        the linear width of the quadrats with which to cut up the geometry (in
        the units the geometry is in)
    min_num : int
        the minimum number of linear quadrat lines (e.g., min_num=3 would
        produce a quadrat grid of 4 squares)
    buffer_amount : numeric
        buffer the quadrat grid lines by quadrat_width times buffer_amount

    Returns
    -------
    shapely MultiPolygon
    """

    # create n evenly spaced points between the min and max x and y bounds
    west, south, east, north = geometry.bounds
    x_num = math.ceil((east-west) / quadrat_width) + 1
    y_num = math.ceil((north-south) / quadrat_width) + 1
    x_points = np.linspace(west, east, num=max(x_num, min_num))
    y_points = np.linspace(south, north, num=max(y_num, min_num))

    # create a quadrat grid of lines at each of the evenly spaced points
    vertical_lines = [LineString([(x, y_points[0]), (x, y_points[-1])]) for x in x_points]
    horizont_lines = [LineString([(x_points[0], y), (x_points[-1], y)]) for y in y_points]
    lines = vertical_lines + horizont_lines

    # buffer each line to distance of the quadrat width divided by 1 billion,
    # take their union, then cut geometry into pieces by these quadrats
    buffer_size = quadrat_width * buffer_amount
    lines_buffered = [line.buffer(buffer_size) for line in lines]
    quadrats = unary_union(lines_buffered)
    multipoly = geometry.difference(quadrats)

    return multipoly


def intersect_index_quadrats(gdf, geometry, quadrat_width=0.05, min_num=3, buffer_amount=1e-9):
    """
    Intersect points with a polygon, using an r-tree spatial index and cutting
    the polygon up into smaller sub-polygons for r-tree acceleration.

    Parameters
    ----------
    gdf : GeoDataFrame
        the set of points to intersect
    geometry : shapely Polygon or MultiPolygon
        the geometry to intersect with the points
    quadrat_width : numeric
        the linear length (in degrees) of the quadrats with which to cut up the
        geometry (default = 0.05, approx 4km at NYC's latitude)
    min_num : int
        the minimum number of linear quadrat lines (e.g., min_num=3 would
        produce a quadrat grid of 4 squares)
    buffer_amount : numeric
        buffer the quadrat grid lines by quadrat_width times buffer_amount

    Returns
    -------
    GeoDataFrame
    """

    # create an empty dataframe to append matches to
    points_within_geometry = pd.DataFrame()

    # cut the geometry into chunks for r-tree spatial index intersecting
    multipoly = quadrat_cut_geometry(geometry, quadrat_width=quadrat_width, buffer_amount=buffer_amount, min_num=min_num)

    # create an r-tree spatial index for the nodes (ie, points)
    start_time = time.time()
    sindex = gdf['geometry'].sindex
    log('Created r-tree spatial index for {:,} points in {:,.2f} seconds'.format(len(gdf), time.time()-start_time))

    # loop through each chunk of the geometry to find approximate and then
    # precisely intersecting points
    start_time = time.time()
    for poly in multipoly:

        # buffer by the tiny distance to account for any space lost in the
        # quadrat cutting, otherwise may miss point(s) that lay directly on
        # quadrat line
        buffer_size = quadrat_width * buffer_amount
        poly = poly.buffer(buffer_size).buffer(0)

        # find approximate matches with r-tree, then precise matches from those
        # approximate ones
        if poly.is_valid and poly.area > 0:
            possible_matches_index = list(sindex.intersection(poly.bounds))
            possible_matches = gdf.iloc[possible_matches_index]
            precise_matches = possible_matches[possible_matches.intersects(poly)]
            points_within_geometry = points_within_geometry.append(precise_matches)

    if len(points_within_geometry) > 0:
        # drop duplicate points, if buffered poly caused an overlap on point(s)
        # that lay directly on a quadrat line
        points_within_geometry = points_within_geometry.drop_duplicates(subset='node')
    else:
        # after simplifying the graph, and given the requested network type,
        # there are no nodes inside the polygon - can't create graph from that
        # so throw error
        raise Exception('There are no nodes within the requested geometry')

    log('Identified {:,} nodes inside polygon in {:,.2f} seconds'.format(len(points_within_geometry), time.time()-start_time))
    return points_within_geometry


def truncate_graph_polygon(G, polygon, retain_all=False, truncate_by_edge=False, quadrat_width=0.05, min_num=3, buffer_amount=1e-9):
    """
    Remove every node in graph that falls outside some shapely Polygon or
    MultiPolygon.

    Parameters
    ----------
    G : networkx multidigraph
    polygon : Polygon or MultiPolygon
        only retain nodes in graph that lie within this geometry
    retain_all : bool
        if True, return the entire graph even if it is not connected
    truncate_by_edge : bool
        if True retain node if it's outside polygon but at least one of node's
        neighbors are within polygon (NOT CURRENTLY IMPLEMENTED)
    quadrat_width : numeric
        passed on to intersect_index_quadrats: the linear length (in degrees) of
        the quadrats with which to cut up the geometry (default = 0.05, approx
        4km at NYC's latitude)
    min_num : int
        passed on to intersect_index_quadrats: the minimum number of linear
        quadrat lines (e.g., min_num=3 would produce a quadrat grid of 4
        squares)
    buffer_amount : numeric
        passed on to intersect_index_quadrats: buffer the quadrat grid lines by
        quadrat_width times buffer_amount

    Returns
    -------
    networkx multidigraph
    """

    start_time = time.time()
    G = G.copy()
    log('Identifying all nodes that lie outside the polygon...')

    # get a GeoDataFrame of all the nodes, for spatial analysis
    node_geom = [Point(data['x'], data['y']) for _, data in G.nodes(data=True)]
    gdf_nodes = gpd.GeoDataFrame({'node':pd.Series(G.nodes()), 'geometry':node_geom})
    gdf_nodes.crs = G.graph['crs']

    # find all the nodes in the graph that lie outside the polygon
    points_within_geometry = intersect_index_quadrats(gdf_nodes, polygon, quadrat_width=quadrat_width, min_num=min_num, buffer_amount=buffer_amount)
    nodes_outside_polygon = gdf_nodes[~gdf_nodes.index.isin(points_within_geometry.index)]

    # now remove from the graph all those nodes that lie outside the place
    # polygon
    start_time = time.time()
    G.remove_nodes_from(nodes_outside_polygon['node'])
    log('Removed {:,} nodes outside polygon in {:,.2f} seconds'.format(len(nodes_outside_polygon), time.time()-start_time))

    # remove any isolated nodes and retain only the largest component (if retain_all is False)
    if not retain_all:
        G = remove_isolated_nodes(G)
        G = get_largest_component(G)

    return G


def add_edge_lengths(G):
    """
    Add length (meters) attribute to each edge by great circle distance between
    nodes u and v.

    Parameters
    ----------
    G : networkx multidigraph

    Returns
    -------
    G : networkx multidigraph
    """

    start_time = time.time()

    # first load all the edges' origin and destination coordinates as a
    # dataframe indexed by u, v, key
    coords = np.array([[u, v, k, G.nodes[u]['y'], G.nodes[u]['x'], G.nodes[v]['y'], G.nodes[v]['x']] for u, v, k in G.edges(keys=True)])
    df_coords = pd.DataFrame(coords, columns=['u', 'v', 'k', 'u_y', 'u_x', 'v_y', 'v_x'])
    df_coords[['u', 'v', 'k']] = df_coords[['u', 'v', 'k']].astype(np.int64)
    df_coords = df_coords.set_index(['u', 'v', 'k'])

    # then calculate the great circle distance with the vectorized function
    gc_distances = great_circle_vec(lat1=df_coords['u_y'],
                                    lng1=df_coords['u_x'],
                                    lat2=df_coords['v_y'],
                                    lng2=df_coords['v_x'])

    # fill nulls with zeros and round to the millimeter
    gc_distances = gc_distances.fillna(value=0).round(3)
    nx.set_edge_attributes(G, name='length', values=gc_distances.to_dict())

    log('Added edge lengths to graph in {:,.2f} seconds'.format(time.time()-start_time))
    return G


def add_path(G, data, one_way):
    """
    Add a path to the graph.

    Parameters
    ----------
    G : networkx multidigraph
    data : dict
        the attributes of the path
    one_way : bool
        if this path is one-way or if it is bi-directional

    Returns
    -------
    None
    """

    # extract the ordered list of nodes from this path element, then delete it
    # so we don't add it as an attribute to the edge later
    path_nodes = data['nodes']
    del data['nodes']

    # set the oneway attribute to the passed-in value, to make it consistent
    # True/False values
    data['oneway'] = one_way

    # zip together the path nodes so you get tuples like (0,1), (1,2), (2,3)
    # and so on
    path_edges = list(zip(path_nodes[:-1], path_nodes[1:]))
    G.add_edges_from(path_edges, **data)

    # if the path is NOT one-way
    if not one_way:
        # reverse the direction of each edge and add this path going the
        # opposite direction
        path_edges_opposite_direction = [(v, u) for u, v in path_edges]
        G.add_edges_from(path_edges_opposite_direction, **data)


def add_paths(G, paths, bidirectional=False):
    """
    Add a collection of paths to the graph.

    Parameters
    ----------
    G : networkx multidigraph
    paths : dict
        the paths from OSM
    bidirectional : bool
        if True, create bidirectional edges for one-way streets


    Returns
    -------
    None
    """

    # the list of values OSM uses in its 'oneway' tag to denote True
    osm_oneway_values = ['yes', 'true', '1', '-1']

    for data in paths.values():

        # if this path is tagged as one-way and if it is not a walking network,
        # then we'll add the path in one direction only
        if ('oneway' in data and data['oneway'] in osm_oneway_values) and not bidirectional:
            if data['oneway'] == '-1':
                # paths with a one-way value of -1 are one-way, but in the
                # reverse direction of the nodes' order, see osm documentation
                data['nodes'] = list(reversed(data['nodes']))
            # add this path (in only one direction) to the graph
            add_path(G, data, one_way=True)

        elif ('junction' in data and data['junction'] == 'roundabout') and not bidirectional:
            # roundabout are also oneway but not tagged as is
            add_path(G, data, one_way=True)

        # else, this path is not tagged as one-way or it is a walking network
        # (you can walk both directions on a one-way street)
        else:
            # add this path (in both directions) to the graph and set its
            # 'oneway' attribute to False. if this is a walking network, this
            # may very well be a one-way street (as cars/bikes go), but in a
            # walking-only network it is a bi-directional edge
            add_path(G, data, one_way=False)

    return G


def create_graph(response_jsons, name='unnamed', retain_all=False, bidirectional=False):
    """
    Create a networkx graph from Overpass API HTTP response objects.

    Parameters
    ----------
    response_jsons : list
        list of dicts of JSON responses from from the Overpass API
    name : string
        the name of the graph
    retain_all : bool
        if True, return the entire graph even if it is not connected
    bidirectional : bool
        if True, create bidirectional edges for one-way streets

    Returns
    -------
    networkx multidigraph
    """

    log('Creating networkx graph from downloaded OSM data...')
    start_time = time.time()

    # make sure we got data back from the server requests
    elements = []
    for response_json in response_jsons:
        elements.extend(response_json['elements'])
    if len(elements) < 1:
        raise EmptyOverpassResponse('There are no data elements in the response JSON objects')

    # create the graph as a MultiDiGraph and set the original CRS to default_crs
    G = nx.MultiDiGraph(name=name, crs=settings.default_crs)

    # extract nodes and paths from the downloaded osm data
    nodes = {}
    paths = {}
    for osm_data in response_jsons:
        nodes_temp, paths_temp = parse_osm_nodes_paths(osm_data)
        for key, value in nodes_temp.items():
            nodes[key] = value
        for key, value in paths_temp.items():
            paths[key] = value

    # add each osm node to the graph
    for node, data in nodes.items():
        G.add_node(node, **data)

    # add each osm way (aka, path) to the graph
    G = add_paths(G, paths, bidirectional=bidirectional)

    # retain only the largest connected component, if caller did not
    # set retain_all=True
    if not retain_all:
        G = get_largest_component(G)

    log('Created graph with {:,} nodes and {:,} edges in {:,.2f} seconds'.format(len(list(G.nodes())), len(list(G.edges())), time.time()-start_time))

    # add length (great circle distance between nodes) attribute to each edge to
    # use as weight
    if len(G.edges) > 0:
        G = add_edge_lengths(G)

    return G


def bbox_from_point(point, distance=1000, project_utm=False, return_crs=False):
    """
    Create a bounding box some distance in each direction (north, south, east,
    and west) from some (lat, lng) point.

    Parameters
    ----------
    point : tuple
        the (lat, lon) point to create the bounding box around
    distance : int
        how many meters the north, south, east, and west sides of the box should
        each be from the point
    project_utm : bool
        if True return bbox as UTM coordinates
    return_crs : bool
        if True and project_utm=True, return the projected CRS

    Returns
    -------
    north, south, east, west : tuple, if return_crs=False
    north, south, east, west, crs_proj : tuple, if return_crs=True
    """

    # reverse the order of the (lat,lng) point so it is (x,y) for shapely, then
    # project to UTM and buffer in meters
    lat, lng = point
    point_proj, crs_proj = project_geometry(Point((lng, lat)))
    buffer_proj = point_proj.buffer(distance)

    if project_utm:
        west, south, east, north = buffer_proj.bounds
        log('Created bounding box {} meters in each direction from {} and projected it: {},{},{},{}'.format(distance, point, north, south, east, west))
    else:
        # if project_utm is False, project back to lat-long then get the
        # bounding coordinates
        buffer_latlong, _ = project_geometry(buffer_proj, crs=crs_proj, to_latlong=True)
        west, south, east, north = buffer_latlong.bounds
        log('Created bounding box {} meters in each direction from {}: {},{},{},{}'.format(distance, point, north, south, east, west))

    if return_crs:
        return north, south, east, west, crs_proj
    else:
        return north, south, east, west


def graph_from_bbox(north, south, east, west, network_type='all_private',
                    simplify=True, retain_all=False, truncate_by_edge=False,
                    name='unnamed', timeout=180, memory=None,
                    max_query_area_size=50*1000*50*1000, clean_periphery=True,
                    infrastructure='way["highway"]', custom_filter=None):
    """
    Create a networkx graph from OSM data within some bounding box.

    Parameters
    ----------
    north : float
        northern latitude of bounding box
    south : float
        southern latitude of bounding box
    east : float
        eastern longitude of bounding box
    west : float
        western longitude of bounding box
    network_type : string
        what type of street network to get
    simplify : bool
        if true, simplify the graph topology
    retain_all : bool
        if True, return the entire graph even if it is not connected
    truncate_by_edge : bool
        if True retain node if it's outside bbox but at least one of node's
        neighbors are within bbox
    name : string
        the name of the graph
    timeout : int
        the timeout interval for requests and to pass to API
    memory : int
        server memory allocation size for the query, in bytes. If none, server
        will use its default allocation size
    max_query_area_size : float
        max size for any part of the geometry, in square degrees: any polygon
        bigger will get divided up for multiple queries to API
    clean_periphery : bool
        if True (and simplify=True), buffer 0.5km to get a graph larger than
        requested, then simplify, then truncate it to requested spatial extent
    infrastructure : string
        download infrastructure of given type (default is streets (ie, 'way["highway"]') but other
        infrastructures may be selected like power grids (ie, 'way["power"~"line"]'))
    custom_filter : string
        a custom network filter to be used instead of the network_type presets

    Returns
    -------
    networkx multidigraph
    """

    if clean_periphery and simplify:
        # create a new buffered bbox 0.5km around the desired one
        buffer_dist = 500
        polygon = Polygon([(west, north), (west, south), (east, south), (east, north)])
        polygon_utm, crs_utm = project_geometry(geometry=polygon)
        polygon_proj_buff = polygon_utm.buffer(buffer_dist)
        polygon_buff, _ = project_geometry(geometry=polygon_proj_buff, crs=crs_utm, to_latlong=True)
        west_buffered, south_buffered, east_buffered, north_buffered = polygon_buff.bounds

        # get the network data from OSM then create the graph
        response_jsons = osm_net_download(north=north_buffered, south=south_buffered,
                                          east=east_buffered, west=west_buffered,
                                          network_type=network_type, timeout=timeout,
                                          memory=memory, max_query_area_size=max_query_area_size,
                                          infrastructure=infrastructure, custom_filter=custom_filter)
        G_buffered = create_graph(response_jsons, name=name, retain_all=retain_all,
                                  bidirectional=network_type in settings.bidirectional_network_types)
        G = truncate_graph_bbox(G_buffered, north, south, east, west, retain_all=True, truncate_by_edge=truncate_by_edge)

        # simplify the graph topology
        G_buffered = simplify_graph(G_buffered)

        # truncate graph by desired bbox to return the graph within the bbox
        # caller wants
        G = truncate_graph_bbox(G_buffered, north, south, east, west, retain_all=retain_all, truncate_by_edge=truncate_by_edge)

        # count how many street segments in buffered graph emanate from each
        # intersection in un-buffered graph, to retain true counts for each
        # intersection, even if some of its neighbors are outside the bbox
        G.graph['streets_per_node'] = count_streets_per_node(G_buffered, nodes=G.nodes())

    else:
        # get the network data from OSM
        response_jsons = osm_net_download(north=north, south=south, east=east,
                                          west=west, network_type=network_type,
                                          timeout=timeout, memory=memory,
                                          max_query_area_size=max_query_area_size,
                                          infrastructure=infrastructure, custom_filter=custom_filter)

        # create the graph, then truncate to the bounding box
        G = create_graph(response_jsons, name=name, retain_all=retain_all,
                         bidirectional=network_type in settings.bidirectional_network_types)
        G = truncate_graph_bbox(G, north, south, east, west, retain_all=retain_all, truncate_by_edge=truncate_by_edge)

        # simplify the graph topology as the last step. don't truncate after
        # simplifying or you may have simplified out to an endpoint
        # beyond the truncation distance, in which case you will then strip out
        # your entire edge
        if simplify:
            G = simplify_graph(G)

    log('graph_from_bbox() returning graph with {:,} nodes and {:,} edges'.format(len(list(G.nodes())), len(list(G.edges()))))
    return  G


def graph_from_point(center_point, distance=1000, distance_type='bbox',
                     network_type='all_private', simplify=True, retain_all=False,
                     truncate_by_edge=False, name='unnamed', timeout=180,
                     memory=None, max_query_area_size=50*1000*50*1000,
                     clean_periphery=True, infrastructure='way["highway"]',
                     custom_filter=None):
    """
    Create a networkx graph from OSM data within some distance of some (lat,
    lon) center point.

    Parameters
    ----------
    center_point : tuple
        the (lat, lon) central point around which to construct the graph
    distance : int
        retain only those nodes within this many meters of the center of the
        graph, with distance determined according to distance_type argument
    distance_type : string
        {'network', 'bbox'} if 'bbox', retain only those nodes within a bounding
        box of the distance parameter. if 'network', retain only those nodes
        within some network distance from the center-most node.
    network_type : string
        what type of street network to get
    simplify : bool
        if true, simplify the graph topology
    retain_all : bool
        if True, return the entire graph even if it is not connected
    truncate_by_edge : bool
        if True retain node if it's outside bbox but at least one of node's
        neighbors are within bbox
    name : string
        the name of the graph
    timeout : int
        the timeout interval for requests and to pass to API
    memory : int
        server memory allocation size for the query, in bytes. If none, server
        will use its default allocation size
    max_query_area_size : float
        max size for any part of the geometry, in square degrees: any polygon
        bigger will get divided up for multiple queries to API
    clean_periphery : bool,
        if True (and simplify=True), buffer 0.5km to get a graph larger than
        requested, then simplify, then truncate it to requested spatial extent
    infrastructure : string
        download infrastructure of given type (default is streets (ie, 'way["highway"]') but other
        infrastructures may be selected like power grids (ie, 'way["power"~"line"]'))
    custom_filter : string
        a custom network filter to be used instead of the network_type presets

    Returns
    -------
    networkx multidigraph
    """

    if distance_type not in ['bbox', 'network']:
        raise InvalidDistanceType('distance_type must be "bbox" or "network"')

    # create a bounding box from the center point and the distance in each
    # direction
    north, south, east, west = bbox_from_point(center_point, distance)

    # create a graph from the bounding box
    G = graph_from_bbox(north, south, east, west, network_type=network_type, simplify=simplify,
                        retain_all=retain_all, truncate_by_edge=truncate_by_edge, name=name,
                        timeout=timeout, memory=memory, max_query_area_size=max_query_area_size,
                        clean_periphery=clean_periphery, infrastructure=infrastructure,
                        custom_filter=custom_filter)

    # if the network distance_type is network, find the node in the graph
    # nearest to the center point, and truncate the graph by network distance
    # from this node
    if distance_type == 'network':
        centermost_node = get_nearest_node(G, center_point)
        G = truncate_graph_dist(G, centermost_node, max_distance=distance)

    log('graph_from_point() returning graph with {:,} nodes and {:,} edges'.format(len(list(G.nodes())), len(list(G.edges()))))
    return G


def graph_from_address(address, distance=1000, distance_type='bbox',
                       network_type='all_private', simplify=True, retain_all=False,
                       truncate_by_edge=False, return_coords=False,
                       name='unnamed', timeout=180, memory=None,
                       max_query_area_size=50*1000*50*1000,
                       clean_periphery=True, infrastructure='way["highway"]',
                       custom_filter=None):
    """
    Create a networkx graph from OSM data within some distance of some address.

    Parameters
    ----------
    address : string
        the address to geocode and use as the central point around which to
        construct the graph
    distance : int
        retain only those nodes within this many meters of the center of the
        graph
    distance_type : string
        {'network', 'bbox'} if 'bbox', retain only those nodes within a bounding
        box of the distance parameter.
        if 'network', retain only those nodes within some network distance from
        the center-most node.
    network_type : string
        what type of street network to get
    simplify : bool
        if true, simplify the graph topology
    retain_all : bool
        if True, return the entire graph even if it is not connected
    truncate_by_edge : bool
        if True retain node if it's outside bbox but at least one of node's
        neighbors are within bbox
    return_coords : bool
        optionally also return the geocoded coordinates of the address
    name : string
        the name of the graph
    timeout : int
        the timeout interval for requests and to pass to API
    memory : int
        server memory allocation size for the query, in bytes. If none, server
        will use its default allocation size
    max_query_area_size
        float, max size for any part of the geometry, in square degrees: any
        polygon bigger will get divided up for multiple queries to API
    clean_periphery : bool,
        if True (and simplify=True), buffer 0.5km to get a graph larger than
        requested, then simplify, then truncate it to requested spatial extent
    infrastructure : string
        download infrastructure of given type (default is streets (ie, 'way["highway"]') but other
        infrastructures may be selected like power grids (ie, 'way["power"~"line"]'))
    custom_filter : string
        a custom network filter to be used instead of the network_type presets

    Returns
    -------
    networkx multidigraph or tuple
        multidigraph or optionally (multidigraph, tuple)
    """

    # geocode the address string to a (lat, lon) point
    point = geocode(query=address)

    # then create a graph from this point
    G = graph_from_point(point, distance, distance_type, network_type=network_type,
                         simplify=simplify, retain_all=retain_all, truncate_by_edge=truncate_by_edge,
                         name=name, timeout=timeout, memory=memory,
                         max_query_area_size=max_query_area_size,
                         clean_periphery=clean_periphery, infrastructure=infrastructure,
                         custom_filter=custom_filter)
    log('graph_from_address() returning graph with {:,} nodes and {:,} edges'.format(len(list(G.nodes())), len(list(G.edges()))))

    if return_coords:
        return G, point
    else:
        return G


def graph_from_polygon(polygon, network_type='all_private', simplify=True,
                       retain_all=False, truncate_by_edge=False, name='unnamed',
                       timeout=180, memory=None,
                       max_query_area_size=50*1000*50*1000,
                       clean_periphery=True, infrastructure='way["highway"]',
                       custom_filter=None):
    """
    Create a networkx graph from OSM data within the spatial boundaries of the
    passed-in shapely polygon.

    Parameters
    ----------
    polygon : shapely Polygon or MultiPolygon
        the shape to get network data within. coordinates should be in units of
        latitude-longitude degrees.
    network_type : string
        what type of street network to get
    simplify : bool
        if true, simplify the graph topology
    retain_all : bool
        if True, return the entire graph even if it is not connected
    truncate_by_edge : bool
        if True retain node if it's outside bbox but at least one of node's
        neighbors are within bbox
    name : string
        the name of the graph
    timeout : int
        the timeout interval for requests and to pass to API
    memory : int
        server memory allocation size for the query, in bytes. If none, server
        will use its default allocation size
    max_query_area_size : float
        max size for any part of the geometry, in square degrees: any polygon
        bigger will get divided up for multiple queries to API
    clean_periphery : bool
        if True (and simplify=True), buffer 0.5km to get a graph larger than
        requested, then simplify, then truncate it to requested spatial extent
    infrastructure : string
        download infrastructure of given type (default is streets
        (ie, 'way["highway"]') but other infrastructures may be selected
        like power grids (ie, 'way["power"~"line"]'))
    custom_filter : string
        a custom network filter to be used instead of the network_type presets

    Returns
    -------
    networkx multidigraph
    """

    # verify that the geometry is valid and is a shapely Polygon/MultiPolygon
    # before proceeding
    if not polygon.is_valid:
        raise TypeError('Shape does not have a valid geometry')
    if not isinstance(polygon, (Polygon, MultiPolygon)):
        raise TypeError('Geometry must be a shapely Polygon or MultiPolygon. If you requested '
                         'graph from place name or address, make sure your query resolves to a '
                         'Polygon or MultiPolygon, and not some other geometry, like a Point. '
                         'See OSMnx documentation for details.')

    if clean_periphery and simplify:
        # create a new buffered polygon 0.5km around the desired one
        buffer_dist = 500
        polygon_utm, crs_utm = project_geometry(geometry=polygon)
        polygon_proj_buff = polygon_utm.buffer(buffer_dist)
        polygon_buffered, _ = project_geometry(geometry=polygon_proj_buff, crs=crs_utm, to_latlong=True)

        # get the network data from OSM,  create the buffered graph, then
        # truncate it to the buffered polygon
        response_jsons = osm_net_download(polygon=polygon_buffered, network_type=network_type,
                                          timeout=timeout, memory=memory,
                                          max_query_area_size=max_query_area_size,
                                          infrastructure=infrastructure, custom_filter=custom_filter)
        G_buffered = create_graph(response_jsons, name=name, retain_all=True,
                                  bidirectional=network_type in settings.bidirectional_network_types)
        G_buffered = truncate_graph_polygon(G_buffered, polygon_buffered, retain_all=True, truncate_by_edge=truncate_by_edge)

        # simplify the graph topology
        G_buffered = simplify_graph(G_buffered)

        # truncate graph by polygon to return the graph within the polygon that
        # caller wants. don't simplify again - this allows us to retain
        # intersections along the street that may now only connect 2 street
        # segments in the network, but in reality also connect to an
        # intersection just outside the polygon
        G = truncate_graph_polygon(G_buffered, polygon, retain_all=retain_all, truncate_by_edge=truncate_by_edge)

        # count how many street segments in buffered graph emanate from each
        # intersection in un-buffered graph, to retain true counts for each
        # intersection, even if some of its neighbors are outside the polygon
        G.graph['streets_per_node'] = count_streets_per_node(G_buffered, nodes=G.nodes())

    else:
        # download a list of API responses for the polygon/multipolygon
        response_jsons = osm_net_download(polygon=polygon, network_type=network_type,
                                          timeout=timeout, memory=memory,
                                          max_query_area_size=max_query_area_size,
                                          infrastructure=infrastructure, custom_filter=custom_filter)

        # create the graph from the downloaded data
        G = create_graph(response_jsons, name=name, retain_all=True,
                         bidirectional=network_type in settings.bidirectional_network_types)

        # truncate the graph to the extent of the polygon
        G = truncate_graph_polygon(G, polygon, retain_all=retain_all, truncate_by_edge=truncate_by_edge)

        # simplify the graph topology as the last step. don't truncate after
        # simplifying or you may have simplified out to an endpoint beyond the
        # truncation distance, in which case you will then strip out your entire
        # edge
        if simplify:
            G = simplify_graph(G)

    log('graph_from_polygon() returning graph with {:,} nodes and {:,} edges'.format(len(list(G.nodes())), len(list(G.edges()))))
    return G


def graph_from_place(query, network_type='all_private', simplify=True,
                     retain_all=False, truncate_by_edge=False, name='unnamed',
                     which_result=1, buffer_dist=None, timeout=180, memory=None,
                     max_query_area_size=50*1000*50*1000, clean_periphery=True,
                     infrastructure='way["highway"]', custom_filter=None):
    """
    Create a networkx graph from OSM data within the spatial boundaries of some
    geocodable place(s).

    The query must be geocodable and OSM must have polygon boundaries for the
    geocode result. If OSM does not have a polygon for this place, you can
    instead get its street network using the graph_from_address function, which
    geocodes the place name to a point and gets the network within some distance
    of that point. Alternatively, you might try to vary the which_result
    parameter to use a different geocode result. For example, the first geocode
    result (ie, the default) might resolve to a point geometry, but the second
    geocode result for this query might resolve to a polygon, in which case you
    can use graph_from_place with which_result=2.

    Parameters
    ----------
    query : string or dict or list
        the place(s) to geocode/download data for
    network_type : string
        what type of street network to get
    simplify : bool
        if true, simplify the graph topology
    retain_all : bool
        if True, return the entire graph even if it is not connected
    truncate_by_edge : bool
        if True retain node if it's outside bbox but at least one of node's
        neighbors are within bbox
    name : string
        the name of the graph
    which_result : int
        max number of results to return and which to process upon receipt
    buffer_dist : float
        distance to buffer around the place geometry, in meters
    timeout : int
        the timeout interval for requests and to pass to API
    memory : int
        server memory allocation size for the query, in bytes. If none, server
        will use its default allocation size
    max_query_area_size : float
        max size for any part of the geometry, in square degrees: any polygon
        bigger will get divided up for multiple queries to API
    clean_periphery : bool
        if True (and simplify=True), buffer 0.5km to get a graph larger than
        requested, then simplify, then truncate it to requested spatial extent
    infrastructure : string
        download infrastructure of given type (default is streets (ie, 'way["highway"]') but other
        infrastructures may be selected like power grids (ie, 'way["power"~"line"]'))
    custom_filter : string
        a custom network filter to be used instead of the network_type presets

    Returns
    -------
    networkx multidigraph
    """

    # create a GeoDataFrame with the spatial boundaries of the place(s)
    if isinstance(query, str) or isinstance(query, dict):
        # if it is a string (place name) or dict (structured place query), then
        # it is a single place
        gdf_place = gdf_from_place(query, which_result=which_result, buffer_dist=buffer_dist)
        name = query
    elif isinstance(query, list):
        # if it is a list, it contains multiple places to get
        gdf_place = gdf_from_places(query, buffer_dist=buffer_dist)
    else:
        raise TypeError('query must be a string or a list of query strings')

    # extract the geometry from the GeoDataFrame to use in API query
    polygon = gdf_place['geometry'].unary_union
    log('Constructed place geometry polygon(s) to query API')

    # create graph using this polygon(s) geometry
    G = graph_from_polygon(polygon, network_type=network_type, simplify=simplify,
                           retain_all=retain_all, truncate_by_edge=truncate_by_edge,
                           name=name, timeout=timeout, memory=memory,
                           max_query_area_size=max_query_area_size,
                           clean_periphery=clean_periphery, infrastructure=infrastructure,
                           custom_filter=custom_filter)

    log('graph_from_place() returning graph with {:,} nodes and {:,} edges'.format(len(list(G.nodes())), len(list(G.edges()))))
    return G


def graph_from_file(filename, bidirectional=False, simplify=True,
                    retain_all=False, name='unnamed'):
    """
    Create a networkx graph from OSM data in an XML file.

    Parameters
    ----------
    filename : string
        the name of a file containing OSM XML data
    bidirectional : bool
        if True, create bidirectional edges for one-way streets
    simplify : bool
        if True, simplify the graph topology
    retain_all : bool
        if True, return the entire graph even if it is not connected
    name : string
        the name of the graph

    Returns
    -------
    networkx multidigraph
    """
    # transmogrify file of OSM XML data into JSON
    response_jsons = [overpass_json_from_file(filename)]

    # create graph using this response JSON
    G = create_graph(response_jsons, bidirectional=bidirectional,
                     retain_all=retain_all, name=name)

    # simplify the graph topology as the last step.
    if simplify:
        G = simplify_graph(G)

    log('graph_from_file() returning graph with {:,} nodes and {:,} edges'.format(len(list(G.nodes())), len(list(G.edges()))))
    return G
