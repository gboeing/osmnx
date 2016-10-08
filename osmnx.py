# The MIT License (MIT)
# 
# Copyright (c) 2016 Geoff Boeing http://geoffboeing.com
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
###################################################################################################
# Module: OSMnx.py
# Description: Retrieve and construct spatial geometries and street networks from OpenStreetMap
###################################################################################################

import json, math, sys, os, io, ast, unicodedata, hashlib, re, random, time, datetime as dt, logging as lg
from collections import OrderedDict
from dateutil import parser as date_parser
import requests, numpy as np, pandas as pd, geopandas as gpd, networkx as nx, matplotlib.pyplot as plt, matplotlib.cm as cm
from matplotlib.collections import LineCollection
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from shapely import wkt
from shapely.ops import unary_union
from geopy.distance import great_circle, vincenty
from geopy.geocoders import Nominatim


###################################################################################################
# global defaults
# you can edit any of these by passing the value to the config() function
###################################################################################################


# default locations to save data, logs, images, and cache
_data_folder = 'data'
_logs_folder = 'logs'
_imgs_folder = 'images'
_cache_folder = 'cache'

_use_cache = False

# write log to file and/or to console
_log_file = False
_log_console = False
_log_level = lg.INFO
_log_name = 'osmnx'
_log_filename = 'osmnx'

# useful osm tags - note that load_graphml expects a consistent set of tag names for parsing
_useful_tags_node = ['ref', 'highway']
_useful_tags_path = ['bridge', 'tunnel', 'oneway', 'lanes', 'ref', 'name', 'highway', 'maxspeed', 'service', 'access']


###################################################################################################
# define functions
###################################################################################################


def config(data_folder=_data_folder, logs_folder=_logs_folder, imgs_folder=_imgs_folder, 
           cache_folder=_cache_folder, use_cache=_use_cache,
           log_file=_log_file, log_console=_log_console, log_level=_log_level, log_name=_log_name, log_filename=_log_filename,
           useful_tags_node=_useful_tags_node, useful_tags_path=_useful_tags_path):
    """
    Configure osmnx by setting the default global vars to desired values.
    
    Parameters
    ---------
    data_folder : string, where to save and load data files
    logs_folder : string, where to write the log files
    imgs_folder : string, where to save figures
    cache_folder : string, where to save the http response cache
    use_cache : bool, if True, use a local cache to save/retrieve http responses instead of calling API repetitively for the same request URL
    log_file : bool, if true, save log output to a log file in logs_folder
    log_console : bool, if true, print log output to the console
    log_level : int, one of the logger.level constants
    log_name : string, name of the logger
    useful_tags_node : list, a list of useful OSM tags to attempt to save from node elements
    useful_tags_path : list, a list of useful OSM tags to attempt to save from path elements
    
    Returns
    -------
    None
    """
    # use the global keyword to update global variables from within this function
    global _use_cache, _cache_folder, _data_folder, _imgs_folder, _logs_folder, _log_console, _log_file, _log_level, _log_name, _log_filename, _useful_tags_node, _useful_tags_path
    
    # set each global variable to the passed-in parameter value
    _use_cache = use_cache
    _cache_folder = cache_folder
    _data_folder = data_folder
    _imgs_folder = imgs_folder
    _logs_folder = logs_folder
    _log_console = log_console
    _log_file = log_file
    _log_level = log_level
    _log_name = log_name
    _log_filename = log_filename
    _useful_tags_node = useful_tags_node
    _useful_tags_path = useful_tags_path
    
    # if logging is turned on, log that we are configured
    if _log_file or _log_console:
        log('Configured osmnx')


def log(message, level=None, name=None, filename=None):
    """
    Write a message to the log file and/or print to the the console.
    
    Parameters
    ----------
    message : string, the content of the message to log
    level : int, one of the logger.level constants
    name : string, name of the logger
    filename : string, name of the log file
    
    Returns
    -------
    None
    """
    
    if level is None:
        level = _log_level
    if name is None:
        name = _log_name
    if filename is None:
        filename = _log_filename
    
    # if logging to file is turned on
    if _log_file:
        # get the current logger (or create a new one, if none), then log message at requested level
        logger = get_logger(level=level, name=name, filename=filename)
        if level == lg.DEBUG:
            logger.debug(message)
        elif level == lg.INFO:
            logger.info(message)
        elif level == lg.WARNING:
            logger.warning(message)
        elif level == lg.ERROR:
            logger.error(message)
    
    # if logging to console is turned on, convert message to ascii and print to the console
    if _log_console:
        # capture current stdout, then switch it to the console, print the message, then switch back to what had been the stdout
        # this prevents logging to notebook - instead, it goes to console
        standard_out = sys.stdout
        sys.stdout = sys.__stdout__
        
        # convert message to ascii for console display so it doesn't break windows terminals
        message = unicodedata.normalize('NFKD', make_str(message)).encode('ascii', errors='replace').decode()
        print(message)
        sys.stdout = standard_out


def get_logger(level=None, name=None, filename=None):
    """
    Create a logger or return the current one if already instantiated.
    
    Parameters
    ----------
    level : int, one of the logger.level constants
    name : string, name of the logger
    filename : string, name of the log file
    
    Returns
    -------
    logger : logger.logger
    """
    
    if level is None:
        level = _log_level
    if name is None:
        name = _log_name
    if filename is None:
        filename = _log_filename
    
    logger = lg.getLogger(name)
    
    # if a logger with this name is not already set up
    if not getattr(logger, 'handler_set', None):
        
        # get today's date and construct a log filename
        todays_date = dt.datetime.today().strftime('%Y_%m_%d')
        log_filename = '{}/{}_{}.log'.format(_logs_folder, filename, todays_date)
        
        # if the logs folder does not already exist, create it
        if not os.path.exists(_logs_folder):
            os.makedirs(_logs_folder)
            
        # create file handler and log formatter and set them up
        handler = lg.FileHandler(log_filename, encoding='utf-8')
        formatter = lg.Formatter('%(asctime)s %(levelname)s %(name)s %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(level)
        logger.handler_set = True
    
    return logger    


def make_str(value):
    """
    Convert a passed-in value to unicode if Python 2, or string if Python 3.
    
    Parameters
    ----------
    value : any type, the value to convert to unicode/string
    
    Returns
    -------
    value : unicode or string
    """
    try:
        # for python 2.x compatibility, use unicode
        return unicode(value)
    except:
        # python 3.x has no unicode type, so if error, use str type
        return str(value)
    
    
def save_to_cache(url, response_json):
    """
    Save a HTTP response json object to the cache. If the request was sent to server via POST instead of GET, 
    then URL should be a GET-style representation of request. Users should always pass OrderedDicts instead of dicts
    of parameters into request functions, so that the parameters stay in the same order each time, producing the same URL string,
    and thus the same hash. Otherwise the cache will eventually contain multiple saved responses for the same
    request because the URL's parameters appeared in a different order each time.
    
    Parameters
    ----------
    url : string, the url of the request
    response_json : dict, the json response
    
    Returns
    -------
    None
    """
    if _use_cache:
        if response_json is None:
            log('Saved nothing to cache because response_json is None')
        else:        
            # create the folder on the disk if it doesn't already exist
            if not os.path.exists(_cache_folder):
                os.makedirs(_cache_folder)

            # hash the url (to make filename shorter than the often extremely long url) 
            filename = hashlib.md5(url.encode('utf-8')).hexdigest()
            cache_path_filename = '{}/{}.json'.format(_cache_folder, filename)
            
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
    url : string, the url of the request
    
    Returns
    -------
    response_json : dict
    """
    # if the tool is configured to use the cache
    if _use_cache:
        # determine the filename by hashing the url
        filename = hashlib.md5(url.encode('utf-8')).hexdigest()
        cache_path_filename = '{}/{}.json'.format(_cache_folder, filename)
        # open the cache file for this url hash if it already exists, otherwise return None
        if os.path.isfile(cache_path_filename):
            response_json = json.load(io.open(cache_path_filename, encoding='utf-8'))
            log('Retrieved response from cache file "{}" for URL "{}"'.format(cache_path_filename, url))
            return response_json


def get_pause_duration(recursive_delay=5):
    """
    Check the Overpass API status endpoint to determine how long to wait until next slot is available.
    
    Parameters
    ----------
    recursive_delay : int, how long to wait between recursive calls if server is currently running a query
    
    Returns
    -------
    pause_duration : int
    """
    
    response = requests.get('http://overpass-api.de/api/status')
    status = response.text.split('\n')[2]
    status_first_token = status.split(' ')[0]

    try:
        # if first token is numeric, it's how many slots you have available - no wait required
        available_slots = int(status_first_token)
        pause_duration = 0
    except:
        # if first token is 'Slot', it tells you when your slot will be free
        if status_first_token == 'Slot':
            utc_time_str = status.split(' ')[3]
            utc_time = date_parser.parse(utc_time_str).replace(tzinfo=None)
            pause_duration = math.ceil((utc_time - dt.datetime.utcnow()).total_seconds())
            pause_duration = max(pause_duration, 1)
        
        # if first token is 'Currently', it is currently running a query so check back in recursive_delay seconds
        elif status_first_token == 'Currently':
            time.sleep(recursive_delay)
            pause_duration = get_pause_duration()
            
        else:
            raise ValueError('Unrecognized server status: "{}"'.format(status))
    
    return pause_duration


def nominatim_request(params, pause_duration=1, timeout=30, error_pause_duration=180):
    """
    Send a request to the Nominatim API via HTTP GET and return the JSON response.
    
    Parameters
    ----------
    params : dict or OrderedDict, key-value pairs of parameters
    pause_duration : int, how long to pause before requests, in seconds
    timeout : int, the timeout interval for the requests library
    error_pause_duration : int, how long to pause in seconds before re-trying requests if error
    
    Returns
    -------
    response_json : dict
    """
    
    # prepare the Nominatim API URL and see if request already exists in the cache
    url = 'https://nominatim.openstreetmap.org/search'
    prepared_url = requests.Request('GET', url, params=params).prepare().url
    cached_response_json = get_from_cache(prepared_url)
    
    if not cached_response_json is None:
        # found this request in the cache, just return it instead of making a new HTTP call
        return cached_response_json
    
    else:
        # if this URL is not already in the cache, pause, then request it
        log('Pausing {:,.2f} seconds before making API GET request'.format(pause_duration))
        time.sleep(pause_duration)
        start_time = time.time()
        log('Requesting {} with timeout={}'.format(prepared_url, timeout))
        response = requests.get(url, params, timeout=timeout)
        
        # get the response size and the domain, log result
        size_kb = len(response.content) / 1000.
        domain = re.findall(r'//(?s)(.*?)/', url)[0]
        log('Downloaded {:,.1f}KB from {} in {:,.2f} seconds'.format(size_kb, domain, time.time()-start_time))
        
        try:
            response_json = response.json()
            save_to_cache(prepared_url, response_json)
        except:
            #429 is 'too many requests' and 504 is 'gateway timeout' from server overload - handle these errors by recursively calling nominatim_request until we get a valid response
            if response.status_code in [429, 504]:
                # pause for error_pause_duration seconds before re-trying request
                log('Server at {} returned status code {} and no JSON data. Re-trying request in {:.2f} seconds.'.format(domain, response.status_code, error_pause_duration), level=lg.WARNING)
                time.sleep(error_pause_duration)
                response_json = nominatim_request(params=params, pause_duration=pause_duration, timeout=timeout)
            
            # else, this was an unhandled status_code, throw an exception
            else:
                log('Server at {} returned status code {} and no JSON data'.format(domain, response.status_code), level=lg.ERROR)
                raise Exception('Server returned no JSON data.\n{} {}\n{}'.format(response, response.reason, response.text))
        
        return response_json


def overpass_request(data, pause_duration=None, timeout=180, error_pause_duration=None):
    """
    Send a request to the Overpass API via HTTP POST and return the JSON response
    
    Parameters
    ----------
    data : dict or OrderedDict, key-value pairs of parameters to post to the API
    pause_duration : int, how long to pause in seconds before requests, if None, will query API status endpoint to find when next slot is available
    timeout : int, the timeout interval for the requests library
    error_pause_duration : int, how long to pause in seconds before re-trying requests if error
    
    Returns
    -------
    response_json : dict
    """
    
    # define the Overpass API URL, then construct a GET-style URL as a string to hash to look up/save to cache
    url = 'http://www.overpass-api.de/api/interpreter'
    prepared_url = requests.Request('GET', url, params=data).prepare().url
    cached_response_json = get_from_cache(prepared_url)
    
    if not cached_response_json is None:
        # found this request in the cache, just return it instead of making a new HTTP call
        return cached_response_json
    
    else:
        # if this URL is not already in the cache, pause, then request it
        if pause_duration is None:
            this_pause_duration = get_pause_duration()
        log('Pausing {:,.2f} seconds before making API POST request'.format(this_pause_duration))
        time.sleep(this_pause_duration)
        start_time = time.time()
        log('Posting to {} with timeout={}, "{}"'.format(url, timeout, data))
        response = requests.post(url, data=data, timeout=timeout)
        
        # get the response size and the domain, log result
        size_kb = len(response.content) / 1000.
        domain = re.findall(r'//(?s)(.*?)/', url)[0]
        log('Downloaded {:,.1f}KB from {} in {:,.2f} seconds'.format(size_kb, domain, time.time()-start_time))
        
        try:
            response_json = response.json()
            if 'remark' in response_json:
                log('Server remark: "{}"'.format(response_json['remark'], level=lg.WARNING))
            save_to_cache(prepared_url, response_json)
        except:
            #429 is 'too many requests' and 504 is 'gateway timeout' from server overload - handle these errors by recursively calling overpass_request until we get a valid response
            if response.status_code in [429, 504]:
                # pause for error_pause_duration seconds before re-trying request
                if error_pause_duration is None:
                    error_pause_duration = get_pause_duration()
                log('Server at {} returned status code {} and no JSON data. Re-trying request in {:.2f} seconds.'.format(domain, response.status_code, error_pause_duration), level=lg.WARNING)
                time.sleep(error_pause_duration)
                response_json = overpass_request(data=data, pause_duration=pause_duration, timeout=timeout)
            
            # else, this was an unhandled status_code, throw an exception
            else:
                log('Server at {} returned status code {} and no JSON data'.format(domain, response.status_code), level=lg.ERROR)
                raise Exception('Server returned no JSON data.\n{} {}\n{}'.format(response, response.reason, response.text))
        
        return response_json
    
    
def osm_polygon_download(query, limit=1, polygon_geojson=1, pause_duration=1):
    """
    Geocode a place and download its boundary geometry from OSM's Nominatim API.
    
    Parameters
    ----------
    query : string or dict, query string or structured query dict to geocode/download
    limit : int, max number of results to return
    polygon_geojson : int, request the boundary geometry polygon from the API, 0=no, 1=yes
    pause_duration : int, time in seconds to pause before API requests
    
    Returns
    -------
    response_json : dict
    """
    # define the parameters
    params = OrderedDict()
    params['format'] = 'json'
    params['limit'] = limit
    params['dedupe'] = 0 #this prevents OSM from de-duping results so we're guaranteed to get precisely 'limit' number of results
    params['polygon_geojson'] = polygon_geojson
    
    # add the structured query dict (if provided) to params, otherwise query with place name string
    if isinstance(query, str):
        params['q'] = query
    elif isinstance(query, dict):
        # add the query keys in alphabetical order so the URL is the same string each time, for caching purposes
        for key in sorted(list(query.keys())):
            params[key] = query[key]
    else:
        raise ValueError('query must be a dict or a string')
    
    # request the URL, return the JSON
    response_json = nominatim_request(params=params, timeout=30)
    return response_json
    

def gdf_from_place(query, gdf_name=None, which_result=1, buffer_dist=None):
    """
    Create a GeoDataFrame from a single place name query.
    
    Parameters
    ----------
    query : string or dict, query string or structured query dict to geocode/download
    gdf_name : string, name attribute metadata for GeoDataFrame (this is used to save shapefile later)
    which_result : int, max number of results to return and which to process upon receipt
    buffer_dist : float, distance to buffer around the place geometry, in meters
    
    Returns
    -------
    gdf : GeoDataFrame
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
        
        # create the GeoDataFrame, name it, and set its original CRS to lat-long (EPSG 4326)
        gdf = gpd.GeoDataFrame.from_features(features)
        gdf.name = gdf_name
        gdf.crs = {'init':'epsg:4326'}
        
        # if buffer_dist was passed in, project the geometry to UTM, buffer it in meters, then project it back to lat-long
        if not buffer_dist is None:
            gdf_utm = project_gdf(gdf)
            gdf_utm['geometry'] = gdf_utm['geometry'].buffer(buffer_dist)
            gdf = project_gdf(gdf_utm, to_latlong=True)
            log('Buffered the GeoDataFrame "{}" to {} meters'.format(gdf.name, buffer_dist))
        
        # return the gdf
        log('Created GeoDataFrame with {} row for query "{}"'.format(len(gdf), query))
        return gdf
    else:
        # if there was no data returned (or fewer results than which_result specified)
        log('OSM returned no results (or fewer than which_result) for query "{}"'.format(query), level=lg.WARNING)
        gdf = gpd.GeoDataFrame()
        gdf.name = gdf_name
        return gdf
        
        
def gdf_from_places(queries, gdf_name='unnamed', buffer_dist=None):
    """
    Create a GeoDataFrame from a list of place names to query.
    
    Parameters
    ----------
    queries : list, list of query strings or structured query dicts to geocode/download, one at a time
    gdf_name : string, name attribute metadata for GeoDataFrame (this is used to save shapefile later)
    buffer_dist : float, distance to buffer around the place geometry, in meters
    
    Returns
    -------
    gdf : GeoDataFrame
    """
    # create an empty GeoDataFrame then append each result as a new row
    gdf = gpd.GeoDataFrame()
    for query in queries:
        gdf = gdf.append(gdf_from_place(query, buffer_dist=buffer_dist))
        
    # reset the index, name the GeoDataFrame
    gdf = gdf.reset_index().drop(labels='index', axis=1)
    gdf.name = gdf_name
    
    # set the original CRS of the GeoDataFrame to lat-long, and return it
    gdf.crs = {'init':'epsg:4326'}
    log('Finished creating GeoDataFrame with {} rows from {} queries'.format(len(gdf), len(queries)))
    return gdf


def make_shp_filename(place_name):
    """
    Create a filename string in a consistent format from a place name string.
    
    Parameters
    ----------
    place_name : string, place name to convert into a filename
    
    Returns
    -------
    filename : string
    """
    name_pieces = list(reversed(place_name.split(', ')))
    filename = '-'.join(name_pieces).lower().replace(' ','_')
    filename = re.sub('[^0-9a-zA-Z_-]+', '', filename)
    return filename


def save_gdf_shapefile(gdf, filename=None, folder=None):
    """
    Save GeoDataFrame as an ESRI shapefile.
    
    Parameters
    ----------
    gdf : GeoDataFrame, the gdf to be saved
    filename : string, what to call the shapefile (file extensions are added automatically)
    folder : string, where to save the shapefile, if none, then default folder
    
    Returns
    -------
    None
    """
    if folder is None:
        folder = _data_folder
    
    if filename is None:
        filename = make_shp_filename(gdf.name)
    
    # give the save folder a filename subfolder to make the full path to the files
    folder_path = '{}/{}'.format(folder, filename)
    
    # if the save folder does not already exist, create it with a filename subfolder
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    gdf.to_file(folder_path)
    log('Saved the GeoDataFrame "{}" as shapefile "{}"'.format(gdf.name, folder_path))
 

def project_gdf(gdf, to_latlong=False):
    """
    Project a GeoDataFrame to the UTM zone appropriate for its geometries' centroid. The simple calculation
    in this function works well for most latitudes, but won't work for some far northern locations like
    Svalbard and parts of far northern Norway.
    
    Parameters
    ----------
    gdf : GeoDataFrame, the gdf to be projected to UTM
    to_latlong : bool, if True, projects to latlong instead of to UTM
    
    Returns
    -------
    gdf : GeoDataFrame
    """
    assert len(gdf) > 0, 'You cannot project an empty GeoDataFrame.'
    start_time = time.time()
    
    if to_latlong:
        # if to_latlong is True, project the gdf to latlong
        latlong_crs = {'init':'epsg:4326'}
        projected_gdf = gdf.to_crs(latlong_crs)
        log('Projected the GeoDataFrame "{}" to EPSG 4326 in {:,.2f} seconds'.format(gdf.name, time.time()-start_time))
    else:
        # else, project the gdf to UTM
        # if GeoDataFrame is already in UTM, just return it
        if (not gdf.crs is None) and ('proj' in gdf.crs) and (gdf.crs['proj'] == 'utm'):
            return gdf
        
        # calculate the centroid of the union of all the geometries in the GeoDataFrame
        avg_longitude = gdf['geometry'].unary_union.centroid.x
        
        # calculate the UTM zone from this avg longitude and define the UTM CRS to project
        utm_zone = int(math.floor((avg_longitude + 180) / 6.) + 1)
        utm_crs = {'datum': 'NAD83',
                   'ellps': 'GRS80',
                   'proj' : 'utm',
                   'zone' : utm_zone,
                   'units': 'm'}
        
        # project the GeoDataFrame to the UTM CRS
        projected_gdf = gdf.to_crs(utm_crs)
        log('Projected the GeoDataFrame "{}" to UTM-{} in {:,.2f} seconds'.format(gdf.name, utm_zone, time.time()-start_time))
    
    projected_gdf.name = gdf.name
    return projected_gdf


    
##########################################################################################################        
#
# End of functions for getting place boundary geometries.
#
# Below are functions for getting and processing street networks.
#
#########################################################################################################    
   
 
 
def get_osm_filter(network_type):
    """
    Create a filter to query OSM for the specified network type.
    
    Parameters
    ----------
    network_type : string, {'walk', 'bike', 'drive', 'drive_service', 'all', 'all_private'} what type of street network to get
    
    Returns
    -------
    osm_filter : string
    """
    filters = {}

    # driving: filter out un-drivable roads, service roads, private ways, and anything specifying motor=no
    # also filter out any non-service roads that are tagged as providing parking, driveway, private, or emergency-access services
    filters['drive'] = ('["area"!~"yes"]["highway"!~"cycleway|footway|path|pedestrian|steps|track|'
                        'proposed|construction|bridleway|abandoned|platform|raceway|service"]'
                        '["motor_vehicle"!~"no"]["motorcar"!~"no"]["access"!~"private"]'
                        '["service"!~"parking|parking_aisle|driveway|private|emergency_access"]')

    # drive+service: allow ways tagged 'service' but filter out certain types of service ways
    filters['drive_service'] = ('["area"!~"yes"]["highway"!~"cycleway|footway|path|pedestrian|steps|track|'
                                'proposed|construction|bridleway|abandoned|platform|raceway"]'
                                '["motor_vehicle"!~"no"]["motorcar"!~"no"]["access"!~"private"]'
                                '["service"!~"parking|parking_aisle|private|emergency_access"]')
    
    # walking: filter out cycle ways, motor ways, private ways, and anything specifying foot=no
    # allow service roads, permitting things like parking lot lanes, alleys, etc that you *can* walk on even if they're not exactly nice walks
    filters['walk'] = ('["area"!~"yes"]["highway"!~"cycleway|motor|proposed|construction|abandoned|platform|raceway"]'
                       '["foot"!~"no"]["service"!~"private"]["access"!~"private"]')

    # biking: filter out foot ways, motor ways, private ways, and anything specifying biking=no
    filters['bike'] = ('["area"!~"yes"]["highway"!~"footway|motor|proposed|construction|abandoned|platform|raceway"]'
                       '["bicycle"!~"no"]["service"!~"private"]["access"!~"private"]')

    # to download all ways, just filter out everything not currently in use or that is private-access only
    filters['all'] = ('["area"!~"yes"]["highway"!~"proposed|construction|abandoned|platform|raceway"]'
                      '["service"!~"private"]["access"!~"private"]')
                      
    # to download all ways, including private-access ones, just filter out everything not currently in use
    filters['all_private'] = '["area"!~"yes"]["highway"!~"proposed|construction|abandoned|platform|raceway"]'
    
    if network_type in filters:
        osm_filter = filters[network_type]
    else:
        raise ValueError('unknown network_type "{}"'.format(network_type))
    
    return osm_filter
 
 
def osm_net_download(polygon=None, north=None, south=None, east=None, west=None, network_type='all_private', timeout=180, memory=None, max_query_area_size=0.4):
    """
    Download OSM ways and nodes within some bounding box from the Overpass API.
    
    Parameters
    ----------
    polygon : shapely Polygon or MultiPolygon, geographic shape to fetch the street network within
    north : float, northern latitude of bounding box
    south : float, southern latitude of bounding box
    east : float, eastern longitude of bounding box
    west : float, western longitude of bounding box
    network_type : string, {'walk', 'bike', 'drive', 'drive_service', 'all', 'all_private'} what type of street network to get
    timeout : int, the timeout interval for requests and to pass to API
    memory : int, server memory allocation size for the query, in bytes. If none, server will use its default allocation size
    max_query_area_size : float, max size for any part of the geometry, in square degrees: any polygon bigger will get divided up for multiple queries to API
    
    Returns
    -------
    response_json : dict
    """
    
    # check if we're querying by polygon or by bounding box based on which argument(s) where passed into this function
    by_poly = not polygon is None
    by_bbox = not (north is None or south is None or east is None or west is None)
    if not (by_poly or by_bbox):
        raise ValueError('You must pass a polygon or north, south, east, and west')
    
    # create a filter to exclude certain kinds of routes based on the requested network_type
    osm_filter = get_osm_filter(network_type)
    response_jsons = []
    
    # pass server memory allocation in bytes for the query to the API
    # if None, pass nothing so the server will use its default allocation size
    # otherwise, define the query's maxsize parameter value as whatever the caller passed in
    if memory is None:
        maxsize = ''
    else:
        maxsize = '[maxsize:{}]'.format(memory)
    
    # define the query to send the API
    # specifying way["highway"] means that all ways returned must have a highway key. the {filters} then remove ways by key/value.
    # the '>' makes it recurse so we get ways and way nodes. maxsize is in bytes.
    if by_bbox:
        # turn bbox into a polygon then subdivide it if it exceeds the max area size
        polygon = Polygon([(west, south), (east, south), (east, north), (west, north)])
        geometry = consolidate_subdivide_geometry(polygon, max_query_area_size=max_query_area_size)
        log('Requesting network data within bounding box from API in {:,} request(s)'.format(len(geometry)))
        start_time = time.time()
        
        # loop through each polygon rectangle in the geometry (there will only be one if original bbox didn't exceed max area size)
        for poly in geometry:
            # represent bbox as south,west,north,east and round lat-longs to 8 decimal places (ie, within 1 mm) so URL strings aren't different due to float rounding issues (for consistent caching)
            west, south, east, north = poly.bounds
            query_template = '[out:json][timeout:{timeout}]{maxsize};(way["highway"]{filters}({south:.8f},{west:.8f},{north:.8f},{east:.8f});>;);out;'
            query_str = query_template.format(north=north, south=south, east=east, west=west, filters=osm_filter, timeout=timeout, maxsize=maxsize)
            response_json = overpass_request(data={'data':query_str}, timeout=timeout)
            response_jsons.append(response_json)
        log('Got all network data within bounding box from API in {:,} request(s) and {:,.2f} seconds'.format(len(geometry), time.time()-start_time))
    
    elif by_poly:
        # divide polygon up into sub-polygons if area exceeds a max size, then get a list of polygon(s) exterior coordinates
        geometry = consolidate_subdivide_geometry(polygon, max_query_area_size=max_query_area_size)
        polygon_coord_strs = get_polygons_coordinates(geometry)
        log('Requesting network data within polygon from API in {:,} request(s)'.format(len(polygon_coord_strs)))
        start_time = time.time()
        
        # pass each polygon exterior coordinates in the list to the API, one at a time
        for polygon_coord_str in polygon_coord_strs:
            query_template = '[out:json][timeout:{timeout}]{maxsize};(way["highway"]{filters}(poly:"{polygon}");>;);out;'
            query_str = query_template.format(polygon=polygon_coord_str, filters=osm_filter, timeout=timeout, maxsize=maxsize)
            response_json = overpass_request(data={'data':query_str}, timeout=timeout)
            response_jsons.append(response_json)
        log('Got all network data within polygon from API in {:,} request(s) and {:,.2f} seconds'.format(len(polygon_coord_strs), time.time()-start_time))
        
    return response_jsons


def consolidate_subdivide_geometry(geometry, max_query_area_size=0.4):
    """
    Consolidate a geometry into a convex hull, then subdivide it into smaller sub-polygons if its area exceeds max size.
    
    Parameters
    ----------
    geometry : shapely Polygon or MultiPolygon, the geometry to consolidate and subdivide
    max_query_area_size : float, max size for any part of the geometry, in square degrees: any polygon bigger will get divided up for multiple queries to API
    
    Returns
    -------
    geometry : Polygon or MultiPolygon
    """
    
    # let the linear length of the quadrats (with which to subdivide the geometry) be the square root of max area size
    quadrat_length = math.sqrt(max_query_area_size)
    
    if not isinstance(geometry, (Polygon, MultiPolygon)):
        raise ValueError('Geometry must be a shapely Polygon or MultiPolygon')
    
    # if geometry is a MultiPolygon OR a single Polygon whose area exceeds the max size, get the convex hull around the geometry
    if isinstance(geometry, MultiPolygon) or (isinstance(geometry, Polygon) and geometry.area > max_query_area_size):
        geometry = geometry.convex_hull
        
    # if geometry area exceeds max size, subdivide it into smaller sub-polygons
    if geometry.area > max_query_area_size:
        geometry = quadrat_cut_geometry(geometry, quadrat_length=quadrat_length)
    
    if isinstance(geometry, Polygon):
        geometry = MultiPolygon([geometry])
        
    return geometry
    
    
def get_polygons_coordinates(geometry):
    """
    Extract exterior coordinates from polygon(s) to pass to OSM in a query by polygon.
    
    Parameters
    ----------
    geometry : shapely Polygon or MultiPolygon, the geometry to extract exterior coordinates from
    
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
        raise ValueError('Geometry must be a shapely Polygon or MultiPolygon')
        
    # convert the exterior coordinates of the polygon(s) to the string format the API expects
    polygon_coord_strs = []
    for coords in polygons_coords:
        s = ''
        separator = ' '
        for coord in list(coords):
            # round floating point lats and longs to ~1 nanometer, so we can hash and cache strings consistently
            s = '{}{}{:.14f}{}{:.14f}'.format(s, separator, coord[1], separator, coord[0])
        polygon_coord_strs.append(s.strip(separator))
    
    return polygon_coord_strs
    

def get_node(element):
    """
    Convert an OSM node element into the format for a networkx node.
    
    Parameters
    ----------
    element : dict, an OSM node element
    
    Returns
    -------
    node : dict
    """
    node = {}
    node['y'] = element['lat']
    node['x'] = element['lon']
    node['osmid'] = element['id']
    if 'tags' in element:
        for useful_tag in _useful_tags_node:
            if useful_tag in element['tags']:
                node[useful_tag] = element['tags'][useful_tag]
    return node
    
    
def get_path(element):
    """
    Convert an OSM way element into the format for a networkx graph path.
    
    Parameters
    ----------
    element : dict, an OSM way element
    
    Returns
    -------
    path : dict
    """
    path = {}
    path['osmid'] = element['id']
    path['nodes'] = element['nodes']
    if 'tags' in element:
        for useful_tag in _useful_tags_path:
            if useful_tag in element['tags']:
                path[useful_tag] = element['tags'][useful_tag] 
    return path
    
    
def parse_osm_nodes_paths(osm_data):
    """
    Construct dicts of nodes and paths with key=osmid and value=dict of attributes.
    
    Parameters
    ----------
    osm_data : dict, JSON response from from the Overpass API
    
    Returns
    -------
    nodes, paths : dict, dict
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
    Remove from a graph all the nodes that have no incident edges (ie, node degree = 0)
    
    Parameters
    ----------
    G : graph, the graph from which to remove nodes
    
    Returns
    -------
    G : graph
    """
    isolated_nodes = [node for node, degree in G.degree().items() if degree < 1]
    G.remove_nodes_from(isolated_nodes)
    log('Removed {:,} isolated nodes'.format(len(isolated_nodes)))
    return G
    
   
def get_largest_component(G, strongly=False):
    """
    Return the largest weakly or strongly connected component from a directed graph.
    
    Parameters
    ----------
    G : graph
    strongly : bool, if True, return the largest strongly instead of weakly connected component
    
    Returns
    -------
    G : graph
    """
    
    start_time = time.time()
    original_len = len(G.nodes())
    
    if strongly:
        # if the graph is not connected and caller did not request retain_all, retain only the largest strongly connected component
        if not nx.is_strongly_connected(G):
            G = max(nx.strongly_connected_component_subgraphs(G), key=len)
            log('Graph was not connected, retained only the largest strongly connected component ({:,} of {:,} total nodes) in {:.2f} seconds'.format(len(G.nodes()), original_len, time.time()-start_time))
    else:
        # if the graph is not connected and caller did not request retain_all, retain only the largest weakly connected component
        if not nx.is_weakly_connected(G):
            G = max(nx.weakly_connected_component_subgraphs(G), key=len)
            log('Graph was not connected, retained only the largest weakly connected component ({:,} of {:,} total nodes) in {:.2f} seconds'.format(len(G.nodes()), original_len, time.time()-start_time))

    return G
    

def truncate_graph_dist(G, source_node, max_distance=1000, weight='length', retain_all=False):
    """
    Remove everything further than some network distance from a specified node in graph.
    
    Parameters
    ----------
    G : graph
    source_node : int, the node in the graph from which to measure network distances to other nodes
    max_distance : int, remove every node in the graph greater than this distance from the source_node
    weight : string, how to weight the graph when measuring distance (default 'length' is how many meters long the edge is)
    retain_all : bool, if True, return the entire graph even if it is not connected
    
    Returns
    -------
    G : graph
    """
    
    # get the shortest distance between the node and every other node, then remove every node further than max_distance away
    start_time = time.time()
    G = G.copy()
    distances = nx.shortest_path_length(G, source=source_node, weight=weight)
    distant_nodes = {key:value for key, value in distances.items() if value > max_distance}
    G.remove_nodes_from(distant_nodes.keys())
    log('Truncated graph by weighted network distance in {:,.2f} seconds'.format(time.time()-start_time))
    
    # remove any isolated nodes and retain only the largest component (if retain_all is True)
    if not retain_all:
        G = remove_isolated_nodes(G)
        G = get_largest_component(G)
    
    return G
    
    
def truncate_graph_bbox(G, north, south, east, west, truncate_by_edge=False, retain_all=False):
    """
    Remove every node in graph that falls outside a bounding box. Needed because overpass returns entire ways that also 
    include nodes outside the bbox if the way (that is, a way with a single OSM ID) has a node inside the bbox at some point.
    
    Parameters
    ----------
    G : graph
    north : float, northern latitude of bounding box
    south : float, southern latitude of bounding box
    east : float, eastern longitude of bounding box
    west : float, western longitude of bounding box
    truncate_by_edge : bool, if True retain node if it's outside bbox but at least one of node's neighbors are within bbox
    retain_all : bool, if True, return the entire graph even if it is not connected
    
    Returns
    -------
    G : graph
    """
    start_time = time.time()
    G = G.copy()
    nodes_outside_bbox = []
    
    for node, data in G.nodes(data=True):
        if data['y'] > north or data['y'] < south or data['x'] > east or data['x'] < west:
            # this node is outside the bounding box
            if not truncate_by_edge:
                # if we're not truncating by edge, add node to list of nodes outside the bounding box
                nodes_outside_bbox.append(node)
            else:
                # if we're truncating by edge, see if any of node's neighbors are within bounding box
                any_neighbors_in_bbox = False
                neighbors = G.successors(node) + G.predecessors(node)
                for neighbor in neighbors:
                    x = G.node[neighbor]['x']
                    y = G.node[neighbor]['y']
                    if y < north and y > south and x < east and x > west:
                        any_neighbors_in_bbox = True

                # if none of its neighbors are within the bounding box, add node to list of nodes outside the bounding box
                if not any_neighbors_in_bbox:
                    nodes_outside_bbox.append(node)
    
    G.remove_nodes_from(nodes_outside_bbox)
    log('Truncated graph by bounding box in {:,.2f} seconds'.format(time.time()-start_time))
    
    # remove any isolated nodes and retain only the largest component (if retain_all is True)
    if not retain_all:
        G = remove_isolated_nodes(G)
        G = get_largest_component(G)

    return G
    

def quadrat_cut_geometry(geometry, quadrat_length=0.025, min_num=3, dist=1e-14):
    """
    Split a Polygon or MultiPolygon up into sub-polygons of a specified size, using quadrats.
    
    Parameters
    ----------
    geometry : shapely Polygon or MultiPolygon, the geometry to split up into smaller sub-polygons
    quadrat_length : the linear length in degrees of the quadrats with which to cut up the geometry
    min_num : the minimum number of linear quadrat lines (e.g., min_num=3 would produce a quadrat grid of 4 squares)
    dist : the distance to buffer the quadrat grid in degrees (e.g., 1e-14 degrees is no more than 1 nanometer)
    
    Returns
    -------
    multipoly : shapely MultiPolygon
    """
    
    # create n evenly spaced points between the min and max x and y bounds
    west, south, east, north = geometry.bounds
    x_num = math.ceil((east-west) / quadrat_length) + 1
    y_num = math.ceil((north-south) / quadrat_length) + 1
    x_points = np.linspace(west, east, num=max(x_num, min_num))
    y_points = np.linspace(south, north, num=max(y_num, min_num))

    # create a quadrat grid of lines at each of the evenly spaced points
    vertical_lines = [LineString([(x, y_points[0]), (x, y_points[-1])]) for x in x_points]
    horizont_lines = [LineString([(x_points[0], y), (x_points[-1], y)]) for y in y_points]
    lines = vertical_lines + horizont_lines

    # buffer each line to distance, take their union, then cut geometry into pieces by these quadrats
    lines_buffered = [line.buffer(dist) for line in lines]
    quadrats = unary_union(lines_buffered)
    multipoly = geometry.difference(quadrats)
    
    return multipoly
    
    
def intersect_index_quadrats(gdf, geometry, quadrat_length=0.025, min_num=3, dist=1e-14):
    """
    Intersect points with a polygon, using an r-tree spatial index and cutting the polygon up into
    smaller sub-polygons for r-tree acceleration.
    
    Parameters
    ----------
    gdf : GeoDataFrame, the set of points to intersect
    geometry : shapely Polygon or MultiPolygon, the geometry to intersect with the points
    quadrat_length : the linear length in degrees of the quadrats with which to cut up the geometry
    min_num : the minimum number of linear quadrat lines (e.g., min_num=3 would produce a quadrat grid of 4 squares)
    dist : the distance to buffer the quadrat grid in degrees (e.g., 1e-14 degrees is no more than 1 nanometer)
    
    Returns
    -------
    points_within_geometry : GeoDataFrame
    """
    
    # create an empty dataframe to append matches to
    points_within_geometry = pd.DataFrame()
    
    # cut the geometry into chunks for r-tree spatial index intersecting
    multipoly = quadrat_cut_geometry(geometry, quadrat_length=quadrat_length, dist=dist)
    
    # create an r-tree spatial index for the nodes (ie, points)
    start_time = time.time()
    sindex = gdf['geometry'].sindex
    log('Created r-tree spatial index for {:,} points in {:,.2f} seconds'.format(len(gdf), time.time()-start_time))
    
    # loop through each chunk of the geometry to find approximate and then precisely intersecting points
    start_time = time.time()
    for poly in multipoly:
        
        # buffer by the <1 micron dist to account for any space lost in the quadrat cutting, otherwise may miss point(s) that lay directly on quadrat line
        poly = poly.buffer(dist).buffer(0)
        
        # find approximate matches with r-tree, then precise matches from those approximate ones
        possible_matches_index = list(sindex.intersection(poly.bounds))
        possible_matches = gdf.iloc[possible_matches_index]
        precise_matches = possible_matches[possible_matches.intersects(poly)]
        points_within_geometry = points_within_geometry.append(precise_matches)
    
    # drop duplicate points, if buffered poly caused an overlap on point(s) that lay directly on a quadrat line
    points_within_geometry = points_within_geometry.drop_duplicates(subset='node')
    
    log('Identified {:,} nodes inside polygon in {:,.2f} seconds'.format(len(points_within_geometry), time.time()-start_time))
    return points_within_geometry
    
    
def truncate_graph_polygon(G, polygon, retain_all=False, truncate_by_edge=False):
    """
    Remove every node in graph that falls outside some shapely Polygon or MultiPolygon.
    
    Parameters
    ----------
    G : graph
    polygon : Polygon or MultiPolygon, only retain nodes in graph that lie within this geometry
    retain_all : bool, if True, return the entire graph even if it is not connected
    truncate_by_edge : bool, if True retain node if it's outside polygon but at least one of node's neighbors are within polygon (NOT CURRENTLY IMPLEMENTED)
    
    Returns
    -------
    G : graph
    """
    
    start_time = time.time()
    G = G.copy()
    log('Identifying all nodes that lie outside the polygon...')
    
    # get a GeoDataFrame of all the nodes, for spatial analysis
    node_geom = [Point(data['x'], data['y']) for _, data in G.nodes(data=True)]
    gdf_nodes = gpd.GeoDataFrame({'node':G.nodes(), 'geometry':node_geom})
    gdf_nodes.crs = G.graph['crs']
    
    # find all the nodes in the graph that lie outside the polygon
    points_within_geometry = intersect_index_quadrats(gdf_nodes, polygon)
    nodes_outside_polygon = gdf_nodes[~gdf_nodes.index.isin(points_within_geometry.index)]    
    
    # now remove from the graph all those nodes that lie outside the place polygon
    start_time = time.time()
    G.remove_nodes_from(nodes_outside_polygon['node'])
    log('Removed {:,} nodes outside polygon in {:,.2f} seconds'.format(len(nodes_outside_polygon), time.time()-start_time))
    
    # remove any isolated nodes
    
    
    # remove any isolated nodes and retain only the largest component (if retain_all is True)
    if not retain_all:
        G = remove_isolated_nodes(G)
        G = get_largest_component(G)
    
    return G
    
    
def add_edge_lengths(G):
    """
    Add length (meters) attribute to each edge by great circle distance between nodes u and v.
    
    Parameters
    ----------
    G : graph
    
    Returns
    -------
    G : graph
    """
    
    start_time = time.time()
    
    for u, v, key, data in G.edges(keys=True, data=True):
        
        #geopy points are (lat, lon) so that's (y, x) for this great_circle calculation
        u_point = (G.node[u]['y'], G.node[u]['x'])
        v_point = (G.node[v]['y'], G.node[v]['x'])
        edge_length = great_circle(u_point, v_point).m 
        data['length'] = edge_length
    
    log('Added edge lengths to graph in {:,.2f} seconds'.format(time.time()-start_time))
    return G
    

def get_nearest_node(G, point, return_dist=False):
    """
    Return the graph node nearest to some specified point.
    
    Parameters
    ----------
    G : graph
    point : tuple, the (lat, lon) point for which we will find the nearest node in the graph
    return_dist : bool, optionally also return the distance between the point and the nearest node
    
    Returns
    -------
    G : graph
    distance : float, optional
    """
    start_time = time.time()
    nodes = G.nodes(data=True)
    
    # the nearest node is the one that minimizes great circle distance between its coordinates and the passed-in point
    # geopy points are (lat, lon) so that's (y, x)
    nearest_node = min(nodes, key=lambda node: great_circle((node[1]['y'], node[1]['x']), point).m)
    log('Found nearest node ({}) to point {} in {:,.2f} seconds'.format(nearest_node[0], point, time.time()-start_time))
    
    if return_dist:
        # if caller requested return_dist, calculate the great circle distance between the point and the nearest node and return it as well
        distance = great_circle((nearest_node[1]['y'], nearest_node[1]['x']), point).m
        return nearest_node[0], distance
    else:
        return nearest_node[0]

        
def add_path(G, data, one_way):
    """
    Add a path to the graph.
    
    Parameters
    ----------
    G : graph
    data : dict, the attributes of the path
    one_way : bool, if this path is one-way or if it is bi-directional
    
    Returns
    -------
    None
    """
    # extract the ordered list of nodes from this path element, then delete it so we don't add it as an attribute to the edge later
    path_nodes = data['nodes']
    del data['nodes']
    
    # set the oneway attribute to the passed-in value, to make it consistent True/False values
    data['oneway'] = one_way
    
    # zip together the path nodes so you get tuples like (0,1), (1,2), (2,3) and so on
    path_edges = list(zip(path_nodes[:-1], path_nodes[1:]))
    G.add_edges_from(path_edges, attr_dict=data)
    
    # if the path is NOT one-way
    if not one_way:
        # reverse the direction of each edge and add this path going the opposite direction
        path_edges_opposite_direction = [(v, u) for u, v in path_edges]
        G.add_edges_from(path_edges_opposite_direction, attr_dict=data)


def add_paths(G, paths, network_type):
    """
    Add a collection of paths to the graph.
    
    Parameters
    ----------
    G : graph
    paths : dict, the paths from OSM
    network_type : string, {'all', 'walk', 'drive', etc}, what type of network
    
    Returns
    -------
    None
    """
    # the list of values OSM uses in its 'oneway' tag to denote True
    osm_oneway_values = ['yes', 'true', '1', '-1']
    
    for data in paths.values():
        
        # if this path is tagged as one-way and if it is not a walking network, then we'll add the path in one direction only
        if ('oneway' in data and data['oneway'] in osm_oneway_values) and not network_type=='walk':
            if data['oneway'] == '-1':
                # paths with a one-way value of -1 are one-way, but in the reverse direction of the nodes' order, see osm documentation
                data['nodes'] = list(reversed(data['nodes']))
            # add this path (in only one direction) to the graph
            add_path(G, data, one_way=True)
        
        # else, this path is not tagged as one-way or it is a walking network (you can walk both directions on a one-way street)
        else:
            # add this path (in both directions) to the graph and set its 'oneway' attribute to False
            # if this is a walking network, this may very well be a one-way street (as cars/bikes go), but in a walking-only network it is a bi-directional edge
            add_path(G, data, one_way=False)
    
    return G


def create_graph(response_jsons, name='unnamed', retain_all=False, network_type='all_private'):
    """
    Create a networkx graph from OSM data.
    
    Parameters
    ----------
    response_jsons : list, list of dicts of JSON responses from from the Overpass API
    name : string, the name of the graph
    retain_all : bool, if True, return the entire graph even if it is not connected
    network_type : string, what type of network to create
    
    Returns
    -------
    G : graph
    """
    log('Creating networkx graph from downloaded OSM data...')
    start_time = time.time()
    
    # make sure we got data back from the server requests
    elements = []
    for response_json in response_jsons:
        elements.extend(response_json['elements'])
    if len(elements) < 1:
        raise ValueError('There are no data elements in the response JSON objects')
    
    # create the graph as a MultiDiGraph and set the original CRS to lat-long
    G = nx.MultiDiGraph(name=name, crs={'init':'epsg:4326'})
    
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
        G.add_node(node, attr_dict=data)
    
    # add each osm way (aka, path) to the graph
    G = add_paths(G, paths, network_type)
    
    # retain only the largest connected component, if caller did not set retain_all=True
    if not retain_all:
        G = get_largest_component(G)
    
    log('Created graph with {:,} nodes and {:,} edges in {:,.2f} seconds'.format(len(G.nodes()), len(G.edges()), time.time()-start_time))
    
    # add length (great circle distance between nodes) attribute to each edge to use as weight
    G = add_edge_lengths(G)

    return G
    
    
def bbox_from_point(point, distance=1000, project_utm=False, utm_crs=None):
    """
    Create a bounding box some distance in each direction (north, south, east, and west) from some (lat, lon) point.
    
    Parameters
    ----------
    point : tuple, the (lat, lon) point to create the bounding box around
    distance : int, how many meters the north, south, east, and west sides of the box should each be from the point
    project_utm : bool, if True return bbox as UTM coordinates
    utm_crs : dict, if project_utm is True, use this CRS to project the bbox to UTM
    
    Returns
    -------
    north, south, east, west : float, float, float, float
    """
    north = vincenty(meters=distance).destination(point, bearing=0).latitude
    south = vincenty(meters=distance).destination(point, bearing=180).latitude
    east = vincenty(meters=distance).destination(point, bearing=90).longitude
    west = vincenty(meters=distance).destination(point, bearing=270).longitude
    
    if project_utm:
        # create a gdf with one geometry element: the lat-long bounding box
        gdf = gpd.GeoDataFrame()
        gdf['geometry'] = None
        gdf.loc[0, 'geometry'] = LineString([(west, north), (east, north), (east, south), (west, south), (west, north)])
        
        # set the original crs to lat-long and then project it to the specified CRS
        gdf.crs = {'init':'epsg:4326'}
        gdf = gdf.to_crs(utm_crs)
        
        # extract the projected bounding box from the gdf then extract the max and min extents from it
        x, y = gdf.loc[0, 'geometry'].xy
        north, south, east, west = max(y), min(y), max(x), min(x)
        log('Created bounding box {} meters in each direction from {} and projected it: {},{},{},{}'.format(distance, point, north, south, east, west))
    else:
        log('Created bounding box {} meters in each direction from {}: {},{},{},{}'.format(distance, point, north, south, east, west))
        
    return north, south, east, west
    
    
def graph_from_bbox(north, south, east, west, network_type='all_private', simplify=True, retain_all=False, truncate_by_edge=False, 
                    name='unnamed', timeout=180, memory=None, max_query_area_size=0.4, clean_periphery=True):
    """
    Create a networkx graph from OSM data within some bounding box.
    
    Parameters
    ----------
    north : float, northern latitude of bounding box
    south : float, southern latitude of bounding box
    east : float, eastern longitude of bounding box
    west : float, western longitude of bounding box
    network_type : string, what type of street network to get
    simplify : bool, if true, simplify the graph topology
    retain_all : bool, if True, return the entire graph even if it is not connected
    truncate_by_edge : bool, if True retain node if it's outside bbox but at least one of node's neighbors are within bbox
    name : string, the name of the graph
    timeout : int, the timeout interval for requests and to pass to API
    memory : int, server memory allocation size for the query, in bytes. If none, server will use its default allocation size
    max_query_area_size : float, max size for any part of the geometry, in square degrees: any polygon bigger will get divided up for multiple queries to API
    clean_periphery : if True (and simplify=True), buffer 0.5km to get a graph larger than requested, then simplify, then truncate it to requested spatial extent
    
    Returns
    -------
    G : graph
    """
    
    if clean_periphery and simplify:
        # create a new buffered bbox 0.5km around the desired one
        buffer_dist = 500
        polygon = Polygon([(west, north), (west, south), (east, south), (east, north)])
        gdf = gpd.GeoDataFrame(columns=['geometry'], dtype=object)
        gdf.crs = {'init':'epsg:4326'}
        gdf.name = ''
        gdf.loc[0, 'geometry'] = polygon
        gdf_utm = project_gdf(gdf)
        gdf_utm['geometry'] = gdf_utm['geometry'].buffer(buffer_dist)
        gdf = project_gdf(gdf_utm, to_latlong=True)
        polygon_buffered = gdf['geometry'].iloc[0]
        west_buffered, south_buffered, east_buffered, north_buffered = polygon_buffered.bounds
        
        # get the network data from OSM then create the graph
        response_jsons = osm_net_download(north=north_buffered, south=south_buffered, east=east_buffered, west=west_buffered, 
                                          network_type=network_type, timeout=timeout, memory=memory, max_query_area_size=max_query_area_size)
        G_buffered = create_graph(response_jsons, name=name, retain_all=retain_all, network_type=network_type)
        G = truncate_graph_bbox(G_buffered, north, south, east, west, retain_all=True, truncate_by_edge=truncate_by_edge)
        
        # simplify the graph topology 
        G_buffered = simplify_graph(G_buffered)
        
        # truncate graph by desired bbox to return the graph within the bbox caller wants
        G = truncate_graph_bbox(G_buffered, north, south, east, west, retain_all=retain_all, truncate_by_edge=truncate_by_edge)
        
        # retain the degree of each node in this smaller graph, saved as graph attribute
        # this retains the true degree of each node, even if some of its neighbors are outside the bbox
        start_time = time.time()
        k_buffered = G_buffered.to_undirected(reciprocal=False).degree()
        nodes = set(G.nodes())
        k = {key:value for key, value in k_buffered.items() if key in nodes}
        G.graph['degree_undirected_buffered'] = k
        log('Got undirected node degrees for graph (before removing peripheral edges) in {:,.2f} seconds'.format(time.time()-start_time))
    
    else:
        # get the network data from OSM
        response_jsons = osm_net_download(north=north, south=south, east=east, west=west, network_type=network_type, timeout=timeout, memory=memory, max_query_area_size=max_query_area_size)
        
        # create the graph, then truncate to the bounding box
        G = create_graph(response_jsons, name=name, retain_all=retain_all, network_type=network_type)
        G = truncate_graph_bbox(G, north, south, east, west, retain_all=retain_all, truncate_by_edge=truncate_by_edge)
        
        # simplify the graph topology as the last step. don't truncate after simplifying or you may have simplified out to an endpoint
        # beyond the truncation distance, in which case you will then strip out your entire edge
        if simplify:
            G = simplify_graph(G)
    
    log('graph_from_bbox() returning graph with {:,} nodes and {:,} edges'.format(len(G.nodes()), len(G.edges())))
    return  G
    
    
def graph_from_point(center_point, distance=1000, distance_type='bbox', network_type='all_private', simplify=True, retain_all=False, 
                     truncate_by_edge=False, name='unnamed', timeout=180, memory=None, max_query_area_size=0.4, clean_periphery=True):
    """
    Create a networkx graph from OSM data within some distance of some (lat, lon) center point.
    
    Parameters
    ----------
    center_point : tuple, the (lat, lon) central point around which to construct the graph
    distance : int, retain only those nodes within this many meters of the center of the graph
    distance_type : string, {'network', 'bbox'} if 'bbox', retain only those nodes within a bounding box of the distance parameter. if 'network', retain only those nodes within some network distance from the center-most node.
    network_type : string, what type of street network to get
    simplify : bool, if true, simplify the graph topology
    retain_all : bool, if True, return the entire graph even if it is not connected
    truncate_by_edge : bool, if True retain node if it's outside bbox but at least one of node's neighbors are within bbox
    name : string, the name of the graph
    timeout : int, the timeout interval for requests and to pass to API
    memory : int, server memory allocation size for the query, in bytes. If none, server will use its default allocation size
    max_query_area_size : float, max size for any part of the geometry, in square degrees: any polygon bigger will get divided up for multiple queries to API
    clean_periphery : if True (and simplify=True), buffer 0.5km to get a graph larger than requested, then simplify, then truncate it to requested spatial extent
    
    Returns
    -------
    G : graph
    """
    
    # create a bounding box from the center point and the distance in each direction
    north, south, east, west = bbox_from_point(center_point, distance)
    if distance_type == 'bbox':
        # if the network distance_type is bbox, create a graph from the bounding box
        G = graph_from_bbox(north, south, east, west, network_type=network_type, simplify=simplify, retain_all=retain_all, 
                            truncate_by_edge=truncate_by_edge, name=name, timeout=timeout, memory=memory, max_query_area_size=max_query_area_size, clean_periphery=clean_periphery)
    elif distance_type == 'network':
        # if the network distance_type is network, create a graph from the bounding box but do not simplify it yet (only simplify a graph after all truncation is performed! otherwise you get weird artifacts)
        G = graph_from_bbox(north, south, east, west, network_type=network_type, simplify=False, retain_all=retain_all, 
                            truncate_by_edge=truncate_by_edge, name=name, timeout=timeout, memory=memory, max_query_area_size=max_query_area_size, clean_periphery=clean_periphery)
        
        # next find the node in the graph nearest to the center point, and truncate the graph by network distance from this node
        centermost_node = get_nearest_node(G, center_point)
        G = truncate_graph_dist(G, centermost_node, max_distance=distance)
        
        # simplify the graph topology as the last step. don't truncate after simplifying or you may have simplified out to an endpoint
        # beyond the truncation distance, in which case you will then strip out your entire edge
        # that's why simplify=False above, so we didn't do it before the truncate_graph_dist() call 2 lines after it
        if simplify:
            G = simplify_graph(G)
    else:
        raise ValueError('distance_type must be "bbox" or "network"')
    
    log('graph_from_point() returning graph with {:,} nodes and {:,} edges'.format(len(G.nodes()), len(G.edges())))
    return G
        
        
def graph_from_address(address, distance=1000, distance_type='bbox', network_type='all_private', simplify=True, retain_all=False, 
                       truncate_by_edge=False, return_coords=False, name='unnamed', timeout=180, memory=None, max_query_area_size=0.4, clean_periphery=True):
    """
    Create a networkx graph from OSM data within some distance of some address.
    
    Parameters
    ----------
    address : string, the address to geocode and use as the central point around which to construct the graph
    distance : int, retain only those nodes within this many meters of the center of the graph
    distance_type : string, {'network', 'bbox'} if 'bbox', retain only those nodes within a bounding box of the distance parameter. if 'network', retain only those nodes within some network distance from the center-most node.
    network_type : string, what type of street network to get
    simplify : bool, if true, simplify the graph topology
    retain_all : bool, if True, return the entire graph even if it is not connected
    truncate_by_edge : bool, if True retain node if it's outside bbox but at least one of node's neighbors are within bbox
    return_coords : bool, optionally also return the geocoded coordinates of the address
    name : string, the name of the graph
    timeout : int, the timeout interval for requests and to pass to API
    memory : int, server memory allocation size for the query, in bytes. If none, server will use its default allocation size
    max_query_area_size : float, max size for any part of the geometry, in square degrees: any polygon bigger will get divided up for multiple queries to API
    clean_periphery : if True (and simplify=True), buffer 0.5km to get a graph larger than requested, then simplify, then truncate it to requested spatial extent
    
    Returns
    -------
    G : graph
    point : tuple, optional
    """
    
    # geocode the address string to a (lat, lon) point
    geolocation = Nominatim().geocode(query=address)
    try:
        point = (geolocation.latitude, geolocation.longitude)
    except:
        raise Exception('Geocoding failed for address "{}"'.format(address))
    
    # then create a graph from this point
    G = graph_from_point(point, distance, distance_type, network_type=network_type, simplify=simplify, retain_all=retain_all, 
                         truncate_by_edge=truncate_by_edge, name=name, timeout=timeout, memory=memory, max_query_area_size=max_query_area_size, clean_periphery=clean_periphery)
    log('graph_from_address() returning graph with {:,} nodes and {:,} edges'.format(len(G.nodes()), len(G.edges())))
    
    if return_coords:
        return G, point
    else:
        return G


def graph_from_polygon(polygon, network_type='all_private', simplify=True, retain_all=False, 
                       truncate_by_edge=False, name='unnamed', timeout=180, memory=None, max_query_area_size=0.4, clean_periphery=True):
    """
    Create a networkx graph from OSM data within the spatial boundaries of the passed-in shapely polygon.
    
    Parameters
    ----------
    polygon : shapely Polygon or MultiPolygon, the shape to get network data within
    network_type : string, what type of street network to get
    simplify : bool, if true, simplify the graph topology
    retain_all : bool, if True, return the entire graph even if it is not connected
    truncate_by_edge : bool, if True retain node if it's outside bbox but at least one of node's neighbors are within bbox
    name : string, the name of the graph
    timeout : int, the timeout interval for requests and to pass to API
    memory : int, server memory allocation size for the query, in bytes. If none, server will use its default allocation size
    max_query_area_size : float, max size for any part of the geometry, in square degrees: any polygon bigger will get divided up for multiple queries to API
    clean_periphery : if True (and simplify=True), buffer 0.5km to get a graph larger than requested, then simplify, then truncate it to requested spatial extent
    
    Returns
    -------
    G : graph
    """
    
    # verify that the geometry is valid and is a shapely Polygon/MultiPolygon before proceeding
    if not polygon.is_valid:
        raise ValueError('Shape does not have a valid geometry')
    if not isinstance(polygon, (Polygon, MultiPolygon)):
        raise ValueError('Geometry must be a shapely Polygon or MultiPolygon')
    
    if clean_periphery and simplify:
        # create a new buffered polygon 0.5km around the desired one
        buffer_dist = 500
        gdf = gpd.GeoDataFrame(columns=['geometry'], dtype=object)
        gdf.crs = {'init':'epsg:4326'}
        gdf.name = ''
        gdf.loc[0, 'geometry'] = polygon
        gdf_utm = project_gdf(gdf)
        gdf_utm['geometry'] = gdf_utm['geometry'].buffer(buffer_dist)
        gdf = project_gdf(gdf_utm, to_latlong=True)
        polygon_buffered = gdf['geometry'].iloc[0]
        
        # get the network data from OSM,  create the buffered graph, then truncate it to the buffered polygon
        response_jsons = osm_net_download(polygon=polygon_buffered, network_type=network_type, timeout=timeout, memory=memory, max_query_area_size=max_query_area_size)
        G_buffered = create_graph(response_jsons, name=name, retain_all=True, network_type=network_type)
        G_buffered = truncate_graph_polygon(G_buffered, polygon_buffered, retain_all=True, truncate_by_edge=truncate_by_edge)
        
        # simplify the graph topology
        G_buffered = simplify_graph(G_buffered)
        
        # truncate graph by polygon to return the graph within the polygon caller wants
        G = truncate_graph_polygon(G_buffered, polygon, retain_all=retain_all, truncate_by_edge=truncate_by_edge)
        
        # retain the degree of each node in this smaller graph, saved as graph attribute
        # this retains the true degree of each node, even if some of its neighbors are outside the polygon
        start_time = time.time()
        k_buffered = G_buffered.to_undirected(reciprocal=False).degree()
        nodes = set(G.nodes())
        k = {key:value for key, value in k_buffered.items() if key in nodes}
        G.graph['degree_undirected_buffered'] = k
        log('Got undirected node degrees for graph (before removing peripheral edges) in {:,.2f} seconds'.format(time.time()-start_time))
    
    else:
        # download a list of API responses for the polygon/multipolygon
        response_jsons = osm_net_download(polygon=polygon, network_type=network_type, timeout=timeout, memory=memory, max_query_area_size=max_query_area_size)
        
        # create the graph from the downloaded data
        G = create_graph(response_jsons, name=name, retain_all=True, network_type=network_type)
        
        # truncate the graph to the extent of the polygon
        G = truncate_graph_polygon(G, polygon, retain_all=retain_all, truncate_by_edge=truncate_by_edge)
        
        # simplify the graph topology as the last step. don't truncate after simplifying or you may have simplified out to an endpoint
        # beyond the truncation distance, in which case you will then strip out your entire edge
        if simplify:
            G = simplify_graph(G)
    
    log('graph_from_polygon() returning graph with {:,} nodes and {:,} edges'.format(len(G.nodes()), len(G.edges())))
    return G
    
        
def graph_from_place(query, network_type='all_private', simplify=True, retain_all=False, 
                     truncate_by_edge=False, name='unnamed', which_result=1, buffer_dist=None, timeout=180, memory=None, max_query_area_size=0.4, clean_periphery=True):
    """
    Create a networkx graph from OSM data within the spatial boundaries of some geocodable place(s).
    
    Parameters
    ----------
    query : string or dict or list, the place(s) to geocode/download data for
    network_type : string, what type of street network to get
    simplify : bool, if true, simplify the graph topology
    retain_all : bool, if True, return the entire graph even if it is not connected
    truncate_by_edge : bool, if True retain node if it's outside bbox but at least one of node's neighbors are within bbox
    name : string, the name of the graph
    which_result : int, max number of results to return and which to process upon receipt
    buffer_dist : float, distance to buffer around the place geometry, in meters
    timeout : int, the timeout interval for requests and to pass to API
    memory : int, server memory allocation size for the query, in bytes. If none, server will use its default allocation size
    max_query_area_size : float, max size for any part of the geometry, in square degrees: any polygon bigger will get divided up for multiple queries to API
    clean_periphery : if True (and simplify=True), buffer 0.5km to get a graph larger than requested, then simplify, then truncate it to requested spatial extent
    
    Returns
    -------
    G : graph
    """
    
    # create a GeoDataFrame with the spatial boundaries of the place(s)
    if isinstance(query, str) or isinstance(query, dict):
        # if it is a string (place name) or dict (structured place query), then it is a single place
        gdf_place = gdf_from_place(query, which_result=which_result, buffer_dist=buffer_dist)
        name = query
    elif isinstance(query, list):
        # if it is a list, it contains multiple places to get
        gdf_place = gdf_from_places(query, buffer_dist=buffer_dist)
    else:
        raise ValueError('query must be a string or a list of query strings')
    
    # extract the geometry from the GeoDataFrame to use in API query
    polygon = gdf_place['geometry'].unary_union
    log('Constructed place geometry polygon(s) to query API')
    
    # create graph using this polygon(s) geometry
    G = graph_from_polygon(polygon, network_type=network_type, simplify=simplify, retain_all=retain_all, 
                           truncate_by_edge=truncate_by_edge, name=name, timeout=timeout, memory=memory, max_query_area_size=max_query_area_size, clean_periphery=clean_periphery)
    
    log('graph_from_place() returning graph with {:,} nodes and {:,} edges'.format(len(G.nodes()), len(G.edges())))
    return G
    

def project_graph(G):
    """
    Project a graph from lat-long to the UTM zone appropriate for its geographic location.
    
    Parameters
    ----------
    G : graph, the networkx graph to be projected to UTM
    
    Returns
    -------
    G_proj : graph
    """
    
    G_proj = G.copy()
    start_time = time.time()
    
    # create a GeoDataFrame of the nodes, name it, convert osmid to str
    nodes = {node:data for node, data in G_proj.nodes(data=True)}
    gdf_nodes = gpd.GeoDataFrame(nodes).T
    gdf_nodes.crs = G_proj.graph['crs']
    gdf_nodes.name = '{}_nodes'.format(G_proj.name)
    gdf_nodes['osmid'] = gdf_nodes['osmid'].astype(np.int64).astype(str)
    
    # create new lat/lon columns just to save that data for later, and create a geometry column from x/y
    gdf_nodes['lon'] = gdf_nodes['x']
    gdf_nodes['lat'] = gdf_nodes['y']
    gdf_nodes['geometry'] = gdf_nodes.apply(lambda row: Point(row['x'], row['y']), axis=1)
    log('Created a GeoDataFrame from graph in {:,.2f} seconds'.format(time.time()-start_time))
    
    # project the nodes GeoDataFrame to UTM
    gdf_nodes_utm = project_gdf(gdf_nodes)
    
    # extract data for all edges that have geometry attribute
    edges_with_geom = []
    for u, v, key, data in G_proj.edges(keys=True, data=True):
        if 'geometry' in data:
            edges_with_geom.append({'u':u, 'v':v, 'key':key, 'geometry':data['geometry']})
    
    # create an edges GeoDataFrame and project to UTM, if there were any edges with a geometry attribute
    # geom attr only exists if graph has been simplified, otherwise you don't have to project anything for the edges because the nodes still contain all spatial data
    if len(edges_with_geom) > 0:
        gdf_edges = gpd.GeoDataFrame(edges_with_geom)
        gdf_edges.crs = G_proj.graph['crs']
        gdf_edges.name = '{}_edges'.format(G_proj.name)
        gdf_edges_utm = project_gdf(gdf_edges)
    
    # extract projected x and y values from the nodes' geometry column
    start_time = time.time()
    gdf_nodes_utm['x'] = gdf_nodes_utm['geometry'].map(lambda point: point.x)
    gdf_nodes_utm['y'] = gdf_nodes_utm['geometry'].map(lambda point: point.y)
    gdf_nodes_utm = gdf_nodes_utm.drop('geometry', axis=1)
    log('Extracted projected node geometries from GeoDataFrame in {:,.2f} seconds'.format(time.time()-start_time))
    
    # clear the graph to make it a blank slate for the projected data
    start_time = time.time()
    edges = G_proj.edges(keys=True, data=True)
    graph_name = G_proj.graph['name']
    G_proj.clear()
    
    # add the projected nodes and all their attributes to the graph
    G_proj.add_nodes_from(gdf_nodes_utm.index)
    attributes = gdf_nodes_utm.to_dict()
    for name in gdf_nodes_utm.columns:
        nx.set_node_attributes(G_proj, name, attributes[name])
    
    # add the edges and all their attributes (including reconstructed geometry, when it exists) to the graph
    for u, v, key, attributes in edges:
        if 'geometry' in attributes:
            row = gdf_edges_utm[(gdf_edges_utm['u']==u) & (gdf_edges_utm['v']==v) & (gdf_edges_utm['key']==key)]
            attributes['geometry'] = row['geometry'].iloc[0]
        G_proj.add_edge(u, v, key, attributes)
    
    # set the graph's CRS attribute to the new, projected CRS and return the projected graph
    G_proj.graph['crs'] = gdf_nodes_utm.crs
    G_proj.graph['name'] = '{}_UTM'.format(graph_name)
    log('Rebuilt projected graph in {:,.2f} seconds'.format(time.time()-start_time))
    return G_proj


def get_undirected(G):
    """
    Convert a directed graph to an undirected graph that maintains parallel edges in opposite directions if geometries differ.
    
    Parameters
    ----------
    G : graph
    
    Returns
    -------
    G_undir : Graph
    """
    # set from/to nodes and then make undirected
    G = G.copy()
    for u, v, key in G.edges(keys=True):
        G.edge[u][v][key]['from'] = u
        G.edge[u][v][key]['to'] = v
    
    G_undir = G.to_undirected(reciprocal=False)
    
    # if edges in both directions (u,v) and (v,u) exist in the graph, 
    # attributes for the new undirected edge will be a combination of the attributes of the directed edges.
    # if both edges exist in digraph and their edge data is different, 
    # only one edge is created with an arbitrary choice of which edge data to use.
    # you need to manually retain edges in both directions between nodes if their geometries are different
    # this is necessary to save shapefiles for weird intersections like the one at 41.8958697,-87.6794924
    # find all edges (u,v) that have a parallel edge going the opposite direction (v,u) with a different osmid
    for u, v, key, data in G.edges(keys=True, data=True):
        try:
            # look at each edge going the opposite direction (from v to u)
            for key2 in G.edge[v][u]:
                # if this edge has geometry and its osmid is different from its reverse's
                if 'geometry' in data and not data['osmid'] == G.edge[v][u][key2]['osmid']:
                    # turn the geometry of each edge into lists of x's and y's
                    geom1 = [list(coords) for coords in data['geometry'].xy]
                    geom2 = [list(coords) for coords in G_undir[u][v][key]['geometry'].xy]
                    # reverse the first edge's list of x's and y's to look for a match in either order
                    geom1_r = [list(reversed(list(coords))) for coords in data['geometry'].xy]
                    # if the edge's geometry doesn't match its reverse's geometry in either order
                    if not (geom1 == geom2 or geom1_r == geom2):
                        # add it as a new edge to the graph to be saved (with key equal to the current largest key plus one)
                        new_key = max(G.edge[u][v]) + 1
                        G_undir.add_edge(u, v, new_key, attr_dict=data)
        except:
            pass
    
    return G_undir
    
    
def save_graph_shapefile(G, filename='graph', folder=None):
    """
    Save graph nodes and edges as ESRI shapefiles to disk.
    
    Parameters
    ----------
    G : graph
    filename : string, the name of the shapefiles (not including file extensions)
    folder : string, the folder to contain the shapefiles, if None, use default data folder
    
    Returns
    -------
    None
    """
    
    start_time = time.time()
    if folder is None:
        folder = _data_folder
    
    # convert directed graph G to an undirected graph for saving as a shapefile
    G_save = G.copy()
    G_save = get_undirected(G_save)
    
    # create a GeoDataFrame of the nodes and set CRS
    nodes = {node:data for node, data in G_save.nodes(data=True)}
    gdf_nodes = gpd.GeoDataFrame(nodes).T
    gdf_nodes.crs = G_save.graph['crs']
    
    # create a geometry column then drop the x and y columns
    gdf_nodes['geometry'] = gdf_nodes.apply(lambda row: Point(row['x'], row['y']), axis=1)
    gdf_nodes = gdf_nodes.drop(['x', 'y'], axis=1)

    # make everything but geometry column a string
    gdf_nodes['osmid'] = gdf_nodes['osmid'].astype(np.int64)
    for col in [c for c in gdf_nodes.columns if not c == 'geometry']:
        gdf_nodes[col] = gdf_nodes[col].fillna('').astype(str)
        
    # create a list to hold our edges, then loop through each edge in the graph
    edges = []
    for u, v, key, data in G_save.edges(keys=True, data=True):

        # for each edge, add key and all attributes in data dict to the edge_details
        edge_details = {'key':key}
        for attr_key in data:
            edge_details[attr_key] = data[attr_key]

        # if edge doesn't already have a geometry attribute, create one now
        if not 'geometry' in data:
            point_u = Point((G_save.node[u]['x'], G_save.node[u]['y']))
            point_v = Point((G_save.node[v]['x'], G_save.node[v]['y']))
            edge_details['geometry'] = LineString([point_u, point_v])
        
        edges.append(edge_details)
    
    # create a geodataframe from the list of edges and set the CRS
    gdf_edges = gpd.GeoDataFrame(edges)
    gdf_edges.crs = G_save.graph['crs']

    # make everything but geometry column a string
    for col in [c for c in gdf_edges.columns if not c == 'geometry']:
        gdf_edges[col] = gdf_edges[col].fillna('').astype(str)
        
    # if the save folder does not already exist, create it with a filename subfolder
    folder = '{}/{}'.format(folder, filename)
    if not os.path.exists(folder):
        os.makedirs(folder)
    
    # save the nodes and edges as separate ESRI shapefiles
    gdf_nodes.to_file('{}/nodes'.format(folder))
    gdf_edges.to_file('{}/edges'.format(folder))
    log('Saved graph "{}" to disk as shapefiles at "{}" in {:,.2f} seconds'.format(G_save.name, folder, time.time()-start_time))
    
    
def save_graphml(G, filename='graph.graphml', folder=None):
    """
    Save graph as GraphML file to disk.
    
    Parameters
    ----------
    G : graph
    filename : string, the name of the graphml file (including file extension)
    folder : string, the folder to contain the file, if None, use default data folder
    
    Returns
    -------
    None
    """
    
    start_time = time.time()
    if folder is None:
        folder = _data_folder
    
    # create a copy and convert all the node/edge attribute values to string or it won't save
    G_save = G.copy()
    for dict_key in G_save.graph:
        # convert all the graph attribute values to strings
        G_save.graph[dict_key] = make_str(G_save.graph[dict_key])
    for node, data in G_save.nodes(data=True):
        for dict_key in data:
            # convert all the node attribute values to strings
            data[dict_key] = make_str(data[dict_key])
    for u, v, key, data in G_save.edges(keys=True, data=True):
        for dict_key in data:
            # convert all the edge attribute values to strings
            data[dict_key] = make_str(data[dict_key])
                
    if not os.path.exists(folder):
        os.makedirs(folder)
    
    nx.write_graphml(G_save, '{}/{}'.format(folder, filename))
    log('Saved graph "{}" to disk as GraphML at "{}/{}" in {:,.2f} seconds'.format(G_save.name, folder, filename, time.time()-start_time))
    
    
def load_graphml(filename, folder=None):
    """
    Load a GraphML file from disk and convert the node/edge attributes to correct data types.
    
    Parameters
    ----------
    filename : string, the name of the graphml file (including file extension)
    folder : string, the folder containing the file, if None, use default data folder
    
    Returns
    -------
    G : graph
    """
    start_time = time.time()
    
    # read the graph from disk
    if folder is None:
        folder = _data_folder
    path = '{}/{}'.format(folder, filename)
    G = nx.MultiDiGraph(nx.read_graphml(path, node_type=int))
    
    # convert graph crs attribute from saved string to correct dict data type
    G.graph['crs'] = ast.literal_eval(G.graph['crs'])
     
    if 'degree_undirected_buffered' in G.graph:
        G.graph['degree_undirected_buffered'] = ast.literal_eval(G.graph['degree_undirected_buffered'])
        
    # convert numeric node tags from string to numeric data types
    log('Converting node and edge attribute data types')
    for node, data in G.nodes(data=True):
        data['osmid'] = int(data['osmid'])
        data['x'] = float(data['x'])
        data['y'] = float(data['y'])
    
    # convert numeric, bool, and list node tags from string to correct data types
    for u, v, key, data in G.edges(keys=True, data=True):

        # first parse oneway to bool and length to float - they should always have only 1 value each
        data['oneway'] = bool(data['oneway'])
        data['length'] = float(data['length'])

        # these attributes might have a single value, or a list if edge's topology was simplified
        for attr in ['highway', 'name', 'bridge', 'tunnel', 'lanes', 'ref', 'maxspeed', 'service', 'access']:
            # if this edge has this attribute, and it starts with '[' and ends with ']', then it's a list to be parsed
            if attr in data and data[attr][0] == '[' and data[attr][-1] == ']':
                # convert the string list to a list type, else leave as single-value string
                data[attr] = ast.literal_eval(data[attr])

        # osmid might have a single value or a list, but if single value, then parse int
        if 'osmid' in data:
            if data['osmid'][0] == '[' and data['osmid'][-1] == ']':
                data['osmid'] = ast.literal_eval(data['osmid'])
            else:
                data['osmid'] = int(data['osmid'])

        # if geometry attribute exists, load the string as well-known text to shapely LineString
        if 'geometry' in data:
            data['geometry'] = wkt.loads(data['geometry'])
    
    log('Loaded graph with {:,} nodes and {:,} edges from disk in {:,.2f} seconds'.format(len(G.nodes()),
                                                                                                len(G.edges()),
                                                                                                time.time()-start_time))
    return G
    
    
############################################################################
#
# Functions for simplification of graph topology
#
############################################################################    


def is_endpoint(G, node, strict=True):
    """
    Return True if the node is a "real" endpoint of an edge in the network, otherwise False.
    OSM data includes lots of nodes that exist only as points to help streets bend around curves.
    An end point is a node that either:
        1. is its own neighbor, ie, it self-loops
        2. or, has no incoming edges or no outgoing edges, ie, all its incident edges point inward or all its incident edges point outward
        3. or, it does not have exactly two neighbors and degree of 2 or 4
        4. or, if strict mode is false, if its edges have different OSM IDs
    
    Parameters
    ----------
    G : graph
    node : int, the node to examine
    strict : bool, if False, allow nodes to be end points even if they fail all other rules but have edges with different OSM IDs
    
    Returns
    -------
    bool
    """
    neighbors = set(G.predecessors(node) + G.successors(node))
    n = len(neighbors)
    d = G.degree(node)
    
    if node in neighbors:
        # if the node appears in its list of neighbors, it self-loops. this is always an endpoint.
        return True
        
    # if node has no incoming edges or no outgoing edges, it must be an end point
    elif G.out_degree(node)==0 or G.in_degree(node)==0:
        return True
    
    elif not (n==2 and (d==2 or d==4)):
        # else, if it does NOT have 2 neighbors AND either 2 or 4 directed edges, it is an endpoint
        # either it has 1 or 3+ neighbors, in which case it is a dead-end or an intersection of multiple streets
        # or it has 2 neighbors but 3 degree (indicating a change from oneway to twoway)
        # or more than 4 degree (indicating a parallel edge) and thus is an endpoint
        return True
    
    elif not strict:
        # non-strict mode
        osmids = []

        # add all the edge OSM IDs for incoming edges
        for u in G.predecessors(node):
            for key in G.edge[u][node]:
                osmids.append(G.edge[u][node][key]['osmid'])

        # add all the edge OSM IDs for outgoing edges        
        for v in G.successors(node):
            for key in G.edge[node][v]:
                osmids.append(G.edge[node][v][key]['osmid'])

        # if there is more than 1 OSM ID in the list of edge OSM IDs then it is an endpoint, if not, it isn't
        return len(set(osmids)) > 1
    
    else:
        # if none of the preceding rules returned true, then it is not an endpoint
        return False
            

def build_path(G, node, endpoints, path):
    """
    Recursively build a path of nodes until you hit an endpoint node.
    
    Parameters
    ----------
    G : graph
    node : int, the current node to start from
    endpoints : set, the set of all nodes in the graph that are endpoints
    path : list, the list of nodes in order in the path so far
    
    Returns
    -------
    paths_to_simplify : list
    """
    # for each successor in the passed-in node
    for successor in G.successors(node):
        if not successor in path:
            # if this successor is already in the path, ignore it, otherwise add it to the path
            path.append(successor)
            if not successor in endpoints:
                # if this successor is not an endpoint, recursively call build_path until you find an endpoint
                path = build_path(G, successor, endpoints, path)
            else:
                # if this successor is an endpoint, we've completed the path, so return it
                return path
    
    if (not path[-1] in endpoints) and (path[0] in G.successors(path[-1])):
        # if the end of the path is not actually an endpoint and the path's first node is a successor of the 
        # path's final node, then this is actually a self loop, so add path's first node to end of path to close it
        path.append(path[0])
    
    return path
    
    
def get_paths_to_simplify(G, strict=True):
    """
    Create a list of all the paths to be simplified between endpoint nodes.
    The path is ordered from the first endpoint, through the interstitial nodes, to the second endpoint.
    
    Parameters
    ----------
    G : graph
    strict : bool, if False, allow nodes to be end points even if they fail all other rules but have edges with different OSM IDs
    
    Returns
    -------
    paths_to_simplify : list
    """
    
    # first identify all the nodes that are endpoints
    start_time = time.time()
    endpoints = set([node for node in G.nodes() if is_endpoint(G, node, strict=strict)])
    log('Identified {:,} edge endpoints in {:,.2f} seconds'.format(len(endpoints), time.time()-start_time))
    
    start_time = time.time()
    paths_to_simplify = []
    
    # for each endpoint node, look at each of its successor nodes
    for node in endpoints:
        for successor in G.successors(node):
            if not successor in endpoints:
                # if the successor is not an endpoint, build a path from the endpoint node to the next endpoint node
                try:
                    path = build_path(G, successor, endpoints, path=[node, successor])
                    paths_to_simplify.append(path)
                except RuntimeError:
                    # recursion errors occur if some connected component is a self-contained ring in which all nodes are not end points
                    # handle it by just ignoring that component and letting its topology remain intact (this should be a rare occurrence)
                    # RuntimeError is what Python <3.5 will throw, Py3.5+ throws RecursionError but it is a subtype of RuntimeError so it still gets handled
                    log('Recursion error: exceeded max depth, moving on to next endpoint successor', level=lg.WARNING)
    
    log('Constructed all paths to simplify in {:,.2f} seconds'.format(time.time()-start_time))
    return paths_to_simplify
    
    
def simplify_graph(G_, strict=True):
    """
    Simplify a graph by removing all nodes that are not end points. 
    Create an edge directly between the end points that encapsulate them,
    but retain the geometry of the original edges, saved as attribute in new edge
    
    Parameters
    ----------
    G_ : graph
    strict : bool, if False, allow nodes to be end points even if they fail all other rules but have edges with different OSM IDs
    
    Returns
    -------
    G : graph
    """
    
    log('Begin topologically simplifying the graph...')
    G = G_.copy()
    initial_node_count = len(G.nodes())
    initial_edge_count = len(G.edges())
    all_nodes_to_remove = []
    all_edges_to_add = []
    
    # construct a list of all the paths that need to be simplified
    paths = get_paths_to_simplify(G, strict=strict)
        
    start_time = time.time()
    for path in paths:
    
        # add the interstitial edges we're removing to a list so we can retain their spatial geometry
        edge_attributes = {}
        for u, v in zip(path[:-1], path[1:]):

            # there shouldn't be multiple edges between interstitial nodes
            edges = G.edge[u][v]
            if not len(edges) == 1:
                log('Multiple edges between "{}" and "{}" found when simplifying'.format(u, v), level=lg.WARNING)

            # the only element in this list as long as above assertion is True (MultiGraphs use keys (the 0 here), indexed with ints from 0 and up)
            edge = edges[0]
            for key in edge:
                if key in edge_attributes:
                    # if this key already exists in the dict, append it to the value list
                    edge_attributes[key].append(edge[key])
                else:
                    # if this key doesn't already exist, set the value to a list containing the one value
                    edge_attributes[key] = [edge[key]]

        for key in edge_attributes:
            # don't touch the length attribute, we'll sum it at the end
            if len(set(edge_attributes[key])) == 1 and not key == 'length':
                # if there's only 1 unique value in this attribute list, consolidate it to the single value (the zero-th)
                edge_attributes[key] = edge_attributes[key][0]
            elif not key == 'length':
                # otherwise, if there are multiple values, keep one of each value
                edge_attributes[key] = list(set(edge_attributes[key]))

        # construct the geometry and sum the lengths of the segments
        edge_attributes['geometry'] = LineString([Point((G.node[node]['x'], G.node[node]['y'])) for node in path])
        edge_attributes['length'] = sum(edge_attributes['length'])

        # add the nodes and edges to their lists for processing at the end
        all_nodes_to_remove.extend(path[1:-1])
        all_edges_to_add.append({'origin':path[0], 
                                 'destination':path[-1], 
                                 'attr_dict':edge_attributes})
    
    # for each edge to add in the list we assembled, create a new edge between the origin and destination
    for edge in all_edges_to_add:
        G.add_edge(edge['origin'], edge['destination'], attr_dict=edge['attr_dict'])

    # finally remove all the interstitial nodes between the new edges
    G.remove_nodes_from(set(all_nodes_to_remove))
    
    msg = 'Simplified graph (from {:,} to {:,} nodes and from {:,} to {:,} edges) in {:,.2f} seconds'
    log(msg.format(initial_node_count, len(G.nodes()), initial_edge_count, len(G.edges()), time.time()-start_time))
    return G
    
    
############################################################################
#
# Start of plotting functions
#
############################################################################


def get_edge_colors_by_attr(G, attr, num_bins=5, cmap='spectral', start=0.1, stop=0.9):
    """
    Get a list of edge colors by binning some continuous-variable attribute into quantiles
    
    Parameters
    ----------
    G : graph
    attr : string, the name of the continuous-variable attribute
    num_bins : int, how many quantiles
    cmap : string, name of a colormap
    start : float, where to start in the colorspace
    stop : float, where to end in the colorspace
    
    Returns
    -------
    colors : list
    """
    bin_labels = range(num_bins)
    attr_values = pd.Series([data[attr] for u, v, key, data in G.edges(keys=True, data=True)])
    cats = pd.qcut(x=attr_values, q=num_bins, labels=bin_labels)
    color_list = [cm.get_cmap(cmap)(x) for x in np.linspace(start, stop, num_bins)]
    colors = [color_list[cat] for cat in cats]
    return colors
    
    
def save_and_show(fig, ax, save, show, close, filename, file_format, dpi, axis_off):
    """
    Save a figure to disk and show it, as specified.
    
    Parameters
    ----------
    fig : figure
    ax : axis
    save : bool, whether to save the figure to disk or not
    show : bool, whether to display the figure or not
    close : close the figure (only if show equals False) to prevent display
    filename : string, the name of the file to save
    file_format : string, the format of the file to save (e.g., 'jpg', 'png', 'svg')
    dpi : int, the resolution of the image file if saving
    axis_off : bool, if True matplotlib axis was turned off by plot_graph so constrain the saved figure's extent to the interior of the axis
    
    Returns
    -------
    fig, ax : figure, axis
    """
    # save the figure if specified
    if save:
        start_time = time.time()
        
        # create the save folder if it doesn't already exist
        if not os.path.exists(_imgs_folder):
            os.makedirs(_imgs_folder)
        path_filename = '{}/{}.{}'.format(_imgs_folder, filename, file_format)
        
        if file_format == 'svg':
            # if the file_format is svg, prep the fig/ax a bit for saving
            ax.axis('off')
            ax.set_position([0, 0, 1, 1])
            ax.patch.set_alpha(0.)
            fig.patch.set_alpha(0.)
            fig.savefig(path_filename, bbox_inches=0, format=file_format, facecolor=fig.get_facecolor(), transparent=True)
        else:
            if axis_off:
                # if axis is turned off, constrain the saved figure's extent to the interior of the axis
                extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            else:
                extent = 'tight'
            fig.savefig(path_filename, dpi=dpi, bbox_inches=extent, format=file_format, facecolor=fig.get_facecolor(), transparent=True)
        log('Saved the figure to disk in {:,.2f} seconds'.format(time.time()-start_time))
    
    # show the figure if specified
    if show:
        start_time = time.time()
        plt.show()
        log('Showed the plot in {:,.2f} seconds'.format(time.time()-start_time))
    # if show=False, close the figure if close=True to prevent display
    elif close:
        plt.close()
        
    return fig, ax
    
    
def plot_graph(G, bbox=None, fig_height=6, fig_width=None, margin=0.02, axis_off=True, bgcolor='w',
               show=True, save=False, close=True, file_format='jpg', filename='temp', dpi=300, annotate=False,
               node_color='#66ccff', node_size=15, node_alpha=1, node_edgecolor='none', node_zorder=1,
               edge_color='#999999', edge_linewidth=1, edge_alpha=1, use_geom=True):
    """
    Plot a networkx spatial graph.
    
    Parameters
    ----------
    G : graph
    bbox : tuple, bounding box as north,south,east,west - if None will calculate from spatial extents of data
    fig_height : int, matplotlib figure height in inches
    fig_width : int, matplotlib figure width in inches
    margin : float, relative margin around the figure
    axis_off : bool, if True turn off the matplotlib axis
    bgcolor : string, the background color of the figure and axis
    show : bool, if True, show the figure
    save : bool, if True, save the figure as an image file to disk
    close : close the figure (only if show equals False) to prevent display
    file_format : string, the format of the file to save (e.g., 'jpg', 'png', 'svg')
    filename : string, the name of the file if saving
    dpi : int, the resolution of the image file if saving
    annotate : bool, if True, annotate the nodes in the figure
    node_color : string, the color of the nodes
    node_size : int, the size of the nodes
    node_alpha : float, the opacity of the nodes
    node_edgecolor : string, the color of the node's marker's border
    node_zorder : int, zorder to plot nodes, edges are always 2, so make node_zorder 1 to plot nodes beneath them or 3 to plot nodes atop them
    edge_color : string, the color of the edges' lines
    edge_linewidth : float, the width of the edges' lines
    edge_alpha : float, the opacity of the edges' lines
    use_geom : bool, if True, use the spatial geometry attribute of the edges to draw geographically accurate edges, rather than just lines straight from node to node
    
    Returns
    -------
    fig, ax : figure, axis    
    """
    
    log('Begin plotting the graph...')
    node_Xs = [float(node['x']) for node in G.node.values()]
    node_Ys = [float(node['y']) for node in G.node.values()]
    
    # get north, south, east, west values either from bbox parameter or from min/max node coordinate values
    if bbox is None:
        north = max(node_Ys)
        south = min(node_Ys)
        east = max(node_Xs)
        west = min(node_Xs)
    else:
        north, south, east, west = bbox
    
    # if caller did not pass in a fig_width, calculate it proportionately from the fig_height and bounding box aspect ratio
    bbox_aspect_ratio = (north-south)/(east-west)
    if fig_width is None:
        fig_width = fig_height / bbox_aspect_ratio
    
    # create the figure and axis
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), facecolor=bgcolor)
    ax.set_axis_bgcolor(bgcolor)
    
    # draw the edges as lines from node to node
    start_time = time.time()
    lines = []
    for u, v, key, data in G.edges(keys=True, data=True):
        if 'geometry' in data and use_geom:
            # if it has a geometry attribute (a list of line segments), add them to the list of lines to plot
            xs, ys = data['geometry'].xy
            lines.append(list(zip(xs, ys)))
        else:
            # if it doesn't have a geometry attribute, the edge is a straight line from node to node
            x1 = G.node[u]['x']
            y1 = G.node[u]['y']
            x2 = G.node[v]['x']
            y2 = G.node[v]['y']
            line = [(x1, y1), (x2, y2)]
            lines.append(line)
    
    # add the lines to the axis as a linecollection
    lc = LineCollection(lines, colors=edge_color, linewidths=edge_linewidth, alpha=edge_alpha, zorder=2)
    ax.add_collection(lc)
    log('Drew the graph edges in {:,.2f} seconds'.format(time.time()-start_time))
    
    # scatter plot the nodes
    ax.scatter(node_Xs, node_Ys, s=node_size, c=node_color, alpha=node_alpha, edgecolor=node_edgecolor, zorder=node_zorder)
    
    # set the extent of the figure
    margin_ns = (north - south) * margin
    margin_ew = (east - west) * margin
    ax.set_ylim((south - margin_ns, north + margin_ns))
    ax.set_xlim((west - margin_ew, east + margin_ew))
    
    # configure axis appearance
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    if axis_off:
        ax.axis('off')
    
    # annotate the axis with node IDs if annotate=True
    if annotate:
        for node, data in G.nodes(data=True):
            ax.annotate(node, xy=(data['x'], data['y']))
            
    # save and show the figure as specified
    fig, ax = save_and_show(fig, ax, save, show, close, filename, file_format, dpi, axis_off)
    return fig, ax


def plot_graph_route(G, route, bbox=None, fig_height=6, fig_width=None, margin=0.02, bgcolor='w',
                     axis_off=True, show=True, save=False, close=True, file_format='jpg', filename='temp', dpi=300, annotate=False,
                     node_color='#999999', node_size=15, node_alpha=1, node_edgecolor='none', node_zorder=1,
                     edge_color='#999999', edge_linewidth=1, edge_alpha=1, use_geom=True,
                     origin_point=None, destination_point=None,
                     route_color='r', route_linewidth=4, route_alpha=0.5, orig_dest_node_alpha=0.5,
                     orig_dest_node_size=100, orig_dest_node_color='r', orig_dest_point_color='b'):
    """
    Plot a route along a networkx spatial graph.
    
    Parameters
    ----------
    G : graph
    route : list, the route as a list of nodes
    bbox : tuple, bounding box as north,south,east,west - if None will calculate from spatial extents of data
    fig_height : int, matplotlib figure height in inches
    fig_width : int, matplotlib figure width in inches
    margin : float, relative margin around the figure
    axis_off : bool, if True turn off the matplotlib axis
    bgcolor : string, the background color of the figure and axis
    show : bool, if True, show the figure
    save : bool, if True, save the figure as an image file to disk
    close : close the figure (only if show equals False) to prevent display
    file_format : string, the format of the file to save (e.g., 'jpg', 'png', 'svg')
    filename : string, the name of the file if saving
    dpi : int, the resolution of the image file if saving
    annotate : bool, if True, annotate the nodes in the figure
    node_color : string, the color of the nodes
    node_size : int, the size of the nodes
    node_alpha : float, the opacity of the nodes
    node_edgecolor : string, the color of the node's marker's border
    node_zorder : int, zorder to plot nodes, edges are always 2, so make node_zorder 1 to plot nodes beneath them or 3 to plot nodes atop them
    edge_color : string, the color of the edges' lines
    edge_linewidth : float, the width of the edges' lines
    edge_alpha : float, the opacity of the edges' lines
    use_geom : bool, if True, use the spatial geometry attribute of the edges to draw geographically accurate edges, rather than just lines straight from node to node
    origin_point : tuple, optional, an origin (lat, lon) point to plot instead of the origin node
    destination_point : tuple, optional, a destination (lat, lon) point to plot instead of the destination node
    route_color : string, the color of the route
    route_linewidth : int, the width of the route line
    route_alpha : float, the opacity of the route line
    orig_dest_node_alpha : float, the opacity of the origin and destination nodes
    orig_dest_node_size : int, the size of the origin and destination nodes
    orig_dest_node_color : string, the color of the origin and destination nodes
    orig_dest_point_color : string, the color of the origin and destination points if being plotted instead of nodes
    
    Returns
    -------
    fig, ax : figure, axis    
    """
    
    # plot the graph but not the route
    fig, ax = plot_graph(G, bbox=bbox, fig_height=fig_height, fig_width=fig_width, margin=margin, axis_off=axis_off, bgcolor=bgcolor,
                         show=False, save=False, close=False, filename=filename, dpi=dpi, annotate=annotate,
                         node_color=node_color, node_size=node_size, node_alpha=node_alpha, node_edgecolor=node_edgecolor, node_zorder=node_zorder,
                         edge_color=edge_color, edge_linewidth=edge_linewidth, edge_alpha=edge_alpha, use_geom=use_geom)
    
    # the origin and destination nodes are the first and last nodes in the route
    origin_node = route[0]
    destination_node = route[-1]
        
    if origin_point is None or destination_point is None:
        # if caller didn't pass points, use the first and last node in route as origin/destination    
        origin_destination_lats = (G.node[origin_node]['y'], G.node[destination_node]['y'])
        origin_destination_lons = (G.node[origin_node]['x'], G.node[destination_node]['x'])
    else:
        # otherwise, use the passed points as origin/destination
        origin_destination_lats = (origin_point[0], destination_point[0])
        origin_destination_lons = (origin_point[1], destination_point[1])
        orig_dest_node_color = orig_dest_point_color
    
    # scatter the origin and destination points
    ax.scatter(origin_destination_lons, origin_destination_lats, s=orig_dest_node_size, 
               c=orig_dest_node_color, alpha=orig_dest_node_alpha, edgecolor=node_edgecolor, zorder=4)
    
    # plot the route lines
    edge_nodes = list(zip(route[:-1], route[1:]))
    lines = []
    for u, v in edge_nodes:
        # if there are parallel edges, select the shortest in length
        data = min([data for data in G.edge[u][v].values()], key=lambda x: x['length'])
        
        # if it has a geometry attribute (ie, a list of line segments)
        if 'geometry' in data and use_geom:
            # add them to the list of lines to plot
            xs, ys = data['geometry'].xy
            lines.append(list(zip(xs, ys)))
        else:
            # if it doesn't have a geometry attribute, the edge is a straight line from node to node
            x1 = G.node[u]['x']
            y1 = G.node[u]['y']
            x2 = G.node[v]['x']
            y2 = G.node[v]['y']
            line = [(x1, y1), (x2, y2)]
            lines.append(line)
    
    # add the lines to the axis as a linecollection    
    lc = LineCollection(lines, colors=route_color, linewidths=route_linewidth, alpha=route_alpha, zorder=3)
    ax.add_collection(lc)
    
    # save and show the figure as specified
    fig, ax = save_and_show(fig, ax, save, show, close, filename, file_format, dpi, axis_off)
    return fig, ax
    
    