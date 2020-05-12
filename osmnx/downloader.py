import io
import json
import hashlib
import math
import requests
import time
import re
import datetime as dt
import os
import logging as lg
from collections import OrderedDict
from dateutil import parser as date_parser
from .errors import *
from .utils import make_str, log

from . import settings


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
                        'elevator|escalator|proposed|construction|bridleway|abandoned|platform|raceway|service"]'
                        '["motor_vehicle"!~"no"]["motorcar"!~"no"]{}'
                        '["service"!~"parking|parking_aisle|driveway|private|emergency_access"]').format(
        settings.default_access)

    # drive+service: allow ways tagged 'service' but filter out certain types of
    # service ways
    filters['drive_service'] = ('["area"!~"yes"]["highway"!~"cycleway|footway|path|pedestrian|steps|track|corridor|'
                                'elevator|escalator|proposed|construction|bridleway|abandoned|platform|raceway"]'
                                '["motor_vehicle"!~"no"]["motorcar"!~"no"]{}'
                                '["service"!~"parking|parking_aisle|private|emergency_access"]').format(
        settings.default_access)

    # walking: filter out cycle ways, motor ways, private ways, and anything
    # specifying foot=no. allow service roads, permitting things like parking
    # lot lanes, alleys, etc that you *can* walk on even if they're not exactly
    # pleasant walks. some cycleways may allow pedestrians, but this filter ignores
    # such cycleways.
    filters['walk'] = ('["area"!~"yes"]["highway"!~"cycleway|motor|proposed|construction|abandoned|platform|raceway"]'
                       '["foot"!~"no"]["service"!~"private"]{}').format(settings.default_access)

    # biking: filter out foot ways, motor ways, private ways, and anything
    # specifying biking=no
    filters['bike'] = ('["area"!~"yes"]["highway"!~"footway|steps|corridor|elevator|escalator|motor|proposed|'
                       'construction|abandoned|platform|raceway"]'
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
            log('Did not save to cache because response_json is None')
        else:
            # create the folder on the disk if it doesn't already exist
            if not os.path.exists(settings.cache_folder):
                os.makedirs(settings.cache_folder)

            # hash the url to make the filename succinct but unique
            filename = hashlib.md5(url.encode('utf-8')).hexdigest()
            cache_filepath = os.path.join(settings.cache_folder, os.extsep.join([filename, 'json']))

            # dump to json, and save to file
            json_str = make_str(json.dumps(response_json))
            with io.open(cache_filepath, 'w', encoding='utf-8') as cache_file:
                cache_file.write(json_str)

            log('Saved response to cache file "{}"'.format(cache_filepath))



def url_in_cache(url):
    """
    Determine if a URL's response exists in the cache.

    Parameters
    ----------
    url : string
        the url to look for in the cache

    Returns
    -------
    filepath : string
        path to cached response for url if it exists in the cache,
        otherwise None
    """

    # hash the url to generate the cache filename
    filename = hashlib.md5(url.encode('utf-8')).hexdigest()
    filepath = os.path.join(settings.cache_folder, os.extsep.join([filename, 'json']))
    
    # if this file exists in the cache, return its full path
    if os.path.isfile(filepath):
        return filepath



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
        cached response for url if it exists in the cache, otherwise None
    """

    # if the tool is configured to use the cache
    if settings.use_cache:
        
        # return cached response for this url if it exists, otherwise return None
        cache_filepath = url_in_cache(url)
        if cache_filepath is not None:
            with io.open(cache_filepath, encoding='utf-8') as cache_file:
                response_json = json.load(cache_file)
            log('Retrieved response from cache file "{}" for URL "{}"'.format(cache_filepath, url))
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
        response = requests.get(settings.overpass_endpoint.rstrip('/') + '/status', headers=get_http_headers())
        status = response.text.split('\n')[3]
        status_first_token = status.split(' ')[0]
    except Exception:
        # if we cannot reach the status endpoint or parse its output, log an
        # error and return default duration
        log('Unable to query {}/status'.format(settings.overpass_endpoint.rstrip('/')), level=lg.ERROR)
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
    params['dedupe'] = 0  # prevent OSM from deduping results so we get precisely 'limit' # of results
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


def nominatim_request(params, type="search", pause_duration=1, timeout=30, error_pause_duration=180):
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
        raise ValueError(
            "The type of Nominatim request is invalid. Please choose one of {{'search', 'reverse', 'lookup'}}")

    # prepare the Nominatim API URL and see if request already exists in the
    # cache
    url = settings.nominatim_endpoint.rstrip('/') + '/{}'.format(type)
    prepared_url = requests.Request('GET', url, params=params).prepare().url
    cached_response_json = get_from_cache(prepared_url)

    if settings.nominatim_key:
        params['key'] = settings.nominatim_key

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
        log('Downloaded {:,.1f}KB from {} in {:,.2f} seconds'.format(size_kb, domain, time.time() - start_time))

        try:
            response_json = response.json()
            save_to_cache(prepared_url, response_json)
        except Exception:
            # 429 is 'too many requests' and 504 is 'gateway timeout' from server
            # overload - handle these errors by recursively calling
            # nominatim_request until we get a valid response
            if response.status_code in [429, 504]:
                # pause for error_pause_duration seconds before re-trying request
                log(
                    'Server at {} returned status code {} and no JSON data. Re-trying request in {:.2f} seconds.'.format(
                        domain,
                        response.status_code,
                        error_pause_duration),
                    level=lg.WARNING)
                time.sleep(error_pause_duration)
                response_json = nominatim_request(params=params, pause_duration=pause_duration, timeout=timeout)

            # else, this was an unhandled status_code, throw an exception
            else:
                log('Server at {} returned status code {} and no JSON data'.format(domain, response.status_code),
                    level=lg.ERROR)
                raise Exception(
                    'Server returned no JSON data.\n{} {}\n{}'.format(response, response.reason, response.text))

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
    url = settings.overpass_endpoint.rstrip('/') + '/interpreter'
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
        log('Downloaded {:,.1f}KB from {} in {:,.2f} seconds'.format(size_kb, domain, time.time() - start_time))

        try:
            response_json = response.json()
            if 'remark' in response_json:
                log('Server remark: "{}"'.format(response_json['remark'], level=lg.WARNING))
            save_to_cache(prepared_url, response_json)
        except Exception:
            # 429 is 'too many requests' and 504 is 'gateway timeout' from server
            # overload - handle these errors by recursively calling
            # overpass_request until we get a valid response
            if response.status_code in [429, 504]:
                # pause for error_pause_duration seconds before re-trying request
                if error_pause_duration is None:
                    error_pause_duration = get_pause_duration()
                log(
                    'Server at {} returned status code {} and no JSON data. Re-trying request in {:.2f} seconds.'.format(
                        domain,
                        response.status_code,
                        error_pause_duration),
                    level=lg.WARNING)
                time.sleep(error_pause_duration)
                response_json = overpass_request(data=data, pause_duration=pause_duration, timeout=timeout)

            # else, this was an unhandled status_code, throw an exception
            else:
                log('Server at {} returned status code {} and no JSON data'.format(domain, response.status_code),
                    level=lg.ERROR)
                raise Exception(
                    'Server returned no JSON data.\n{} {}\n{}'.format(response, response.reason, response.text))

        return response_json



