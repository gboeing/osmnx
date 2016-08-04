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
# Module: osmnx.py
# Description: Retrieve and construct spatial geometries and street networks from OpenStreetMap
###################################################################################################

from __future__ import unicode_literals
import json, math, os, io, hashlib, re, time, datetime as dt, logging as lg
from collections import OrderedDict
import requests, numpy as np, pandas as pd, geopandas as gpd, networkx as nx, matplotlib.pyplot as plt, matplotlib.cm as cm
from matplotlib import collections as mc
from shapely.geometry import Point, LineString
from geopy.distance import great_circle, vincenty
from geopy.geocoders import Nominatim

# default locations to save data, logs, images, and cache
_data_folder = 'osmnx_data'
_logs_folder = 'osmnx_logs'
_imgs_folder = 'osmnx_images'
_cache_folder = 'osmnx_cache'

_use_cache = False

# write log to file and/or to console
_file_log = False
_print_log = False

# useful osm tags
_useful_tags_node = ['ref', 'highway']
_useful_tags_path = ['bridge', 'tunnel', 'oneway', 'lanes', 'ref', 'name', 'highway', 'maxspeed']


def init(data_folder=_data_folder, logs_folder=_logs_folder, imgs_folder=_imgs_folder, 
         cache_folder=_cache_folder, use_cache=_use_cache,
         file_log=_file_log, print_log=_print_log,
         useful_tags_node=_useful_tags_node, useful_tags_path=_useful_tags_path):
    """
    Initialize osmnx and set the default global vars to desired values.
    
    Parameters
    ---------
    data_folder : string, where to save and load data files
    logs_folder : string, where to write the log files
    imgs_folder : string, where to save figures
    cache_folder : string, where to save the http response cache
    use_cache : bool, if True, use a local cache to save/retrieve http responses instead of calling API repetitively for the same request URL
    file_log : bool, if true, save log output to a log file in logs_folder
    print_log : bool, if true, print log output to the console
    useful_tags_node : list, a list of useful OSM tags to attempt to save from node elements
    useful_tags_path : list, a list of useful OSM tags to attempt to save from path elements
    
    Returns
    -------
    None
    """
    
    global _use_cache, _cache_folder, _data_folder, _imgs_folder, _logs_folder, _print_log, _file_log, _useful_tags_node, _useful_tags_path
    _use_cache = use_cache
    _cache_folder = cache_folder
    _data_folder = data_folder
    _imgs_folder = imgs_folder
    _logs_folder = logs_folder
    _print_log = print_log
    _file_log = file_log
    _useful_tags_node = useful_tags_node
    _useful_tags_path = useful_tags_path
    if _file_log:
        log('Initialized osmnx')


def log(message, level=lg.INFO):
    """
    Write a message to the log file and/or print to the the console.
    
    Parameters
    ----------
    message : string, the content of the message to log
    level : int, one of the logger.level constants
    
    Returns
    -------
    None
    """
    if _file_log:
        logger = get_logger()
        if level == lg.DEBUG:
            logger.debug(message)
        elif level == lg.INFO:
            logger.info(message)
        elif level == lg.WARNING:
            logger.warning(message)
        elif level == lg.ERROR:
            logger.error(message)
            
    if _print_log:
        print(message)


def get_logger(name='osmnx', level=lg.INFO):
    """
    Create a logger or return the current one if already instantiated.
    
    Parameters
    ----------
    name : string, name of the logger
    level : int, one of the logger.level constants
    
    Returns
    -------
    logger : logger.logger
    """
    logger = lg.getLogger(name)
    if not getattr(logger, 'handler_set', None):
        todays_date = dt.datetime.today().strftime('%Y_%m_%d')
        log_filename = '{}/{}_{}.log'.format(_logs_folder, name, todays_date)
        if not os.path.exists(_logs_folder):
            os.makedirs(_logs_folder)
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
    Save a http response json object to the cache. Users should always pass OrderedDicts instead of dicts
    of parameters in, so that the parameters stay in the same order each time, producing the same URL string,
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

            # hash the url (to make filename shorter than long url) 
            filename = hashlib.md5(url.encode('utf-8')).hexdigest()
            cache_path_filename = '{}/{}.json'.format(_cache_folder, filename)
            
            # dump to json, and save to file
            json_str = make_str(json.dumps(response_json))
            with io.open(cache_path_filename, 'w', encoding='utf-8') as cache_file:
                cache_file.write(json_str)
            log('Saved response from URL "{}" to cache file "{}"'.format(url, cache_path_filename))
        

def get_from_cache(url):
    """
    Retrieve a http response json object from the cache.
    
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
        
        
def make_request(url, params=None, pause_duration=1, timeout=30):
    """
    Request a URL and return the JSON response
    
    Parameters
    ----------
    url : string, the url of the request
    params : dict or OrderedDict, key-value pairs of parameters
    pause_duration : int, how long to pause before requests, in seconds
    timeout : int, the timeout interval for requests and to pass to API when possible
    
    Returns
    -------
    response_json : dict
    """
    # you have to pass a timeout to the overpass api for longer queries. add this here,
    # making it match the timeout that requests is using (the value passed in the func arg)
    # we recursively pass params_original to make_request each time the request times out 
    # so we can increase the timeout interval for both requests and reformat the string for 
    # the overpass server parameter
    params_original = params.copy()
    if isinstance(params, dict) and 'data' in params and '[timeout:{timeout}]' in params['data']:
        params['data'] = params['data'].format(timeout=timeout)
    
    # prepare the URL and see if it exists in the cache - if so, just use that response instead of making a new HTTP call
    prepared_url = requests.Request('GET', url, params=params).prepare().url
    cached_response_json = get_from_cache(prepared_url)
    if not cached_response_json is None:
        return cached_response_json
    else:
        # if this URL is not already in the cache, go fetch the response
        log('Pausing {:,.2f} seconds before making API request'.format(pause_duration))
        time.sleep(pause_duration)
        start_time = time.time()
        log('Requesting {} with timeout={}'.format(prepared_url, timeout))
        try:
            response = requests.get(url, params, timeout=timeout)
            size_kb = len(response.content) / 1000.0
            domain = re.findall(r'//(?s)(.*?)/', url)[0]
            log('Downloaded {:,.1f}KB from {} in {:,.2f} seconds'.format(size_kb, domain, time.time()-start_time))
            response_json = response.json()
            save_to_cache(prepared_url, response_json)
        except requests.exceptions.Timeout:
            # if the request timed out, call it again recursively, doubling the timeout interval each time
            log('Request timed out after {:,.2f} seconds. Increasing timeout interval and re-trying.'.format(time.time()-start_time), level=lg.WARNING)
            response = make_request(url=url, params=params_original, pause_duration=pause_duration, timeout=timeout*2)
            response_json = response.json()
        
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
    # define the Nominatim API endpoint and parameters
    url = 'https://nominatim.openstreetmap.org/search'
    params = OrderedDict()
    params['format'] = 'json'
    params['limit'] = limit
    params['polygon_geojson'] = polygon_geojson
    
    # add the structured query dict (if provided) to params, otherwise query with place name string
    if isinstance(query, str):
        params['q'] = query
    elif isinstance(query, dict):
        for key in query:
            params[key] = query[key]
    else:
        raise ValueError('query must be a dict or a string')
    
    # request the URL, return the JSON
    response_json = make_request(url, params)
    return response_json
    

def gdf_from_place(query, gdf_name=None, which_result=1):
    """
    Create a GeoDataFrame from a single place name query.
    
    Parameters
    ----------
    query : string or dict, query string or structured query dict to geocode/download
    gdf_name : string, name attribute metadata for GeoDataFrame (this is used to save shapefile later)
    which_result : int, max number of results to return and which to process upon receipt
    
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
        bbox_south, bbox_north, bbox_west, bbox_east = [float(x) for x in data[which_result - 1]['boundingbox']]
        geometry = data[which_result - 1]['geojson']
        place = data[which_result - 1]['display_name']
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
        
        # create the GeoDataFrame, name it
        gdf = gpd.GeoDataFrame.from_features(features)
        gdf.name = gdf_name
        
        # set the original CRS of the GeoDataFrame to lat-long, and return it
        gdf.crs = {'init':'epsg:4326'}
        log('Created GeoDataFrame with {} row for query "{}"'.format(len(gdf), query))
        return gdf
    else:
        # if there was no data returned (or fewer results than which_result specified)
        log('OSM returned no results (or fewer than which_result) for query "{}"'.format(query), level=lg.WARNING)
        gdf = gpd.GeoDataFrame()
        gdf.name = gdf_name
        return gdf
        

def gdf_from_places(queries, gdf_name='unnamed'):
    """
    Create a GeoDataFrame from a  list of place names to query.
    
    Parameters
    ----------
    query : string or dict, query string or structured query dict to geocode/download
    gdf_name : string, name attribute metadata for GeoDataFrame (this is used to save shapefile later)
    
    Returns
    -------
    gdf : GeoDataFrame
    """
    # create an empty GeoDataFrame then append each result as a new row
    gdf = gpd.GeoDataFrame()
    for query in queries:
        gdf = gdf.append(gdf_from_place(query))
        
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
    return '{}.shp'.format(filename)


def save_shapefile(gdf):
    """
    Save GeoDataFrame as an ESRI shapefile.
    
    Parameters
    ----------
    gdf : GeoDataFrame, the gdf to be saved
    
    Returns
    -------
    None
    """
    filename = make_shp_filename(gdf.name)
    folder_path = '{}/{}'.format(_data_folder, filename)
    file_path_name = '{}/{}'.format(folder_path, filename)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    gdf.to_file(file_path_name)
    log('Saved the GeoDataFrame "{}" as shapefile "{}"'.format(gdf.name, file_path_name))
 

def project_gdf(gdf):
    """
    Project a GeoDataFrame to the UTM zone appropriate for its geometries' centroid. The simple calculation
    in this function works well for most latitudes, but won't work for some far northern locations like
    Svalbard and parts of far northern Norway.
    
    Parameters
    ----------
    gdf : GeoDataFrame, the gdf to be projected to UTM
    
    Returns
    -------
    gdf : GeoDataFrame
    """
    assert len(gdf) > 0, 'You cannot project an empty GeoDataFrame.'
    start_time = time.time()
    
    # if GeoDataFrame is already in UTM, just return it
    if (not gdf.crs is None) and ('proj' in gdf.crs) and (gdf.crs['proj'] == 'utm'):
        return gdf
    
    # calculate the centroid of the union of all the geometries in the GeoDataFrame
    avg_longitude = gdf['geometry'].unary_union.centroid.x
    
    # calculate the UTM zone from this avg longitude and define the UTM CRS to project
    utm_zone = int(math.floor((avg_longitude + 180) / 6.0) + 1)
    utm_crs = {'datum': 'NAD83',
               'ellps': 'GRS80',
               'proj' : 'utm',
               'zone' : utm_zone,
               'units': 'm'}
    
    # project the GeoDataFrame to the UTM CRS
    projected_gdf = gdf.to_crs(utm_crs)
    projected_gdf.name = gdf.name
    log('Projected the GeoDataFrame "{}" to UTM {} in {:,.2f} seconds'.format(projected_gdf.name, utm_zone, time.time()-start_time))
    return projected_gdf


    
##########################################################################################################        
#
# End of functions for getting place boundary geometries.
#
# Below are functions for getting and processing street networks.
#
#########################################################################################################    
   
 
 
def osm_net_download(north, south, east, west, network_type='all', pause_duration=1, timeout=180):
    """
    Download OSM ways and nodes within some bounding box from the Overpass API.
    
    Parameters
    ----------
    north : float, northern latitude of bounding box
    south : float, southern latitude of bounding box
    east : float, eastern longitude of bounding box
    west : float, western longitude of bounding box
    network_type : string, what type of street network to get
    pause_duration : int, time to pause in seconds between API requests
    timeout : int, the timeout interval for requests and to pass to API when possible
    
    Returns
    -------
    response_json : dict
    """
    url = 'http://www.overpass-api.de/api/interpreter'
    
    # define the query to send the API. put timeout in double brackets so it remains unformatted until it gets to the next function, make_request() (see comments below)
    # represent bbox as south,west,north,east (and the '>' makes it recurse so we get ways and way nodes)
    data = '[out:json][timeout:{{timeout}}];( way ["highway"] {filters} ({south},{west},{north},{east}); >;);out;' 
    
    # create a filter to exclude certain kinds of routes based on the requested network_type
    if network_type == 'foot':
        filters = '["highway"!~"motor"]' #for walking/biking accessible routes, exclude ways that are motor-only
    elif network_type == 'motor':
        filters = '["highway"!~"foot|cycle|steps|path|pedestrian|track|service"]' #for car accessible routes, exclude ways that are non-motorized-only
    elif network_type == 'all':
        filters = '' #for all routes, filter out nothing
    else:
        raise ValueError('unknown network_type "{}"'.format(network_type))
    
    # format everything but timeout here. we pass the timeout int along to make_request() where it adds it 
    # to the params OrderedDict so that make_request() can call itself recursively, increasing the timeout interval each
    # time, if the query is really big and causing the request to timeout.
    data = data.format(north=north, south=south, east=east, west=west, filters=filters)
    
    # request the URL, return the JSON
    params = OrderedDict()
    params['data'] = data
    response_json = make_request(url, params, timeout=timeout)
    return response_json
    

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
    
    
def remove_orphan_nodes(G):
    """
    Remove from a graph all the nodes that have no edges.
    
    Parameters
    ----------
    G : graph, the graph from which to remove nodes
    
    Returns
    -------
    G : graph
    """
    degree = G.degree()
    orphaned_nodes = [node for node in degree.keys() if degree[node] < 1]
    G.remove_nodes_from(orphaned_nodes)
    log('Removed {:,} orphaned nodes'.format(len(orphaned_nodes)))
    return G
    
   
def get_largest_subgraph(G, retain_all=False):
    """
    Return the largest connected subgraph from a graph.
    
    Parameters
    ----------
    G : graph
    retain_all : bool, if True, return the entire graph even if it is not connected
    
    Returns
    -------
    G : graph
    """
    # if the graph is not connected and caller did not request retain_all, retain only the largest connected subgraph
    if (not retain_all) and (not nx.is_weakly_connected(G)):
        original_len = len(G.nodes())
        G = max(nx.weakly_connected_component_subgraphs(G), key=len)
        log('Graph was not connected, retained only the largest connected subgraph ({:,} of {:,} total nodes)'.format(len(G.nodes()), original_len))
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
    start_time = time.time()
    distances = nx.shortest_path_length(G, source=source_node, weight=weight)
    distant_nodes = {key:value for key, value in distances.items() if value > max_distance}
    G.remove_nodes_from(distant_nodes.keys())
    log('Truncated graph by weighted network distance in {:,.2f} seconds'.format(time.time()-start_time))
    
    # remove any orphaned nodes, keep only the largest subgraph (if retain_all is True), and return G
    G = remove_orphan_nodes(G)
    G = get_largest_subgraph(G, retain_all)
    return G
    
    
def truncate_graph_bbox(G, north, south, east, west, retain_all=False):
    """
    Remove every node in graph that falls outside a bounding box. Needed because osm seems to return entire ways that also 
    include nodes outside the bbox if the way (that is, a way with a single OSM ID) has a node inside the bbox at some point.
    
    Parameters
    ----------
    G : graph
    north : float, northern latitude of bounding box
    south : float, southern latitude of bounding box
    east : float, eastern longitude of bounding box
    west : float, western longitude of bounding box
    retain_all : bool, if True, return the entire graph even if it is not connected
    
    Returns
    -------
    G : graph
    """
    start_time = time.time()
    nodes_outside_bbox = []
    for node, data in G.nodes(data=True):
        if data['y'] > north or data['y'] < south or data['x'] > east or data['x'] < west:
            nodes_outside_bbox.append(node)
    
    G.remove_nodes_from(nodes_outside_bbox)
    log('Truncated graph by bounding box in {:,.2f} seconds'.format(time.time()-start_time))
    
    # remove any orphaned nodes, keep only the largest subgraph (if retain_all is True), and return G
    G = remove_orphan_nodes(G)
    G = get_largest_subgraph(G, retain_all)

    return G
    

def truncate_graph_polygon(G, polygon, retain_all=False):
    """
    Remove every node in graph that falls outside some shapely Polygon or MultiPolygon.
    
    Parameters
    ----------
    G : graph
    polygon : Polygon or MultiPolygon, only retain nodes in graph that lie within this geometry
    retain_all : bool, if True, return the entire graph even if it is not connected
    
    Returns
    -------
    G : graph
    """
    # find all the nodes in the graph that lie outside the polygon
    start_time = time.time()
    log('Identifying all nodes that lie outside the polygon')
    geometry = [Point(data['x'], data['y']) for _, data in G.nodes(data=True)]
    gdf_nodes = gpd.GeoDataFrame({'node':G.nodes(), 'geometry':geometry})
    gdf_nodes.crs = G.graph['crs']
    nodes_outside_polygon = gdf_nodes[~gdf_nodes.intersects(polygon)]
    log('Found {:,} nodes outside polygon in {:,.2f} seconds'.format(len(nodes_outside_polygon), time.time()-start_time))
    
    # now remove from the graph all those nodes that lie outside the place polygon
    start_time = time.time()
    G.remove_nodes_from(nodes_outside_polygon['node'])
    log('Truncated graph by polygon in {:,.2f} seconds'.format(time.time()-start_time))
    
    # remove any orphaned nodes
    G = remove_orphan_nodes(G)
    
    # keep only the largest subgraph (if retain_all is True), then return G
    G = get_largest_subgraph(G, retain_all)
    
    return G
    
    
def add_edge_lengths(G):
    """
    Add length (meters) attribute to each edge by great circle distance between vertices u and v.
    
    Parameters
    ----------
    G : graph
    
    Returns
    -------
    G : graph
    """
    for u, v, key, data in G.edges(keys=True, data=True):
        #geopy points are (lat, lon) so that's (y, x)
        u_point = (G.node[u]['y'], G.node[u]['x'])
        v_point = (G.node[v]['y'], G.node[v]['x'])
        edge_length = great_circle(u_point, v_point).m 
        data['length'] = edge_length
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
    nearest_node = min(nodes, key=lambda node: great_circle((node[1]['y'], node[1]['x']), point).m)
    log('Found nearest node ({}) to point {} in {:,.2f} seconds'.format(nearest_node[0], point, time.time()-start_time))
    
    if return_dist:
        distance = great_circle((nearest_node[1]['y'], nearest_node[1]['x']), point).m #geopy points are (lat, lon) so that's (y, x)
        return nearest_node[0], distance
    else:
        return nearest_node[0]

        
def add_path(G, data, one_way):
    path_nodes = data['nodes']
    del data['nodes']
    data['oneway'] = one_way
    
    # zip together the path nodes so you get tuples like (0,1), (1,2), (2,3) and so on
    path_edges = list(zip(path_nodes[:-1], path_nodes[1:]))
    G.add_edges_from(path_edges, attr_dict=data)
    if not one_way:
        # if it's not one-way, reverse the direction of each edge and add this path going the opposite direction
        path_edges_opposite_direction = [(v, u) for u, v in path_edges]
        G.add_edges_from(path_edges_opposite_direction, attr_dict=data)

        
def add_paths(G, paths):
    
    osm_oneway_values = ['yes', 'true', '1', '-1']
    
    for path, data in paths.items():
        
        if 'oneway' in data and data['oneway'] in osm_oneway_values:
            # path is one-way
            if data['oneway'] == '-1':
                # this one-way needs to have its nodes' order reversed, see osm wiki documentation
                data['nodes'] = list(reversed(data['nodes']))
            add_path(G, data, one_way=True)
        
        else:
            # path is not one-way
            add_path(G, data, one_way=False)
    
    return G


def create_graph(osm_data, name='unnamed', retain_all=False, directed=False):
    """
    Create a networkx graph from OSM data.
    
    Parameters
    ----------
    osm_data : dict, JSON response from from the Overpass API
    name : string, the name of the graph
    retain_all : bool, if True, return the entire graph even if it is not connected
    directed : bool, if True, create a directed graph
    
    Returns
    -------
    G : graph
    """
    log('Creating networkx graph from downloaded OSM data')
    start_time = time.time()
    
    # create the graph as a MultiDiGraph and set the original CRS to lat-long
    G = nx.MultiDiGraph(name=name, crs={'init':'epsg:4326'})
    
    # extract nodes and paths from the downloaded osm data
    nodes, paths = parse_osm_nodes_paths(osm_data)    
    
    # add each osm node to the graph
    for node, data in nodes.items():
        G.add_node(node, attr_dict=data)
    
    # add each osm way (aka, path) to the graph by unpacking the data dict as keyword args (can't pass attribute dict to this function)
    G = add_paths(G, paths)
    
    # retain only the largest connected subgraph, if caller did not set retain_all=True
    G = get_largest_subgraph(G, retain_all=retain_all)
    
    # add length (great circle distance between vertices) attribute to each edge to use as weight
    G = add_edge_lengths(G)
    
    # change the node labels from osm ids to the standard sequential integers
    #G = nx.convert_node_labels_to_integers(G)
    log('Created graph with {:,} nodes and {:,} edges in {:,.2f} seconds'.format(len(G.nodes()), len(G.edges()), time.time()-start_time))

    return G
    
    
def bbox_from_point(point, distance=1000):
    """
    Create a bounding box some distance in each direction (north, south, east, and west) from some (lat, lon) point.
    
    Parameters
    ----------
    point : tuple, the (lat, lon) point to create the bounding box around
    distance : int, how many meters the north, south, east, and west sides of the box should each be from the point
    
    Returns
    -------
    north, south, east, west : float, float, float, float
    """
    north = vincenty(meters=distance).destination(point, bearing=0).latitude
    south = vincenty(meters=distance).destination(point, bearing=180).latitude
    east = vincenty(meters=distance).destination(point, bearing=90).longitude
    west = vincenty(meters=distance).destination(point, bearing=270).longitude
    log('Created bounding box {} meters in each direction from {}: {},{},{},{}'.format(distance, point, north, south, east, west))
    return north, south, east, west
    
    
def graph_from_bbox(north, south, east, west, network_type='all', simplify=True, retain_all=False, directed=False, name='unnamed'):
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
    directed : bool, if True, create a directed graph
    name : string, the name of the graph
    
    Returns
    -------
    G : graph
    """
    
    # get the network data from OSM, create the graph, then truncate to the bounding box
    osm_data = osm_net_download(north, south, east, west, network_type=network_type)
    G = create_graph(osm_data, name=name, retain_all=retain_all, directed=directed)
    G = truncate_graph_bbox(G, north, south, east, west)

    # simplify the graph topology as the last step. don't truncate after simplifying or you may have simplified out to an endpoint
    # beyond the truncation distance, in which case you will then strip out your entire edge
    if simplify:
        G = simplify_graph(G)
    
    log('graph_from_bbox() returning graph with {:,} nodes and {:,} edges'.format(len(G.nodes()), len(G.edges())))
    return  G
    
    
def graph_from_point(center_point, distance=1000, distance_type='bbox', network_type='all', simplify=True, retain_all=False, directed=False, name='unnamed'):
    """
    Create a networkx graph from OSM data within some distance of some (lat, lon) center point.
    
    Parameters
    ----------
    center_point : tuple, the (lat, lon) central point around which to construct the graph
    distance : int, retain only those nodes within this many meters of the center of the graph
    distance_type : string, if 'bbox', retain only those nodes within a bounding box of the distance parameter
                            if 'network', retain only those nodes within some network distance from the center-most node
    network_type : string, what type of street network to get
    simplify : bool, if true, simplify the graph topology
    retain_all : bool, if True, return the entire graph even if it is not connected
    directed : bool, if True, create a directed graph
    name : string, the name of the graph
    
    Returns
    -------
    G : graph
    """
    # create a bounding box from the center point and the distance in each direction
    north, south, east, west = bbox_from_point(center_point, distance)
    if distance_type == 'bbox':
        # if the network distance_type is bbox, create a graph from the bounding box
        G = graph_from_bbox(north, south, east, west, network_type=network_type, simplify=simplify, retain_all=retain_all, directed=directed, name=name)
    elif distance_type == 'network':
        # if the network distance_type is network, create a graph from the bounding box but do not simplify it yet (only simplify a graph after all truncation is performed! otherwise you get weird artifacts)
        G = graph_from_bbox(north, south, east, west, network_type=network_type, simplify=False, retain_all=retain_all, directed=directed, name=name)
        
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
        
        
def graph_from_address(address, distance=1000, distance_type='bbox', network_type='all', simplify=True, retain_all=False, directed=False, return_coords=False, name='unnamed', geocoder_timeout=30):
    """
    Create a networkx graph from OSM data within some distance of some address.
    
    Parameters
    ----------
    address : string, the address to geocode and use as the central point around which to construct the graph
    distance : int, retain only those nodes within this many meters of the center of the graph
    distance_type : string, if 'bbox', retain only those nodes within a bounding box of the distance parameter
                            if 'network', retain only those nodes within some network distance from the center-most node
    network_type : string, what type of street network to get
    simplify : bool, if true, simplify the graph topology
    retain_all : bool, if True, return the entire graph even if it is not connected
    directed : bool, if True, create a directed graph
    return_coords : bool, optionally also return the geocoded coordinates of the address
    name : string, the name of the graph
    geocoder_timeout : int, how many seconds to wait for server response before the geocoder times-out
    
    Returns
    -------
    G : graph
    point : tuple, optional
    """
    
    # geocode the address string to a (lat, lon) point
    geolocation = Nominatim().geocode(query=address, timeout=geocoder_timeout)
    point = (geolocation.latitude, geolocation.longitude)
    
    # then create a graph from this point
    G = graph_from_point(point, distance, distance_type, network_type=network_type, simplify=simplify, retain_all=retain_all, directed=directed, name=name)
    log('graph_from_address() returning graph with {:,} nodes and {:,} edges'.format(len(G.nodes()), len(G.edges())))
    
    if return_coords:
        return G, point
    else:
        return G
        
        
def graph_from_place(query, network_type='all', simplify=True, retain_all=False, directed=False, name='unnamed', which_result=1):
    """
    Create a networkx graph from OSM data within the spatial boundaries of some geocodable place(s).
    
    Parameters
    ----------
    query : string or list, the places to geocode/download data for
    network_type : string, what type of street network to get
    simplify : bool, if true, simplify the graph topology
    retain_all : bool, if True, return the entire graph even if it is not connected
    directed : bool, if True, create a directed graph
    name : string, the name of the graph
    which_result : int, max number of results to return and which to process upon receipt
    
    Returns
    -------
    G : graph
    """
    # create a GeoDataFrame with the spatial boundaries of the place(s)
    if isinstance(query, str):
        gdf_place = gdf_from_place(query, which_result=1)
        name = query
    elif isinstance(query, list):
        gdf_place = gdf_from_places(query)
    else:
        raise ValueError('query must be a string or a list of query strings')
    
    # get the bounding box containing the place(s) then get the graph within that bounding box
    north = gdf_place['bbox_north'].max()
    south = gdf_place['bbox_south'].min()
    east = gdf_place['bbox_east'].max()
    west = gdf_place['bbox_west'].min()
    
    # retain_all is true and truncate is false here because we'll handle that in truncate_graph_polygon() later
    G = graph_from_bbox(north, south, east, west, network_type=network_type, simplify=False, retain_all=True, directed=directed, name=name)
    
    # truncate the graph to the shape of the place(s) polygon then return it
    polygon = gdf_place['geometry'].unary_union
    G = truncate_graph_polygon(G, polygon, retain_all=retain_all)
    
    # simplify the graph topology as the last step. don't truncate after simplifying or you may have simplified out to an endpoint
    # beyond the truncation distance, in which case you will then strip out your entire edge
    # that's why simplify=False above, so we didn't do it before the truncate_graph_polygon() call 2 lines after it
    if simplify:
        G = simplify_graph(G)
        
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
    gdf_nodes['osmid'] = gdf_nodes['osmid'].astype(str)
    
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

    
def save_graph(G, filename='graph'):
    """
    Save graph as file to disk.
    
    Parameters
    ----------
    G : graph
    filename : string
    
    Returns
    -------
    None
    """
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
                
    if not os.path.exists(_data_folder):
        os.makedirs(_data_folder)
    
    nx.write_graphml(G_save, '{}/{}.graphml'.format(_data_folder, filename))
    
    
############################################################################
#
# Functions for simplification of graph topology
#
############################################################################    


def is_endpoint(G, node, strict=True):
    """
    Return True if the node is a "real" endpoint of an edge in the network, otherwise False.
    OSM data includes lots of nodes that exist only as points to help streets bend around curves.
    An endpoint is a node that either has degree 1 (so it's a dead end), or has degree > 2
    (so it's the intersection of multiple streets), or has degree = 2 but its two adjacent
    edges have different OSM IDs, so it is the intersection of different paths.
    
    Parameters
    ----------
    G : graph
    node : int, the node to examine
    strict : bool, if True, only consider intersections or dead-ends (ie, nodes with degree not = 2)
                   if False, also consider nodes with degree = 2 if its two adjacent edges have different OSM IDs
    
    Returns
    -------
    bool
    """
    neighbors = list(set(G.predecessors(node) + G.successors(node)))
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
        # or it has 2 neighbors but 3 degrees (indicating a change from oneway to twoway) or more than 4
        # degrees (indicating a parallel edge) and thus is an endpoint
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

        # if there is more than 1 OSM ID in the list of IDs then it is an endpoint, if not, it isn't
        return len(set(osmids)) > 1
    
    else:
        # if none of the preceding rules returned true, then it is not an endpoint
        return False
            

def find_end_points(G, node, end_points=[], nodes_to_remove=[], previous_node=-1, strict=True):
    """
    Given a node in a graph, if it has degree=2, find the first node in both directions with degree not 2.
    This function works recursively to identify all the nodes of degree 2 between two nodes of degreee not 2.
    
    Parameters
    ----------
    G : graph
    node : int, the node to start with
    end_points : list, a list of the first nodes found in both directions with degree not 2016
    nodes_to_remove : list, a list of all nodes found with degree=2 before finding an endpoint
    previous_node : int, the last node processed
    
    Returns
    -------
    end_points, nodes_to_remove : list, list
    """
    if not is_endpoint(G, node, strict):
        # if this is not an endpoint/intersection, we will remove this node
        nodes_to_remove.append(node)
    
    neighbors = list(set(G.predecessors(node) + G.successors(node)))
    for neighbor in neighbors:
        # look at each neighbor of this node
        if not neighbor == previous_node:
            # if this neighbor is the previous node we just looked at, just ignore it
            if is_endpoint(G, neighbor, strict):
                # the neighbor is an endpoint/intersection, at to the list
                end_points.append(neighbor)
            else:
                # otherwise, if this neighbor is not an endpoint or intersection, recursively call this function
                find_end_points(G, neighbor, end_points=end_points, nodes_to_remove=nodes_to_remove, previous_node=node, strict=strict)
    
    return end_points, nodes_to_remove
    
    
def build_path(G, origin, destination, nodes_to_remove, current_node, path=[], previous_node=-1):

    # get all node's successors that appear in nodes_to_remove until you hit the other end points
    for successor in G.successors(current_node):
        
        # if this successor is the destination 
        # and we are not on the first iteration 
        # and this successor is not the previous node we just visited
        if (successor==destination) and (not previous_node==-1) and (not successor==previous_node):
            # then this is the final destination, add it to the path and return it
            path.append(successor)
            return path
        
        # if this neighbor is not the node we processed in the previous iteration and it is in nodes_to_remove
        if (not successor==previous_node) and (successor in nodes_to_remove):
            # add it to the path, then call this function recursively to find the next step in the path
            path.append(successor)
            path = build_path(G, origin, destination, nodes_to_remove, successor, path, current_node)
            return path
        
        
def build_paths(G, end_points, nodes_to_remove):
    
    # try starting at the first endpoint as the origin
    origin, destination = end_points
    
    # if none of its successors are in nodes_to_remove, this is a one-way at you've started at the end 
    if set(G.successors(origin)).isdisjoint(set(nodes_to_remove)):
        # so, swap origin/destination
        origin, destination = destination, origin

    # add path as a single element to the list of paths (we might also add its reverse momentarily)
    path = build_path(G, origin, destination, nodes_to_remove, current_node=origin, path=[origin])
    paths = [path]
    
    # assert that all edges in path are either all one-way or all two-way
    oneway_values = []
    for u, v in list(zip(path[:-1], path[1:])):
        edges = G.edge[u][v]
        assert len(edges) == 1
        oneway_values.append(G.edge[u][v][0]['oneway'])
    assert len(set(oneway_values))==1
    
    # if not one-way, reverse the path and add as a second element to list paths
    if not oneway_values[0]:
        paths.append(list(reversed(path)))
        
    return paths


def simplify_graph(G_, strict=True):
    """
    Simplify a graph by removing all nodes where degree=2. 
    Create an edge directly between the nodes outside them where degree != 2, 
    but retain the geometry of the original edges, saved as attr in new edge
    
    Parameters
    ----------
    G_ : graph
    
    Returns
    -------
    G : graph
    """
    start_time = time.time()
    log('Begin topologically simplifying the graph')
    
    # make a copy to not alter the original
    G = G_.copy()
    initial_node_count = len(G.nodes())
    initial_edge_count = len(G.edges())
    
    all_nodes_to_remove = []
    all_edges_to_add = []
    
    for node in G.nodes():
    # look at each node in the graph
        
        if node not in all_nodes_to_remove:
        # if this node is in this list, ignore it: it's not an endpoint and we already marked it for removal
            
            if not is_endpoint(G, node, strict):
            # this node is not an intersection or an endpoint
            # therefore it is not a network node, just a point to help a street bend around a curve
                
                try:
                    # find the 'end points' on either side of this node
                    # nodes_to_remove is all the other nodes (that are neither endpoints nor intersections) between the end_points
                    end_points, nodes_to_remove = find_end_points(G, node, end_points=[], nodes_to_remove=[], strict=strict)
                    
                    # build a path of the nodes in sequence between the origin and destination,
                    # making sure each node in the path is in nodes_to_remove so we don't create
                    # a path from some alternate route between these two end point nodes
                    paths = build_paths(G, end_points, nodes_to_remove)
                    
                    # there will be one element in paths if it was one-way, or two if it was two-way
                    for path in paths:
                        
                        # add the interstitial edges we're removing to a list so we can draw with spatial accuracy later
                        edge_attributes = {}
                        for u, v in list(zip(path[:-1], path[1:])):
                            
                            # there should never be multiple edges between interstitial nodes
                            edges = G.edge[u][v]
                            assert len(edges) == 1 
                            
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
                        points = [Point((G.node[node]['x'], G.node[node]['y'])) for node in path]
                        edge_attributes['geometry'] = LineString(points)
                        edge_attributes['length'] = sum(edge_attributes['length'])

                        # add the nodes and edges to their lists for processing at the end
                        all_nodes_to_remove.extend(nodes_to_remove)
                        all_edges_to_add.append({'origin':path[0], 
                                                 'destination':path[-1], 
                                                 'attr_dict':edge_attributes})
                    
                except RuntimeError:
                    # recursion errors occur if some subgraph is a self-contained ring in which all nodes have degree 2
                    # handle it by just ignoring that subgraph and letting its topology remain intact (this should be a rare occurrence)
                    # RuntimeError is what Python <3.5 will throw, Py3.5+ throws RecursionError but it is a subtype of RuntimeError so it still gets handled
                    log('Recursion error: encountered subgraph where all nodes have degree=2', level=lg.WARNING)
    
    # for each edge to add in the list we assembled, create a new edge between the origin and destination
    for edge in all_edges_to_add:
        G.add_edge(edge['origin'], edge['destination'], attr_dict=edge['attr_dict'])

    # finally remove all the interstitial nodes between the new edges
    G.remove_nodes_from(list(set(all_nodes_to_remove)))
    
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
    
    
def plot_graph(G, bbox=None, fig_height=6, fig_width=None, margin=0.02, axis_off=True,
               show=True, save=False, filename='temp.jpg', dpi=300, annotate=False,
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
    show : bool, if True, show the figure
    save : bool, if True, save the figure as an image file to disk
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
    
    log('Begin plotting the graph')
    node_Xs = [float(node['x']) for node in G.node.values()]
    node_Ys = [float(node['y']) for node in G.node.values()]
    
    if bbox is None:
        north = max(node_Ys)
        south = min(node_Ys)
        east = max(node_Xs)
        west = min(node_Xs)
    else:
        north, south, east, west = bbox
    bbox_aspect_ratio = (north-south)/(east-west)
    
    fig, ax = plt.subplots(figsize=(fig_height / bbox_aspect_ratio, fig_height))
    
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
            
    lc = mc.LineCollection(lines, colors=edge_color, linewidths=edge_linewidth, alpha=edge_alpha, zorder=2)
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
    
    if annotate:
        for node, data in G.nodes(data=True):
            ax.annotate(node, xy=(data['x'], data['y']))
            
    # save the figure if specified
    if save:
        start_time = time.time()
        if not os.path.exists(_imgs_folder):
            os.makedirs(_imgs_folder)
        path_filename = '{}/{}'.format(_imgs_folder, filename)
        fig.savefig(path_filename, dpi=dpi, bbox_inches='tight')
        log('Saved the figure to disk in {:,.2f} seconds'.format(time.time()-start_time))
    
    if show:
        start_time = time.time()
        plt.show()
        log('Showed the plot in {:,.2f} seconds'.format(time.time()-start_time))
    
    
    
    return fig, ax


def plot_graph_route(G, route, bbox=None, fig_height=6, fig_width=None, margin=0.02,
                     axis_off=True, show=True, save=False, filename='temp.jpg', dpi=300, annotate=False,
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
    show : bool, if True, show the figure
    save : bool, if True, save the figure as an image file to disk
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
    fig, ax = plot_graph(G, bbox=bbox, fig_height=fig_height, fig_width=fig_width, margin=margin, axis_off=axis_off,
                         show=False, save=save, filename=filename, dpi=dpi, annotate=annotate,
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
                
    lc = mc.LineCollection(lines, colors=route_color, linewidths=route_linewidth, alpha=route_alpha, zorder=3)
    ax.add_collection(lc)
    
    if show:
        start_time = time.time()
        plt.show()
        log('Showed the plot in {:,.2f} seconds'.format(time.time()-start_time))
        
    return fig, ax
    
    