################################################################################
# Module: utils.py
# Description: Utility functions for configuration, logging, geocoding,
#              geospatial analysis, etc.
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

import os
import sys
import time
import unicodedata
import math
import warnings
import logging as lg
import datetime as dt
import networkx as nx
import numpy as np
import pandas as pd
import bz2
import xml.sax
import io
import requests
from itertools import chain
from collections import Counter

from . import settings


def config(data_folder=settings.data_folder,
           logs_folder=settings.logs_folder,
           imgs_folder=settings.imgs_folder,
           cache_folder=settings.cache_folder,
           use_cache=settings.use_cache,
           log_file=settings.log_file,
           log_console=settings.log_console,
           log_level=settings.log_level,
           log_name=settings.log_name,
           log_filename=settings.log_filename,
           useful_tags_node=settings.useful_tags_node,
           useful_tags_path=settings.useful_tags_path):
    """
    Configure osmnx by setting the default global vars to desired values.

    Parameters
    ---------
    data_folder : string
        where to save and load data files
    logs_folder : string
        where to write the log files
    imgs_folder : string
        where to save figures
    cache_folder : string
        where to save the http response cache
    use_cache : bool
        if True, use a local cache to save/retrieve http responses instead of
        calling API repetitively for the same request URL
    log_file : bool
        if true, save log output to a log file in logs_folder
    log_console : bool
        if true, print log output to the console
    log_level : int
        one of the logger.level constants
    log_name : string
        name of the logger
    useful_tags_node : list
        a list of useful OSM tags to attempt to save from node elements
    useful_tags_path : list
        a list of useful OSM tags to attempt to save from path elements

    Returns
    -------
    None
    """

    # set each global variable to the passed-in parameter value
    settings.use_cache = use_cache
    settings.cache_folder = cache_folder
    settings.data_folder = data_folder
    settings.imgs_folder = imgs_folder
    settings.logs_folder = logs_folder
    settings.log_console = log_console
    settings.log_file = log_file
    settings.log_level = log_level
    settings.log_name = log_name
    settings.log_filename = log_filename
    settings.useful_tags_node = useful_tags_node
    settings.useful_tags_path = useful_tags_path

    # if logging is turned on, log that we are configured
    if settings.log_file or settings.log_console:
        log('Configured osmnx')


def log(message, level=None, name=None, filename=None):
    """
    Write a message to the log file and/or print to the the console.

    Parameters
    ----------
    message : string
        the content of the message to log
    level : int
        one of the logger.level constants
    name : string
        name of the logger
    filename : string
        name of the log file

    Returns
    -------
    None
    """

    if level is None:
        level = settings.log_level
    if name is None:
        name = settings.log_name
    if filename is None:
        filename = settings.log_filename

    # if logging to file is turned on
    if settings.log_file:
        # get the current logger (or create a new one, if none), then log
        # message at requested level
        logger = get_logger(level=level, name=name, filename=filename)
        if level == lg.DEBUG:
            logger.debug(message)
        elif level == lg.INFO:
            logger.info(message)
        elif level == lg.WARNING:
            logger.warning(message)
        elif level == lg.ERROR:
            logger.error(message)

    # if logging to console is turned on, convert message to ascii and print to
    # the console
    if settings.log_console:
        # capture current stdout, then switch it to the console, print the
        # message, then switch back to what had been the stdout. this prevents
        # logging to notebook - instead, it goes to console
        standard_out = sys.stdout
        sys.stdout = sys.__stdout__

        # convert message to ascii for console display so it doesn't break
        # windows terminals
        message = unicodedata.normalize('NFKD', make_str(message)).encode('ascii', errors='replace').decode()
        print(message)
        sys.stdout = standard_out


def get_logger(level=None, name=None, filename=None):
    """
    Create a logger or return the current one if already instantiated.

    Parameters
    ----------
    level : int
        one of the logger.level constants
    name : string
        name of the logger
    filename : string
        name of the log file

    Returns
    -------
    logger.logger
    """

    if level is None:
        level = settings.log_level
    if name is None:
        name = settings.log_name
    if filename is None:
        filename = settings.log_filename

    logger = lg.getLogger(name)

    # if a logger with this name is not already set up
    if not getattr(logger, 'handler_set', None):

        # get today's date and construct a log filename
        todays_date = dt.datetime.today().strftime('%Y_%m_%d')
        log_filename = '{}/{}_{}.log'.format(settings.logs_folder, filename, todays_date)

        # if the logs folder does not already exist, create it
        if not os.path.exists(settings.logs_folder):
            os.makedirs(settings.logs_folder)

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
    value : any
        the value to convert to unicode/string

    Returns
    -------
    unicode or string
    """
    try:
        # for python 2.x compatibility, use unicode
        return unicode(value)
    except NameError:
        # python 3.x has no unicode type, so if error, use str type
        return str(value)


def get_largest_component(G, strongly=False):
    """
    Return the largest weakly or strongly connected component from a directed
    graph.

    Parameters
    ----------
    G : networkx multidigraph
    strongly : bool
        if True, return the largest strongly instead of weakly connected
        component

    Returns
    -------
    networkx multidigraph
    """

    start_time = time.time()
    original_len = len(list(G.nodes()))

    if strongly:
        # if the graph is not connected and caller did not request retain_all,
        # retain only the largest strongly connected component
        if not nx.is_strongly_connected(G):
            G = max(nx.strongly_connected_component_subgraphs(G), key=len)
            msg = ('Graph was not connected, retained only the largest strongly '
                   'connected component ({:,} of {:,} total nodes) in {:.2f} seconds')
            log(msg.format(len(list(G.nodes())), original_len, time.time()-start_time))
    else:
        # if the graph is not connected and caller did not request retain_all,
        # retain only the largest weakly connected component
        if not nx.is_weakly_connected(G):
            G = max(nx.weakly_connected_component_subgraphs(G), key=len)
            msg = ('Graph was not connected, retained only the largest weakly '
                   'connected component ({:,} of {:,} total nodes) in {:.2f} seconds')
            log(msg.format(len(list(G.nodes())), original_len, time.time()-start_time))

    return G


def great_circle_vec(lat1, lng1, lat2, lng2, earth_radius=6371009):
    """
    Vectorized function to calculate the great-circle distance between two
    points or between vectors of points.

    Parameters
    ----------
    lat1 : float or array of float
    lng1 : float or array of float
    lat2 : float or array of float
    lng2 : float or array of float
    earth_radius : numeric
        radius of earth in units in which distance will be returned (default is
        meters)

    Returns
    -------
    distance : float or array of float
        distance or vector of distances from (lat1, lng1) to (lat2, lng2) in
        units of earth_radius
    """

    phi1 = np.deg2rad(90 - lat1)
    phi2 = np.deg2rad(90 - lat2)

    theta1 = np.deg2rad(lng1)
    theta2 = np.deg2rad(lng2)

    cos = (np.sin(phi1) * np.sin(phi2) * np.cos(theta1 - theta2) + np.cos(phi1) * np.cos(phi2))

    # ignore warnings during this calculation because numpy warns it cannot
    # calculate arccos for self-loops since u==v
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        arc = np.arccos(cos)

    # return distance in units of earth_radius
    distance = arc * earth_radius
    return distance


def euclidean_dist_vec(y1, x1, y2, x2):
    """
    Vectorized function to calculate the euclidean distance between two points
    or between vectors of points.

    Parameters
    ----------
    y1 : float or array of float
    x1 : float or array of float
    y2 : float or array of float
    x2 : float or array of float

    Returns
    -------
    distance : float or array of float
        distance or vector of distances from (x1, y1) to (x2, y2) in graph units
    """

    # euclid's formula
    distance = ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5
    return distance


def get_nearest_node(G, point, method='greatcircle', return_dist=False):
    """
    Return the graph node nearest to some specified (lat, lng) or (y, x) point,
    and optionally the distance between the node and the point. This function
    can use either a great circle or euclidean distance calculator.

    Parameters
    ----------
    G : networkx multidigraph
    point : tuple
        The (lat, lng) or (y, x) point for which we will find the nearest node
        in the graph
    method : str {'greatcircle', 'euclidean'}
        Which method to use for calculating distances to find nearest node.
        If 'greatcircle', graph nodes' coordinates must be in units of decimal
        degrees. If 'euclidean', graph nodes' coordinates must be projected.
    return_dist : bool
        Optionally also return the distance (in meters if great circle, or graph
        node coordinate units if euclidean) between the point and the nearest
        node.

    Returns
    -------
    int or tuple of (int, float)
        Nearest node ID or optionally a tuple of (node ID, dist), where dist is
        the distance (in meters if great circle, or graph node coordinate units
        if euclidean) between the point and nearest node
    """
    start_time = time.time()
	
    if not G or (G.number_of_nodes() == 0):
        raise ValueError('G argument must be not be empty or should contain at least one node')

    # dump graph node coordinates into a pandas dataframe indexed by node id
    # with x and y columns
    coords = np.array([[node, data['x'], data['y']] for node, data in G.nodes(data=True)])
    df = pd.DataFrame(coords, columns=['node', 'x', 'y']).set_index('node')

    # add columns to the dataframe representing the (constant) coordinates of
    # the reference point
    df['reference_y'] = point[0]
    df['reference_x'] = point[1]

    # calculate the distance between each node and the reference point
    if method == 'greatcircle':
        # calculate distance vector using great circle distances (ie, for
        # spherical lat-long geometries)
        distances = great_circle_vec(lat1=df['reference_y'],
                                     lng1=df['reference_x'],
                                     lat2=df['y'],
                                     lng2=df['x'])

    elif method == 'euclidean':
        # calculate distance vector using euclidean distances (ie, for projected
        # planar geometries)
        distances = euclidean_dist_vec(y1=df['reference_y'],
                                       x1=df['reference_x'],
                                       y2=df['y'],
                                       x2=df['x'])

    else:
        raise ValueError('method argument must be either "greatcircle" or "euclidean"')

    # nearest node's ID is the index label of the minimum distance
    nearest_node = int(distances.idxmin())
    log('Found nearest node ({}) to point {} in {:,.2f} seconds'.format(nearest_node, point, time.time()-start_time))

    # if caller requested return_dist, return distance between the point and the
    # nearest node as well
    if return_dist:
        return nearest_node, distances.loc[nearest_node]
    else:
        return nearest_node


def get_bearing(origin_point, destination_point):
    """
    Calculate the bearing between two lat-long points. Each tuple should
    represent (lat, lng) as decimal degrees.

    Parameters
    ----------
    origin_point : tuple
    destination_point : tuple

    Returns
    -------
    bearing : float
        the compass bearing in decimal degrees from the origin point
        to the destination point
    """

    if not (isinstance(origin_point, tuple) and isinstance(destination_point, tuple)):
        raise TypeError('origin_point and destination_point must be (lat, lng) tuples')

    # get latitudes and the difference in longitude, as radians
    lat1 = math.radians(origin_point[0])
    lat2 = math.radians(destination_point[0])
    diff_lng = math.radians(destination_point[1] - origin_point[1])

    # calculate initial bearing from -180 degrees to +180 degrees
    x = math.sin(diff_lng) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - (math.sin(lat1) * math.cos(lat2) * math.cos(diff_lng))
    initial_bearing = math.atan2(x, y)

    # normalize initial bearing to 0 degrees to 360 degrees to get compass bearing
    initial_bearing = math.degrees(initial_bearing)
    bearing = (initial_bearing + 360) % 360

    return bearing


def add_edge_bearings(G):
    """
    Calculate the compass bearing from origin node to destination node for each
    edge in the directed graph then add each bearing as a new edge attribute.

    Parameters
    ----------
    G : networkx multidigraph

    Returns
    -------
    G : networkx multidigraph
    """

    for u, v, data in G.edges(keys=False, data=True):
        origin_point = (G.nodes[u]['y'], G.nodes[u]['x'])
        destination_point = (G.nodes[v]['y'], G.nodes[v]['x'])
        bearing = get_bearing(origin_point, destination_point)
        data['bearing'] = bearing

    return G


def geocode(query):
    """
    Geocode a query string to (lat, lon) with the Nominatim geocoder.

    Parameters
    ----------
    query : string
        the query string to geocode

    Returns
    -------
    point : tuple
        the (lat, lon) coordinates returned by the geocoder
    """

    # send the query to the nominatim geocoder and parse the json response
    url_template = 'https://nominatim.openstreetmap.org/search?format=json&limit=1&q={}'
    url = url_template.format(query)
    response = requests.get(url, timeout=60)
    results = response.json()

    # if results were returned, parse lat and long out of the result
    if len(results) > 0 and 'lat' in results[0] and 'lon' in results[0]:
        lat = float(results[0]['lat'])
        lon = float(results[0]['lon'])
        point = (lat, lon)
        log('Geocoded "{}" to {}'.format(query, point))
        return point
    else:
        raise Exception('Nominatim geocoder returned no results for query "{}"'.format(query))


def get_route_edge_attributes(G, route, attribute, minimize_key='length'):
    """
    Get a list of attribute values for each edge in a path.

    Parameters
    ----------
    G : networkx multidigraph
    route : list
        list of nodes in the path
    attribute : string
        the name of the attribute to get the value of for each edge
    minimize_key : string
        if there are parallel edges between two nodes, select the one with the
        lowest value of minimize_key

    Returns
    -------
    attribute_values : list
        list of edge attribute values
    """

    attribute_values = []
    for u, v in zip(route[:-1], route[1:]):
        # if there are parallel edges between two nodes, select the one with the
        # lowest value of minimize_key
        data = min(G.get_edge_data(u, v).values(), key=lambda x: x[minimize_key])
        attribute_values.append(data[attribute])
    return attribute_values


def count_streets_per_node(G, nodes=None):
    """
    Count how many street segments emanate from each node (i.e., intersections and dead-ends) in this graph.

    If nodes is passed, then only count the nodes in the graph with those IDs.

    Parameters
    ----------
    G : networkx multidigraph
    nodes : iterable
        the set of node IDs to get counts for

    Returns
    ----------
    streets_per_node : dict
        counts of how many streets emanate from each node with keys=node id and values=count
    """

    start_time = time.time()

    # to calculate the counts, get undirected representation of the graph. for
    # each node, get the list of the set of unique u,v,key edges, including
    # parallel edges but excluding self-loop parallel edges (this is necessary
    # because bi-directional self-loops will appear twice in the undirected
    # graph as you have u,v,key0 and u,v,key1 where u==v when you convert from
    # MultiDiGraph to MultiGraph - BUT, one-way self-loops will appear only
    # once. to get consistent accurate counts of physical streets, ignoring
    # directionality, we need the list of the set of unique edges...). then,
    # count how many times the node appears in the u,v tuples in the list. this
    # is the count of how many street segments emanate from this node. finally,
    # create a dict of node id:count
    G_undir = G.to_undirected(reciprocal=False)
    all_edges = G_undir.edges(keys=False)
    if nodes is None:
        nodes = G_undir.nodes()

    # get all unique edges - this throws away any parallel edges (including
    # those in self-loops)
    all_unique_edges = set(all_edges)

    # get all edges (including parallel edges) that are not self-loops
    non_self_loop_edges = [e for e in all_edges if not e[0]==e[1]]

    # get a single copy of each self-loop edge (ie, if it's bi-directional, we
    # ignore the parallel edge going the reverse direction and keep only one
    # copy)
    set_non_self_loop_edges = set(non_self_loop_edges)
    self_loop_edges = [e for e in all_unique_edges if e not in set_non_self_loop_edges]

    # final list contains all unique edges, including each parallel edge, unless
    # the parallel edge is a self-loop, in which case it doesn't double-count
    # the self-loop
    edges = non_self_loop_edges + self_loop_edges

    # flatten the list of (u,v) tuples
    edges_flat = list(chain.from_iterable(edges))

    # count how often each node appears in the list of flattened edge endpoints
    counts = Counter(edges_flat)
    streets_per_node = {node:counts[node] for node in nodes}
    msg = ('Got the counts of undirected street segments incident to each node '
           '(before removing peripheral edges) in {:,.2f} seconds')
    log(msg.format(time.time()-start_time))
    return streets_per_node


class OSMContentHandler (xml.sax.handler.ContentHandler):
    ''' SAX content handler for OSM XML.
    
        Used to build an Overpass-like response JSON object in self.object. For format
        notes, see http://wiki.openstreetmap.org/wiki/OSM_XML#OSM_XML_file_format_notes
        and http://overpass-api.de/output_formats.html#json
    '''
    def __init__(self):
        self._element = None
        self.object = {'elements': []}

    def startElement(self, name, attrs):
        if name == 'osm':
            self.object.update({k: attrs[k] for k in attrs.keys()
                if k in ('version', 'generator')})
        
        elif name in ('node', 'way'):
            self._element = dict(type=name, tags={}, nodes=[], **attrs)
            self._element.update({k: float(attrs[k]) for k in attrs.keys()
                if k in ('lat', 'lon')})
            self._element.update({k: int(attrs[k]) for k in attrs.keys()
                if k in ('id', 'uid', 'version', 'changeset')})
        
        elif name == 'tag':
            self._element['tags'].update({attrs['k']: attrs['v']})
        
        elif name == 'nd':
            self._element['nodes'].append(int(attrs['ref']))
        
        elif name == 'relation':
            # Placeholder for future relation support.
            # Look for nested members and tags.
            pass

    def endElement(self, name):
        if name in ('node', 'way'):
            self.object['elements'].append(self._element)


def overpass_json_from_file(filename):
    ''' Read OSM XML from input filename and return Overpass-like JSON.
    '''
    _, ext = os.path.splitext(filename)
    
    if ext == '.bz2':
        # Use Python 2/3 compatible BZ2File()
        opener = lambda fn: bz2.BZ2File(fn)
    else:
        # Assume an unrecognized file extension is just XML
        opener = lambda fn: open(fn, mode='rb')
    
    with opener(filename) as file:
        handler = OSMContentHandler()
        xml.sax.parse(file, handler)
        
        #import pprint; pprint.pprint(handler.object)
        
        return handler.object