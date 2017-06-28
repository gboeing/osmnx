###################################################################################################
# Module: utils.py
# Description: Utility functions for configuration, logging, etc
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
###################################################################################################

import os
import sys
import time
import unicodedata
import logging as lg
import datetime as dt
import networkx as nx
import numpy as np
import pandas as pd
import requests

from . import globals


def config(data_folder=globals.data_folder,
           logs_folder=globals.logs_folder,
           imgs_folder=globals.imgs_folder,
           cache_folder=globals.cache_folder,
           use_cache=globals.use_cache,
           log_file=globals.log_file,
           log_console=globals.log_console,
           log_level=globals.log_level,
           log_name=globals.log_name,
           log_filename=globals.log_filename,
           useful_tags_node=globals.useful_tags_node,
           useful_tags_path=globals.useful_tags_path):
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
        if True, use a local cache to save/retrieve http responses instead of calling API repetitively for the same request URL
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
    globals.use_cache = use_cache
    globals.cache_folder = cache_folder
    globals.data_folder = data_folder
    globals.imgs_folder = imgs_folder
    globals.logs_folder = logs_folder
    globals.log_console = log_console
    globals.log_file = log_file
    globals.log_level = log_level
    globals.log_name = log_name
    globals.log_filename = log_filename
    globals.useful_tags_node = useful_tags_node
    globals.useful_tags_path = useful_tags_path

    # if logging is turned on, log that we are configured
    if globals.log_file or globals.log_console:
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
        level = globals.log_level
    if name is None:
        name = globals.log_name
    if filename is None:
        filename = globals.log_filename

    # if logging to file is turned on
    if globals.log_file:
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
    if globals.log_console:
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
        level = globals.log_level
    if name is None:
        name = globals.log_name
    if filename is None:
        filename = globals.log_filename

    logger = lg.getLogger(name)

    # if a logger with this name is not already set up
    if not getattr(logger, 'handler_set', None):

        # get today's date and construct a log filename
        todays_date = dt.datetime.today().strftime('%Y_%m_%d')
        log_filename = '{}/{}_{}.log'.format(globals.logs_folder, filename, todays_date)

        # if the logs folder does not already exist, create it
        if not os.path.exists(globals.logs_folder):
            os.makedirs(globals.logs_folder)

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
    except:
        # python 3.x has no unicode type, so if error, use str type
        return str(value)


def get_largest_component(G, strongly=False):
    """
    Return the largest weakly or strongly connected component from a directed graph.

    Parameters
    ----------
    G : networkx multidigraph
    strongly : bool
        if True, return the largest strongly instead of weakly connected component

    Returns
    -------
    networkx multidigraph
    """

    start_time = time.time()
    original_len = len(list(G.nodes()))

    if strongly:
        # if the graph is not connected and caller did not request retain_all, retain only the largest strongly connected component
        if not nx.is_strongly_connected(G):
            G = max(nx.strongly_connected_component_subgraphs(G), key=len)
            log('Graph was not connected, retained only the largest strongly connected component ({:,} of {:,} total nodes) in {:.2f} seconds'.format(len(list(G.nodes())), original_len, time.time()-start_time))
    else:
        # if the graph is not connected and caller did not request retain_all, retain only the largest weakly connected component
        if not nx.is_weakly_connected(G):
            G = max(nx.weakly_connected_component_subgraphs(G), key=len)
            log('Graph was not connected, retained only the largest weakly connected component ({:,} of {:,} total nodes) in {:.2f} seconds'.format(len(list(G.nodes())), original_len, time.time()-start_time))

    return G


def great_circle_vec(lat1, lng1, lat2, lng2, earth_radius=6371009):
    """
    Vectorized function to calculate the great-circle distance between two points or between vectors of points.

    Parameters
    ----------
    lat1 : float or array of float
    lng1 : float or array of float
    lat2 : float or array of float
    lng2 : float or array of float
    earth_radius : numeric
        radius of earth in units in which distance will be returned (default is meters)

    Returns
    -------
    distance : float or array of float
        distance or vector of distances from (lat1, lng1) to (lat2, lng2) in units of earth_radius
    """

    phi1 = np.deg2rad(90 - lat1)
    phi2 = np.deg2rad(90 - lat2)

    theta1 = np.deg2rad(lng1)
    theta2 = np.deg2rad(lng2)

    cos = (np.sin(phi1) * np.sin(phi2) * np.cos(theta1 - theta2) + np.cos(phi1) * np.cos(phi2))
    arc = np.arccos(cos)

    # return distance in units of earth_radius
    distance = arc * earth_radius
    return distance


def get_nearest_node(G, point, return_dist=False):
    """
    Return the graph node nearest to some specified point.

    Parameters
    ----------
    G : networkx multidigraph
    point : tuple
        the (lat, lon) point for which we will find the nearest node in the graph
    return_dist : bool
        optionally also return the distance between the point and the nearest node

    Returns
    -------
    int or tuple of (int, float)
        nearest node ID or optionally (node ID, dist), where dist is the distance between
        the point and nearest node
    """
    start_time = time.time()

    # dump graph node coordinates into a pandas dataframe indexed by node id with x and y columns
    coords = np.array([[node, data['x'], data['y']] for node, data in G.nodes(data=True)])
    df = pd.DataFrame(coords, columns=['node', 'x', 'y']).set_index('node')

    # add columns to the dataframe representing the (constant) coordinates of the reference point
    df['reference_y'] = point[0]
    df['reference_x'] = point[1]

    # calculate the distance between each node and the reference point
    distances = great_circle_vec(lat1=df['reference_y'],
                                 lng1=df['reference_x'],
                                 lat2=df['y'],
                                 lng2=df['x'])

    # nearest node's ID is the index label of the minimum distance
    nearest_node = int(distances.idxmin())
    log('Found nearest node ({}) to point {} in {:,.2f} seconds'.format(nearest_node, point, time.time()-start_time))

    # if caller requested return_dist, return distance between the point and the nearest node as well
    if return_dist:
        return nearest_node, distances.loc[nearest_node]
    else:
        return nearest_node


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
        # if there are parallel edges between two nodes, select the one with the lowest value of minimize_key
        data = min([data for data in G.edge[u][v].values()], key=lambda x: x[minimize_key])
        attribute_values.append(data[attribute])
    return attribute_values
