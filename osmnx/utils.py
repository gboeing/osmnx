################################################################################
# Module: utils.py
# Description: Utility functions for configuration, logging, geocoding,
#              geospatial analysis, etc.
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

import bz2
import datetime as dt
import io
import logging as lg
import math
import os
import networkx as nx
import numpy as np
import pandas as pd
import requests
import sys
import time
import unicodedata
import warnings
import xml.sax
from collections import Counter
from itertools import chain
from shapely.geometry import Point
from shapely.geometry import MultiPoint
from shapely.geometry import LineString
from shapely.geometry import MultiLineString
from shapely.geometry import MultiPolygon
from shapely.geometry import Polygon

from . import settings

# scipy and sklearn are optional dependencies for faster nearest node search
try:
    from scipy.spatial import cKDTree
except ImportError as e:
    cKDTree = None
try:
    from sklearn.neighbors import BallTree
except ImportError as e:
    BallTree = None

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
           useful_tags_path=settings.useful_tags_path,
           osm_xml_node_attrs=settings.osm_xml_node_attrs,
           osm_xml_node_tags=settings.osm_xml_node_tags,
           osm_xml_way_attrs=settings.osm_xml_way_attrs,
           osm_xml_way_tags=settings.osm_xml_way_tags,
           default_access=settings.default_access,
           default_crs=settings.default_crs,
           default_user_agent=settings.default_user_agent,
           default_referer=settings.default_referer,
           default_accept_language=settings.default_accept_language):
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
    default_access : string
        default filter for OSM "access" key
    default_crs : string
        default CRS to set when creating graphs
    default_user_agent : string
        HTTP header user-agent
    default_referer : string
        HTTP header referer
    default_accept_language : string
        HTTP header accept-language

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
    settings.useful_tags_node = list(set(
        useful_tags_node + osm_xml_node_attrs + osm_xml_node_tags))
    settings.useful_tags_path = list(set(
        useful_tags_path + osm_xml_way_attrs + osm_xml_way_tags))
    settings.osm_xml_node_attrs = osm_xml_node_attrs
    settings.osm_xml_node_tags = osm_xml_node_tags
    settings.osm_xml_way_attrs = osm_xml_way_attrs
    settings.osm_xml_way_tags = osm_xml_way_tags
    settings.default_access = default_access
    settings.default_crs = default_crs
    settings.default_user_agent = default_user_agent
    settings.default_referer = default_referer
    settings.default_accept_language = default_accept_language

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
        log_filename = os.path.join(settings.logs_folder, '{}_{}.log'.format(filename, todays_date))

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



def induce_subgraph(G, node_subset):
    """
    Induce a subgraph of G.

    Parameters
    ----------
    G : networkx multidigraph
    node_subset : list-like
        the subset of nodes to induce a subgraph of G

    Returns
    -------
    G2 : networkx multidigraph
        the subgraph of G induced by node_subset
    """

    node_subset = set(node_subset)

    # copy nodes into new graph
    G2 = G.__class__()
    G2.add_nodes_from((n, G.nodes[n]) for n in node_subset)

    # copy edges to new graph, including parallel edges
    if G2.is_multigraph:
        G2.add_edges_from((n, nbr, key, d)
            for n, nbrs in G.adj.items() if n in node_subset
            for nbr, keydict in nbrs.items() if nbr in node_subset
            for key, d in keydict.items())
    else:
        G2.add_edges_from((n, nbr, d)
            for n, nbrs in G.adj.items() if n in node_subset
            for nbr, d in nbrs.items() if nbr in node_subset)

    # update graph attribute dict, and return graph
    G2.graph.update(G.graph)
    return G2



def get_largest_component(G, strongly=False):
    """
    Return a subgraph of the largest weakly or strongly connected component
    from a directed graph.

    Parameters
    ----------
    G : networkx multidigraph
    strongly : bool
        if True, return the largest strongly instead of weakly connected
        component

    Returns
    -------
    G : networkx multidigraph
        the largest connected component subgraph from the original graph
    """

    start_time = time.time()
    original_len = len(list(G.nodes()))

    if strongly:
        # if the graph is not connected retain only the largest strongly connected component
        if not nx.is_strongly_connected(G):

            # get all the strongly connected components in graph then identify the largest
            sccs = nx.strongly_connected_components(G)
            largest_scc = max(sccs, key=len)
            G = induce_subgraph(G, largest_scc)

            msg = ('Graph was not connected, retained only the largest strongly '
                   'connected component ({:,} of {:,} total nodes) in {:.2f} seconds')
            log(msg.format(len(list(G.nodes())), original_len, time.time()-start_time))
    else:
        # if the graph is not connected retain only the largest weakly connected component
        if not nx.is_weakly_connected(G):

            # get all the weakly connected components in graph then identify the largest
            wccs = nx.weakly_connected_components(G)
            largest_wcc = max(wccs, key=len)
            G = induce_subgraph(G, largest_wcc)

            msg = ('Graph was not connected, retained only the largest weakly '
                   'connected component ({:,} of {:,} total nodes) in {:.2f} seconds')
            log(msg.format(len(list(G.nodes())), original_len, time.time()-start_time))

    return G



def great_circle_vec(lat1, lng1, lat2, lng2, earth_radius=6371009):
    """
    Vectorized function to calculate the great-circle distance between two
    points or between vectors of points, using haversine.

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
    distance : float or vector of floats
        distance or vector of distances from (lat1, lng1) to (lat2, lng2) in
        units of earth_radius
    """

    phi1 = np.deg2rad(lat1)
    phi2 = np.deg2rad(lat2)
    d_phi = phi2 - phi1

    theta1 = np.deg2rad(lng1)
    theta2 = np.deg2rad(lng2)
    d_theta = theta2 - theta1

    h = np.sin(d_phi / 2) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(d_theta / 2) ** 2
    h = np.minimum(1.0, h) # protect against floating point errors

    arc = 2 * np.arcsin(np.sqrt(h))

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



def get_nearest_node(G, point, method='haversine', return_dist=False):
    """
    Return the graph node nearest to some specified (lat, lng) or (y, x) point,
    and optionally the distance between the node and the point. This function
    can use either a haversine or euclidean distance calculator.

    Parameters
    ----------
    G : networkx multidigraph
    point : tuple
        The (lat, lng) or (y, x) point for which we will find the nearest node
        in the graph
    method : str {'haversine', 'euclidean'}
        Which method to use for calculating distances to find nearest node.
        If 'haversine', graph nodes' coordinates must be in units of decimal
        degrees. If 'euclidean', graph nodes' coordinates must be projected.
    return_dist : bool
        Optionally also return the distance (in meters if haversine, or graph
        node coordinate units if euclidean) between the point and the nearest
        node.

    Returns
    -------
    int or tuple of (int, float)
        Nearest node ID or optionally a tuple of (node ID, dist), where dist is
        the distance (in meters if haversine, or graph node coordinate units
        if euclidean) between the point and nearest node
    """
    start_time = time.time()

    if not G or (G.number_of_nodes() == 0):
        raise ValueError('G argument must be not be empty or should contain at least one node')

    # dump graph node coordinates into a pandas dataframe indexed by node id
    # with x and y columns
    coords = [[node, data['x'], data['y']] for node, data in G.nodes(data=True)]
    df = pd.DataFrame(coords, columns=['node', 'x', 'y']).set_index('node')

    # add columns to the dataframe representing the (constant) coordinates of
    # the reference point
    df['reference_y'] = point[0]
    df['reference_x'] = point[1]

    # calculate the distance between each node and the reference point
    if method == 'haversine':
        # calculate distance vector using haversine (ie, for
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
        raise ValueError('method argument must be either "haversine" or "euclidean"')

    # nearest node's ID is the index label of the minimum distance
    nearest_node = distances.idxmin()
    log('Found nearest node ({}) to point {} in {:,.2f} seconds'.format(nearest_node, point, time.time()-start_time))

    # if caller requested return_dist, return distance between the point and the
    # nearest node as well
    if return_dist:
        return nearest_node, distances.loc[nearest_node]
    else:
        return nearest_node


def get_nearest_edge(G, point):
    """
    Return the nearest edge to a pair of coordinates. Pass in a graph and a tuple
    with the coordinates. We first get all the edges in the graph. Secondly we compute
    the euclidean distance from the coordinates to the segments determined by each edge.
    The last step is to sort the edge segments in ascending order based on the distance
    from the coordinates to the edge. In the end, the first element in the list of edges
    will be the closest edge that we will return as a tuple containing the shapely
    geometry and the u, v nodes.

    Parameters
    ----------
    G : networkx multidigraph
    point : tuple
        The (lat, lng) or (y, x) point for which we will find the nearest edge
        in the graph

    Returns
    -------
    closest_edge_to_point : tuple (shapely.geometry, u, v)
        A geometry object representing the segment and the coordinates of the two
        nodes that determine the edge section, u and v, the OSM ids of the nodes.
    """
    start_time = time.time()

    gdf = graph_to_gdfs(G, nodes=False, fill_edge_geometry=True)
    graph_edges = gdf[["geometry", "u", "v"]].values.tolist()

    edges_with_distances = [
        (
            graph_edge,
            Point(tuple(reversed(point))).distance(graph_edge[0])
        )
        for graph_edge in graph_edges
    ]

    edges_with_distances = sorted(edges_with_distances, key=lambda x: x[1])
    closest_edge_to_point = edges_with_distances[0][0]

    geometry, u, v = closest_edge_to_point

    log('Found nearest edge ({}) to point {} in {:,.2f} seconds'.format((u, v), point, time.time() - start_time))

    return geometry, u, v


def get_nearest_nodes(G, X, Y, method=None):
    """
    Return the graph nodes nearest to a list of points. Pass in points
    as separate vectors of X and Y coordinates. The 'kdtree' method
    is by far the fastest with large data sets, but only finds approximate
    nearest nodes if working in unprojected coordinates like lat-lng (it
    precisely finds the nearest node if working in projected coordinates).
    The 'balltree' method is second fastest with large data sets, but it
    is precise if working in unprojected coordinates like lat-lng.

    Parameters
    ----------
    G : networkx multidigraph
    X : list-like
        The vector of longitudes or x's for which we will find the nearest
        node in the graph
    Y : list-like
        The vector of latitudes or y's for which we will find the nearest
        node in the graph
    method : str {None, 'kdtree', 'balltree'}
        Which method to use for finding nearest node to each point.
        If None, we manually find each node one at a time using
        osmnx.utils.get_nearest_node and haversine. If 'kdtree' we use
        scipy.spatial.cKDTree for very fast euclidean search. If
        'balltree', we use sklearn.neighbors.BallTree for fast
        haversine search.

    Returns
    -------
    nn : array
        list of nearest node IDs
    """

    start_time = time.time()

    if method is None:

        # calculate nearest node one at a time for each point
        nn = [get_nearest_node(G, (y, x), method='haversine') for x, y in zip(X, Y)]

    elif method == 'kdtree':

        # check if we were able to import scipy.spatial.cKDTree successfully
        if not cKDTree:
            raise ImportError('The scipy package must be installed to use this optional feature.')

        # build a k-d tree for euclidean nearest node search
        nodes = pd.DataFrame({'x':nx.get_node_attributes(G, 'x'),
                              'y':nx.get_node_attributes(G, 'y')})
        tree = cKDTree(data=nodes[['x', 'y']], compact_nodes=True, balanced_tree=True)

        # query the tree for nearest node to each point
        points = np.array([X, Y]).T
        dist, idx = tree.query(points, k=1)
        nn = nodes.iloc[idx].index

    elif method == 'balltree':

        # check if we were able to import sklearn.neighbors.BallTree successfully
        if not BallTree:
            raise ImportError('The scikit-learn package must be installed to use this optional feature.')

        # haversine requires data in form of [lat, lng] and inputs/outputs in units of radians
        nodes = pd.DataFrame({'x':nx.get_node_attributes(G, 'x'),
                              'y':nx.get_node_attributes(G, 'y')})
        nodes_rad = np.deg2rad(nodes[['y', 'x']].astype(np.float))
        points = np.array([Y.astype(np.float), X.astype(np.float)]).T
        points_rad = np.deg2rad(points)

        # build a ball tree for haversine nearest node search
        tree = BallTree(nodes_rad, metric='haversine')

        # query the tree for nearest node to each point
        idx = tree.query(points_rad, k=1, return_distance=False)
        nn = nodes.iloc[idx[:,0]].index

    else:
        raise ValueError('You must pass a valid method name, or None.')

    log('Found nearest nodes to {:,} points in {:,.2f} seconds'.format(len(X), time.time()-start_time))

    return np.array(nn)


def get_nearest_edges(G, X, Y, method=None, dist=0.0001):
    """
    Return the graph edges nearest to a list of points. Pass in points
    as separate vectors of X and Y coordinates. The 'kdtree' method
    is by far the fastest with large data sets, but only finds approximate
    nearest edges if working in unprojected coordinates like lat-lng (it
    precisely finds the nearest edge if working in projected coordinates).
    The 'balltree' method is second fastest with large data sets, but it
    is precise if working in unprojected coordinates like lat-lng.

    Parameters
    ----------
    G : networkx multidigraph
    X : list-like
        The vector of longitudes or x's for which we will find the nearest
        edge in the graph. For projected graphs, use the projected coordinates,
        usually in meters.
    Y : list-like
        The vector of latitudes or y's for which we will find the nearest
        edge in the graph. For projected graphs, use the projected coordinates,
        usually in meters.
    method : str {None, 'kdtree', 'balltree'}
        Which method to use for finding nearest edge to each point.
        If None, we manually find each edge one at a time using
        osmnx.utils.get_nearest_edge. If 'kdtree' we use
        scipy.spatial.cKDTree for very fast euclidean search. Recommended for
        projected graphs. If 'balltree', we use sklearn.neighbors.BallTree for
        fast haversine search. Recommended for unprojected graphs.

    dist : float
        spacing length along edges. Units are the same as the geom; Degrees for
        unprojected geometries and meters for projected geometries. The smaller
        the value, the more points are created.

    Returns
    -------
    ne : ndarray
        array of nearest edges represented by their startpoint and endpoint ids,
        u and v, the OSM ids of the nodes.

    Info
    ----
    The method creates equally distanced points along the edges of the network.
    Then, these points are used in a kdTree or BallTree search to identify which
    is nearest.Note that this method will not give the exact perpendicular point
    along the edge, but the smaller the *dist* parameter, the closer the solution
    will be.

    Code is adapted from an answer by JHuw from this original question:
    https://gis.stackexchange.com/questions/222315/geopandas-find-nearest-point
    -in-other-dataframe
    """
    start_time = time.time()

    if method is None:
        # calculate nearest edge one at a time for each point
        ne = [get_nearest_edge(G, (x, y)) for x, y in zip(X, Y)]
        ne = [(u, v) for _, u, v in ne]

    elif method == 'kdtree':

        # check if we were able to import scipy.spatial.cKDTree successfully
        if not cKDTree:
            raise ImportError('The scipy package must be installed to use this optional feature.')

        # transform graph into DataFrame
        edges = graph_to_gdfs(G, nodes=False, fill_edge_geometry=True)

        # transform edges into evenly spaced points
        edges['points'] = edges.apply(lambda x: redistribute_vertices(x.geometry, dist), axis=1)

        # develop edges data for each created points
        extended = edges['points'].apply([pd.Series]).stack().reset_index(level=1, drop=True).join(edges).reset_index()

        # Prepare btree arrays
        nbdata = np.array(list(zip(extended['Series'].apply(lambda x: x.x),
                                   extended['Series'].apply(lambda x: x.y))))

        # build a k-d tree for euclidean nearest node search
        btree = cKDTree(data=nbdata, compact_nodes=True, balanced_tree=True)

        # query the tree for nearest node to each point
        points = np.array([X, Y]).T
        dist, idx = btree.query(points, k=1)  # Returns ids of closest point
        eidx = extended.loc[idx, 'index']
        ne = edges.loc[eidx, ['u', 'v']]

    elif method == 'balltree':

        # check if we were able to import sklearn.neighbors.BallTree successfully
        if not BallTree:
            raise ImportError('The scikit-learn package must be installed to use this optional feature.')

        # transform graph into DataFrame
        edges = graph_to_gdfs(G, nodes=False, fill_edge_geometry=True)

        # transform edges into evenly spaced points
        edges['points'] = edges.apply(lambda x: redistribute_vertices(x.geometry, dist), axis=1)

        # develop edges data for each created points
        extended = edges['points'].apply([pd.Series]).stack().reset_index(level=1, drop=True).join(edges).reset_index()

        # haversine requires data in form of [lat, lng] and inputs/outputs in units of radians
        nodes = pd.DataFrame({'x': extended['Series'].apply(lambda x: x.x),
                              'y': extended['Series'].apply(lambda x: x.y)})
        nodes_rad = np.deg2rad(nodes[['y', 'x']].values.astype(np.float))
        points = np.array([Y, X]).T
        points_rad = np.deg2rad(points)

        # build a ball tree for haversine nearest node search
        tree = BallTree(nodes_rad, metric='haversine')

        # query the tree for nearest node to each point
        idx = tree.query(points_rad, k=1, return_distance=False)
        eidx = extended.loc[idx[:, 0], 'index']
        ne = edges.loc[eidx, ['u', 'v']]

    else:
        raise ValueError('You must pass a valid method name, or None.')

    log('Found nearest edges to {:,} points in {:,.2f} seconds'.format(len(X), time.time() - start_time))

    return np.array(ne)


def redistribute_vertices(geom, dist):
    """
    Redistribute the vertices on a projected LineString or MultiLineString. The distance
    argument is only approximate since the total distance of the linestring may not be
    a multiple of the preferred distance. This function works on only [Multi]LineString
    geometry types.

    This code is adapted from an answer by Mike T from this original question:
    https://stackoverflow.com/questions/34906124/interpolating-every-x-distance-along-multiline-in-shapely

    Parameters
    ----------
    geom : LineString or MultiLineString
        a Shapely geometry
    dist : float
        spacing length along edges. Units are the same as the geom; Degrees for unprojected geometries and meters
        for projected geometries. The smaller the value, the more points are created.

    Returns
    -------
        list of Point geometries : list
    """
    if geom.geom_type == 'LineString':
        num_vert = int(round(geom.length / dist))
        if num_vert == 0:
            num_vert = 1
        return [geom.interpolate(float(n) / num_vert, normalized=True)
                for n in range(num_vert + 1)]
    elif geom.geom_type == 'MultiLineString':
        parts = [redistribute_vertices(part, dist)
                 for part in geom]
        return type(geom)([p for p in parts if not p.is_empty])
    else:
        raise ValueError('unhandled geometry {}'.format(geom.geom_type))


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

        if u == v:
            # a self-loop has an undefined compass bearing
            data['bearing'] = np.nan

        else:
            # calculate bearing from edge's origin to its destination
            origin_point = (G.nodes[u]['y'], G.nodes[u]['x'])
            destination_point = (G.nodes[v]['y'], G.nodes[v]['x'])
            bearing = get_bearing(origin_point, destination_point)

            # round to thousandth of a degree
            data['bearing'] = round(bearing, 3)

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



def get_route_edge_attributes(G, route, attribute=None, minimize_key='length', retrieve_default=None):
    """
    Get a list of attribute values for each edge in a path.

    Parameters
    ----------
    G : networkx multidigraph
    route : list
        list of nodes in the path
    attribute : string
        the name of the attribute to get the value of for each edge.
        If not specified, the complete data dict is returned for each edge.
    minimize_key : string
        if there are parallel edges between two nodes, select the one with the
        lowest value of minimize_key
    retrieve_default : Callable[Tuple[Any, Any], Any]
        Function called with the edge nodes as parameters to retrieve a default value, if the edge does not
        contain the given attribute. Per default, a `KeyError` is raised
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
        if attribute is None:
            attribute_value = data
        elif retrieve_default is not None:
            attribute_value = data.get(attribute, retrieve_default(u, v))
        else:
            attribute_value = data[attribute]
        attribute_values.append(attribute_value)
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



def round_polygon_coords(p, precision):
    """
    Round the coordinates of a shapely Polygon to some decimal precision.

    Parameters
    ----------
    p : shapely Polygon
        the polygon to round the coordinates of
    precision : int
        decimal precision to round coordinates to

    Returns
    -------
    new_poly : shapely Polygon
        the polygon with rounded coordinates
    """

    # round the coordinates of the Polygon exterior
    new_exterior = [[round(x, precision) for x in c] for c in p.exterior.coords]

    # round the coordinates of the (possibly multiple, possibly none) Polygon interior(s)
    new_interiors = []
    for interior in p.interiors:
        new_interiors.append([[round(x, precision) for x in c] for c in interior.coords])

    # construct a new Polygon with the rounded coordinates
    # buffer by zero to clean self-touching or self-crossing polygons
    new_poly = Polygon(shell=new_exterior, holes=new_interiors).buffer(0)
    return new_poly



def round_multipolygon_coords(mp, precision):
    """
    Round the coordinates of a shapely MultiPolygon to some decimal precision.

    Parameters
    ----------
    mp : shapely MultiPolygon
        the MultiPolygon to round the coordinates of
    precision : int
        decimal precision to round coordinates to

    Returns
    -------
    MultiPolygon
    """

    return MultiPolygon([round_polygon_coords(p, precision) for p in mp])



def round_point_coords(pt, precision):
    """
    Round the coordinates of a shapely Point to some decimal precision.

    Parameters
    ----------
    pt : shapely Point
        the Point to round the coordinates of
    precision : int
        decimal precision to round coordinates to

    Returns
    -------
    Point
    """

    return Point([round(x, precision) for x in pt.coords[0]])



def round_multipoint_coords(mpt, precision):
    """
    Round the coordinates of a shapely MultiPoint to some decimal precision.

    Parameters
    ----------
    mpt : shapely MultiPoint
        the MultiPoint to round the coordinates of
    precision : int
        decimal precision to round coordinates to

    Returns
    -------
    MultiPoint
    """

    return MultiPoint([round_point_coords(pt, precision) for pt in mpt])



def round_linestring_coords(ls, precision):
    """
    Round the coordinates of a shapely LineString to some decimal precision.

    Parameters
    ----------
    ls : shapely LineString
        the LineString to round the coordinates of
    precision : int
        decimal precision to round coordinates to

    Returns
    -------
    LineString
    """

    return LineString([[round(x, precision) for x in c] for c in ls.coords])



def round_multilinestring_coords(mls, precision):
    """
    Round the coordinates of a shapely MultiLineString to some decimal precision.

    Parameters
    ----------
    mls : shapely MultiLineString
        the MultiLineString to round the coordinates of
    precision : int
        decimal precision to round coordinates to

    Returns
    -------
    MultiLineString
    """

    return MultiLineString([round_linestring_coords(ls, precision) for ls in mls])



def round_shape_coords(shape, precision):
    """
    Round the coordinates of a shapely geometry to some decimal precision.

    Parameters
    ----------
    shape : shapely geometry, one of Point, MultiPoint, LineString,
            MultiLineString, Polygon, or MultiPolygon
        the geometry to round the coordinates of
    precision : int
        decimal precision to round coordinates to

    Returns
    -------
    shapely geometry
    """

    if isinstance(shape, Point):
        return round_point_coords(shape, precision)

    elif isinstance(shape, MultiPoint):
        return round_multipoint_coords(shape, precision)

    elif isinstance(shape, LineString):
        return round_linestring_coords(shape, precision)

    elif isinstance(shape, MultiLineString):
        return round_multilinestring_coords(shape, precision)

    elif isinstance(shape, Polygon):
        return round_polygon_coords(shape, precision)

    elif isinstance(shape, MultiPolygon):
        return round_multipolygon_coords(shape, precision)

    else:
        raise TypeError('cannot round coordinates of unhandled geometry type: {}'.format(type(shape)))



class OSMContentHandler (xml.sax.handler.ContentHandler):
    """
    SAX content handler for OSM XML.

    Used to build an Overpass-like response JSON object in self.object. For format
    notes, see http://wiki.openstreetmap.org/wiki/OSM_XML#OSM_XML_file_format_notes
    and http://overpass-api.de/output_formats.html#json
    """

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
    """
    Read OSM XML from input filename and return Overpass-like JSON.

    Parameters
    ----------
    filename : string
        name of file containing OSM XML data

    Returns
    -------
    OSMContentHandler object
    """

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
        return handler.object



def bbox_to_poly(north, south, east, west):
    """
    Convenience function to parse bbox -> poly
    """

    return Polygon([(west, south), (east, south), (east, north), (west, north)])



def citation():
    """
    Print the OSMnx package's citation information.

    Boeing, G. 2017. OSMnx: New Methods for Acquiring, Constructing, Analyzing,
    and Visualizing Complex Street Networks. Computers, Environment and Urban
    Systems, 65(126-139). doi:10.1016/j.compenvurbsys.2017.05.004
    """

    cite = ("To cite OSMnx, use:\n\n"
            "Boeing, G. 2017. OSMnx: New Methods for Acquiring, Constructing, Analyzing, "
            "and Visualizing Complex Street Networks. Computers, Environment and Urban "
            "Systems, 65(126-139). doi:10.1016/j.compenvurbsys.2017.05.004"
            "\n\n"
            "BibTeX entry for LaTeX users:\n\n"

            "@article{boeing_osmnx_2017,\n"
            "    title = {{OSMnx}: {New} {Methods} for {Acquiring}, {Constructing}, {Analyzing}, and {Visualizing} {Complex} {Street} {Networks}},\n"
            "    volume = {65},\n"
            "    doi = {10.1016/j.compenvurbsys.2017.05.004},\n"
            "    number = {126-139},\n"
            "    journal = {Computers, Environment and Urban Systems},\n"
            "    author = {Boeing, G.},\n"
            "    year = {2017}\n"
            "}")

    print(cite)

from .save_load import graph_to_gdfs
