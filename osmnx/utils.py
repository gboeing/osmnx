import os
import json
import io
import sys
import time
import hashlib
import unicodedata
import logging as lg
import datetime as dt
import networkx as nx
from . import globals


def config(data_folder=globals.global_data_folder, 
           logs_folder=globals.global_logs_folder, 
           imgs_folder=globals.global_imgs_folder, 
           cache_folder=globals.global_cache_folder, 
           use_cache=globals.global_use_cache,
           log_file=globals.global_log_file, 
           log_console=globals.global_log_console, 
           log_level=globals.global_log_level, 
           log_name=globals.global_log_name, 
           log_filename=globals.global_log_filename,
           useful_tags_node=globals.global_useful_tags_node, 
           useful_tags_path=globals.global_useful_tags_path):
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
    
    # set each global variable to the passed-in parameter value
    globals.global_use_cache = use_cache
    globals.global_cache_folder = cache_folder
    globals.global_data_folder = data_folder
    globals.global_imgs_folder = imgs_folder
    globals.global_logs_folder = logs_folder
    globals.global_log_console = log_console
    globals.global_log_file = log_file
    globals.global_log_level = log_level
    globals.global_log_name = log_name
    globals.global_log_filename = log_filename
    globals.global_useful_tags_node = useful_tags_node
    globals.global_useful_tags_path = useful_tags_path
    
    # if logging is turned on, log that we are configured
    if globals.global_log_file or globals.global_log_console:
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
        level = globals.global_log_level
    if name is None:
        name = globals.global_log_name
    if filename is None:
        filename = globals.global_log_filename
    
    # if logging to file is turned on
    if globals.global_log_file:
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
    if globals.global_log_console:
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
        level = globals.global_log_level
    if name is None:
        name = globals.global_log_name
    if filename is None:
        filename = globals.global_log_filename
    
    logger = lg.getLogger(name)
    
    # if a logger with this name is not already set up
    if not getattr(logger, 'handler_set', None):
        
        # get today's date and construct a log filename
        todays_date = dt.datetime.today().strftime('%Y_%m_%d')
        log_filename = '{}/{}_{}.log'.format(globals.global_logs_folder, filename, todays_date)
        
        # if the logs folder does not already exist, create it
        if not os.path.exists(globals.global_logs_folder):
            os.makedirs(globals.global_logs_folder)
            
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
    if globals.global_use_cache:
        if response_json is None:
            log('Saved nothing to cache because response_json is None')
        else:        
            # create the folder on the disk if it doesn't already exist
            if not os.path.exists(globals.global_cache_folder):
                os.makedirs(globals.global_cache_folder)

            # hash the url (to make filename shorter than the often extremely long url) 
            filename = hashlib.md5(url.encode('utf-8')).hexdigest()
            cache_path_filename = '{}/{}.json'.format(globals.global_cache_folder, filename)
            
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
    if globals.global_use_cache:
        # determine the filename by hashing the url
        filename = hashlib.md5(url.encode('utf-8')).hexdigest()
        cache_path_filename = '{}/{}.json'.format(globals.global_cache_folder, filename)
        # open the cache file for this url hash if it already exists, otherwise return None
        if os.path.isfile(cache_path_filename):
            response_json = json.load(io.open(cache_path_filename, encoding='utf-8'))
            log('Retrieved response from cache file "{}" for URL "{}"'.format(cache_path_filename, url))
            return response_json        
            
            
            
            
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
    
    
    