import datetime as dt
import logging as lg
import networkx as nx
import numpy as np
import os
import sys
import unicodedata
from . import settings


def citation():
    """
    Print the OSMnx package's citation information.

    Boeing, G. 2017. OSMnx: New Methods for Acquiring, Constructing, Analyzing,
    and Visualizing Complex Street Networks. Computers, Environment and Urban
    Systems, 65(126-139). https://doi.org/10.1016/j.compenvurbsys.2017.05.004
    """

    cite = ("To cite OSMnx, use:\n\n"
            "Boeing, G. 2017. OSMnx: New Methods for Acquiring, Constructing, Analyzing, "
            "and Visualizing Complex Street Networks. Computers, Environment and Urban "
            "Systems, 65(126-139). https://doi.org/10.1016/j.compenvurbsys.2017.05.004"
            "\n\n"
            "BibTeX entry for LaTeX users:\n\n"

            "@article{boeing_osmnx_2017,\n"
            "    title = {{OSMnx}: {New} {Methods} for {Acquiring}, {Constructing}, {Analyzing}, and {Visualizing} {Complex} {Street} {Networks}},\n"
            "    volume = {65},\n"
            "    doi = {10.1016/j.compenvurbsys.2017.05.004},\n"
            "    number = {126-139},\n"
            "    journal = {Computers, Environment and Urban Systems},\n"
            "    author = {Boeing, Geoff},\n"
            "    year = {2017}\n"
            "}")

    print(cite)



def ts(style='datetime', template=None):
    """
    Get current timestamp as string

    Parameters
    ----------
    style : string
        format the timestamp with this built-in template. must be one
        of {'datetime', 'date', 'time'}
    template : string
        if not None, format the timestamp with this template

    Returns
    -------
    ts : string
        the string timestamp
    """

    if template is None:
        if style == 'datetime':
            template = '{:%Y-%m-%d %H:%M:%S}'
        elif style == 'date':
            template = '{:%Y-%m-%d}'
        elif style == 'time':
            template = '{:%H:%M:%S}'
        else:
            raise ValueError('unknown timestamp style "{}"'.format(style))

    ts = template.format(dt.datetime.now())
    return ts



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
    if (sys.version_info > (3, 0)):
        # python 3.x has no unicode type, so if error, use str type
        return str(value)
    else:
        # for python 2.x compatibility, use unicode
        return unicode(value)


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
           default_accept_language=settings.default_accept_language,
           nominatim_endpoint=settings.nominatim_endpoint,
           nominatim_key=settings.nominatim_key,
           overpass_endpoint=settings.overpass_endpoint,
           all_oneway=settings.all_oneway):
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
    nominatim_endpoint : string
        which API endpoint to use for nominatim queries
    nominatim_key : string
        your API key, if you are using an endpoint that requires one
    overpass_endpoint : string
        which API endpoint to use for overpass queries
    all_oneway : boolean
        if True, forces all paths to be loaded as oneway ways, preserving
        the original order of nodes stored in the OSM way XML.

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
    settings.useful_tags_node = list(set(useful_tags_node + osm_xml_node_attrs + osm_xml_node_tags))
    settings.useful_tags_path = list(set(useful_tags_path + osm_xml_way_attrs + osm_xml_way_tags))
    settings.osm_xml_node_attrs = osm_xml_node_attrs
    settings.osm_xml_node_tags = osm_xml_node_tags
    settings.osm_xml_way_attrs = osm_xml_way_attrs
    settings.osm_xml_way_tags = osm_xml_way_tags
    settings.default_access = default_access
    settings.default_crs = default_crs
    settings.default_user_agent = default_user_agent
    settings.default_referer = default_referer
    settings.default_accept_language = default_accept_language
    settings.nominatim_endpoint = nominatim_endpoint
    settings.nominatim_key = nominatim_key
    settings.overpass_endpoint = overpass_endpoint
    settings.all_oneway = all_oneway

    # if logging is turned on, log that we are configured
    if settings.log_file or settings.log_console:
        log('Configured osmnx')


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
    h = np.minimum(1.0, h)  # protect against floating point errors

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
        log_filename = os.path.join(settings.logs_folder, '{}_{}.log'.format(filename, ts(style='date')))

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


def get_unique_nodes_ordered_from_way(way_edges_df):
    """
    Function to recover the original order of nodes from a dataframe
    of edges associated with a single OSM way.

    Parameters
    ----------
    way_edges_df : pandas.DataFrame()
        Dataframe containing columns 'u' and 'v' corresponding to
        origin/desitination nodes.

    Returns
    -------
    unique_ordered_nodes : list
        An ordered list of unique node IDs

    NOTE: If the edges do not all connect (e.g. [(1, 2), (2,3),
    (10, 11), (11, 12), (12, 13)]), then this method will return
    only those nodes associated with the largest component of
    connected edges, even if subsequent connected chunks are contain
    more total nodes. This is done to ensure a proper topological
    representation of nodes in the XML way records because if there
    are unconnected components, the sorting algorithm cannot recover
    their original order. I don't believe that we would ever encounter
    this kind of disconnected structure of nodes within a given way,
    but as best I could tell it is not explicitly forbidden in the
    OSM XML design schema. I'm using a print statement right now to
    tell the user whether or not any nodes have been dropped and
    how many.
    """

    G = nx.MultiDiGraph()
    all_nodes = list(way_edges_df['u'].values) + \
        list(way_edges_df['v'].values)

    G.add_nodes_from(all_nodes)
    G.add_edges_from(way_edges_df[['u', 'v']].values)
    wccs = nx.weakly_connected_components(G)
    largest_wcc = max(wccs, key=len)
    node_subset = set(largest_wcc)

    # NOTE: this code (L387-403) is copied from geo_utils.py
    # which cannot be imported here without triggering a
    # circular import error. This should be fixed next time the
    # code base is refactored

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

    unique_ordered_nodes = list(nx.topological_sort(G2))
    num_unique_nodes = len(np.unique(all_nodes))

    if len(unique_ordered_nodes) < num_unique_nodes:
        print('Recovered order for {0} of {1} nodes'.format(
            len(unique_ordered_nodes), num_unique_nodes))

    return unique_ordered_nodes