import sys
import os
import datetime as dt
import unicodedata
import numpy as np
import logging as lg
from . import settings


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
           overpass_endpoint=settings.overpass_endpoint):
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
