"""General utility functions."""

import datetime as dt
import logging as lg
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

    Returns
    -------
    None
    """
    cite = (
        "Citation:\n\n"
        "Boeing, G. 2017. OSMnx: New Methods for Acquiring, "
        "Constructing, Analyzing, and Visualizing Complex Street "
        "Networks. Computers, Environment and Urban Systems, 65(126-139). "
        "https://doi.org/10.1016/j.compenvurbsys.2017.05.004\n\n"
        "BibTeX entry for LaTeX users:\n\n"
        "@article{boeing_osmnx_2017,\n"
        "    title = {{OSMnx}: {New} {Methods} for {Acquiring}, "
        "{Constructing}, {Analyzing}, and {Visualizing} {Complex} "
        "{Street} {Networks}},\n"
        "    volume = {65},\n"
        "    doi = {10.1016/j.compenvurbsys.2017.05.004},\n"
        "    number = {126-139},\n"
        "    journal = {Computers, Environment and Urban Systems},\n"
        "    author = {Boeing, Geoff},\n"
        "    year = {2017}\n"
        "}"
    )

    print(cite)


def ts(style="datetime", template=None):
    """
    Get current timestamp as string.

    Parameters
    ----------
    style : string
        format the timestamp with this built-in template. must be one
        of {'datetime', 'date', 'time'}
    template : string
        if not None, format the timestamp with this template instead of
        one of the built-in styles

    Returns
    -------
    ts : string
        the string timestamp
    """
    if template is None:
        if style == "datetime":
            template = "{:%Y-%m-%d %H:%M:%S}"
        elif style == "date":
            template = "{:%Y-%m-%d}"
        elif style == "time":
            template = "{:%H:%M:%S}"
        else:
            raise ValueError(f'unrecognized timestamp style "{style}"')

    ts = template.format(dt.datetime.now())
    return ts


def config(
    data_folder=settings.data_folder,
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
    all_oneway=settings.all_oneway,
):
    """
    Configure OSMnx by setting the default global settings' values.

    Note that any parameters not passed-in by the caller are set to their
    default values.

    Parameters
    ----------
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
    log_filename : string
        filename of the log
    useful_tags_node : list
        a list of useful OSM tags to attempt to save from node elements
    useful_tags_path : list
        a list of useful OSM tags to attempt to save from path elements
    osm_xml_node_attrs : list
        list of node attributes for .osm xml files
    osm_xml_node_tags : list
        list of node tags for .osm xml files
    osm_xml_way_attrs : list
        list of edge attributes for .osm xml files
    osm_xml_way_tags : list
        list of edge tags for .osm xml files
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
        the original order of nodes stored in the OSM way XML. Only use if
        specifically saving to .osm xml file with save_graph_xml function.

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
        log("Configured osmnx")


def log(message, level=None, name=None, filename=None):
    """
    Write a message to the logger.

    This logs to file and/or prints to the console, depending on the current
    configuration of settings.log_file and settings.log_console.

    Parameters
    ----------
    message : string
        the message to log
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
        logger = _get_logger(level=level, name=name, filename=filename)
        if level == lg.DEBUG:
            logger.debug(message)
        elif level == lg.INFO:
            logger.info(message)
        elif level == lg.WARNING:
            logger.warning(message)
        elif level == lg.ERROR:
            logger.error(message)

    # if logging to console, convert message to ascii and print to console
    if settings.log_console:
        # capture current stdout, then switch it to the console, print the
        # message, then switch back to what had been the stdout. this prevents
        # logging to notebook in jupyter: instead, it goes to terminal
        standard_out = sys.stdout
        sys.stdout = sys.__stdout__

        # prepend timestamp
        message = f"{ts()} {message}"

        # convert to ascii so it doesn't break windows terminals
        message = (
            unicodedata.normalize("NFKD", str(message)).encode("ascii", errors="replace").decode()
        )
        print(message)
        sys.stdout = standard_out


def _get_logger(level=None, name=None, filename=None):
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
    logger : logging.logger
    """
    if level is None:
        level = settings.log_level
    if name is None:
        name = settings.log_name
    if filename is None:
        filename = settings.log_filename

    logger = lg.getLogger(name)

    # if a logger with this name is not already set up
    if not getattr(logger, "handler_set", None):

        # get today's date and construct a log filename
        log_filename = os.path.join(settings.logs_folder, f'{filename}_{ts(style="date")}.log')

        # if the logs folder does not already exist, create it
        if not os.path.exists(settings.logs_folder):
            os.makedirs(settings.logs_folder)

        # create file handler and log formatter and set them up
        handler = lg.FileHandler(log_filename, encoding="utf-8")
        formatter = lg.Formatter("%(asctime)s %(levelname)s %(name)s %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(level)
        logger.handler_set = True

    return logger
