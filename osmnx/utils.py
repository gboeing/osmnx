"""General utility functions."""

import datetime as dt
import logging as lg
import sys
import unicodedata
from pathlib import Path

from . import _version
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
    style : string {"datetime", "date", "time"}
        format the timestamp with this built-in template
    template : string
        if not None, format the timestamp with this template instead of one of
        the built-in styles

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
    all_oneway=settings.all_oneway,
    bidirectional_network_types=settings.bidirectional_network_types,
    cache_folder=settings.cache_folder,
    cache_only_mode=settings.cache_only_mode,
    data_folder=settings.data_folder,
    default_accept_language=settings.default_accept_language,
    default_access=settings.default_access,
    default_crs=settings.default_crs,
    default_referer=settings.default_referer,
    default_user_agent=settings.default_user_agent,
    elevation_provider=settings.elevation_provider,
    imgs_folder=settings.imgs_folder,
    log_console=settings.log_console,
    log_file=settings.log_file,
    log_filename=settings.log_filename,
    log_level=settings.log_level,
    log_name=settings.log_name,
    logs_folder=settings.logs_folder,
    max_query_area_size=settings.max_query_area_size,
    memory=settings.memory,
    nominatim_endpoint=settings.nominatim_endpoint,
    nominatim_key=settings.nominatim_key,
    osm_xml_node_attrs=settings.osm_xml_node_attrs,
    osm_xml_node_tags=settings.osm_xml_node_tags,
    osm_xml_way_attrs=settings.osm_xml_way_attrs,
    osm_xml_way_tags=settings.osm_xml_way_tags,
    overpass_endpoint=settings.overpass_endpoint,
    overpass_settings=settings.overpass_settings,
    timeout=settings.timeout,
    use_cache=settings.use_cache,
    useful_tags_node=settings.useful_tags_node,
    useful_tags_way=settings.useful_tags_way,
):
    """
    Configure OSMnx by setting the default global settings' values.

    Any parameters not passed by the caller are (re-)set to their original
    default values.

    Parameters
    ----------
    all_oneway : bool
        Only use if specifically saving to .osm XML file with save_graph_xml
        function. if True, forces all ways to be loaded as oneway ways,
        preserving the original order of nodes stored in the OSM way XML.
    bidirectional_network_types : list
        network types for which a fully bidirectional graph will be created
    cache_folder : string or pathlib.Path
        path to folder in which to save/load HTTP response cache
    data_folder : string or pathlib.Path
        path to folder in which to save/load graph files by default
    cache_only_mode : bool
        If True, download network data from Overpass then raise a
        CacheOnlyModeInterrupt error for user to catch. This prevents graph
        building from taking place and instead just saves OSM response data to
        cache. Useful for sequentially caching lots of raw data (as you can
        only query Overpass one request at a time) then using the cache to
        quickly build many graphs simultaneously with multiprocessing.
    default_accept_language : string
        HTTP header accept-language
    default_access : string
        default filter for OSM "access" key
    default_crs : string
        default coordinate reference system to set when creating graphs
    default_referer : string
        HTTP header referer
    default_user_agent : string
        HTTP header user-agent
    elevation_provider : string {"google", "airmap"}
        the API provider to use for adding node elevations
    imgs_folder : string or pathlib.Path
        path to folder in which to save plot images by default
    log_file : bool
        if True, save log output to a file in logs_folder
    log_filename : string
        name of the log file, without file extension
    log_console : bool
        if True, print log output to the console (terminal window)
    log_level : int
        one of Python's logger.level constants
    log_name : string
        name of the logger
    logs_folder : string or pathlib.Path
        path to folder in which to save log files
    max_query_area_size : int
        maximum area for any part of the geometry in meters: any polygon
        bigger than this will get divided up for multiple queries to API
        (default 50km x 50km)
    memory : int
        Overpass server memory allocation size for the query, in bytes. If
        None, server will use its default allocation size. Use with caution.
    nominatim_endpoint : string
        the API endpoint to use for nominatim queries
    nominatim_key : string
        your API key, if you are using an endpoint that requires one
    osm_xml_node_attrs : list
        node attributes for saving .osm XML files with save_graph_xml function
    osm_xml_node_tags : list
        node tags for saving .osm XML files with save_graph_xml function
    osm_xml_way_attrs : list
        edge attributes for saving .osm XML files with save_graph_xml function
    osm_xml_way_tags : list
        edge tags for for saving .osm XML files with save_graph_xml function
    overpass_endpoint : string
        the API endpoint to use for overpass queries
    overpass_settings : string
        Settings string for overpass queries. For example, to query historical
        OSM data as of a certain date:
        `'[out:json][timeout:90][date:"2019-10-28T19:20:00Z"]'`.
        Use with caution.
    timeout : int
        the timeout interval for the HTTP request and for API to use while
        running the query
    use_cache : bool
        if True, cache HTTP responses locally instead of calling API
        repeatedly for the same request
    useful_tags_node : list
        OSM "node" tags to add as graph node attributes, when present
    useful_tags_way : list
        OSM "way" tags to add as graph edge attributes, when present

    Returns
    -------
    None
    """
    # set each global setting to the argument value
    settings.all_oneway = all_oneway
    settings.bidirectional_network_types = bidirectional_network_types
    settings.cache_folder = cache_folder
    settings.cache_only_mode = cache_only_mode
    settings.data_folder = data_folder
    settings.default_accept_language = default_accept_language
    settings.default_access = default_access
    settings.default_crs = default_crs
    settings.default_referer = default_referer
    settings.default_user_agent = default_user_agent
    settings.elevation_provider = elevation_provider
    settings.imgs_folder = imgs_folder
    settings.log_console = log_console
    settings.log_file = log_file
    settings.log_filename = log_filename
    settings.log_level = log_level
    settings.log_name = log_name
    settings.logs_folder = logs_folder
    settings.max_query_area_size = max_query_area_size
    settings.memory = memory
    settings.nominatim_endpoint = nominatim_endpoint
    settings.nominatim_key = nominatim_key
    settings.osm_xml_node_attrs = osm_xml_node_attrs
    settings.osm_xml_node_tags = osm_xml_node_tags
    settings.osm_xml_way_attrs = osm_xml_way_attrs
    settings.osm_xml_way_tags = osm_xml_way_tags
    settings.overpass_endpoint = overpass_endpoint
    settings.overpass_settings = overpass_settings
    settings.timeout = timeout
    settings.use_cache = use_cache
    settings.useful_tags_node = useful_tags_node
    settings.useful_tags_way = useful_tags_way

    log(f"Configured OSMnx {_version.__version__}")
    log(f"HTTP response caching is {'on' if settings.use_cache else 'off'}")


def log(message, level=None, name=None, filename=None):
    """
    Write a message to the logger.

    This logs to file and/or prints to the console (terminal), depending on
    the current configuration of settings.log_file and settings.log_console.

    Parameters
    ----------
    message : string
        the message to log
    level : int
        one of Python's logger.level constants
    name : string
        name of the logger
    filename : string
        name of the log file, without file extension

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
        one of Python's logger.level constants
    name : string
        name of the logger
    filename : string
        name of the log file, without file extension

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
        log_filename = Path(settings.logs_folder) / f'{filename}_{ts(style="date")}.log'

        # if the logs folder does not already exist, create it
        log_filename.parent.mkdir(parents=True, exist_ok=True)

        # create file handler and log formatter and set them up
        handler = lg.FileHandler(log_filename, encoding="utf-8")
        formatter = lg.Formatter("%(asctime)s %(levelname)s %(name)s %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(level)
        logger.handler_set = True

    return logger
