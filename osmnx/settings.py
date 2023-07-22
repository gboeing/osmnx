"""
Global settings that can be configured by the user.

all_oneway : bool
    Only use if specifically saving to .osm XML file with the `save_graph_xml`
    function. If True, forces all ways to be loaded as oneway ways, preserving
    the original order of nodes stored in the OSM way XML. This also retains
    original OSM string values for oneway attribute values, rather than
    converting them to a True/False bool. Default is `False`.
bidirectional_network_types : list
    Network types for which a fully bidirectional graph will be created.
    Default is `["walk"]`.
cache_folder : string or pathlib.Path
    Path to folder in which to save/load HTTP response cache, if the
    `use_cache` setting equals `True`. Default is `"./cache"`.
cache_only_mode : bool
    If True, download network data from Overpass then raise a
    `CacheOnlyModeInterrupt` error for user to catch. This prevents graph
    building from taking place and instead just saves OSM response data to
    cache. Useful for sequentially caching lots of raw data (as you can
    only query Overpass one request at a time) then using the local cache to
    quickly build many graphs simultaneously with multiprocessing. Default is
    `False`.
data_folder : string or pathlib.Path
    Path to folder in which to save/load graph files by default. Default is
    `"./data"`.
default_accept_language : string
    HTTP header accept-language. Default is `"en"`.
default_access : string
    Default filter for OSM "access" key. Default is `'["access"!~"private"]'`.
    Note that also filtering out "access=no" ways prevents including
    transit-only bridges (e.g., Tilikum Crossing) from appearing in drivable
    road network (e.g., `'["access"!~"private|no"]'`). However, some drivable
    tollroads have "access=no" plus a "access:conditional" key to clarify when
    it is accessible, so we can't filter out all "access=no" ways by default.
    Best to be permissive here then remove complicated combinations of tags
    programatically after the full graph is downloaded and constructed.
default_crs : string
    Default coordinate reference system to set when creating graphs. Default
    is `"epsg:4326"`.
default_referer : string
    HTTP header referer. Default is
    `"OSMnx Python package (https://github.com/gboeing/osmnx)"`.
default_user_agent : string
    HTTP header user-agent. Default is
    `"OSMnx Python package (https://github.com/gboeing/osmnx)"`.
doh_url_template : string
    Endpoint to resolve DNS-over-HTTPS if local DNS resolution fails. Set to
    None to disable DoH, but see `downloader._config_dns` documentation for
    caveats. Default is: `"https://8.8.8.8/resolve?name={hostname}"`
imgs_folder : string or pathlib.Path
    Path to folder in which to save plotted images by default. Default is
    `"./images"`.
log_file : bool
    If True, save log output to a file in logs_folder. Default is `False`.
log_filename : string
    Name of the log file, without file extension. Default is `"osmnx"`.
log_console : bool
    If True, print log output to the console (terminal window). Default is
    `False`.
log_level : int
    One of Python's logger.level constants. Default is `logging.INFO`.
log_name : string
    Name of the logger. Default is `"OSMnx"`.
logs_folder : string or pathlib.Path
    Path to folder in which to save log files. Default is `"./logs"`.
max_query_area_size : int
    Maximum area for any part of the geometry in meters: any polygon bigger
    than this will get divided up for multiple queries to the API. Default is
    `2500000000`.
memory : int
    Overpass server memory allocation size for the query, in bytes. If
    None, server will use its default allocation size. Use with caution.
    Default is `None`.
nominatim_endpoint : string
    The base API url to use for Nominatim queries. Default is
    `"https://nominatim.openstreetmap.org/"`.
nominatim_key : string
    Your Nominatim API key, if you are using an API instance that requires
    one. Default is `None`.
osm_xml_node_attrs : list
    Node attributes for saving .osm XML files with `save_graph_xml` function.
    Default is `["id", "timestamp", "uid", "user", "version", "changeset",
    "lat", "lon"]`.
osm_xml_node_tags : list
    Node tags for saving .osm XML files with `save_graph_xml` function.
    Default is `["highway"]`.
osm_xml_way_attrs : list
    Edge attributes for saving .osm XML files with `save_graph_xml` function.
    Default is `["id", "timestamp", "uid", "user", "version", "changeset"]`.
osm_xml_way_tags : list
    Edge tags for for saving .osm XML files with `save_graph_xml` function.
    Default is `["highway", "lanes", "maxspeed", "name", "oneway"]`.
overpass_endpoint : string
    The base API url to use for overpass queries. Default is
    `"https://overpass-api.de/api"`.
overpass_rate_limit : bool
    If True, check the Overpass server status endpoint for how long to
    pause before making request. Necessary if server uses slot management,
    but can be set to False if you are running your own overpass instance
    without rate limiting. Default is `True`.
overpass_settings : string
    Settings string for Overpass queries. Default is
    `"[out:json][timeout:{timeout}]{maxsize}"`. By default, the {timeout} and
    {maxsize} values are set dynamically by OSMnx when used.
    To query, for example, historical OSM data as of a certain date:
    `'[out:json][timeout:90][date:"2019-10-28T19:20:00Z"]'`. Use with caution.
requests_kwargs : dict
    Optional keyword args to pass to the requests package when connecting
    to APIs, for example to configure authentication or provide a path to
    a local certificate file. More info on options such as auth, cert,
    verify, and proxies can be found in the requests package advanced docs.
    Default is `{}`.
timeout : int
    The timeout interval in seconds for HTTP requests, and (when applicable)
    for API to use while running the query. Default is `180`.
use_cache : bool
    If True, cache HTTP responses locally instead of calling API repeatedly
    for the same request. Default is `True`.
useful_tags_node : list
    OSM "node" tags to add as graph node attributes, when present in the data
    retrieved from OSM. Default is `["ref", "highway"]`.
useful_tags_way : list
    OSM "way" tags to add as graph edge attributes, when present in the data
    retrieved from OSM. Default is `["bridge", "tunnel", "oneway", "lanes",
    "ref", "name", "highway", "maxspeed", "service", "access", "area",
    "landuse", "width", "est_width", "junction"]`.
"""

import logging as lg

all_oneway = False
bidirectional_network_types = ["walk"]
cache_folder = "./cache"
cache_only_mode = False
data_folder = "./data"
default_accept_language = "en"
default_access = '["access"!~"private"]'
default_crs = "epsg:4326"
default_referer = "OSMnx Python package (https://github.com/gboeing/osmnx)"
default_user_agent = "OSMnx Python package (https://github.com/gboeing/osmnx)"
doh_url_template = "https://8.8.8.8/resolve?name={hostname}"
imgs_folder = "./images"
log_console = False
log_file = False
log_filename = "osmnx"
log_level = lg.INFO
log_name = "OSMnx"
logs_folder = "./logs"
max_query_area_size = 50 * 1000 * 50 * 1000
memory = None
nominatim_endpoint = "https://nominatim.openstreetmap.org/"
nominatim_key = None
osm_xml_node_attrs = ["id", "timestamp", "uid", "user", "version", "changeset", "lat", "lon"]
osm_xml_node_tags = ["highway"]
osm_xml_way_attrs = ["id", "timestamp", "uid", "user", "version", "changeset"]
osm_xml_way_tags = ["highway", "lanes", "maxspeed", "name", "oneway"]
overpass_endpoint = "https://overpass-api.de/api"
overpass_rate_limit = True
overpass_settings = "[out:json][timeout:{timeout}]{maxsize}"
requests_kwargs = {}
timeout = 180
use_cache = True
useful_tags_node = ["ref", "highway"]
useful_tags_way = [
    "bridge",
    "tunnel",
    "oneway",
    "lanes",
    "ref",
    "name",
    "highway",
    "maxspeed",
    "service",
    "access",
    "area",
    "landuse",
    "width",
    "est_width",
    "junction",
]
