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
cache_folder : str or pathlib.Path
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
data_folder : str or pathlib.Path
    Path to folder in which to save/load graph files by default. Default is
    `"./data"`.
default_accept_language : str
    Do not use, deprecated. Use `http_accept_language` instead.
default_access : str
    Default filter for OSM "access" key. Default is `'["access"!~"private"]'`.
    Note that also filtering out "access=no" ways prevents including
    transit-only bridges (e.g., Tilikum Crossing) from appearing in drivable
    road network (e.g., `'["access"!~"private|no"]'`). However, some drivable
    tollroads have "access=no" plus a "access:conditional" key to clarify when
    it is accessible, so we can't filter out all "access=no" ways by default.
    Best to be permissive here then remove complicated combinations of tags
    programatically after the full graph is downloaded and constructed.
default_crs : str
    Default coordinate reference system to set when creating graphs. Default
    is `"epsg:4326"`.
default_referer : str
    Do not use, deprecated. Use `http_referer` instead.
default_user_agent : str
    Do not use, deprecated. Use `http_user_agent` instead.
doh_url_template : str
    Endpoint to resolve DNS-over-HTTPS if local DNS resolution fails. Set to
    None to disable DoH, but see `downloader._config_dns` documentation for
    caveats. Default is: `"https://8.8.8.8/resolve?name={hostname}"`
elevation_url_template : str
    Endpoint of the Google Maps Elevation API (or equivalent), containing
    exactly two parameters: `locations` and `key`. Default is:
    `"https://maps.googleapis.com/maps/api/elevation/json?locations={locations}&key={key}"`
    One example of an alternative equivalent would be Open Topo Data:
    `"https://api.opentopodata.org/v1/aster30m?locations={locations}&key={key}"`
http_accept_language : str
    HTTP header accept-language. Default is `"en"`. Note that Nominatim's
    default language is "en" and it can sort result importance scores
    differently if a different language is specified.
http_referer : str
    HTTP header referer. Default is
    `"OSMnx Python package (https://github.com/gboeing/osmnx)"`.
http_user_agent : str
    HTTP header user-agent. Default is
    `"OSMnx Python package (https://github.com/gboeing/osmnx)"`.
imgs_folder : str or pathlib.Path
    Path to folder in which to save plotted images by default. Default is
    `"./images"`.
log_file : bool
    If True, save log output to a file in logs_folder. Default is `False`.
log_filename : str
    Name of the log file, without file extension. Default is `"osmnx"`.
log_console : bool
    If True, print log output to the console (terminal window). Default is
    `False`.
log_level : int
    One of Python's logger.level constants. Default is `logging.INFO`.
log_name : str
    Name of the logger. Default is `"OSMnx"`.
logs_folder : str or pathlib.Path
    Path to folder in which to save log files. Default is `"./logs"`.
max_query_area_size : int
    Maximum area for any part of the geometry in meters: any polygon bigger
    than this will get divided up for multiple queries to the API. Default is
    `2500000000`.
memory : int
    Do not use, deprecated. Use `overpass_memory` instead.
nominatim_endpoint : str
    Do not use, deprecated. Use `nominatim_url` instead.
nominatim_key : str
    Your Nominatim API key, if you are using an API instance that requires
    one. Default is `None`.
nominatim_url : str
    The base API url to use for Nominatim queries. Default is
    `"https://nominatim.openstreetmap.org/"`.
osm_xml_node_attrs : list
    Do not use, deprecated.
osm_xml_node_tags : list
    Do not use, deprecated.
osm_xml_way_attrs : list
    Do not use, deprecated.
osm_xml_way_tags : list
    Do not use, deprecated.
overpass_endpoint : str
    Do not use, deprecated. Use `overpass_url` instead.
overpass_memory : int | None
    Overpass server memory allocation size for the query, in bytes. If
    None, server will choose its default allocation size. Use with caution.
    Default is `None`.
overpass_rate_limit : bool
    If True, check the Overpass server status endpoint for how long to
    pause before making request. Necessary if server uses slot management,
    but can be set to False if you are running your own overpass instance
    without rate limiting. Default is `True`.
overpass_settings : str
    Settings string for Overpass queries. Default is
    `"[out:json][timeout:{timeout}]{maxsize}"`. By default, the {timeout} and
    {maxsize} values are set dynamically by OSMnx when used.
    To query, for example, historical OSM data as of a certain date:
    `'[out:json][timeout:90][date:"2019-10-28T19:20:00Z"]'`. Use with caution.
overpass_url : str
    The base API url to use for Overpass queries. Default is
    `"https://overpass-api.de/api"`.
requests_kwargs : dict
    Optional keyword args to pass to the requests package when connecting
    to APIs, for example to configure authentication or provide a path to
    a local certificate file. More info on options such as auth, cert,
    verify, and proxies can be found in the requests package advanced docs.
    Default is `{}`.
requests_timeout : int
    The timeout interval in seconds for HTTP requests, and (when applicable)
    for Overpass server to use for executing the query. Default is `180`.
timeout : int
    Do not use, deprecated. Use `requests_timeout` instead.
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
default_accept_language = None
default_access = '["access"!~"private"]'
default_crs = "epsg:4326"
default_referer = None
default_user_agent = None
doh_url_template = "https://8.8.8.8/resolve?name={hostname}"
elevation_url_template = (
    "https://maps.googleapis.com/maps/api/elevation/json?locations={locations}&key={key}"
)
http_accept_language = "en"
http_referer = "OSMnx Python package (https://github.com/gboeing/osmnx)"
http_user_agent = "OSMnx Python package (https://github.com/gboeing/osmnx)"
imgs_folder = "./images"
log_console = False
log_file = False
log_filename = "osmnx"
log_level = lg.INFO
log_name = "OSMnx"
logs_folder = "./logs"
max_query_area_size = 50 * 1000 * 50 * 1000
memory = None
nominatim_endpoint = None
nominatim_key = None
nominatim_url = "https://nominatim.openstreetmap.org/"
osm_xml_node_attrs = None
osm_xml_node_tags = None
osm_xml_way_attrs = None
osm_xml_way_tags = None
overpass_endpoint = None
overpass_memory = None
overpass_rate_limit = True
overpass_settings = "[out:json][timeout:{timeout}]{maxsize}"
overpass_url = "https://overpass-api.de/api"
requests_kwargs: dict = {}
requests_timeout = 180
timeout = None
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
