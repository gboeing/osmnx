"""Global settings that can be configured by user with utils.config()."""

import logging as lg

import requests

# default locations to save data, logs, images, and cache
data_folder = "./data"
logs_folder = "./logs"
imgs_folder = "./images"
cache_folder = "./cache"

# cache server responses
use_cache = True

# if True, download net data from Overpass then raise CacheOnlyModeInterrupt.
# useful for sequentially caching lots of raw data then using that cache to
# quickly build many graphs simultaneously with multiprocessing
cache_only_mode = False

# write log to file and/or to console
log_file = False
log_console = False
log_level = lg.INFO
log_name = "OSMnx"
log_filename = "osmnx"

# OSM node/way tags to add as graph node/edge attributes when these tags are
# present in the data retrieved from OSM
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

# tags and attributes for generating an OSM XML file
osm_xml_node_attrs = ["id", "timestamp", "uid", "user", "version", "changeset", "lat", "lon"]
osm_xml_node_tags = ["highway"]
osm_xml_way_attrs = ["id", "timestamp", "uid", "user", "version", "changeset"]
osm_xml_way_tags = ["highway", "lanes", "maxspeed", "name", "oneway"]

# default settings string for overpass queries: timeout and maxsize need to
# be dynamically set where used
overpass_settings = "[out:json][timeout:{timeout}]{maxsize}"

# the timeout interval for the HTTP request and for API to use while running
# the query
timeout = 180

# overpass server memory allocation size for the query, in bytes. If None,
# server will use its default allocation size
memory = None

# maximum area for any part of the geometry in meters: any polygon bigger than
# this will get divided up for multiple queries to API (default 50km x 50km)
max_query_area_size = 50 * 1000 * 50 * 1000

# default filter for OSM "access" key. filtering out "access=no" ways prevents
# including transit-only bridges like tilikum crossing from appearing in drivable
# road network (e.g., '["access"!~"private|no"]'). however, some drivable
# tollroads have "access=no" plus a "access:conditional" key to clarify when
# it is accessible, so we can't filter out all "access=no" ways by default.
# best to be permissive here then remove complicated combinations of tags
# programatically after the full graph is downloaded and constructed.
default_access = '["access"!~"private"]'

# The network types for which a bidirectional graph will be created
bidirectional_network_types = ["walk"]

# all one-way mode to maintain original OSM node order when creating graphs
# specifically to save to .osm xml file with save_graph_xml function.
# this also retains original OSM string values for oneway attribute rather
# than converting it to a True/False bool
all_oneway = False

# default CRS to set when creating graphs
default_crs = "epsg:4326"

# default HTTP request headers
default_user_agent = "OSMnx Python package (https://github.com/gboeing/osmnx)"
default_referer = "OSMnx Python package (https://github.com/gboeing/osmnx)"
default_accept_language = "en"

# base API endpoint to use for nominatim queries
# and your nominatim API key, if you are using a host that requires one
nominatim_endpoint = "https://nominatim.openstreetmap.org/"
nominatim_key = None

# base API endpoint to use for overpass queries
overpass_endpoint = "https://overpass-api.de/api"

# whether to check the overpass server status endpoint for how long to pause
# between requests: necessary if server uses slot management, but can be set
# to False if you are running your own overpass instance without rate limiting
overpass_rate_limit = True

# which API provider to use for adding node elevations. default is "google"
# for Google Maps Elevation API but also accepts "airmap"
# this setting is deprecated and will be removed in a future release
elevation_provider = "google"

# requests session that can be overridden based on user's needs
# to include providing certificates for client/server authentication
session = requests.Session()
