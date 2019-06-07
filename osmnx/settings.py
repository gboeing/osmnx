################################################################################
# Module: settings.py
# Description: Global settings, can be configured by user by passing values to
#              utils.config()
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

import logging as lg


# locations to save data, logs, images, and cache
data_folder = 'data'
logs_folder = 'logs'
imgs_folder = 'images'
cache_folder = 'cache'

# cache server responses
use_cache = False

# write log to file and/or to console
log_file = False
log_console = False
log_level = lg.INFO
log_name = 'osmnx'
log_filename = 'osmnx'

# useful osm tags - note that load_graphml expects a consistent set of tag names
# for parsing
useful_tags_node = ['ref', 'highway']
useful_tags_path = ['bridge', 'tunnel', 'oneway', 'lanes', 'ref', 'name',
                    'highway', 'maxspeed', 'service', 'access', 'area',
                    'landuse', 'width', 'est_width', 'junction']

# tags and attributes for generating an OSM XML file
osm_xml_node_attrs = [
    'id', 'timestamp', 'uid', 'user', 'version', 'changeset', 'lat', 'lon']
osm_xml_node_tags = ['highway']
osm_xml_way_attrs = ['id', 'timestamp', 'uid', 'user', 'version', 'changeset']
osm_xml_way_tags = ['highway', 'lanes', 'maxspeed', 'name', 'oneway']


# default filter for OSM "access" key. filtering out "access=no" ways prevents
# including transit-only bridges like tilikum crossing from appearing in drivable
# road network (e.g., '["access"!~"private|no"]'). however, some drivable
# tollroads have "access=no" plus a "access:conditional" key to clarify when
# it is accessible, so we can't filter out all "access=no" ways by default.
# best to be permissive here then remove complicated combinations of tags in
# python after the full graph is downloaded and constructed.
default_access = '["access"!~"private"]'

# The network types for which a bidirectional graph will be created
bidirectional_network_types = ['walk']

# default CRS to set when creating graphs
default_crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

# default HTTP request headers
default_user_agent = 'Python OSMnx package (https://github.com/gboeing/osmnx)'
default_referer = 'Python OSMnx package (https://github.com/gboeing/osmnx)'
default_accept_language = 'en'

# which API endpoint to use for nominatim queries
# and your API key, if you are using a commercial endpoint that requires it
nominatim_endpoint = "https://nominatim.openstreetmap.org/"
nominatim_key = None
