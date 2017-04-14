###################################################################################################
# Module: globals.py
# Description: Global defaults, can be configured by user by passing values to utils.config()
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
###################################################################################################

import logging as lg


# default locations to save data, logs, images, and cache
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

# useful osm tags - note that load_graphml expects a consistent set of tag names for parsing
useful_tags_node = ['ref', 'highway']
useful_tags_path = ['bridge', 'tunnel', 'oneway', 'lanes', 'ref', 'name', 'highway', 'maxspeed', 'service', 'access', 'area', 'landuse', 'width', 'est_width']
