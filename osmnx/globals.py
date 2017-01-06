import logging as lg


###################################################################################################
# global defaults
# you can edit any of these by passing the value to the config() function
###################################################################################################


# default locations to save data, logs, images, and cache
global_data_folder = 'data'
global_logs_folder = 'logs'
global_imgs_folder = 'images'
global_cache_folder = 'cache'

global_use_cache = False

# write log to file and/or to console
global_log_file = False
global_log_console = False
global_log_level = lg.INFO
global_log_name = 'osmnx'
global_log_filename = 'osmnx'

# useful osm tags - note that load_graphml expects a consistent set of tag names for parsing
global_useful_tags_node = ['ref', 'highway']
global_useful_tags_path = ['bridge', 'tunnel', 'oneway', 'lanes', 'ref', 'name', 'highway', 'maxspeed', 'service', 'access', 'area', 'landuse', 'width', 'est_width']



