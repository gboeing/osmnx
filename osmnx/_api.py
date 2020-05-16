################################################################################
# Module: api.py
# Description: Expose the OSMnx API
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

from .core import gdf_from_place
from .core import gdf_from_places
from .core import graph_from_address
from .core import graph_from_bbox
from .core import graph_from_file
from .core import graph_from_place
from .core import graph_from_point
from .core import graph_from_polygon

from .elevation import add_edge_grades
from .elevation import add_node_elevations

from .footprints import footprints_from_address
from .footprints import footprints_from_place
from .footprints import footprints_from_point
from .footprints import footprints_from_polygon
from .footprints import plot_footprints

from .plot import plot_figure_ground
from .plot import plot_graph
from .plot import plot_graph_folium
from .plot import plot_graph_route
from .plot import plot_graph_routes
from .plot import plot_route_folium
from .plot import plot_shape

from .pois import pois_from_address
from .pois import pois_from_place
from .pois import pois_from_point
from .pois import pois_from_polygon

from .projection import project_graph

from .save_load import load_graphml
from .save_load import save_graph_geopackage
from .save_load import save_graph_osm
from .save_load import save_graph_shapefile
from .save_load import save_graphml

from .simplification import clean_intersections
from .simplification import consolidate_intersections
from .simplification import simplify_graph

from .stats import basic_stats
from .stats import extended_stats

from .utils import citation
from .utils import config
from .utils import log
from .utils import ts

from .utils_geo import add_edge_bearings
from .utils_geo import geocode
from .utils_geo import get_nearest_edge
from .utils_geo import get_nearest_edges
from .utils_geo import get_nearest_node
from .utils_geo import get_nearest_nodes

from .utils_graph import gdfs_to_graph
from .utils_graph import get_undirected
from .utils_graph import graph_to_gdfs
