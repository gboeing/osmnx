# ruff: noqa: PLC0414
"""
Expose the old v1 API for backwards compatibility.

This allows common functionality to be accessed directly via the
ox.function_name() shortcut by exposing these functions directly in the
package's namespace.
"""

from .bearing import add_edge_bearings as add_edge_bearings
from .bearing import orientation_entropy as orientation_entropy
from .convert import graph_from_gdfs as graph_from_gdfs
from .convert import graph_to_gdfs as graph_to_gdfs
from .distance import nearest_edges as nearest_edges
from .distance import nearest_nodes as nearest_nodes
from .elevation import add_edge_grades as add_edge_grades
from .elevation import add_node_elevations_google as add_node_elevations_google
from .elevation import add_node_elevations_raster as add_node_elevations_raster
from .features import features_from_address as features_from_address
from .features import features_from_bbox as features_from_bbox
from .features import features_from_place as features_from_place
from .features import features_from_point as features_from_point
from .features import features_from_polygon as features_from_polygon
from .features import features_from_xml as features_from_xml
from .geocoder import geocode as geocode
from .geocoder import geocode_to_gdf as geocode_to_gdf
from .graph import graph_from_address as graph_from_address
from .graph import graph_from_bbox as graph_from_bbox
from .graph import graph_from_place as graph_from_place
from .graph import graph_from_point as graph_from_point
from .graph import graph_from_polygon as graph_from_polygon
from .graph import graph_from_xml as graph_from_xml
from .io import load_graphml as load_graphml
from .io import save_graph_geopackage as save_graph_geopackage
from .io import save_graph_xml as save_graph_xml
from .io import save_graphml as save_graphml
from .plot import plot_figure_ground as plot_figure_ground
from .plot import plot_footprints as plot_footprints
from .plot import plot_graph as plot_graph
from .plot import plot_graph_route as plot_graph_route
from .plot import plot_graph_routes as plot_graph_routes
from .plot import plot_orientation as plot_orientation
from .projection import project_graph as project_graph
from .routing import add_edge_speeds as add_edge_speeds
from .routing import add_edge_travel_times as add_edge_travel_times
from .routing import k_shortest_paths as k_shortest_paths
from .routing import shortest_path as shortest_path
from .simplification import consolidate_intersections as consolidate_intersections
from .simplification import simplify_graph as simplify_graph
from .stats import basic_stats as basic_stats
from .utils import citation as citation
from .utils import log as log
from .utils import ts as ts
