"""Truncate graph by distance, bounding box, or polygon."""

import geopandas as gpd
import networkx as nx
import pandas as pd
from shapely.geometry import Point

from . import utils
from . import utils_geo
from . import utils_graph


def truncate_graph_dist(G, source_node, max_dist=1000, weight="length", retain_all=False):
    """
    Remove everything farther than some network distance from specified node.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    source_node : int
        the node in the graph from which to measure network distances to other
        nodes
    max_dist : int
        remove every node in the graph greater than this distance from the
        source_node (along the network)
    weight : string
        how to weight the graph when measuring distance (default 'length' is
        how many meters long the edge is)
    retain_all : bool
        if True, return the entire graph even if it is not connected

    Returns
    -------
    G : networkx.MultiDiGraph
        the truncated graph
    """
    # get the shortest distance between the node and every other node, then
    # remove every node further than max_dist away
    G = G.copy()
    distances = nx.shortest_path_length(G, source=source_node, weight=weight)
    distant_nodes = {key: value for key, value in dict(distances).items() if value > max_dist}
    G.remove_nodes_from(distant_nodes.keys())
    utils.log("Truncated graph by weighted network distance")

    # remove any isolated nodes and retain only the largest component (if
    # retain_all is True)
    if not retain_all:
        G = utils_graph.remove_isolated_nodes(G)
        G = utils_graph.get_largest_component(G)

    return G


def truncate_graph_bbox(G, north, south, east, west, truncate_by_edge=False, retain_all=False):
    """
    Remove every node in graph that falls outside a bounding box.

    Needed because overpass returns entire ways that also include nodes outside
    the bbox if the way (that is, a way with a single OSM ID) has a node inside
    the bbox at some point.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    north : float
        northern latitude of bounding box
    south : float
        southern latitude of bounding box
    east : float
        eastern longitude of bounding box
    west : float
        western longitude of bounding box
    truncate_by_edge : bool
        if True retain node if it's outside bbox but at least one of node's
        neighbors are within bbox
    retain_all : bool
        if True, return the entire graph even if it is not connected

    Returns
    -------
    G : networkx.MultiDiGraph
        the truncated graph
    """
    G = G.copy()
    nodes_outside_bbox = []

    for node, data in G.nodes(data=True):
        if data["y"] > north or data["y"] < south or data["x"] > east or data["x"] < west:
            # this node is outside the bounding box
            if not truncate_by_edge:
                # if we're not truncating by edge, add node to list of nodes
                # outside the bounding box
                nodes_outside_bbox.append(node)
            else:
                # if we're truncating by edge, see if any of node's neighbors
                # are within bounding box
                any_neighbors_in_bbox = False
                neighbors = list(G.successors(node)) + list(G.predecessors(node))
                for neighbor in neighbors:
                    x = G.nodes[neighbor]["x"]
                    y = G.nodes[neighbor]["y"]
                    if y < north and y > south and x < east and x > west:
                        any_neighbors_in_bbox = True
                        break

                # if none of its neighbors are within the bounding box, add node
                # to list of nodes outside the bounding box
                if not any_neighbors_in_bbox:
                    nodes_outside_bbox.append(node)

    G.remove_nodes_from(nodes_outside_bbox)
    utils.log("Truncated graph by bounding box")

    # remove any isolated nodes and retain only the largest component (if
    # retain_all is True)
    if not retain_all:
        G = utils_graph.remove_isolated_nodes(G)
        G = utils_graph.get_largest_component(G)

    return G


def truncate_graph_polygon(
    G, polygon, retain_all=False, truncate_by_edge=False, quadrat_width=0.05, min_num=3
):
    """
    Remove every node in graph that outside a shapely (Multi)Polygon.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        only retain nodes in graph that lie within this geometry
    retain_all : bool
        if True, return the entire graph even if it is not connected
    truncate_by_edge : bool
        if True retain node if it's outside polygon but at least one of node's
        neighbors are within polygon
    quadrat_width : numeric
        passed on to intersect_index_quadrats: the linear length (in degrees) of
        the quadrats with which to cut up the geometry (default = 0.05, approx
        4km at NYC's latitude)
    min_num : int
        passed on to intersect_index_quadrats: the minimum number of linear
        quadrat lines (e.g., min_num=3 would produce a quadrat grid of 4
        squares)

    Returns
    -------
    G : networkx.MultiDiGraph
        the truncated graph
    """
    G = G.copy()
    utils.log("Identifying all nodes that lie outside the polygon...")

    # get a GeoDataFrame of all the nodes
    node_geom = [Point(data["x"], data["y"]) for _, data in G.nodes(data=True)]
    gdf_nodes = gpd.GeoDataFrame({"node": list(G.nodes()), "geometry": node_geom})
    gdf_nodes.crs = G.graph["crs"]

    # find all the nodes in the graph that lie outside the polygon
    points_within_geometry = utils_geo._intersect_index_quadrats(
        gdf_nodes, polygon, quadrat_width=quadrat_width, min_num=min_num
    )
    nodes_outside_polygon = gdf_nodes[~gdf_nodes.index.isin(points_within_geometry.index)]

    if truncate_by_edge:
        nodes_to_remove = []
        for node in nodes_outside_polygon["node"]:
            neighbors = pd.Series(list(G.successors(node)) + list(G.predecessors(node)))
            # check if all the neighbors of this node also lie outside polygon
            if neighbors.isin(nodes_outside_polygon["node"]).all():
                nodes_to_remove.append(node)
    else:
        nodes_to_remove = nodes_outside_polygon["node"]

    # now remove from the graph all those nodes that lie outside the place
    # polygon
    G.remove_nodes_from(nodes_to_remove)
    utils.log(f"Removed {len(nodes_outside_polygon)} nodes outside polygon")

    # remove any isolated nodes and retain only the largest component (if retain_all is False)
    if not retain_all:
        G = utils_graph.remove_isolated_nodes(G)
        G = utils_graph.get_largest_component(G)

    return G
