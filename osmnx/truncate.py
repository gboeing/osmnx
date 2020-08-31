"""Truncate graph by distance, bounding box, or polygon."""

import networkx as nx

from . import utils
from . import utils_geo
from . import utils_graph


def truncate_graph_dist(G, source_node, max_dist=1000, weight="length", retain_all=False):
    """
    Remove every node farther than some network distance from source_node.

    This function can be slow for large graphs, as it must calculate shortest
    path distances between source_node and every other graph node.

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
        if True, return the entire graph even if it is not connected.
        otherwise, retain only the largest weakly connected component.

    Returns
    -------
    G : networkx.MultiDiGraph
        the truncated graph
    """
    # get the shortest distance between the node and every other node
    distances = nx.shortest_path_length(G, source=source_node, weight=weight)

    # then identify every node further than max_dist away
    distant_nodes = {k: v for k, v in distances.items() if v > max_dist}

    # make a copy to not edit the original graph object the caller passed in
    G = G.copy()
    G.remove_nodes_from(distant_nodes)

    # remove any isolated nodes and retain only the largest component (if
    # retain_all is True)
    if not retain_all:
        G = utils_graph.remove_isolated_nodes(G)
        G = utils_graph.get_largest_component(G)

    utils.log(f"Truncated graph by {weight}-weighted network distance")
    return G


def truncate_graph_bbox(
    G,
    north,
    south,
    east,
    west,
    truncate_by_edge=False,
    retain_all=False,
    quadrat_width=0.05,
    min_num=3,
):
    """
    Remove every node in graph that falls outside a bounding box.

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
        if True, retain nodes outside bounding box if at least one of node's
        neighbors is within the bounding box
    retain_all : bool
        if True, return the entire graph even if it is not connected.
        otherwise, retain only the largest weakly connected component.
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
    # convert bounding box to a polygon, then truncate
    polygon = utils_geo.bbox_to_poly(north, south, east, west)
    G = truncate_graph_polygon(
        G,
        polygon,
        retain_all=retain_all,
        truncate_by_edge=truncate_by_edge,
        quadrat_width=quadrat_width,
        min_num=min_num,
    )

    utils.log("Truncated graph by bounding box")
    return G


def truncate_graph_polygon(
    G, polygon, retain_all=False, truncate_by_edge=False, quadrat_width=0.05, min_num=3
):
    """
    Remove every node in graph that falls outside a (Multi)Polygon.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        only retain nodes in graph that lie within this geometry
    retain_all : bool
        if True, return the entire graph even if it is not connected.
        otherwise, retain only the largest weakly connected component.
    truncate_by_edge : bool
        if True, retain nodes outside boundary polygon if at least one of
        node's neighbors is within the polygon
    quadrat_width : numeric
        passed on to intersect_index_quadrats: the linear length (in degrees)
        of the quadrats with which to cut up the geometry (default = 0.05,
        approx 4km at NYC's latitude)
    min_num : int
        passed on to intersect_index_quadrats: the minimum number of linear
        quadrat lines (e.g., min_num=3 would produce a quadrat grid of 4
        squares)

    Returns
    -------
    G : networkx.MultiDiGraph
        the truncated graph
    """
    utils.log("Identifying all nodes that lie outside the polygon...")

    # first identify all nodes whose point geometries lie within the polygon
    gs_nodes = utils_graph.graph_to_gdfs(G, edges=False)[["geometry"]]
    to_keep = utils_geo._intersect_index_quadrats(gs_nodes, polygon, quadrat_width, min_num)

    if len(to_keep) == 0:
        # no graph nodes within the polygon: can't create a graph from that
        raise ValueError("Found no graph nodes within the requested polygon")

    # now identify all nodes whose point geometries lie outside the polygon
    gs_nodes_outside_poly = gs_nodes[~gs_nodes.index.isin(to_keep)]
    nodes_outside_poly = set(gs_nodes_outside_poly.index)

    if truncate_by_edge:
        # retain nodes outside boundary polygon if at least one of node's
        # neighbors is within the polygon
        nodes_to_remove = set()
        for node in nodes_outside_poly:
            # if all the neighbors of this node also lie outside polygon, then
            # mark this node for removal
            neighbors = set(G.successors(node)) | set(G.predecessors(node))
            if neighbors.issubset(nodes_outside_poly):
                nodes_to_remove.add(node)
    else:
        nodes_to_remove = nodes_outside_poly

    # now remove from the graph all those nodes that lie outside the polygon
    # make a copy to not edit the original graph object the caller passed in
    G = G.copy()
    G.remove_nodes_from(nodes_to_remove)
    utils.log(f"Removed {len(nodes_to_remove)} nodes outside polygon")

    if not retain_all:
        # remove any isolated nodes and retain only the largest component
        G = utils_graph.remove_isolated_nodes(G)
        G = utils_graph.get_largest_component(G)

    utils.log("Truncated graph by polygon")
    return G
