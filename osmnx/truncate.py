"""Truncate graph by distance, bounding box, or polygon."""

from warnings import warn

import networkx as nx

from . import convert
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
        node in graph from which to measure network distances to other nodes
    max_dist : float
        remove every node in the graph that is greater than this distance (in
        same units as `weight` attribute) along the network from `source_node`
    weight : string
        graph edge attribute to use to measure distance
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
    distant_nodes = {k for k, v in distances.items() if v > max_dist}
    unreachable_nodes = G.nodes - distances.keys()

    # make a copy to not mutate original graph object caller passed in
    G = G.copy()
    G.remove_nodes_from(distant_nodes | unreachable_nodes)

    # remove any isolated nodes and retain only the largest component (if
    # retain_all is True)
    if not retain_all:
        G = utils_graph.remove_isolated_nodes(G, warn=False)
        G = largest_component(G)

    utils.log(f"Truncated graph by {weight}-weighted network distance")
    return G


def truncate_graph_bbox(
    G,
    north=None,
    south=None,
    east=None,
    west=None,
    bbox=None,
    truncate_by_edge=False,
    retain_all=False,
    quadrat_width=None,
    min_num=None,
):
    """
    Remove every node in graph that falls outside a bounding box.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    north : float
        deprecated, do not use
    south : float
        deprecated, do not use
    east : float
        deprecated, do not use
    west : float
        deprecated, do not use
    bbox : tuple of floats
        bounding box as (north, south, east, west)
    truncate_by_edge : bool
        if True, retain nodes outside bounding box if at least one of node's
        neighbors is within the bounding box
    retain_all : bool
        if True, return the entire graph even if it is not connected.
        otherwise, retain only the largest weakly connected component.
    quadrat_width : float
        deprecated, do not use
    min_num : int
        deprecated, do not use

    Returns
    -------
    G : networkx.MultiDiGraph
        the truncated graph
    """
    if not (north is None and south is None and east is None and west is None):
        msg = (
            "The `north`, `south`, `east`, and `west` parameters are deprecated and "
            "will be removed in the v2.0.0 release. Use the `bbox` parameter instead. "
            "Note that the expected order of coordinates in `bbox` will change in the "
            "v2.0.0 release to `(left, bottom, right, top)`. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
        )
        warn(msg, FutureWarning, stacklevel=2)
        bbox = (north, south, east, west)

    # convert bounding box to a polygon, then truncate
    polygon = utils_geo.bbox_to_poly(bbox=bbox)
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
    G, polygon, retain_all=False, truncate_by_edge=False, quadrat_width=None, min_num=None
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
    quadrat_width : float
        deprecated, do not use
    min_num : int
        deprecated, do not use

    Returns
    -------
    G : networkx.MultiDiGraph
        the truncated graph
    """
    if quadrat_width is not None or min_num is not None:
        warn(
            "The `quadrat_width` and `min_num` parameters are deprecated and "
            "will be removed in the v2.0.0 release. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123",
            FutureWarning,
            stacklevel=2,
        )

    utils.log("Identifying all nodes that lie outside the polygon...")

    # first identify all nodes whose point geometries lie within the polygon
    gs_nodes = convert.graph_to_gdfs(G, edges=False)[["geometry"]]
    to_keep = utils_geo._intersect_index_quadrats(gs_nodes, polygon)

    if not to_keep:
        # no graph nodes within the polygon: can't create a graph from that
        msg = "Found no graph nodes within the requested polygon"
        raise ValueError(msg)

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
    # make a copy to not mutate original graph object caller passed in
    G = G.copy()
    G.remove_nodes_from(nodes_to_remove)
    utils.log(f"Removed {len(nodes_to_remove):,} nodes outside polygon")

    if not retain_all:
        # remove any isolated nodes and retain only the largest component
        G = utils_graph.remove_isolated_nodes(G, warn=False)
        G = largest_component(G)

    utils.log("Truncated graph by polygon")
    return G


def largest_component(G, strongly=False):
    """
    Get subgraph of G's largest weakly/strongly connected component.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    strongly : bool
        if True, return the largest strongly instead of weakly connected
        component

    Returns
    -------
    G : networkx.MultiDiGraph
        the largest connected component subgraph of the original graph
    """
    if strongly:
        kind = "strongly"
        is_connected = nx.is_strongly_connected
        connected_components = nx.strongly_connected_components
    else:
        kind = "weakly"
        is_connected = nx.is_weakly_connected
        connected_components = nx.weakly_connected_components

    if not is_connected(G):
        # get all the connected components in graph then identify the largest
        largest_cc = max(connected_components(G), key=len)
        n = len(G)

        # induce (frozen) subgraph then unfreeze it by making new MultiDiGraph
        G = nx.MultiDiGraph(G.subgraph(largest_cc))
        utils.log(f"Got largest {kind} connected component ({len(G):,} of {n:,} total nodes)")

    return G
