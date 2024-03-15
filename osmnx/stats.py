"""
Calculate geometric and topological network measures.

This module defines streets as the edges in an undirected representation of
the graph. Using undirected graph edges prevents double-counting bidirectional
edges of a two-way street, but may double-count a divided road's separate
centerlines with different end point nodes. If `clean_periphery=True` when the
graph was created (which is the default parameterization), then you will get
accurate node degrees (and in turn streets-per-node counts) even at the
periphery of the graph.

You can use NetworkX directly for additional topological network measures.
"""

import itertools
import logging as lg
from collections import Counter

import networkx as nx
import numpy as np

from . import convert
from . import distance
from . import projection
from . import simplification
from . import utils


def streets_per_node(G):
    """
    Count streets (undirected edges) incident on each node.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph

    Returns
    -------
    spn : dict
        dictionary with node ID keys and street count values
    """
    spn = dict(nx.get_node_attributes(G, "street_count"))
    if set(spn) != set(G.nodes):
        utils.log("Graph nodes changed since `street_count`s were calculated", level=lg.WARNING)
    return spn


def streets_per_node_avg(G):
    """
    Calculate graph's average count of streets per node.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph

    Returns
    -------
    spna : float
        average count of streets per node
    """
    spn_vals = streets_per_node(G).values()
    return sum(spn_vals) / len(G.nodes)


def streets_per_node_counts(G):
    """
    Calculate streets-per-node counts.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph

    Returns
    -------
    spnc : dict
        dictionary keyed by count of streets incident on each node, and with
        values of how many nodes in the graph have this count
    """
    spn_vals = list(streets_per_node(G).values())
    return {i: spn_vals.count(i) for i in range(int(max(spn_vals)) + 1)}


def streets_per_node_proportions(G):
    """
    Calculate streets-per-node proportions.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph

    Returns
    -------
    spnp : dict
        dictionary keyed by count of streets incident on each node, and with
        values of what proportion of nodes in the graph have this count
    """
    n = len(G.nodes)
    spnc = streets_per_node_counts(G)
    return {i: count / n for i, count in spnc.items()}


def intersection_count(G=None, min_streets=2):
    """
    Count the intersections in a graph.

    Intersections are defined as nodes with at least `min_streets` number of
    streets incident on them.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    min_streets : int
        a node must have at least `min_streets` incident on them to count as
        an intersection

    Returns
    -------
    count : int
        count of intersections in graph
    """
    spn = streets_per_node(G)
    node_ids = set(G.nodes)
    return sum(count >= min_streets and node in node_ids for node, count in spn.items())


def street_segment_count(Gu):
    """
    Count the street segments in a graph.

    Parameters
    ----------
    Gu : networkx.MultiGraph
        undirected input graph

    Returns
    -------
    count : int
        count of street segments in graph
    """
    if nx.is_directed(Gu):  # pragma: no cover
        msg = "`Gu` must be undirected"
        raise ValueError(msg)
    return len(Gu.edges)


def street_length_total(Gu):
    """
    Calculate graph's total street segment length.

    Parameters
    ----------
    Gu : networkx.MultiGraph
        undirected input graph

    Returns
    -------
    length : float
        total length (meters) of streets in graph
    """
    if nx.is_directed(Gu):  # pragma: no cover
        msg = "`Gu` must be undirected"
        raise ValueError(msg)
    return sum(d["length"] for u, v, d in Gu.edges(data=True))


def edge_length_total(G):
    """
    Calculate graph's total edge length.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph

    Returns
    -------
    length : float
        total length (meters) of edges in graph
    """
    return sum(d["length"] for u, v, d in G.edges(data=True))


def self_loop_proportion(Gu):
    """
    Calculate percent of edges that are self-loops in a graph.

    A self-loop is defined as an edge from node `u` to node `v` where `u==v`.

    Parameters
    ----------
    Gu : networkx.MultiGraph
        undirected input graph

    Returns
    -------
    proportion : float
        proportion of graph edges that are self-loops
    """
    if nx.is_directed(Gu):  # pragma: no cover
        msg = "`Gu` must be undirected"
        raise ValueError(msg)
    return sum(u == v for u, v, k in Gu.edges) / len(Gu.edges)


def circuity_avg(Gu):
    """
    Calculate average street circuity using edges of undirected graph.

    Circuity is the sum of edge lengths divided by the sum of straight-line
    distances between edge endpoints. Calculates straight-line distance as
    euclidean distance if projected or great-circle distance if unprojected.

    Parameters
    ----------
    Gu : networkx.MultiGraph
        undirected input graph

    Returns
    -------
    circuity_avg : float
        the graph's average undirected edge circuity
    """
    if nx.is_directed(Gu):  # pragma: no cover
        msg = "`Gu` must be undirected"
        raise ValueError(msg)

    # extract the edges' endpoint nodes' coordinates
    coords = np.array(
        [
            (Gu.nodes[u]["y"], Gu.nodes[u]["x"], Gu.nodes[v]["y"], Gu.nodes[v]["x"])
            for u, v, _ in Gu.edges
        ]
    )
    y1 = coords[:, 0]
    x1 = coords[:, 1]
    y2 = coords[:, 2]
    x2 = coords[:, 3]

    # calculate straight-line distances as euclidean distances if projected or
    # great-circle distances if unprojected
    if projection.is_projected(Gu.graph["crs"]):
        sl_dists = distance.euclidean(y1=y1, x1=x1, y2=y2, x2=x2)
    else:
        sl_dists = distance.great_circle(lat1=y1, lon1=x1, lat2=y2, lon2=x2)

    # return the ratio, handling possible division by zero
    sl_dists_total = sl_dists[~np.isnan(sl_dists)].sum()
    try:
        return edge_length_total(Gu) / sl_dists_total
    except ZeroDivisionError:
        return None


def count_streets_per_node(G, nodes=None):
    """
    Count how many physical street segments connect to each node in a graph.

    This function uses an undirected representation of the graph and special
    handling of self-loops to accurately count physical streets rather than
    directed edges. Note: this function is automatically run by all the
    `graph.graph_from_x` functions prior to truncating the graph to the
    requested boundaries, to add accurate `street_count` attributes to each
    node even if some of its neighbors are outside the requested graph
    boundaries.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    nodes : list
        which node IDs to get counts for. if None, use all graph nodes,
        otherwise calculate counts only for these node IDs

    Returns
    -------
    streets_per_node : dict
        counts of how many physical streets connect to each node, with keys =
        node ids and values = counts
    """
    if nodes is None:
        nodes = G.nodes

    # get one copy of each self-loop edge, because bi-directional self-loops
    # appear twice in the undirected graph (u,v,0 and u,v,1 where u=v), but
    # one-way self-loops will appear only once
    Gu = G.to_undirected(reciprocal=False, as_view=True)
    self_loop_edges = set(nx.selfloop_edges(Gu))

    # get all non-self-loop undirected edges, including parallel edges
    non_self_loop_edges = [e for e in Gu.edges(keys=False) if e not in self_loop_edges]

    # make list of all unique edges including each parallel edge unless the
    # parallel edge is a self-loop, in which case we don't double-count it
    all_unique_edges = non_self_loop_edges + list(self_loop_edges)

    # flatten list of (u, v) edge tuples to count how often each node appears
    edges_flat = itertools.chain.from_iterable(all_unique_edges)
    counts = Counter(edges_flat)
    streets_per_node = {node: counts[node] for node in nodes}

    utils.log("Counted undirected street segments incident on each node")
    return streets_per_node


def basic_stats(G, area=None, clean_int_tol=None):
    """
    Calculate basic descriptive geometric and topological measures of a graph.

    Density measures are only calculated if `area` is provided and clean
    intersection measures are only calculated if `clean_int_tol` is provided.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    area : float
        if not None, calculate density measures and use this value (in square
        meters) as the denominator
    clean_int_tol : float
        if not None, calculate consolidated intersections count (and density,
        if `area` is also provided) and use this tolerance value; refer to the
        `simplification.consolidate_intersections` function documentation for
        details

    Returns
    -------
    stats : dict
        dictionary containing the following keys
          - `circuity_avg` - see `circuity_avg` function documentation
          - `clean_intersection_count` - see `clean_intersection_count` function documentation
          - `clean_intersection_density_km` - `clean_intersection_count` per sq km
          - `edge_density_km` - `edge_length_total` per sq km
          - `edge_length_avg` - `edge_length_total / m`
          - `edge_length_total` - see `edge_length_total` function documentation
          - `intersection_count` - see `intersection_count` function documentation
          - `intersection_density_km` - `intersection_count` per sq km
          - `k_avg` - graph's average node degree (in-degree and out-degree)
          - `m` - count of edges in graph
          - `n` - count of nodes in graph
          - `node_density_km` - `n` per sq km
          - `self_loop_proportion` - see `self_loop_proportion` function documentation
          - `street_density_km` - `street_length_total` per sq km
          - `street_length_avg` - `street_length_total / street_segment_count`
          - `street_length_total` - see `street_length_total` function documentation
          - `street_segment_count` - see `street_segment_count` function documentation
          - `streets_per_node_avg` - see `streets_per_node_avg` function documentation
          - `streets_per_node_counts` - see `streets_per_node_counts` function documentation
          - `streets_per_node_proportions` - see `streets_per_node_proportions` function documentation
    """
    Gu = convert.to_undirected(G)
    stats = {}

    stats["n"] = len(G.nodes)
    stats["m"] = len(G.edges)
    stats["k_avg"] = 2 * stats["m"] / stats["n"]
    stats["edge_length_total"] = edge_length_total(G)
    stats["edge_length_avg"] = stats["edge_length_total"] / stats["m"]
    stats["streets_per_node_avg"] = streets_per_node_avg(G)
    stats["streets_per_node_counts"] = streets_per_node_counts(G)
    stats["streets_per_node_proportions"] = streets_per_node_proportions(G)
    stats["intersection_count"] = intersection_count(G)
    stats["street_length_total"] = street_length_total(Gu)
    stats["street_segment_count"] = street_segment_count(Gu)
    stats["street_length_avg"] = stats["street_length_total"] / stats["street_segment_count"]
    stats["circuity_avg"] = circuity_avg(Gu)
    stats["self_loop_proportion"] = self_loop_proportion(Gu)

    # calculate clean intersection counts if requested
    if clean_int_tol:
        stats["clean_intersection_count"] = len(
            simplification.consolidate_intersections(
                G, tolerance=clean_int_tol, rebuild_graph=False, dead_ends=False
            )
        )

    # can only calculate density measures if area was provided
    if area is not None:
        area_km = area / 1_000_000  # convert m^2 to km^2
        stats["node_density_km"] = stats["n"] / area_km
        stats["intersection_density_km"] = stats["intersection_count"] / area_km
        stats["edge_density_km"] = stats["edge_length_total"] / area_km
        stats["street_density_km"] = stats["street_length_total"] / area_km
        if clean_int_tol:
            stats["clean_intersection_density_km"] = stats["clean_intersection_count"] / area_km

    return stats
