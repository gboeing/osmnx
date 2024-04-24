"""
Calculate geometric and topological network measures.

This module defines streets as the edges in an undirected representation of
the graph. Using undirected graph edges prevents double-counting bidirectional
edges of a two-way street, but may double-count a divided road's separate
centerlines with different end point nodes. Due to OSMnx's periphery cleaning
when the graph was created, you will get accurate node degrees (and in turn
streets-per-node counts) even at the periphery of the graph.

You can use NetworkX directly for additional topological network measures.
"""

from __future__ import annotations

import logging as lg
from collections import Counter
from itertools import chain
from typing import TYPE_CHECKING
from typing import Any

import networkx as nx
import numpy as np

from . import convert
from . import distance
from . import projection
from . import simplification
from . import utils

if TYPE_CHECKING:
    from collections.abc import Iterable


def streets_per_node(G: nx.MultiDiGraph) -> dict[int, int]:
    """
    Retrieve nodes' `street_count` attribute values.

    See also the `count_streets_per_node` function for the calculation.

    Parameters
    ----------
    G
        Input graph.

    Returns
    -------
    spn
        Dictionary with node ID keys and street count values.
    """
    # ensure each count value has type int (otherwise could be type np.int64)
    # if user has projected the graph bc GeoDataFrames use np.int64 for ints
    spn = {k: int(v) for k, v in nx.get_node_attributes(G, "street_count").items()}
    if set(spn) != set(G.nodes):
        msg = "Graph nodes changed since `street_count`s were calculated"
        utils.log(msg, level=lg.WARNING)
    return spn


def streets_per_node_avg(G: nx.MultiDiGraph) -> float:
    """
    Calculate graph's average count of streets per node.

    Parameters
    ----------
    G
        Input graph.

    Returns
    -------
    spna
        Average count of streets per node.
    """
    spn_vals = streets_per_node(G).values()
    return float(sum(spn_vals) / len(G.nodes))


def streets_per_node_counts(G: nx.MultiDiGraph) -> dict[int, int]:
    """
    Calculate streets-per-node counts.

    Parameters
    ----------
    G
        Input graph.

    Returns
    -------
    spnc
        Dictionary keyed by count of streets incident on each node, and with
        values of how many nodes in the graph have this count.
    """
    spn_vals = list(streets_per_node(G).values())
    return {i: spn_vals.count(i) for i in range(int(max(spn_vals)) + 1)}


def streets_per_node_proportions(G: nx.MultiDiGraph) -> dict[int, float]:
    """
    Calculate streets-per-node proportions.

    Parameters
    ----------
    G
        Input graph.

    Returns
    -------
    spnp
        Dictionary keyed by count of streets incident on each node, and with
        values of what proportion of nodes in the graph have this count.
    """
    n = len(G.nodes)
    spnc = streets_per_node_counts(G)
    return {i: count / n for i, count in spnc.items()}


def intersection_count(G: nx.MultiDiGraph, *, min_streets: int = 2) -> int:
    """
    Count the intersections in a graph.

    Intersections are defined as nodes with at least `min_streets` number of
    streets incident on them.

    Parameters
    ----------
    G
        Input graph.
    min_streets
        A node must have at least `min_streets` incident on them to count as
        an intersection.

    Returns
    -------
    count
        Count of intersections in graph.
    """
    spn = streets_per_node(G)
    node_ids = set(G.nodes)
    count = sum(c >= min_streets and n in node_ids for n, c in spn.items())

    # ensure count value has type int (otherwise could be type np.int64) if
    # user has projected the graph bc GeoDataFrames use np.int64 for ints
    return int(count)


def street_segment_count(Gu: nx.MultiGraph) -> int:
    """
    Count the street segments in a graph.

    Parameters
    ----------
    Gu
        Undirected input graph.

    Returns
    -------
    count
        Count of street segments in graph.
    """
    if nx.is_directed(Gu):  # pragma: no cover
        msg = "`Gu` must be undirected."
        raise ValueError(msg)
    return len(Gu.edges)


def street_length_total(Gu: nx.MultiGraph) -> float:
    """
    Calculate graph's total street segment length.

    Parameters
    ----------
    Gu
        Undirected input graph.

    Returns
    -------
    length
        Total length (meters) of streets in graph.
    """
    if nx.is_directed(Gu):  # pragma: no cover
        msg = "`Gu` must be undirected."
        raise ValueError(msg)
    return float(sum(d["length"] for u, v, d in Gu.edges(data=True)))


def edge_length_total(G: nx.MultiGraph) -> float:
    """
    Calculate graph's total edge length.

    Parameters
    ----------
    G
        Input graph.

    Returns
    -------
    length
        Total length (meters) of edges in graph.
    """
    return float(sum(d["length"] for u, v, d in G.edges(data=True)))


def self_loop_proportion(Gu: nx.MultiGraph) -> float:
    """
    Calculate percent of edges that are self-loops in a graph.

    A self-loop is defined as an edge from node `u` to node `v` where `u==v`.

    Parameters
    ----------
    Gu
        Undirected input graph.

    Returns
    -------
    proportion
        Proportion of graph edges that are self-loops.
    """
    if nx.is_directed(Gu):  # pragma: no cover
        msg = "`Gu` must be undirected."
        raise ValueError(msg)
    return float(sum(u == v for u, v, k in Gu.edges) / len(Gu.edges))


def circuity_avg(Gu: nx.MultiGraph) -> float | None:
    """
    Calculate average street circuity using edges of undirected graph.

    Circuity is the sum of edge lengths divided by the sum of straight-line
    distances between edge endpoints. Calculates straight-line distance as
    euclidean distance if projected or great-circle distance if unprojected.
    Returns None if the edge lengths sum to zero.

    Parameters
    ----------
    Gu
        Undirected input graph.

    Returns
    -------
    circuity_avg
        The graph's average undirected edge circuity.
    """
    if nx.is_directed(Gu):  # pragma: no cover
        msg = "`Gu` must be undirected."
        raise ValueError(msg)

    # extract the edges' endpoint nodes' coordinates
    n = Gu.nodes
    coords = np.array([(n[u]["y"], n[u]["x"], n[v]["y"], n[v]["x"]) for u, v, _ in Gu.edges])
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
        return float(edge_length_total(Gu) / sl_dists_total)
    except ZeroDivisionError:
        return None


def count_streets_per_node(
    G: nx.MultiDiGraph,
    *,
    nodes: Iterable[int] | None = None,
) -> dict[int, int]:
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
    G
        Input graph.
    nodes
        Which node IDs to get counts for. If None, use all graph nodes.
        Otherwise calculate counts only for these node IDs.

    Returns
    -------
    streets_per_node
        Counts of how many physical streets connect to each node, with keys =
        node ids and values = counts.
    """
    if nodes is None:
        nodes = G.nodes

    # get one copy of each self-loop edge, because bi-directional self-loops
    # appear twice in the undirected graph (u,v,0 and u,v,1 where u=v), but
    # one-way self-loops will appear only once
    Gu = G.to_undirected(reciprocal=False, as_view=True)
    self_loop_edges = set(nx.selfloop_edges(Gu, keys=False))

    # get all non-self-loop undirected edges, including parallel edges
    non_self_loop_edges = [e for e in Gu.edges(keys=False) if e not in self_loop_edges]

    # make list of all unique edges including each parallel edge unless the
    # parallel edge is a self-loop, in which case we don't double-count it
    all_unique_edges = non_self_loop_edges + list(self_loop_edges)

    # flatten list of (u, v) edge tuples to count how often each node appears
    edges_flat = chain.from_iterable(all_unique_edges)
    counts = Counter(edges_flat)
    streets_per_node = {node: counts[node] for node in nodes}

    msg = "Counted undirected street segments incident on each node"
    utils.log(msg, level=lg.INFO)
    return streets_per_node


def basic_stats(
    G: nx.MultiDiGraph,
    *,
    area: float | None = None,
    clean_int_tol: float | None = None,
) -> dict[str, Any]:
    """
    Calculate basic descriptive geometric and topological measures of a graph.

    Density measures are only calculated if `area` is provided and clean
    intersection measures are only calculated if `clean_int_tol` is provided.

    Parameters
    ----------
    G
        Input graph.
    area
        If not None, calculate density measures and use `area` (in square
        meters) as the denominator.
    clean_int_tol
        If not None, calculate consolidated intersections count (and density,
        if `area` is also provided) and use this tolerance value. Refer to the
        `simplification.consolidate_intersections` function documentation for
        details.

    Returns
    -------
    stats
        Dictionary containing the following keys:
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
    stats: dict[str, Any] = {}

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
                G,
                tolerance=clean_int_tol,
                rebuild_graph=False,
                dead_ends=False,
            ),
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
