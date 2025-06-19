"""Validate that graph objects satisfy OSMnx expectations."""

from __future__ import annotations

import logging as lg
from numbers import Real

import geopandas as gpd
import networkx as nx

from ._errors import GraphValidationError
from .utils import log


def _validate_nodes(G: nx.MultiDiGraph, strict: bool) -> tuple[bool, str, str]:  # noqa: FBT001
    """
    Validate that a graph's nodes satisfy OSMnx expectations.

    Parameters
    ----------
    G
        The input graph.
    strict
        If `True`, enforce "should" rules in addition to "must" rules.

    Returns
    -------
    is_valid, err_msg, warn_msg
        Whether validation passed, plus any error or warning messages.
    """
    # assume nodes are valid but try to falsify that through a series of tests
    is_valid = True
    err_msg = ""
    warn_msg = ""

    # ERR: must have at least 1 node
    if not len(G.nodes) > 0:
        err_msg += "G must have at least 1 node. "
        is_valid = False

    # otherwise, it has at least 1 node, so validate the node attributes
    else:
        # ERR: nodes must have "x" and "y" data attributes
        if not all("x" in d and "y" in d for d in dict(G.nodes(data=True)).values()):
            err_msg += "Nodes must have 'x' and 'y' data attributes. "
            is_valid = False

        # WARN: nodes' "x" and "y" data attributes should be type Real
        valid_xs = all(isinstance(x, Real) for x in nx.get_node_attributes(G, name="x").values())
        valid_ys = all(isinstance(y, Real) for y in nx.get_node_attributes(G, name="y").values())
        if not (valid_xs and valid_ys):
            warn_msg += "Node 'x' and 'y' data attributes should be numeric. "
            if strict:
                is_valid = False

        # WARN: nodes should have "street_count" data attributes
        if not all("street_count" in d for d in dict(G.nodes(data=True)).values()):
            warn_msg += "Nodes should have 'street_count' data attributes. "
            if strict:
                is_valid = False

        # WARN: nodes' "x" and "y" data attributes should be type Real
        valid_xs = all(isinstance(x, Real) for x in nx.get_node_attributes(G, name="x").values())
        valid_ys = all(isinstance(y, Real) for y in nx.get_node_attributes(G, name="y").values())
        if not (valid_xs and valid_ys):
            warn_msg += "Node 'x' and 'y' data attributes should be numeric. "
            if strict:
                is_valid = False

        # WARN: node IDs should be type int
        if not all(isinstance(n, int) for n in G.nodes):
            warn_msg += "Node IDs should be type int. "
            if strict:
                is_valid = False

    return is_valid, err_msg, warn_msg


def _validate_edges(G: nx.MultiDiGraph, strict: bool) -> tuple[bool, str, str]:  # noqa: FBT001
    """
    Validate that a graph's edges satisfy OSMnx expectations.

    Parameters
    ----------
    G
        The input graph.
    strict
        If `True`, enforce optional rules in addition to required rules. These
        optional rules primarily enforce expected attribute data types.

    Returns
    -------
    is_valid, err_msg, warn_msg
        Whether validation passed, plus any error or warning messages.
    """
    # assume edges are valid but try to falsify that through a series of tests
    is_valid = True
    err_msg = ""
    warn_msg = ""

    # ERR: must have at least 1 edge
    if not len(G.edges) > 0:
        err_msg += "G must have at least 1 edge. "
        is_valid = False

    # otherwise, it has at least 1 edge, so validate the edge attributes
    else:
        # ERR: edges must have "osmid" data attributes
        edge_osmids = nx.get_edge_attributes(G, name="osmid")
        if set(edge_osmids) != set(G.edges):
            err_msg += "Edges must have 'osmid' data attributes. "
            is_valid = False

        # WARN: edge "osmid" data attributes should be type int or list[int]
        if not all(isinstance(x, (int, list)) for x in edge_osmids.values()):
            warn_msg += "Edge 'osmid' data attributes should be type `int` or `list[int]`. "
            if strict:
                is_valid = False

        # ERR: edges must have "length" data attributes
        edge_lengths = nx.get_edge_attributes(G, name="length")
        if set(edge_lengths) != set(G.edges):
            err_msg += "Edges must have 'length' data attributes. "
            is_valid = False

        # WARN: edge "length" data attributes should be numeric
        if not all(isinstance(x, Real) for x in edge_lengths.values()):
            warn_msg += "Edge 'length' data attributes should be numeric. "
            if strict:
                is_valid = False

    return is_valid, err_msg, warn_msg


def _validate_graph_attrs(G: nx.MultiDiGraph) -> tuple[bool, str, str]:
    """
    Validate that a graph's attributes satisfy OSMnx expectations.

    Parameters
    ----------
    G
        The input graph.

    Returns
    -------
    is_valid, err_msg, warn_msg
        Whether validation passed, plus any error or warning messages.
    """
    # assume G is valid but try to falsify that through a series of tests
    is_valid = True
    err_msg = ""
    warn_msg = ""

    # ERR: must be a NetworkX MultiDiGraph
    if not isinstance(G, nx.MultiDiGraph):
        err_msg += "G must be a NetworkX MultiDiGraph. "
        is_valid = False

    # ERR: must have top-level graph, nodes, and edges attributes
    if not (hasattr(G, "graph") and hasattr(G, "nodes") and hasattr(G, "edges")):
        err_msg += "G must have top-level graph, nodes, and edges attributes. "
        is_valid = False

    # ERR: graph attr dict must have a "crs" key defining its CRS
    crs = getattr(G, "graph", {}).get("crs")
    if crs is None:
        err_msg += "G.graph must have a 'crs' data attribute. "
        is_valid = False

    # ERR: graph attr dict "crs" value must be a valid pyproj CRS
    else:
        try:
            _ = gpd.GeoSeries(crs=crs).crs
        except RuntimeError:  # RuntimeError is parent of pyproj CRSError
            err_msg += "G.graph['crs'] must be a valid CRS. "
            is_valid = False

    return is_valid, err_msg, warn_msg


def _validate_graph(G: nx.MultiDiGraph, strict: bool) -> None:  # noqa: FBT001
    """
    Validate that a graph object satisfies OSMnx expectations.

    Raises `ox._errors.GraphValidationError` if validation fails.

    Parameters
    ----------
    G
        The input graph.
    strict
        If `True`, enforce "should" rules in addition to "must" rules.
    """
    # validate graph, nodes, and edges
    is_valid_graph, err_msg_graph, warn_msg_graph = _validate_graph_attrs(G)
    is_valid_nodes, err_msg_nodes, warn_msg_nodes = _validate_nodes(G, strict)
    is_valid_edges, err_msg_edges, warn_msg_edges = _validate_edges(G, strict)

    # assemble results
    is_valid = is_valid_graph and is_valid_nodes and is_valid_edges
    err_msg = err_msg_graph + err_msg_nodes + err_msg_edges
    warn_msg = warn_msg_graph + warn_msg_nodes + warn_msg_edges

    # all done, report any errors/warnings
    if is_valid:
        msg = "Successfully validated graph."
        log(msg, level=lg.INFO)
    else:
        log_level = lg.ERROR if len(err_msg) > 0 else lg.WARNING
        log(err_msg + warn_msg, level=log_level)
        raise GraphValidationError(err_msg + warn_msg)
