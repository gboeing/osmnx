"""Validate that graphs and GeoDataFrames satisfy OSMnx expectations."""

from __future__ import annotations

import logging as lg
from numbers import Real
from warnings import warn

import geopandas as gpd
import networkx as nx
import numpy as np

from ._errors import ValidationError
from .utils import log


def _verify_numeric_edge_attribute(G: nx.MultiDiGraph, attr: str, *, strict: bool = True) -> None:
    """
    Verify attribute values are numeric and non-null across graph edges.

    Raises a ValidationError if this attribute contains non-numeric
    values, and issues a UserWarning if this attribute is missing or null on
    any edges.

    Parameters
    ----------
    G
        Input graph.
    attr
        Name of the edge attribute to verify.
    strict
        If `True`, elevate warnings to errors.
    """
    is_valid = True
    valid_msg = "Verified {attr!r} values are numeric and non-null across graph edges."
    warn_msg = ""
    err_msg = ""

    try:
        values_float = (np.array(tuple(G.edges(data=attr)))[:, 2]).astype(float)
        if np.isnan(values_float).any():
            warn_msg += f"The attribute {attr!r} is missing or null on some edges."
            if strict:
                is_valid = False
    except ValueError:
        err_msg += f"The edge attribute {attr!r} contains non-numeric values."
        is_valid = False

    _report_validation(is_valid, valid_msg, warn_msg, err_msg)


def _validate_features_gdf(gdf: gpd.GeoDataFrame) -> None:
    """
    Validate that features GeoDataFrame satisfies OSMnx expectations.

    Raises a `ValidationError` if validation fails.

    Parameters
    ----------
    gdf
        GeoDataFrame of features uniquely multi-indexed by
        `(element_type, osmid)`.
    """
    is_valid = True
    valid_msg = "Validated features GeoDataFrame."
    warn_msg = ""
    err_msg = ""

    # ensure gdf is uniquely indexed
    if not gdf.index.is_unique:
        err_msg += "`gdf` must be uniquely indexed. "
        is_valid = False

    # ensure gdf is multi-indexed with 2 levels (element_type and osmid) and
    # that the element types are all either node, way, or relation
    features_index_levels = 2
    check1 = gdf.index.nlevels == features_index_levels
    element_types = set(gdf.index.get_level_values(0))
    check2 = element_types.issubset({"node", "way", "relation"})
    if not (check1 and check2):
        err_msg += "`gdf` must be multi-indexed by `(element_type, osmid)`. "
        is_valid = False

    # ensure gdf has an active geometry column with all valid non-null geoms
    if (gdf.active_geometry_name is None) or (
        gdf.geometry.isna() | gdf.geometry.is_empty | ~gdf.geometry.is_valid
    ).any():
        err_msg += "`gdf` must contain valid, non-null geometries`. "
        is_valid = False

    _report_validation(is_valid, valid_msg, warn_msg, err_msg)


def _validate_node_edge_gdfs(
    gdf_nodes: gpd.GeoDataFrame,
    gdf_edges: gpd.GeoDataFrame,
    *,
    strict: bool = True,
) -> None:
    """
    Validate that node/edge GeoDataFrames can be converted to a MultiDiGraph.

    Raises a `ValidationError` if validation fails.

    Parameters
    ----------
    gdf_nodes
        GeoDataFrame of graph nodes uniquely indexed by `osmid`.
    gdf_edges
        GeoDataFrame of graph edges uniquely multi-indexed by `(u, v, key)`.
    strict
        If `True`, elevate warnings to errors.
    """
    is_valid = True
    valid_msg = "Validated that node/edge GeoDataFrames can be converted to a MultiDiGraph."
    warn_msg = ""
    err_msg = ""

    # ensure type is GeoDataFrame
    if not (isinstance(gdf_nodes, gpd.GeoDataFrame) and isinstance(gdf_edges, gpd.GeoDataFrame)):
        # if they are not both GeoDataFrames
        err_msg += "`gdf_nodes` and `gdf_edges` must be GeoDataFrames. "
        is_valid = False
    # if they are both GeoDataFrames...
    # warn user if geometry values differ from coordinates in x/y columns,
    # because we ignore the geometry column
    elif gdf_nodes.active_geometry_name is not None:
        msg = (
            "Will ignore the `gdf_nodes` 'geometry' column, though its values "
            "differ from the coordinates in the 'x' and 'y' columns. "
        )
        try:
            all_x_match = (gdf_nodes.geometry.x == gdf_nodes["x"]).all()
            all_y_match = (gdf_nodes.geometry.y == gdf_nodes["y"]).all()
            if not (all_x_match and all_y_match):
                # warn if x/y coords don't match geometry column
                warn_msg += msg
                if strict:
                    is_valid = False
        except ValueError:
            # warn if geometry column contains non-point geometry types
            warn_msg += msg
            if strict:
                is_valid = False

    # ensure gdf_nodes has x and y columns representing node geometries
    if not ("x" in gdf_nodes.columns and "y" in gdf_nodes.columns):
        err_msg += "`gdf_nodes` must have 'x' and 'y' columns. "
        is_valid = False

    # ensure gdf_nodes and gdf_edges are uniquely indexed
    if not (gdf_nodes.index.is_unique and gdf_edges.index.is_unique):
        err_msg += "`gdf_nodes` and `gdf_edges` must each be uniquely indexed. "
        is_valid = False

    # ensure 1) gdf_edges are multi-indexed with 3 levels and 2) that its u
    # and v values (first two index levels) all appear among gdf_nodes index
    edges_index_levels = 3
    check1 = gdf_edges.index.nlevels == edges_index_levels
    try:
        uv = set(gdf_edges.index.get_level_values(0)) | set(gdf_edges.index.get_level_values(1))
        check2 = uv.issubset(set(gdf_nodes.index))
    except IndexError:
        check2 = False
    if not (check1 and check2):
        err_msg += "`gdf_edges` must be multi-indexed by `(u, v, key)`. "
        is_valid = False

    _report_validation(is_valid, valid_msg, warn_msg, err_msg)


def _validate_nodes(G: nx.MultiDiGraph, strict: bool) -> tuple[bool, str, str]:  # noqa: FBT001
    """
    Validate that a graph's nodes satisfy OSMnx expectations.

    Parameters
    ----------
    G
        The input graph.
    strict
        If `True`, elevate warnings to errors.

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
        If `True`, elevate warnings to errors.

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


def _validate_graph(G: nx.MultiDiGraph, *, strict: bool = True) -> None:
    """
    Validate that a graph object satisfies OSMnx expectations.

    Raises `ox._errors.ValidationError` if validation fails.

    Parameters
    ----------
    G
        The input graph.
    strict
        If `True`, elevate warnings to errors.
    """
    # validate graph, nodes, and edges
    is_valid_graph, err_msg_graph, warn_msg_graph = _validate_graph_attrs(G)
    is_valid_nodes, err_msg_nodes, warn_msg_nodes = _validate_nodes(G, strict)
    is_valid_edges, err_msg_edges, warn_msg_edges = _validate_edges(G, strict)

    # report results
    is_valid = is_valid_graph and is_valid_nodes and is_valid_edges
    err_msg = err_msg_graph + err_msg_nodes + err_msg_edges
    warn_msg = warn_msg_graph + warn_msg_nodes + warn_msg_edges
    valid_msg = "Successfully validated graph."
    _report_validation(is_valid, valid_msg, warn_msg, err_msg)


def _report_validation(is_valid: bool, valid_msg: str, warn_msg: str, err_msg: str) -> None:  # noqa: FBT001
    """
    Report validation results by logging, warning, or raising an exception.

    Parameters
    ----------
    is_valid
        Whether or not the validation succeeded.
    valid_msg
        The message to log if validation succeeded.
    warn_msg
        Any warning messages to log and either issue a warning or include in
        error message.
    err_msg
        Any error messages to include when raising exception if validation
        failed.
    """
    if is_valid:
        log(valid_msg, level=lg.INFO)
        if warn_msg != "":
            log(warn_msg, level=lg.WARNING)
            warn(warn_msg, category=UserWarning, stacklevel=2)
    else:
        log(err_msg + warn_msg, level=lg.ERROR)
        raise ValidationError(err_msg + warn_msg)
