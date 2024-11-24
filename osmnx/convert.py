"""Convert spatial graphs to/from different data types."""

from __future__ import annotations

import itertools
import logging as lg
from typing import Any
from typing import Literal
from typing import overload
from warnings import warn

import geopandas as gpd
import networkx as nx
import pandas as pd
from shapely import LineString
from shapely import Point

from . import utils


# nodes and edges are both missing (therefore both default true)
@overload
def graph_to_gdfs(
    G: nx.MultiGraph | nx.MultiDiGraph,
    *,
    node_geometry: bool = True,
    fill_edge_geometry: bool = True,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]: ...


# both present/True
@overload
def graph_to_gdfs(
    G: nx.MultiGraph | nx.MultiDiGraph,
    *,
    nodes: Literal[True],
    edges: Literal[True],
    node_geometry: bool = True,
    fill_edge_geometry: bool = True,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]: ...


# both present, nodes true, edges false
@overload
def graph_to_gdfs(
    G: nx.MultiGraph | nx.MultiDiGraph,
    *,
    nodes: Literal[True],
    edges: Literal[False],
    node_geometry: bool = True,
    fill_edge_geometry: bool = True,
) -> gpd.GeoDataFrame: ...


# both present, nodes false, edges true
@overload
def graph_to_gdfs(
    G: nx.MultiGraph | nx.MultiDiGraph,
    *,
    nodes: Literal[False],
    edges: Literal[True],
    node_geometry: bool = True,
    fill_edge_geometry: bool = True,
) -> gpd.GeoDataFrame: ...


# nodes missing (therefore default true), edges present/true
@overload
def graph_to_gdfs(
    G: nx.MultiGraph | nx.MultiDiGraph,
    *,
    edges: Literal[True],
    node_geometry: bool = True,
    fill_edge_geometry: bool = True,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]: ...


# nodes missing (therefore default true), edges present/false
@overload
def graph_to_gdfs(
    G: nx.MultiGraph | nx.MultiDiGraph,
    *,
    edges: Literal[False],
    node_geometry: bool = True,
    fill_edge_geometry: bool = True,
) -> gpd.GeoDataFrame: ...


# nodes present/true, edges missing (therefore default true)
@overload
def graph_to_gdfs(
    G: nx.MultiGraph | nx.MultiDiGraph,
    *,
    nodes: Literal[True],
    edges: bool = True,
    node_geometry: bool = True,
    fill_edge_geometry: bool = True,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]: ...


# nodes present/false, edges missing (therefore default true)
@overload
def graph_to_gdfs(
    G: nx.MultiGraph | nx.MultiDiGraph,
    *,
    nodes: Literal[False],
    edges: bool = True,
    node_geometry: bool = True,
    fill_edge_geometry: bool = True,
) -> gpd.GeoDataFrame: ...


def graph_to_gdfs(
    G: nx.MultiGraph | nx.MultiDiGraph,
    *,
    nodes: bool = True,
    edges: bool = True,
    node_geometry: bool = True,
    fill_edge_geometry: bool = True,
) -> gpd.GeoDataFrame | tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Convert a MultiGraph or MultiDiGraph to node and/or edge GeoDataFrames.

    This function is the inverse of `graph_from_gdfs`.

    Parameters
    ----------
    G
        Input graph.
    nodes
        If True, convert graph nodes to a GeoDataFrame and return it.
    edges
        If True, convert graph edges to a GeoDataFrame and return it.
    node_geometry
        If True, create a geometry column from node "x" and "y" attributes.
    fill_edge_geometry
        If True, fill missing edge geometry fields using endpoint nodes'
        coordinates to create a LineString.

    Returns
    -------
    gdf_nodes or gdf_edges or (gdf_nodes, gdf_edges)
        `gdf_nodes` is indexed by `osmid` and `gdf_edges` is multi-indexed by
        `(u, v, key)` following normal MultiGraph/MultiDiGraph structure.
    """
    crs = G.graph["crs"]

    if nodes:
        if len(G.nodes) == 0:  # pragma: no cover
            msg = "Graph contains no nodes."
            raise ValueError(msg)

        uvk, data = zip(*G.nodes(data=True))

        if node_geometry:
            # convert node x/y attributes to Points for geometry column
            node_geoms = (Point(d["x"], d["y"]) for d in data)
            gdf_nodes = gpd.GeoDataFrame(data, index=uvk, crs=crs, geometry=list(node_geoms))
        else:
            gdf_nodes = gpd.GeoDataFrame(data, index=uvk)

        gdf_nodes.index = gdf_nodes.index.rename("osmid")
        msg = "Created nodes GeoDataFrame from graph"
        utils.log(msg, level=lg.INFO)

    if edges:
        if len(G.edges) == 0:  # pragma: no cover
            msg = "Graph contains no edges."
            raise ValueError(msg)

        u, v, k, data = zip(*G.edges(keys=True, data=True))

        if fill_edge_geometry:
            node_coords = {n: (G.nodes[n]["x"], G.nodes[n]["y"]) for n in G}
            edge_geoms = (
                d.get("geometry", LineString((node_coords[u], node_coords[v])))
                for u, v, _, d in G.edges(keys=True, data=True)
            )
            gdf_edges = gpd.GeoDataFrame(data, crs=crs, geometry=list(edge_geoms))

        else:
            gdf_edges = gpd.GeoDataFrame(data)
            if "geometry" not in gdf_edges.columns:
                # if no edges have a geometry attribute, create null column
                gdf_edges = gdf_edges.set_geometry([None] * len(gdf_edges))
            gdf_edges = gdf_edges.set_crs(crs)

        # add u, v, key attributes as index
        gdf_edges["u"] = u
        gdf_edges["v"] = v
        gdf_edges["key"] = k
        gdf_edges = gdf_edges.set_index(["u", "v", "key"])

        msg = "Created edges GeoDataFrame from graph"
        utils.log(msg, level=lg.INFO)

    if nodes and edges:
        return gdf_nodes, gdf_edges

    if nodes:
        return gdf_nodes

    if edges:
        return gdf_edges

    # otherwise
    msg = "You must request nodes or edges or both."
    raise ValueError(msg)


def _validate_node_edge_gdfs(
    gdf_nodes: gpd.GeoDataFrame,
    gdf_edges: gpd.GeoDataFrame,
) -> None:
    """
    Validate that node/edge GeoDataFrames can be converted to a MultiDiGraph.

    Raises a `ValueError` if validation fails.

    Parameters
    ----------
    gdf_nodes
        GeoDataFrame of graph nodes uniquely indexed by `osmid`.
    gdf_edges
        GeoDataFrame of graph edges uniquely multi-indexed by `(u, v, key)`.
    graph_attrs

    Returns
    -------
    None
    """
    # ensure gdf_nodes contains x and y columns representing node geometries
    if not ("x" in gdf_nodes.columns and "y" in gdf_nodes.columns):  # pragma: no cover
        msg = "`gdf_nodes` must contain 'x' and 'y' columns."
        raise ValueError(msg)

    # ensure gdf_nodes and gdf_edges are uniquely indexed
    if not (gdf_nodes.index.is_unique and gdf_edges.index.is_unique):  # pragma: no cover
        msg = "`gdf_nodes` and `gdf_edges` must each be uniquely indexed."
        raise ValueError(msg)

    # ensure 1) gdf_edges are multi-indexed with 3 levels and 2) that its u
    # and v values (first two index levels) all appear among gdf_nodes index
    edges_index_levels = 3
    check1 = gdf_edges.index.nlevels == edges_index_levels
    uv = set(gdf_edges.index.get_level_values(0)) | set(gdf_edges.index.get_level_values(1))
    check2 = uv.issubset(set(gdf_nodes.index))
    if not (check1 and check2):  # pragma: no cover
        msg = "`gdf_edges` must be multi-indexed by `(u, v, key)`."
        raise ValueError(msg)

    # warn user if geometry values differ from coordinates in x/y columns,
    # because we discard the geometry column
    if gdf_nodes.active_geometry_name is not None:  # pragma: no cover
        msg = (
            "Discarding the `gdf_nodes` 'geometry' column, though its values "
            "differ from the coordinates in the 'x' and 'y' columns."
        )
        try:
            all_x_match = (gdf_nodes.geometry.x == gdf_nodes["x"]).all()
            all_y_match = (gdf_nodes.geometry.y == gdf_nodes["y"]).all()
            if not (all_x_match and all_y_match):
                # warn if x/y coords don't match geometry column
                warn(msg, category=UserWarning, stacklevel=2)
        except ValueError:  # pragma: no cover
            # warn if geometry column contains non-point geometry types
            warn(msg, category=UserWarning, stacklevel=2)


def graph_from_gdfs(
    gdf_nodes: gpd.GeoDataFrame,
    gdf_edges: gpd.GeoDataFrame,
    *,
    graph_attrs: dict[str, Any] | None = None,
) -> nx.MultiDiGraph:
    """
    Convert node and edge GeoDataFrames to a MultiDiGraph.

    This function is the inverse of `graph_to_gdfs` and is designed to work in
    conjunction with it. However, you can convert arbitrary node and edge
    GeoDataFrames as long as 1) `gdf_nodes` is uniquely indexed by `osmid`, 2)
    `gdf_nodes` contains `x` and `y` coordinate columns representing node
    geometries, 3) `gdf_edges` is uniquely multi-indexed by `(u, v, key)`
    (following normal MultiDiGraph structure). This allows you to load any
    node/edge Shapefiles or GeoPackage layers as GeoDataFrames then convert
    them to a MultiDiGraph for network analysis.

    Note that any `geometry` attribute on `gdf_nodes` is discarded, since `x`
    and `y` provide the necessary node geometry information instead.

    Parameters
    ----------
    gdf_nodes
        GeoDataFrame of graph nodes uniquely indexed by `osmid`.
    gdf_edges
        GeoDataFrame of graph edges uniquely multi-indexed by `(u, v, key)`.
    graph_attrs
        The new `G.graph` attribute dictionary. If None, use `gdf_edges`'s CRS
        as the only graph-level attribute (`gdf_edges` must have its `crs`
        attribute set).

    Returns
    -------
    G
    """
    _validate_node_edge_gdfs(gdf_nodes, gdf_edges)

    # drop geometry column from gdf_nodes (since we use x and y for geometry
    # information), but warn the user if the geometry values differ from the
    # coordinates in the x and y columns. this results in a df instead of gdf.
    if gdf_nodes.active_geometry_name is None:  # pragma: no cover
        df_nodes = pd.DataFrame(gdf_nodes)
    else:
        df_nodes = gdf_nodes.drop(columns=gdf_nodes.active_geometry_name)

    # create graph and add graph-level attribute dict
    if graph_attrs is None:
        graph_attrs = {"crs": gdf_edges.crs}
    G = nx.MultiDiGraph(**graph_attrs)

    # add edges and their attributes to graph, but filter out null attribute
    # values so that edges only get attributes with non-null values
    attr_names = gdf_edges.columns.to_list()
    for (u, v, k), attr_vals in zip(gdf_edges.index, gdf_edges.to_numpy()):
        data_all = zip(attr_names, attr_vals)
        data = {name: val for name, val in data_all if isinstance(val, list) or pd.notna(val)}
        G.add_edge(u, v, key=k, **data)

    # add any nodes with no incident edges, since they wouldn't be added above
    G.add_nodes_from(set(df_nodes.index) - set(G.nodes))

    # now all nodes are added, so set nodes' attributes
    for col in df_nodes.columns:
        nx.set_node_attributes(G, name=col, values=df_nodes[col].dropna())

    msg = "Created graph from node/edge GeoDataFrames"
    utils.log(msg, level=lg.INFO)
    return G


def to_digraph(G: nx.MultiDiGraph, *, weight: str = "length") -> nx.DiGraph:
    """
    Convert MultiDiGraph to DiGraph.

    Chooses between parallel edges by minimizing `weight` attribute value. See
    also `to_undirected` to convert MultiDiGraph to MultiGraph.

    Parameters
    ----------
    G
        Input graph.
    weight
        Attribute value to minimize when choosing between parallel edges.

    Returns
    -------
    G
    """
    # make a copy to not mutate original graph object caller passed in
    G = G.copy()
    to_remove: list[tuple[int, int, int]] = []

    # identify all the parallel edges in the MultiDiGraph
    parallels = ((u, v) for u, v in G.edges(keys=False) if G.number_of_edges(u, v) > 1)

    # among all sets of parallel edges, remove all except the one with the
    # minimum "weight" attribute value
    for u, v in set(parallels):
        k_min, _ = min(G.get_edge_data(u, v).items(), key=lambda x: x[1][weight])
        to_remove.extend((u, v, k) for k in G[u][v] if k != k_min)

    G.remove_edges_from(to_remove)
    msg = "Converted MultiDiGraph to DiGraph"
    utils.log(msg, level=lg.INFO)

    return nx.DiGraph(G)


def to_undirected(G: nx.MultiDiGraph) -> nx.MultiGraph:
    """
    Convert MultiDiGraph to undirected MultiGraph.

    Maintains parallel edges only if their geometries differ. See also
    `to_digraph` to convert MultiDiGraph to DiGraph.

    Parameters
    ----------
    G
        Input graph.

    Returns
    -------
    Gu
    """
    # make a copy to not mutate original graph object caller passed in
    G = G.copy()

    # set from/to nodes before making graph undirected
    for u, v, d in G.edges(data=True):
        d["from"] = u
        d["to"] = v

        # add geometry if missing, to compare parallel edges' geometries
        if "geometry" not in d:
            point_u = (G.nodes[u]["x"], G.nodes[u]["y"])
            point_v = (G.nodes[v]["x"], G.nodes[v]["y"])
            d["geometry"] = LineString([point_u, point_v])

    # increment parallel edges' keys so we don't retain only one edge of sets
    # of true parallel edges when we convert from MultiDiGraph to MultiGraph
    G = _update_edge_keys(G)

    # convert MultiDiGraph to MultiGraph, retaining edges in both directions
    # of parallel edges and self-loops for now
    Gu = nx.MultiGraph(**G.graph)
    Gu.add_nodes_from(G.nodes(data=True))
    Gu.add_edges_from(G.edges(keys=True, data=True))

    # the previous operation added all directed edges from G as undirected
    # edges in Gu. we now have duplicate edges for each bidirectional parallel
    # edge or self-loop. so, look through the edges and remove any duplicates.
    duplicate_edges = set()
    for u1, v1, key1, data1 in Gu.edges(keys=True, data=True):
        # if we haven't already flagged this edge as a duplicate
        if (u1, v1, key1) not in duplicate_edges:
            # look at every other edge between u and v, one at a time
            for key2 in Gu[u1][v1]:
                # don't compare this edge to itself
                if key1 != key2:
                    # compare the first edge's data to the second's
                    # if they match up, flag the duplicate for removal
                    data2 = Gu.edges[u1, v1, key2]
                    if _is_duplicate_edge(data1, data2):
                        duplicate_edges.add((u1, v1, key2))

    Gu.remove_edges_from(duplicate_edges)
    msg = "Converted MultiDiGraph to undirected MultiGraph"
    utils.log(msg, level=lg.INFO)

    return Gu


def _is_duplicate_edge(data1: dict[str, Any], data2: dict[str, Any]) -> bool:
    """
    Check if two graph edge data dicts have the same `osmid` and `geometry`.

    Parameters
    ----------
    data1
        The first edge's attribute data.
    data2
        The second edge's attribute data.

    Returns
    -------
    is_dupe
    """
    is_dupe = False

    # if either edge's osmid contains multiple values (due to simplification)
    # compare them as sets to see if they contain the same values
    osmid1 = set(data1["osmid"]) if isinstance(data1["osmid"], list) else data1["osmid"]
    osmid2 = set(data2["osmid"]) if isinstance(data2["osmid"], list) else data2["osmid"]

    # if they contain the same osmid or set of osmids (due to simplification)
    if osmid1 == osmid2:
        # if both edges have geometry attributes and they match each other
        if ("geometry" in data1) and ("geometry" in data2):
            if _is_same_geometry(data1["geometry"], data2["geometry"]):
                is_dupe = True

        # if neither edge has a geometry attribute
        elif ("geometry" not in data1) and ("geometry" not in data2):
            is_dupe = True

        # if one edge has geometry attribute but the other doesn't: not dupes
        else:
            pass

    return is_dupe


def _is_same_geometry(ls1: LineString, ls2: LineString) -> bool:
    """
    Determine if two LineString geometries are the same (in either direction).

    Check both the normal and reversed orders of their constituent points.

    Parameters
    ----------
    ls1
        The first LineString geometry.
    ls2
        The second LineString geometry.

    Returns
    -------
    is_same
    """
    # extract coordinates from each LineString geometry
    geom1 = [tuple(coords) for coords in ls1.xy]
    geom2 = [tuple(coords) for coords in ls2.xy]

    # reverse the first LineString's coordinates' direction
    geom1_r = [tuple(reversed(coords)) for coords in ls1.xy]

    # if second geometry matches first in either direction, return True
    return geom2 in (geom1, geom1_r)


def _update_edge_keys(G: nx.MultiDiGraph) -> nx.MultiDiGraph:
    """
    Increment key of one edge of parallel edges that differ in geometry.

    For example, two streets from `u` to `v` that bow away from each other as
    separate streets, rather than opposite direction edges of a single street.
    Increment one of these edge's keys so that they do not match across
    `(u, v, k)` or `(v, u, k)` so we can add both to an undirected MultiGraph.

    Parameters
    ----------
    G
        Input graph.

    Returns
    -------
    G
    """
    # identify all the edges that are duplicates based on a sorted combination
    # of their origin, destination, and key. that is, edge uv will match edge vu
    # as a duplicate, but only if they have the same key
    edges = graph_to_gdfs(G, nodes=False, fill_edge_geometry=False)
    edges["uvk"] = ["_".join([*sorted([str(u), str(v)]), str(k)]) for u, v, k in edges.index]
    mask = edges["uvk"].duplicated(keep=False)
    dupes = edges[mask].dropna(subset=["geometry"])

    different_streets = []
    groups = dupes[["geometry", "uvk"]].groupby("uvk")

    # for each group of duplicate edges
    for _, group in groups:
        # for each pair of edges within this group
        for geom1, geom2 in itertools.combinations(group["geometry"], 2):
            # if they don't have the same geometry, flag them as different
            # streets: flag edge uvk, but not edge vuk, otherwise we would
            # increment both their keys and they'll still duplicate each other
            if not _is_same_geometry(geom1, geom2):
                different_streets.append(group.index[0])

    # for each unique different street, increment its key to make it unique
    for u, v, k in set(different_streets):
        new_key = max(list(G[u][v]) + list(G[v][u])) + 1
        G.add_edge(u, v, key=new_key, **G.get_edge_data(u, v, k))
        G.remove_edge(u, v, key=k)

    return G
