"""Simplify, correct, and consolidate spatial graph nodes and edges."""

from __future__ import annotations

import logging as lg
from typing import TYPE_CHECKING
from typing import Any

import geopandas as gpd
import networkx as nx
import pandas as pd
from shapely import LineString
from shapely import Point

from . import convert
from . import stats
from . import utils
from ._errors import GraphSimplificationError

if TYPE_CHECKING:
    from collections.abc import Iterable
    from collections.abc import Iterator


def _is_endpoint(
    G: nx.MultiDiGraph,
    node: int,
    node_attrs_include: Iterable[str] | None,
    edge_attrs_differ: Iterable[str] | None,
) -> bool:
    """
    Determine if a node is a true endpoint of an edge.

    Return True if the node is a "true" endpoint of an edge in the network,
    otherwise False. OpenStreetMap data includes many nodes that exist only as
    geometric vertices to allow ways to curve. `node` is a true edge endpoint
    if it satisfies at least 1 of the following 5 rules:

    1) It is its own neighbor (ie, it self-loops).

    2) Or, it has no incoming edges or no outgoing edges (ie, all its incident
    edges are inbound or all its incident edges are outbound).

    3) Or, it does not have exactly two neighbors and degree of 2 or 4.

    4) Or, if `node_attrs_include` is not None and it has one or more of the
    attributes in `node_attrs_include`.

    5) Or, if `edge_attrs_differ` is not None and its incident edges have
    different values than each other for any of the edge attributes in
    `edge_attrs_differ`.

    Parameters
    ----------
    G
        Input graph.
    node
        The ID of the node.
    node_attrs_include
        Node attribute names for relaxing the strictness of endpoint
        determination. A node is always an endpoint if it possesses one or
        more of the attributes in `node_attrs_include`.
    edge_attrs_differ
        Edge attribute names for relaxing the strictness of endpoint
        determination. A node is always an endpoint if its incident edges have
        different values than each other for any attribute in
        `edge_attrs_differ`.

    Returns
    -------
    endpoint
    """
    neighbors = set(list(G.predecessors(node)) + list(G.successors(node)))
    n = len(neighbors)
    d = G.degree(node)

    # RULE 1
    # if the node appears in its list of neighbors, it self-loops: this is
    # always an endpoint
    if node in neighbors:
        return True

    # RULE 2
    # if node has no incoming edges or no outgoing edges, it is an endpoint
    if G.out_degree(node) == 0 or G.in_degree(node) == 0:
        return True

    # RULE 3
    # else, if it does NOT have 2 neighbors AND either 2 or 4 directed edges,
    # it is an endpoint. either it has 1 or 3+ neighbors, in which case it is
    # a dead-end or an intersection of multiple streets or it has 2 neighbors
    # but 3 degree (indicating a change from oneway to twoway) or more than 4
    # degree (indicating a parallel edge) and thus is an endpoint
    if not ((n == 2) and (d in {2, 4})):  # noqa: PLR2004
        return True

    # RULE 4
    # non-strict mode: does it contain an attr denoting that it is an endpoint
    if node_attrs_include is not None and len(set(node_attrs_include) & G.nodes[node].keys()) > 0:
        return True

    # RULE 5
    # non-strict mode: do its incident edges have different attr values? for
    # each attribute to check, collect the attribute's values in all inbound
    # and outbound edges. if there is more than 1 unique value then this node
    # is an endpoint
    if edge_attrs_differ is not None:
        for attr in edge_attrs_differ:
            in_values = {v for _, _, v in G.in_edges(node, data=attr, keys=False)}
            out_values = {v for _, _, v in G.out_edges(node, data=attr, keys=False)}
            if len(in_values | out_values) > 1:
                return True

    # if none of the preceding rules passed, then it is not an endpoint
    return False


def _build_path(
    G: nx.MultiDiGraph,
    endpoint: int,
    endpoint_successor: int,
    endpoints: set[int],
) -> list[int]:
    """
    Build a path of nodes from one endpoint node to next endpoint node.

    Parameters
    ----------
    G
        Input graph.
    endpoint
        Ehe endpoint node from which to start the path.
    endpoint_successor
        The successor of endpoint through which the path to the next endpoint
        will be built.
    endpoints
        The set of all nodes in the graph that are endpoints.

    Returns
    -------
    path
        The first and last items in the resulting path list are endpoint
        nodes, and all other items are interstitial nodes that can be removed
        subsequently.
    """
    # start building path from endpoint node through its successor
    path = [endpoint, endpoint_successor]

    # for each successor of the endpoint's successor
    for this_successor in G.successors(endpoint_successor):
        successor = this_successor
        if successor not in path:
            # if this successor is already in the path, ignore it, otherwise add
            # it to the path
            path.append(successor)
            while successor not in endpoints:
                # find successors (of current successor) not in path
                successors = [n for n in G.successors(successor) if n not in path]

                # 99%+ of the time there will be only 1 successor: add to path
                if len(successors) == 1:
                    successor = successors[0]
                    path.append(successor)

                # handle relatively rare cases or OSM digitization quirks
                elif len(successors) == 0:
                    if endpoint in G.successors(successor):
                        # we have come to the end of a self-looping edge, so
                        # add first node to end of path to close it and return
                        return [*path, endpoint]

                    # otherwise, this can happen due to OSM digitization error
                    # where a one-way street turns into a two-way here, but
                    # duplicate incoming one-way edges are present
                    msg = f"Unexpected simplify pattern handled near {successor}"
                    utils.log(msg, level=lg.WARNING)
                    return path
                else:  # pragma: no cover
                    # if successor has >1 successors, then successor must have
                    # been an endpoint because you can go in 2 new directions.
                    # this should never occur in practice
                    msg = f"Impossible simplify pattern failed near {successor}."
                    raise GraphSimplificationError(msg)

            # if this successor is an endpoint, we've completed the path
            return path

    # if endpoint_successor has no successors not already in the path, return
    # the current path: this is usually due to a digitization quirk on OSM
    return path


def _get_paths_to_simplify(
    G: nx.MultiDiGraph,
    node_attrs_include: Iterable[str] | None,
    edge_attrs_differ: Iterable[str] | None,
) -> Iterator[list[int]]:
    """
    Generate all the paths to be simplified between endpoint nodes.

    The path is ordered from the first endpoint, through the interstitial nodes,
    to the second endpoint.

    Parameters
    ----------
    G
        Input graph.
    node_attrs_include
        Node attribute names for relaxing the strictness of endpoint
        determination. A node is always an endpoint if it possesses one or
        more of the attributes in `node_attrs_include`.
    edge_attrs_differ
        Edge attribute names for relaxing the strictness of endpoint
        determination. A node is always an endpoint if its incident edges have
        different values than each other for any attribute in
        `edge_attrs_differ`.

    Yields
    ------
    path_to_simplify
    """
    # first identify all the nodes that are endpoints
    endpoints = {n for n in G.nodes if _is_endpoint(G, n, node_attrs_include, edge_attrs_differ)}
    msg = f"Identified {len(endpoints):,} edge endpoints"
    utils.log(msg, level=lg.INFO)

    # for each endpoint node, look at each of its successor nodes
    for endpoint in endpoints:
        for successor in G.successors(endpoint):
            if successor not in endpoints:
                # if endpoint node's successor is not an endpoint, build path
                # from the endpoint node, through the successor, and on to the
                # next endpoint node
                yield _build_path(G, endpoint, successor, endpoints)


def _remove_rings(
    G: nx.MultiDiGraph,
    node_attrs_include: Iterable[str] | None,
    edge_attrs_differ: Iterable[str] | None,
) -> nx.MultiDiGraph:
    """
    Remove all graph components that consist only of a single chordless cycle.

    This identifies all connected components in the graph that consist only of
    a single isolated self-contained ring, and removes them from the graph.

    Parameters
    ----------
    G
        Input graph.
    node_attrs_include
        Node attribute names for relaxing the strictness of endpoint
        determination. A node is always an endpoint if it possesses one or
        more of the attributes in `node_attrs_include`.
    edge_attrs_differ
        Edge attribute names for relaxing the strictness of endpoint
        determination. A node is always an endpoint if its incident edges have
        different values than each other for any attribute in
        `edge_attrs_differ`.

    Returns
    -------
    G
        Graph with all chordless cycle components removed.
    """
    to_remove = set()
    for wcc in nx.weakly_connected_components(G):
        if not any(_is_endpoint(G, n, node_attrs_include, edge_attrs_differ) for n in wcc):
            to_remove.update(wcc)
    G.remove_nodes_from(to_remove)
    return G


def simplify_graph(  # noqa: C901, PLR0912
    G: nx.MultiDiGraph,
    *,
    node_attrs_include: Iterable[str] | None = None,
    edge_attrs_differ: Iterable[str] | None = None,
    remove_rings: bool = True,
    track_merged: bool = False,
    edge_attr_aggs: dict[str, Any] | None = None,
) -> nx.MultiDiGraph:
    """
    Simplify a graph's topology by removing interstitial nodes.

    This simplifies the graph's topology by removing all nodes that are not
    intersections or dead-ends, by creating an edge directly between the end
    points that encapsulate them while retaining the full geometry of the
    original edges, saved as a new `geometry` attribute on the new edge.

    Note that only simplified edges receive a `geometry` attribute. Some of
    the resulting consolidated edges may comprise multiple OSM ways, and if
    so, their unique attribute values are stored as a list. Optionally, the
    simplified edges can receive a `merged_edges` attribute that contains a
    list of all the `(u, v)` node pairs that were merged together.

    Use the `node_attrs_include` or `edge_attrs_differ` parameters to relax
    simplification strictness. For example, `edge_attrs_differ=["osmid"]` will
    retain every node whose incident edges have different OSM IDs. This lets
    you keep nodes at elbow two-way intersections (but be aware that sometimes
    individual blocks have multiple OSM IDs within them too). You could also
    use this parameter to retain nodes where sidewalks or bike lanes begin/end
    in the middle of a block. Or for example, `node_attrs_include=["highway"]`
    will retain every node with a "highway" attribute (regardless of its
    value), even if it does not represent a street junction.

    Parameters
    ----------
    G
        Input graph.
    node_attrs_include
        Node attribute names for relaxing the strictness of endpoint
        determination. A node is always an endpoint if it possesses one or
        more of the attributes in `node_attrs_include`.
    edge_attrs_differ
        Edge attribute names for relaxing the strictness of endpoint
        determination. A node is always an endpoint if its incident edges have
        different values than each other for any attribute in
        `edge_attrs_differ`.
    remove_rings
        If True, remove any graph components that consist only of a single
        chordless cycle (i.e., an isolated self-contained ring).
    track_merged
        If True, add `merged_edges` attribute on simplified edges, containing
        a list of all the `(u, v)` node pairs that were merged together.
    edge_attr_aggs
        Allows user to aggregate edge segment attributes when simplifying an
        edge. Keys are edge attribute names and values are aggregation
        functions to apply to these attributes when they exist for a set of
        edges being merged. Edge attributes not in `edge_attr_aggs` will
        contain the unique values across the merged edge segments. If None,
        defaults to `{"length": sum, "travel_time": sum}`.

    Returns
    -------
    G
        Topologically simplified graph, with a new `geometry` attribute on
        each simplified edge.
    """
    if G.graph.get("simplified"):  # pragma: no cover
        msg = "This graph has already been simplified, cannot simplify it again."
        raise GraphSimplificationError(msg)

    msg = "Begin topologically simplifying the graph..."
    utils.log(msg, level=lg.INFO)

    # default edge segment attributes to aggregate upon simplification
    if edge_attr_aggs is None:
        edge_attr_aggs = {"length": sum, "travel_time": sum}

    # make a copy to not mutate original graph object caller passed in
    G = G.copy()
    initial_node_count = len(G)
    initial_edge_count = len(G.edges)
    all_nodes_to_remove = []
    all_edges_to_add = []

    # generate each path that needs to be simplified
    for path in _get_paths_to_simplify(G, node_attrs_include, edge_attrs_differ):
        # add the interstitial edges we're removing to a list so we can retain
        # their spatial geometry
        merged_edges = []
        path_attributes: dict[str, Any] = {}
        for u, v in zip(path[:-1], path[1:]):
            if track_merged:
                # keep track of the edges that were merged
                merged_edges.append((u, v))

            # there should rarely be multiple edges between interstitial nodes
            # usually happens if OSM has duplicate ways digitized for just one
            # street... we will keep only one of the edges (see below)
            edge_count = G.number_of_edges(u, v)
            if edge_count != 1:
                msg = f"Found {edge_count} edges between {u} and {v} when simplifying"
                utils.log(msg, level=lg.WARNING)

            # get edge between these nodes: if multiple edges exist between
            # them (see above), we retain only one in the simplified graph
            # We can't assume that there exists an edge from u to v
            # with key=0, so we get a list of all edges from u to v
            # and just take the first one.
            edge_data = next(iter(G.get_edge_data(u, v).values()))
            for attr in edge_data:
                if attr in path_attributes:
                    # if this key already exists in the dict, append it to the
                    # value list
                    path_attributes[attr].append(edge_data[attr])
                else:
                    # if this key doesn't already exist, set the value to a list
                    # containing the one value
                    path_attributes[attr] = [edge_data[attr]]

        # consolidate the path's edge segments' attribute values
        for attr_name, attr_values in path_attributes.items():
            if attr_name in edge_attr_aggs:
                # if this attribute's values must be aggregated, do so now
                agg_func = edge_attr_aggs[attr_name]
                path_attributes[attr_name] = agg_func(attr_values)
            elif len(set(attr_values)) == 1:
                # if there's only 1 unique value, keep that single value
                path_attributes[attr_name] = attr_values[0]
            else:
                # otherwise, if there are multiple uniques, keep one of each
                path_attributes[attr_name] = list(set(attr_values))

        # construct the new consolidated edge's geometry for this path
        path_attributes["geometry"] = LineString(
            [Point((G.nodes[node]["x"], G.nodes[node]["y"])) for node in path],
        )

        if track_merged:
            # add the merged edges as a new attribute of the simplified edge
            path_attributes["merged_edges"] = merged_edges

        # add the nodes and edge to their lists for processing at the end
        all_nodes_to_remove.extend(path[1:-1])
        all_edges_to_add.append(
            {"origin": path[0], "destination": path[-1], "attr_dict": path_attributes},
        )

    # for each edge to add in the list we assembled, create a new edge between
    # the origin and destination
    for edge in all_edges_to_add:
        G.add_edge(edge["origin"], edge["destination"], **edge["attr_dict"])

    # finally remove all the interstitial nodes between the new edges
    G.remove_nodes_from(set(all_nodes_to_remove))

    if remove_rings:
        G = _remove_rings(G, node_attrs_include, edge_attrs_differ)

    # mark the graph as having been simplified
    G.graph["simplified"] = True

    msg = (
        f"Simplified graph: {initial_node_count:,} to {len(G):,} nodes, "
        f"{initial_edge_count:,} to {len(G.edges):,} edges"
    )
    utils.log(msg, level=lg.INFO)
    return G


def consolidate_intersections(
    G: nx.MultiDiGraph,
    *,
    tolerance: float | dict[int, float] = 10,
    rebuild_graph: bool = True,
    dead_ends: bool = False,
    reconnect_edges: bool = True,
    node_attr_aggs: dict[str, Any] | None = None,
) -> nx.MultiDiGraph | gpd.GeoSeries:
    """
    Consolidate intersections comprising clusters of nearby nodes.

    Merges nearby nodes and returns either their centroids or a rebuilt graph
    with consolidated intersections and reconnected edge geometries. The
    `tolerance` argument can be a single value applied to all nodes or
    individual per-node values. It should be adjusted to approximately match
    street design standards in the specific street network, and you should use
    a projected graph to work in meaningful and consistent units like meters.
    Note: `tolerance` represents a per-node buffering radius. For example, to
    consolidate nodes within 10 meters of each other, use `tolerance=5`.

    When `rebuild_graph` is False, it uses a purely geometric (and relatively
    fast) algorithm to identify "geometrically close" nodes, merge them, and
    return the merged intersections' centroids. When `rebuild_graph` is True,
    it uses a topological (and slower but more accurate) algorithm to identify
    "topologically close" nodes, merge them, then rebuild/return the graph.
    Returned graph's node IDs represent clusters rather than "osmid" values.
    Refer to nodes' "osmid_original" attributes for original "osmid" values.
    If multiple nodes were merged together, the "osmid_original" attribute is
    a list of merged nodes' "osmid" values.

    Divided roads are often represented by separate centerline edges. The
    intersection of two divided roads thus creates 4 nodes, representing where
    each edge intersects a perpendicular edge. These 4 nodes represent a
    single intersection in the real world. A similar situation occurs with
    roundabouts and traffic circles. This function consolidates nearby nodes
    by buffering them to an arbitrary distance, merging overlapping buffers,
    and taking their centroid.

    Parameters
    ----------
    G
        A projected graph.
    tolerance
        Nodes are buffered to this distance (in graph's geometry's units) and
        subsequent overlaps are dissolved into a single node. If scalar, then
        that single value will be used for all nodes. If dict (mapping node
        IDs to individual values), then those values will be used per node and
        any missing node IDs will not be buffered.
    rebuild_graph
        If True, consolidate the nodes topologically, rebuild the graph, and
        return as MultiDiGraph. Otherwise, consolidate the nodes geometrically
        and return the consolidated node points as GeoSeries.
    dead_ends
        If False, discard dead-end nodes to return only street-intersection
        points.
    reconnect_edges
        If True, reconnect edges (and their geometries) to the consolidated
        nodes in rebuilt graph, and update the edge length attributes. If
        False, the returned graph has no edges (which is faster if you just
        need topologically consolidated intersection counts). Ignored if
        `rebuild_graph` is not True.
    node_attr_aggs
        Allows user to aggregate node attributes values when merging nodes.
        Keys are node attribute names and values are aggregation functions
        (anything accepted as an argument by `pandas.agg`). Node attributes
        not in `node_attr_aggs` will contain the unique values across the
        merged nodes. If None, defaults to `{"elevation": numpy.mean}`.

    Returns
    -------
    G or gs
        If `rebuild_graph=True`, returns MultiDiGraph with consolidated
        intersections and (optionally) reconnected edge geometries. If
        `rebuild_graph=False`, returns GeoSeries of Points representing the
        centroids of street intersections.
    """
    # if dead_ends is False, discard dead-ends to retain only intersections
    if not dead_ends:
        spn = stats.streets_per_node(G)
        dead_end_nodes = [node for node, count in spn.items() if count <= 1]

        # make a copy to not mutate original graph object caller passed in
        G = G.copy()
        G.remove_nodes_from(dead_end_nodes)

    if rebuild_graph:
        if len(G.nodes) == 0 or len(G.edges) == 0:
            # cannot rebuild a graph with no nodes or no edges, just return it
            return G

        # otherwise
        return _consolidate_intersections_rebuild_graph(
            G,
            tolerance,
            reconnect_edges,
            node_attr_aggs,
        )

    # otherwise, if we're not rebuilding the graph
    if len(G) == 0:
        # if graph has no nodes, just return empty GeoSeries
        return gpd.GeoSeries(crs=G.graph["crs"])

    # otherwise, return the centroids of the merged intersection polygons
    return _merge_nodes_geometric(G, tolerance).centroid


def _merge_nodes_geometric(
    G: nx.MultiDiGraph,
    tolerance: float | dict[int, float],
) -> gpd.GeoSeries:
    """
    Geometrically merge nodes within some distance of each other.

    Parameters
    ----------
    G
        A projected graph.
    tolerance
        Nodes are buffered to this distance (in graph's geometry's units) and
        subsequent overlaps are dissolved into a single node. If scalar, then
        that single value will be used for all nodes. If dict (mapping node
        IDs to individual values), then those values will be used per node and
        any missing node IDs will not be buffered.

    Returns
    -------
    merged
        The merged overlapping polygons of the buffered nodes.
    """
    gdf_nodes = convert.graph_to_gdfs(G, edges=False)

    if isinstance(tolerance, dict):
        # create series of tolerances reindexed like nodes, then buffer, then
        # fill nulls (resulting from missing tolerances) with original points,
        # then merge overlapping geometries
        tols = pd.Series(tolerance).reindex(gdf_nodes.index)
        merged = gdf_nodes.buffer(tols).fillna(gdf_nodes["geometry"]).union_all()
    else:
        # buffer nodes then merge overlapping geometries
        merged = gdf_nodes.buffer(tolerance).union_all()

    # extract the member geometries if it's a multi-geometry
    merged = merged.geoms if hasattr(merged, "geoms") else merged
    return gpd.GeoSeries(merged, crs=G.graph["crs"])


def _consolidate_intersections_rebuild_graph(  # noqa: C901,PLR0912,PLR0915
    G: nx.MultiDiGraph,
    tolerance: float | dict[int, float],
    reconnect_edges: bool,  # noqa: FBT001
    node_attr_aggs: dict[str, Any] | None,
) -> nx.MultiDiGraph:
    """
    Consolidate intersections comprising clusters of nearby nodes.

    Merge nodes and return a rebuilt graph with consolidated intersections and
    reconnected edge geometries.

    Parameters
    ----------
    G
        A projected graph.
    tolerance
        Nodes are buffered to this distance (in graph's geometry's units) and
        subsequent overlaps are dissolved into a single node. If scalar, then
        that single value will be used for all nodes. If dict (mapping node
        IDs to individual values), then those values will be used per node and
        any missing node IDs will not be buffered.
    reconnect_edges
        If True, reconnect edges (and their geometries) to the consolidated
        nodes in rebuilt graph, and update the edge length attributes. If
        False, the returned graph has no edges (which is faster if you just
        need topologically consolidated intersection counts).
    node_attr_aggs
        Allows user to aggregate node attributes values when merging nodes.
        Keys are node attribute names and values are aggregation functions
        (anything accepted as an argument by `pandas.agg`). Node attributes
        not in `node_attr_aggs` will contain the unique values across the
        merged nodes. If None, defaults to `{"elevation": "mean"}`.

    Returns
    -------
    Gc
        A rebuilt graph with consolidated intersections and (optionally)
        reconnected edge geometries.
    """
    # default node attributes to aggregate upon consolidation
    if node_attr_aggs is None:
        node_attr_aggs = {"elevation": "mean"}

    # STEP 1
    # buffer nodes to passed-in distance and merge overlaps. turn merged nodes
    # into gdf and get centroids of each cluster as x, y
    node_clusters = gpd.GeoDataFrame(geometry=_merge_nodes_geometric(G, tolerance))
    centroids = node_clusters.centroid
    node_clusters["x"] = centroids.x
    node_clusters["y"] = centroids.y

    # STEP 2
    # attach each node to its cluster of merged nodes. first get the original
    # graph's node points then spatial join to give each node the label of
    # cluster it's within. make cluster labels type string.
    node_points = convert.graph_to_gdfs(G, edges=False).drop(columns=["x", "y"])
    gdf = gpd.sjoin(node_points, node_clusters, how="left", predicate="within")
    gdf = gdf.drop(columns="geometry").rename(columns={"index_right": "cluster"})
    gdf["cluster"] = gdf["cluster"].astype(str)

    # STEP 3
    # if a cluster contains multiple components (i.e., it's not connected)
    # move each component to its own cluster (otherwise you will connect
    # nodes together that are not truly connected, e.g., nearby deadends or
    # surface streets with bridge).
    for cluster_label, nodes_subset in gdf.groupby("cluster"):
        if len(nodes_subset) > 1:
            # identify all the (weakly connected) component in cluster
            wccs = list(nx.weakly_connected_components(G.subgraph(nodes_subset.index)))
            if len(wccs) > 1:
                # if there are multiple components in this cluster
                for suffix, wcc in enumerate(wccs):
                    # set subcluster xy to the centroid of just these nodes
                    idx = list(wcc)
                    subcluster_centroid = node_points.loc[idx].union_all().centroid
                    gdf.loc[idx, "x"] = subcluster_centroid.x
                    gdf.loc[idx, "y"] = subcluster_centroid.y
                    # move to subcluster by appending suffix to cluster label
                    gdf.loc[idx, "cluster"] = f"{cluster_label}-{suffix}"

    # give nodes unique integer IDs (subclusters with suffixes are strings)
    gdf["cluster"] = gdf["cluster"].factorize()[0]

    # STEP 4
    # create new empty graph and copy over misc graph data
    Gc = nx.MultiDiGraph()
    Gc.graph = G.graph

    # STEP 5
    # create a new node for each cluster of merged nodes
    # regroup now that we potentially have new cluster labels from step 3
    groups = gdf.groupby("cluster")
    for cluster_label, nodes_subset in groups:
        osmids = nodes_subset.index.to_list()
        if len(osmids) == 1:
            # if cluster is a single node, add that node to new graph
            osmid = osmids[0]
            Gc.add_node(cluster_label, osmid_original=osmid, **G.nodes[osmid])
        else:
            # if cluster is multiple merged nodes, create one new node with
            # attributes to represent the merged nodes' non-null values
            node_attrs = {
                "osmid_original": osmids,
                "x": nodes_subset["x"].iloc[0],
                "y": nodes_subset["y"].iloc[0],
            }
            for col in set(nodes_subset.columns):
                # get the unique non-null values (we won't add null attrs)
                unique_vals = list(set(nodes_subset[col].dropna()))
                if len(unique_vals) > 0 and col in node_attr_aggs:
                    # if this attribute's values must be aggregated, do so now
                    node_attrs[col] = nodes_subset[col].agg(node_attr_aggs[col])
                elif col == "street_count":
                    # if user doesn't specifically handle street_count with an
                    # agg function, just skip it here then calculate it later
                    continue
                elif len(unique_vals) == 1:
                    # if there's 1 unique value for this attribute, keep it
                    node_attrs[col] = unique_vals[0]
                elif len(unique_vals) > 1:
                    # if there are multiple unique values, keep one of each
                    node_attrs[col] = unique_vals
            Gc.add_node(cluster_label, **node_attrs)

    if len(G.edges) == 0 or not reconnect_edges:
        # if reconnect_edges is False or there are no edges in original graph
        # (after dead-end removed), then skip edges and return new graph as-is
        return Gc

    # STEP 6
    # create new edge from cluster to cluster for each edge in original graph
    gdf_edges = convert.graph_to_gdfs(G, nodes=False)
    for u, v, k, data in G.edges(keys=True, data=True):
        u2 = gdf.loc[u, "cluster"]
        v2 = gdf.loc[v, "cluster"]

        # only create the edge if we're not connecting the cluster
        # to itself, but always add original self-loops
        if (u2 != v2) or (u == v):
            data["u_original"] = u
            data["v_original"] = v
            if "geometry" not in data:
                data["geometry"] = gdf_edges.loc[(u, v, k), "geometry"]
            Gc.add_edge(u2, v2, **data)

    # STEP 7
    # for every group of merged nodes with more than 1 node in it, extend the
    # edge geometries to reach the new node point
    for cluster_label, nodes_subset in groups:
        # but only if there were multiple nodes merged together,
        # otherwise it's the same old edge as in original graph
        if len(nodes_subset) > 1:
            # get coords of merged nodes point centroid to prepend or
            # append to the old edge geom's coords
            x = Gc.nodes[cluster_label]["x"]
            y = Gc.nodes[cluster_label]["y"]
            xy = [(x, y)]

            # for each edge incident on this new merged node, update its
            # geometry to extend to/from the new node's point coords
            in_edges = set(Gc.in_edges(cluster_label, keys=True))
            out_edges = set(Gc.out_edges(cluster_label, keys=True))
            for u, v, k in in_edges | out_edges:
                old_coords = list(Gc.edges[u, v, k]["geometry"].coords)
                new_coords = xy + old_coords if cluster_label == u else old_coords + xy
                new_geom = LineString(new_coords)
                Gc.edges[u, v, k]["geometry"] = new_geom

                # update the edge length attribute, given the new geometry
                Gc.edges[u, v, k]["length"] = new_geom.length

    # calculate street_count attribute for all nodes lacking it
    null_nodes = [n for n, sc in Gc.nodes(data="street_count") if sc is None]
    street_counts = stats.count_streets_per_node(Gc, nodes=null_nodes)
    nx.set_node_attributes(Gc, street_counts, name="street_count")

    return Gc
