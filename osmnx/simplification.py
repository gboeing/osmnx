"""Simplify, correct, and consolidate network topology."""

import logging as lg
from warnings import warn

import geopandas as gpd
import networkx as nx
from shapely.geometry import LineString
from shapely.geometry import MultiPolygon
from shapely.geometry import Point
from shapely.geometry import Polygon

from . import convert
from . import stats
from . import utils
from ._errors import GraphSimplificationError


def _is_endpoint(G, node, endpoint_attrs):
    """
    Determine if a node is a true endpoint of an edge.

    Return True if the node is a "true" endpoint of an edge in the network,
    otherwise False. OpenStreetMap data includes many nodes that exist only as
    geometric vertices to allow ways to curve. A true edge endpoint is a node
    that satisfies at least 1 of the following 4 rules:

    1) It is its own neighbor (ie, it self-loops).

    2) Or, it has no incoming edges or no outgoing edges (ie, all its incident
    edges are inbound or all its incident edges are outbound).

    3) Or, it does not have exactly two neighbors and degree of 2 or 4.

    4) Or, if `endpoint_attrs` is not None, and its incident edges have
    different values than each other for any of the edge attributes in
    `endpoint_attrs`.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    node : int
        the node to examine
    endpoint_attrs : iterable
        An iterable of edge attribute names for relaxing the strictness of
        endpoint determination. If not None, a node is an endpoint if its
        incident edges have different values then each other for any of the
        edge attributes in `endpoint_attrs`.

    Returns
    -------
    bool
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
    # non-strict mode: do its incident edges have different attr values? for
    # each attribute to check, collect the attribute's values in all inbound
    # and outbound edges. if there is more than 1 unique value then then this
    # node is an endpoint
    if endpoint_attrs is not None:
        for attr in endpoint_attrs:
            in_values = {v for _, _, v in G.in_edges(node, data=attr, keys=False)}
            out_values = {v for _, _, v in G.out_edges(node, data=attr, keys=False)}
            if len(in_values | out_values) > 1:
                return True

    # if none of the preceding rules passed, then it is not an endpoint
    return False


def _build_path(G, endpoint, endpoint_successor, endpoints):
    """
    Build a path of nodes from one endpoint node to next endpoint node.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    endpoint : int
        the endpoint node from which to start the path
    endpoint_successor : int
        the successor of endpoint through which the path to the next endpoint
        will be built
    endpoints : set
        the set of all nodes in the graph that are endpoints

    Returns
    -------
    path : list
        the first and last items in the resulting path list are endpoint
        nodes, and all other items are interstitial nodes that can be removed
        subsequently
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
                        return path + [endpoint]

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
                    msg = f"Impossible simplify pattern failed near {successor}"
                    raise GraphSimplificationError(msg)

            # if this successor is an endpoint, we've completed the path
            return path

    # if endpoint_successor has no successors not already in the path, return
    # the current path: this is usually due to a digitization quirk on OSM
    return path


def _get_paths_to_simplify(G, endpoint_attrs):
    """
    Generate all the paths to be simplified between endpoint nodes.

    The path is ordered from the first endpoint, through the interstitial nodes,
    to the second endpoint.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    endpoint_attrs : iterable
        An iterable of edge attribute names for relaxing the strictness of
        endpoint determination. If not None, a node is an endpoint if its
        incident edges have different values then each other for any of the
        edge attributes in `endpoint_attrs`.

    Yields
    ------
    path_to_simplify : list
        a generator of paths to simplify
    """
    # first identify all the nodes that are endpoints
    endpoints = {n for n in G.nodes if _is_endpoint(G, n, endpoint_attrs)}
    utils.log(f"Identified {len(endpoints):,} edge endpoints")

    # for each endpoint node, look at each of its successor nodes
    for endpoint in endpoints:
        for successor in G.successors(endpoint):
            if successor not in endpoints:
                # if endpoint node's successor is not an endpoint, build path
                # from the endpoint node, through the successor, and on to the
                # next endpoint node
                yield _build_path(G, endpoint, successor, endpoints)


def _remove_rings(G, endpoint_attrs):
    """
    Remove all self-contained rings from a graph.

    This identifies any connected components that form a self-contained ring
    without any endpoints, and removes them from the graph.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    endpoint_attrs : iterable
        An iterable of edge attribute names for relaxing the strictness of
        endpoint determination. If not None, a node is an endpoint if its
        incident edges have different values then each other for any of the
        edge attributes in `endpoint_attrs`.

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with self-contained rings removed
    """
    nodes_in_rings = set()
    for wcc in nx.weakly_connected_components(G):
        if not any(_is_endpoint(G, n, endpoint_attrs) for n in wcc):
            nodes_in_rings.update(wcc)
    G.remove_nodes_from(nodes_in_rings)
    return G


def simplify_graph(  # noqa: C901
    G,
    strict=None,
    edge_attrs_differ=None,
    endpoint_attrs=None,
    remove_rings=True,
    track_merged=False,
):
    """
    Simplify a graph's topology by removing interstitial nodes.

    This simplifies graph topology by removing all nodes that are not
    intersections or dead-ends, by creating an edge directly between the end
    points that encapsulate them while retaining the full geometry of the
    original edges, saved as a new `geometry` attribute on the new edge.

    Note that only simplified edges receive a `geometry` attribute. Some of
    the resulting consolidated edges may comprise multiple OSM ways, and if
    so, their multiple attribute values are stored as a list. Optionally, the
    simplified edges can receive a `merged_edges` attribute that contains a
    list of all the (u, v) node pairs that were merged together.

    Use the `edge_attrs_differ` parameter to relax simplification strictness. For
    example, `edge_attrs_differ=['osmid']` will retain every node whose incident
    edges have different OSM IDs. This lets you keep nodes at elbow two-way
    intersections (but be aware that sometimes individual blocks have multiple
    OSM IDs within them too). You could also use this parameter to retain
    nodes where sidewalks or bike lanes begin/end in the middle of a block.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    strict : bool
        deprecated, do not use
    edge_attrs_differ : iterable
        An iterable of edge attribute names for relaxing the strictness of
        endpoint determination. If not None, a node is an endpoint if its
        incident edges have different values then each other for any of the
        edge attributes in `edge_attrs_differ`.
    endpoint_attrs : iterable
        deprecated, do not use
    remove_rings : bool
        if True, remove isolated self-contained rings that have no endpoints
    track_merged : bool
        if True, add `merged_edges` attribute on simplified edges, containing
        a list of all the (u, v) node pairs that were merged together

    Returns
    -------
    G : networkx.MultiDiGraph
        topologically simplified graph, with a new `geometry` attribute on
        each simplified edge
    """
    if endpoint_attrs is not None:
        msg = (
            "The `endpoint_attrs` parameter has been deprecated and will be removed "
            "in the v2.0.0 release. Use the `edge_attrs_differ` parameter instead. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
        )
        warn(msg, FutureWarning, stacklevel=2)
        edge_attrs_differ = endpoint_attrs

    if strict is not None:
        msg = (
            "The `strict` parameter has been deprecated and will be removed in "
            "the v2.0.0 release. Use the `edge_attrs_differ` parameter instead to "
            "relax simplification strictness. For example, `edge_attrs_differ=None` "
            "reproduces the old `strict=True` behvavior and `edge_attrs_differ=['osmid']` "
            "reproduces the old `strict=False` behavior. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
        )
        warn(msg, FutureWarning, stacklevel=2)
        # maintain old behavior if strict is passed during deprecation
        edge_attrs_differ = None if strict else ["osmid"]

    if "simplified" in G.graph and G.graph["simplified"]:  # pragma: no cover
        msg = "This graph has already been simplified, cannot simplify it again."
        raise GraphSimplificationError(msg)

    utils.log("Begin topologically simplifying the graph...")

    # define edge segment attributes to sum upon edge simplification
    attrs_to_sum = {"length", "travel_time"}

    # make a copy to not mutate original graph object caller passed in
    G = G.copy()
    initial_node_count = len(G)
    initial_edge_count = len(G.edges)
    all_nodes_to_remove = []
    all_edges_to_add = []

    # generate each path that needs to be simplified
    for path in _get_paths_to_simplify(G, edge_attrs_differ):
        # add the interstitial edges we're removing to a list so we can retain
        # their spatial geometry
        merged_edges = []
        path_attributes = {}
        for u, v in zip(path[:-1], path[1:]):
            if track_merged:
                # keep track of the edges that were merged
                merged_edges.append((u, v))

            # there should rarely be multiple edges between interstitial nodes
            # usually happens if OSM has duplicate ways digitized for just one
            # street... we will keep only one of the edges (see below)
            edge_count = G.number_of_edges(u, v)
            if edge_count != 1:
                utils.log(f"Found {edge_count} edges between {u} and {v} when simplifying")

            # get edge between these nodes: if multiple edges exist between
            # them (see above), we retain only one in the simplified graph
            # We can't assume that there exists an edge from u to v
            # with key=0, so we get a list of all edges from u to v
            # and just take the first one.
            edge_data = list(G.get_edge_data(u, v).values())[0]
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
        for attr in path_attributes:
            if attr in attrs_to_sum:
                # if this attribute must be summed, sum it now
                path_attributes[attr] = sum(path_attributes[attr])
            elif len(set(path_attributes[attr])) == 1:
                # if there's only 1 unique value in this attribute list,
                # consolidate it to the single value (the zero-th):
                path_attributes[attr] = path_attributes[attr][0]
            else:
                # otherwise, if there are multiple values, keep one of each
                path_attributes[attr] = list(set(path_attributes[attr]))

        # construct the new consolidated edge's geometry for this path
        path_attributes["geometry"] = LineString(
            [Point((G.nodes[node]["x"], G.nodes[node]["y"])) for node in path]
        )

        if track_merged:
            # add the merged edges as a new attribute of the simplified edge
            path_attributes["merged_edges"] = merged_edges

        # add the nodes and edge to their lists for processing at the end
        all_nodes_to_remove.extend(path[1:-1])
        all_edges_to_add.append(
            {"origin": path[0], "destination": path[-1], "attr_dict": path_attributes}
        )

    # for each edge to add in the list we assembled, create a new edge between
    # the origin and destination
    for edge in all_edges_to_add:
        G.add_edge(edge["origin"], edge["destination"], **edge["attr_dict"])

    # finally remove all the interstitial nodes between the new edges
    G.remove_nodes_from(set(all_nodes_to_remove))

    if remove_rings:
        G = _remove_rings(G, edge_attrs_differ)

    # mark the graph as having been simplified
    G.graph["simplified"] = True

    msg = (
        f"Simplified graph: {initial_node_count:,} to {len(G):,} nodes, "
        f"{initial_edge_count:,} to {len(G.edges):,} edges"
    )
    utils.log(msg)
    return G


def consolidate_intersections(
    G, tolerance=10, rebuild_graph=True, dead_ends=False, reconnect_edges=True
):
    """
    Consolidate intersections comprising clusters of nearby nodes.

    Merges nearby nodes and returns either their centroids or a rebuilt graph
    with consolidated intersections and reconnected edge geometries. The
    tolerance argument should be adjusted to approximately match street design
    standards in the specific street network, and you should always use a
    projected graph to work in meaningful and consistent units like meters.
    Note the tolerance represents a per-node buffering radius: for example, to
    consolidate nodes within 10 meters of each other, use tolerance=5.

    When rebuild_graph=False, it uses a purely geometrical (and relatively
    fast) algorithm to identify "geometrically close" nodes, merge them, and
    return just the merged intersections' centroids. When rebuild_graph=True,
    it uses a topological (and slower but more accurate) algorithm to identify
    "topologically close" nodes, merge them, then rebuild/return the graph.
    Returned graph's node IDs represent clusters rather than osmids. Refer to
    nodes' osmid_original attributes for original osmids. If multiple nodes
    were merged together, the osmid_original attribute is a list of merged
    nodes' osmids.

    Divided roads are often represented by separate centerline edges. The
    intersection of two divided roads thus creates 4 nodes, representing where
    each edge intersects a perpendicular edge. These 4 nodes represent a
    single intersection in the real world. A similar situation occurs with
    roundabouts and traffic circles. This function consolidates nearby nodes
    by buffering them to an arbitrary distance, merging overlapping buffers,
    and taking their centroid.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        a projected graph
    tolerance : float
        nodes are buffered to this distance (in graph's geometry's units) and
        subsequent overlaps are dissolved into a single node
    rebuild_graph : bool
        if True, consolidate the nodes topologically, rebuild the graph, and
        return as networkx.MultiDiGraph. if False, consolidate the nodes
        geometrically and return the consolidated node points as
        geopandas.GeoSeries
    dead_ends : bool
        if False, discard dead-end nodes to return only street-intersection
        points
    reconnect_edges : bool
        ignored if rebuild_graph is not True. if True, reconnect edges and
        their geometries in rebuilt graph to the consolidated nodes and update
        edge length attributes; if False, returned graph has no edges (which
        is faster if you just need topologically consolidated intersection
        counts).

    Returns
    -------
    networkx.MultiDiGraph or geopandas.GeoSeries
        if rebuild_graph=True, returns MultiDiGraph with consolidated
        intersections and reconnected edge geometries. if rebuild_graph=False,
        returns GeoSeries of shapely Points representing the centroids of
        street intersections
    """
    # if dead_ends is False, discard dead-ends to retain only intersections
    if not dead_ends:
        spn = stats.streets_per_node(G)
        dead_end_nodes = [node for node, count in spn.items() if count <= 1]

        # make a copy to not mutate original graph object caller passed in
        G = G.copy()
        G.remove_nodes_from(dead_end_nodes)

    if rebuild_graph:
        if not G or not G.edges:
            # cannot rebuild a graph with no nodes or no edges, just return it
            return G

        # otherwise
        return _consolidate_intersections_rebuild_graph(G, tolerance, reconnect_edges)

    # otherwise, if we're not rebuilding the graph
    if not G:
        # if graph has no nodes, just return empty GeoSeries
        return gpd.GeoSeries(crs=G.graph["crs"])

    # otherwise, return the centroids of the merged intersection polygons
    return _merge_nodes_geometric(G, tolerance).centroid


def _merge_nodes_geometric(G, tolerance):
    """
    Geometrically merge nodes within some distance of each other.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        a projected graph
    tolerance : float
        buffer nodes to this distance (in graph's geometry's units) then merge
        overlapping polygons into a single polygon via a unary union operation

    Returns
    -------
    merged : GeoSeries
        the merged overlapping polygons of the buffered nodes
    """
    # buffer nodes GeoSeries then get unary union to merge overlaps
    merged = convert.graph_to_gdfs(G, edges=False)["geometry"].buffer(tolerance).unary_union

    # if only a single node results, make it iterable to convert to GeoSeries
    merged = MultiPolygon([merged]) if isinstance(merged, Polygon) else merged
    return gpd.GeoSeries(merged.geoms, crs=G.graph["crs"])


def _consolidate_intersections_rebuild_graph(G, tolerance=10, reconnect_edges=True):
    """
    Consolidate intersections comprising clusters of nearby nodes.

    Merge nodes and return a rebuilt graph with consolidated intersections and
    reconnected edge geometries.

    The tolerance argument should be adjusted to approximately match street
    design standards in the specific street network, and you should always use
    a projected graph to work in meaningful and consistent units like meters.

    Returned graph's node IDs represent clusters rather than osmids. Refer to
    nodes' osmid_original attributes for original osmids. If multiple nodes
    were merged together, the osmid_original attribute is a list of merged
    nodes' osmids.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        a projected graph
    tolerance : float
        nodes are buffered to this distance (in graph's geometry's units) and
        subsequent overlaps are dissolved into a single node
    reconnect_edges : bool
        ignored if rebuild_graph is not True. if True, reconnect edges and
        their geometries in rebuilt graph to the consolidated nodes and update
        edge length attributes; if False, returned graph has no edges (which
        is faster if you just need topologically consolidated intersection
        counts).

    Returns
    -------
    H : networkx.MultiDiGraph
        a rebuilt graph with consolidated intersections and reconnected
        edge geometries
    """
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
    node_points = convert.graph_to_gdfs(G, edges=False)[["geometry"]]
    gdf = gpd.sjoin(node_points, node_clusters, how="left", predicate="within")
    gdf = gdf.drop(columns="geometry").rename(columns={"index_right": "cluster"})
    gdf["cluster"] = gdf["cluster"].astype(str)

    # STEP 3
    # if a cluster contains multiple components (i.e., it's not connected)
    # move each component to its own cluster (otherwise you will connect
    # nodes together that are not truly connected, e.g., nearby deadends or
    # surface streets with bridge).
    groups = gdf.groupby("cluster")
    for cluster_label, nodes_subset in groups:
        if len(nodes_subset) > 1:
            # identify all the (weakly connected) component in cluster
            wccs = list(nx.weakly_connected_components(G.subgraph(nodes_subset.index)))
            if len(wccs) > 1:
                # if there are multiple components in this cluster
                for suffix, wcc in enumerate(wccs):
                    # set subcluster xy to the centroid of just these nodes
                    idx = list(wcc)
                    subcluster_centroid = node_points.loc[idx].unary_union.centroid
                    gdf.loc[idx, "x"] = subcluster_centroid.x
                    gdf.loc[idx, "y"] = subcluster_centroid.y
                    # move to subcluster by appending suffix to cluster label
                    gdf.loc[idx, "cluster"] = f"{cluster_label}-{suffix}"

    # give nodes unique integer IDs (subclusters with suffixes are strings)
    gdf["cluster"] = gdf["cluster"].factorize()[0]

    # STEP 4
    # create new empty graph and copy over misc graph data
    H = nx.MultiDiGraph()
    H.graph = G.graph

    # STEP 5
    # create a new node for each cluster of merged nodes
    # regroup now that we potentially have new cluster labels from step 3
    groups = gdf.groupby("cluster")
    for cluster_label, nodes_subset in groups:
        osmids = nodes_subset.index.to_list()
        if len(osmids) == 1:
            # if cluster is a single node, add that node to new graph
            osmid = osmids[0]
            H.add_node(cluster_label, osmid_original=osmid, **G.nodes[osmid])
        else:
            # if cluster is multiple merged nodes, create one new node to
            # represent them
            H.add_node(
                cluster_label,
                osmid_original=str(osmids),
                x=nodes_subset["x"].iloc[0],
                y=nodes_subset["y"].iloc[0],
            )

    # calculate street_count attribute for all nodes lacking it
    null_nodes = [n for n, sc in H.nodes(data="street_count") if sc is None]
    street_count = stats.count_streets_per_node(H, nodes=null_nodes)
    nx.set_node_attributes(H, street_count, name="street_count")

    if not G.edges or not reconnect_edges:
        # if reconnect_edges is False or there are no edges in original graph
        # (after dead-end removed), then skip edges and return new graph as-is
        return H

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
            H.add_edge(u2, v2, **data)

    # STEP 7
    # for every group of merged nodes with more than 1 node in it, extend the
    # edge geometries to reach the new node point
    for cluster_label, nodes_subset in groups:
        # but only if there were multiple nodes merged together,
        # otherwise it's the same old edge as in original graph
        if len(nodes_subset) > 1:
            # get coords of merged nodes point centroid to prepend or
            # append to the old edge geom's coords
            x = H.nodes[cluster_label]["x"]
            y = H.nodes[cluster_label]["y"]
            xy = [(x, y)]

            # for each edge incident on this new merged node, update its
            # geometry to extend to/from the new node's point coords
            in_edges = set(H.in_edges(cluster_label, keys=True))
            out_edges = set(H.out_edges(cluster_label, keys=True))
            for u, v, k in in_edges | out_edges:
                old_coords = list(H.edges[u, v, k]["geometry"].coords)
                new_coords = xy + old_coords if cluster_label == u else old_coords + xy
                new_geom = LineString(new_coords)
                H.edges[u, v, k]["geometry"] = new_geom

                # update the edge length attribute, given the new geometry
                H.edges[u, v, k]["length"] = new_geom.length

    return H
