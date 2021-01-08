"""Simplify, correct, and consolidate network topology."""

import logging as lg

import geopandas as gpd
import networkx as nx
import pandas as pd
from shapely.geometry import LineString
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.ops import unary_union

from . import utils
from . import utils_graph


def _is_endpoint(G, node, strict=True):
    """
    Is node a true endpoint of an edge.

    Return True if the node is a "real" endpoint of an edge in the network,
    otherwise False. OSM data includes lots of nodes that exist only as points
    to help streets bend around curves. An end point is a node that either:
    1) is its own neighbor, ie, it self-loops.
    2) or, has no incoming edges or no outgoing edges, ie, all its incident
    edges point inward or all its incident edges point outward.
    3) or, it does not have exactly two neighbors and degree of 2 or 4.
    4) or, if strict mode is false, if its edges have different OSM IDs.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    node : int
        the node to examine
    strict : bool
        if False, allow nodes to be end points even if they fail all other rules
        but have edges with different OSM IDs

    Returns
    -------
    bool
    """
    neighbors = set(list(G.predecessors(node)) + list(G.successors(node)))
    n = len(neighbors)
    d = G.degree(node)

    # rule 1
    if node in neighbors:
        # if the node appears in its list of neighbors, it self-loops
        # this is always an endpoint.
        return True

    # rule 2
    elif G.out_degree(node) == 0 or G.in_degree(node) == 0:
        # if node has no incoming edges or no outgoing edges, it is an endpoint
        return True

    # rule 3
    elif not (n == 2 and (d == 2 or d == 4)):
        # else, if it does NOT have 2 neighbors AND either 2 or 4 directed
        # edges, it is an endpoint. either it has 1 or 3+ neighbors, in which
        # case it is a dead-end or an intersection of multiple streets or it has
        # 2 neighbors but 3 degree (indicating a change from oneway to twoway)
        # or more than 4 degree (indicating a parallel edge) and thus is an
        # endpoint
        return True

    # rule 4
    elif not strict:
        # non-strict mode: do its incident edges have different OSM IDs?
        osmids = []

        # add all the edge OSM IDs for incoming edges
        for u in G.predecessors(node):
            for key in G[u][node]:
                osmids.append(G.edges[u, node, key]["osmid"])

        # add all the edge OSM IDs for outgoing edges
        for v in G.successors(node):
            for key in G[node][v]:
                osmids.append(G.edges[node, v, key]["osmid"])

        # if there is more than 1 OSM ID in the list of edge OSM IDs then it is
        # an endpoint, if not, it isn't
        return len(set(osmids)) > 1

    # if none of the preceding rules returned true, then it is not an endpoint
    else:
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
    for successor in G.successors(endpoint_successor):
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
                    else:  # pragma: no cover
                        # this can happen due to OSM digitization error where
                        # a one-way street turns into a two-way here, but
                        # duplicate incoming one-way edges are present
                        utils.log(
                            f"Unexpected simplify pattern handled near {successor}", level=lg.WARN
                        )
                        return path
                else:  # pragma: no cover
                    # if successor has >1 successors, then successor must have
                    # been an endpoint because you can go in 2 new directions.
                    # this should never occur in practice
                    raise Exception(f"Unexpected simplify pattern failed near {successor}")

            # if this successor is an endpoint, we've completed the path
            return path

    # if endpoint_successor has no successors not already in the path, return
    # the current path: this is usually due to a digitization quirk on OSM
    return path


def _get_paths_to_simplify(G, strict=True):
    """
    Generate all the paths to be simplified between endpoint nodes.

    The path is ordered from the first endpoint, through the interstitial nodes,
    to the second endpoint.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    strict : bool
        if False, allow nodes to be end points even if they fail all other rules
        but have edges with different OSM IDs

    Yields
    ------
    path_to_simplify : list
    """
    # first identify all the nodes that are endpoints
    endpoints = set([n for n in G.nodes if _is_endpoint(G, n, strict=strict)])
    utils.log(f"Identified {len(endpoints)} edge endpoints")

    # for each endpoint node, look at each of its successor nodes
    for endpoint in endpoints:
        for successor in G.successors(endpoint):
            if successor not in endpoints:
                # if endpoint node's successor is not an endpoint, build path
                # from the endpoint node, through the successor, and on to the
                # next endpoint node
                yield _build_path(G, endpoint, successor, endpoints)


def simplify_graph(G, strict=True, remove_rings=True):
    """
    Simplify a graph's topology by removing interstitial nodes.

    Simplifies graph topology by removing all nodes that are not intersections
    or dead-ends. Create an edge directly between the end points that
    encapsulate them, but retain the geometry of the original edges, saved as
    a new `geometry` attribute on the new edge. Note that only simplified
    edges receive a `geometry` attribute. Some of the resulting consolidated
    edges may comprise multiple OSM ways, and if so, their multiple attribute
    values are stored as a list.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    strict : bool
        if False, allow nodes to be end points even if they fail all other
        rules but have incident edges with different OSM IDs. Lets you keep
        nodes at elbow two-way intersections, but sometimes individual blocks
        have multiple OSM IDs within them too.
    remove_rings : bool
        if True, remove isolated self-contained rings that have no endpoints

    Returns
    -------
    G : networkx.MultiDiGraph
        topologically simplified graph, with a new `geometry` attribute on
        each simplified edge
    """
    if "simplified" in G.graph and G.graph["simplified"]:
        raise Exception("This graph has already been simplified, cannot simplify it again.")

    utils.log("Begin topologically simplifying the graph...")

    # make a copy to not mutate original graph object caller passed in
    G = G.copy()
    initial_node_count = len(G)
    initial_edge_count = len(G.edges)
    all_nodes_to_remove = []
    all_edges_to_add = []

    # generate each path that needs to be simplified
    for path in _get_paths_to_simplify(G, strict=strict):

        # add the interstitial edges we're removing to a list so we can retain
        # their spatial geometry
        edge_attributes = dict()
        for u, v in zip(path[:-1], path[1:]):

            # there should rarely be multiple edges between interstitial nodes
            # usually happens if OSM has duplicate ways digitized for just one
            # street... we will keep only one of the edges (see below)
            if G.number_of_edges(u, v) != 1:
                utils.log(f"Found multiple edges between {u} and {v} when simplifying")

            # get edge between these nodes: if multiple edges exist between
            # them (see above), we retain only one in the simplified graph
            edge = G.edges[u, v, 0]
            for key in edge:
                if key in edge_attributes:
                    # if this key already exists in the dict, append it to the
                    # value list
                    edge_attributes[key].append(edge[key])
                else:
                    # if this key doesn't already exist, set the value to a list
                    # containing the one value
                    edge_attributes[key] = [edge[key]]

        for key in edge_attributes:
            # don't touch the length attribute, we'll sum it at the end
            if len(set(edge_attributes[key])) == 1 and key != "length":
                # if there's only 1 unique value in this attribute list,
                # consolidate it to the single value (the zero-th)
                edge_attributes[key] = edge_attributes[key][0]
            elif key != "length":
                # otherwise, if there are multiple values, keep one of each value
                edge_attributes[key] = list(set(edge_attributes[key]))

        # construct the geometry and sum the lengths of the segments
        edge_attributes["geometry"] = LineString(
            [Point((G.nodes[node]["x"], G.nodes[node]["y"])) for node in path]
        )
        edge_attributes["length"] = sum(edge_attributes["length"])

        # add the nodes and edges to their lists for processing at the end
        all_nodes_to_remove.extend(path[1:-1])
        all_edges_to_add.append(
            {"origin": path[0], "destination": path[-1], "attr_dict": edge_attributes}
        )

    # for each edge to add in the list we assembled, create a new edge between
    # the origin and destination
    for edge in all_edges_to_add:
        G.add_edge(edge["origin"], edge["destination"], **edge["attr_dict"])

    # finally remove all the interstitial nodes between the new edges
    G.remove_nodes_from(set(all_nodes_to_remove))

    if remove_rings:
        # remove any connected components that form a self-contained ring
        # without any endpoints
        wccs = nx.weakly_connected_components(G)
        nodes_in_rings = set()
        for wcc in wccs:
            if not any(_is_endpoint(G, n) for n in wcc):
                nodes_in_rings.update(wcc)
        G.remove_nodes_from(nodes_in_rings)

    # mark graph as having been simplified
    G.graph["simplified"] = True

    msg = (
        f"Simplified graph: {initial_node_count} to {len(G)} nodes, "
        f"{initial_edge_count} to {len(G.edges)} edges"
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
        spn = nx.get_node_attributes(G, "street_count")
        if set(spn) != set(G.nodes):
            utils.log("Graph nodes changed since `street_count`s were calculated", level=lg.WARN)
        dead_end_nodes = [node for node, count in spn.items() if count <= 1]

        # make a copy to not mutate original graph object caller passed in
        G = G.copy()
        G.remove_nodes_from(dead_end_nodes)

    if rebuild_graph:
        if not G or not G.edges:
            # cannot rebuild a graph with no nodes or no edges, just return it
            return G
        else:
            return _consolidate_intersections_rebuild_graph(G, tolerance, reconnect_edges)

    else:
        if not G:
            # if graph has no nodes, just return empty GeoSeries
            return gpd.GeoSeries(crs=G.graph["crs"])
        else:
            # return the centroids of the merged intersection polygons
            return _merge_nodes_geometric(G, tolerance).centroid


def _merge_nodes_geometric(G, tolerance, chunk=True):
    """
    Geometrically merge nodes within some distance of each other.

    If chunk=True, it sorts the nodes GeoSeries by geometry x and y values (to
    make unary_union faster), then buffers by tolerance. Next it divides the
    nodes GeoSeries into n-sized chunks, where n = the square root of the
    number of nodes. Then it runs unary_union on each chunk, and then runs
    unary_union on the resulting unary unions. This is much faster on large
    graphs (n>100000) because of how unary_union's runtime scales with vertex
    count. But chunk=False is usually faster on small and medium sized graphs.
    This hacky method will hopefully be made obsolete when shapely becomes
    vectorized by incorporating the pygeos codebase.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        a projected graph
    tolerance : float
        nodes are buffered to this distance (in graph's geometry's units) and
        subsequent overlaps are dissolved into a single polygon
    chunk : bool
        if True, divide nodes into geometrically sorted chunks to improve the
        speed of unary_union operation by running it on each chunk and then
        running it on the results of those runs

    Returns
    -------
    merged_nodes : GeoSeries
        the merged overlapping polygons of the buffered nodes
    """
    # yield n-sized chunks of list one at a time
    def get_chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i : i + n]

    gs_nodes = utils_graph.graph_to_gdfs(G, edges=False)["geometry"]

    if not chunk:
        # buffer nodes GeoSeries then get unary union to merge overlaps
        merged_nodes = gs_nodes.buffer(tolerance).unary_union

    else:
        # sort nodes GeoSeries by geometry x/y vals to make unary_union faster
        idx = pd.DataFrame((gs_nodes.x, gs_nodes.y)).T.sort_values([0, 1]).index

        # buffer nodes by tolerance then generate n-sized chunks of nodes
        n = int(len(G) ** 0.5)
        chunks = get_chunks(gs_nodes.loc[idx].buffer(tolerance).values, n)

        # unary_union each chunk then unary_union those unary unions
        merged_nodes = unary_union([unary_union(chunk) for chunk in chunks])

    # if only a single node results, make it iterable to convert to GeoSeries
    if isinstance(merged_nodes, Polygon):
        merged_nodes = [merged_nodes]

    return gpd.GeoSeries(list(merged_nodes), crs=G.graph["crs"])


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
    # cluster it's within
    node_points = utils_graph.graph_to_gdfs(G, edges=False)[["geometry"]]
    gdf = gpd.sjoin(node_points, node_clusters, how="left", op="within")
    gdf = gdf.drop(columns="geometry").rename(columns={"index_right": "cluster"})

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
                suffix = 0
                for wcc in wccs:
                    # set subcluster xy to the centroid of just these nodes
                    subcluster_centroid = node_points.loc[wcc].unary_union.centroid
                    gdf.loc[wcc, "x"] = subcluster_centroid.x
                    gdf.loc[wcc, "y"] = subcluster_centroid.y
                    # move to subcluster by appending suffix to cluster label
                    gdf.loc[wcc, "cluster"] = f"{cluster_label}-{suffix}"
                    suffix += 1

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

    if not G.edges or not reconnect_edges:
        # if reconnect_edges is False or there are no edges in original graph
        # (after dead-end removed), then skip edges and return new graph as-is
        return H

    # STEP 6
    # create new edge from cluster to cluster for each edge in original graph
    gdf_edges = utils_graph.graph_to_gdfs(G, nodes=False)
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

            # for each edge incident to this new merged node, update its
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
