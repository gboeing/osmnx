"""Graph utility functions."""

import itertools
from collections import Counter

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from shapely.geometry import LineString
from shapely.geometry import Point

from . import distance
from . import utils


def graph_to_gdfs(G, nodes=True, edges=True, node_geometry=True, fill_edge_geometry=True):
    """
    Convert a graph to node and/or edge GeoDataFrames.

    This function is the inverse of `graph_from_gdfs`.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    nodes : bool
        if True, convert graph nodes to a GeoDataFrame and return it
    edges : bool
        if True, convert graph edges to a GeoDataFrame and return it
    node_geometry : bool
        if True, create a geometry column from node x and y data
    fill_edge_geometry : bool
        if True, fill in missing edge geometry fields using nodes u and v

    Returns
    -------
    geopandas.GeoDataFrame or tuple
        gdf_nodes or gdf_edges or tuple of (gdf_nodes, gdf_edges)
    """
    crs = G.graph["crs"]

    if nodes:

        nodes, data = zip(*G.nodes(data=True))

        if node_geometry:
            # convert node x/y attributes to Points for geometry column
            geom = (Point(d["x"], d["y"]) for d in data)
            gdf_nodes = gpd.GeoDataFrame(data, index=nodes, crs=crs, geometry=list(geom))
        else:
            gdf_nodes = gpd.GeoDataFrame(data, index=nodes)

        utils.log("Created nodes GeoDataFrame from graph")

    if edges:

        if len(G.edges) < 1:
            raise ValueError("Graph has no edges, cannot convert to a GeoDataFrame.")

        u, v, k, data = zip(*G.edges(keys=True, data=True))

        if fill_edge_geometry:

            # subroutine to get geometry for every edge: if edge already has
            # geometry return it, otherwise create it using the incident nodes
            x_lookup = nx.get_node_attributes(G, "x")
            y_lookup = nx.get_node_attributes(G, "y")

            def make_geom(u, v, data, x=x_lookup, y=y_lookup):
                if "geometry" in data:
                    return data["geometry"]
                else:
                    return LineString((Point((x[u], y[u])), Point((x[v], y[v]))))

            geom = map(make_geom, u, v, data)
            gdf_edges = gpd.GeoDataFrame(data, crs=crs, geometry=list(geom))

        else:
            gdf_edges = gpd.GeoDataFrame(data)
            if "geometry" not in gdf_edges.columns:
                # if no edges have a geometry attribute, create null column
                gdf_edges["geometry"] = np.nan
            gdf_edges.set_geometry("geometry")
            gdf_edges.crs = crs

        # add u, v, key attributes as columns
        gdf_edges["u"] = u
        gdf_edges["v"] = v
        gdf_edges["key"] = k

        utils.log("Created edges GeoDataFrame from graph")

    if nodes and edges:
        return gdf_nodes, gdf_edges
    elif nodes:
        return gdf_nodes
    elif edges:
        return gdf_edges
    else:
        raise ValueError("You must request nodes or edges or both.")


def graph_from_gdfs(gdf_nodes, gdf_edges, graph_attrs=None):
    """
    Convert node and edge GeoDataFrames to a MultiDiGraph.

    This function is the inverse of `graph_to_gdfs`.

    Parameters
    ----------
    gdf_nodes : geopandas.GeoDataFrame
        GeoDataFrame of graph nodes
    gdf_edges : geopandas.GeoDataFrame
        GeoDataFrame of graph edges, must have crs attribute set
    graph_attrs : dict
        the new G.graph attribute dict; if None, add crs as the only
        graph-level attribute

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    if graph_attrs is None:
        graph_attrs = {"crs": gdf_edges.crs}
    G = nx.MultiDiGraph(**graph_attrs)

    # add the nodes then each node's non-null attributes
    G.add_nodes_from(gdf_nodes.index)
    for col in gdf_nodes.columns:
        nx.set_node_attributes(G, name=col, values=gdf_nodes[col].dropna())

    # add each edge and its non-null attributes
    for (u, v, k), row in gdf_edges.set_index(["u", "v", "key"]).iterrows():
        d = {label: val for label, val in row.items() if isinstance(val, list) or pd.notnull(val)}
        G.add_edge(u, v, k, **d)

    utils.log("Created graph from node/edge GeoDataFrames")
    return G


def induce_subgraph(G, node_subset):
    """
    Induce a subgraph of G: deprecated, do not use.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    node_subset : list-like
        the subset of nodes to induce a subgraph of G

    Returns
    -------
    networkx.MultiDiGraph
        the subgraph of G induced by node_subset
    """
    import warnings

    msg = (
        "The induce_subgraph function has been deprecated and will be "
        "removed in a future release. Use G.subgraph(nodes) instead."
    )
    warnings.warn(msg)

    # induce a (frozen) subgraph then unfreeze it by making new MultiDiGraph
    return nx.MultiDiGraph(G.subgraph(node_subset))


def get_largest_component(G, strongly=False):
    """
    Get subgraph of MultiDiGraph's largest weakly/strongly connected component.

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
        utils.log(f"Got largest {kind} connected component ({len(G)} of {n} total nodes)")

    return G


def get_route_edge_attributes(
    G, route, attribute=None, minimize_key="length", retrieve_default=None
):
    """
    Get a list of attribute values for each edge in a path.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    route : list
        list of nodes IDs constituting the path
    attribute : string
        the name of the attribute to get the value of for each edge. If None,
        the complete data dict is returned for each edge.
    minimize_key : string
        if there are parallel edges between two nodes, select the one with the
        lowest value of minimize_key
    retrieve_default : Callable[Tuple[Any, Any], Any]
        function called with the edge nodes as parameters to retrieve a
        default value, if the edge does not contain the given attribute
        (otherwise a `KeyError` is raised)

    Returns
    -------
    attribute_values : list
        list of edge attribute values
    """
    attribute_values = []
    for u, v in zip(route[:-1], route[1:]):
        # if there are parallel edges between two nodes, select the one with the
        # lowest value of minimize_key
        data = min(G.get_edge_data(u, v).values(), key=lambda x: x[minimize_key])
        if attribute is None:
            attribute_value = data
        elif retrieve_default is not None:
            attribute_value = data.get(attribute, retrieve_default(u, v))
        else:
            attribute_value = data[attribute]
        attribute_values.append(attribute_value)
    return attribute_values


def count_streets_per_node(G, nodes=None):
    """
    Count how many street segments emanate from each node in this graph.

    If nodes is passed, then only count the nodes in the graph with those IDs.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    nodes : iterable
        the set of node IDs to get counts for

    Returns
    -------
    streets_per_node : dict
        counts of how many streets emanate from each node with
        keys=node id and values=count
    """
    # to calculate the counts, get undirected representation of the graph. for
    # each node, get the list of the set of unique u,v,key edges, including
    # parallel edges but excluding self-loop parallel edges (this is necessary
    # because bi-directional self-loops will appear twice in the undirected
    # graph as you have u,v,key0 and u,v,key1 where u==v when you convert from
    # MultiDiGraph to MultiGraph - BUT, one-way self-loops will appear only
    # once. to get consistent accurate counts of physical streets, ignoring
    # directionality, we need the list of the set of unique edges...). then,
    # count how many times the node appears in the u,v tuples in the list. this
    # is the count of how many street segments emanate from this node. finally,
    # create a dict of node id:count
    G_undir = G.to_undirected(reciprocal=False)
    all_edges = G_undir.edges(keys=False)
    if nodes is None:
        nodes = G_undir.nodes()

    # get all unique edges - this throws away any parallel edges (including
    # those in self-loops)
    all_unique_edges = set(all_edges)

    # get all edges (including parallel edges) that are not self-loops
    non_self_loop_edges = [e for e in all_edges if not e[0] == e[1]]

    # get a single copy of each self-loop edge (ie, if it's bi-directional, we
    # ignore the parallel edge going the reverse direction and keep only one
    # copy)
    set_non_self_loop_edges = set(non_self_loop_edges)
    self_loop_edges = [e for e in all_unique_edges if e not in set_non_self_loop_edges]

    # final list contains all unique edges, including each parallel edge, unless
    # the parallel edge is a self-loop, in which case it doesn't double-count
    # the self-loop
    edges = non_self_loop_edges + self_loop_edges

    # flatten the list of (u,v) tuples
    edges_flat = list(itertools.chain.from_iterable(edges))

    # count how often each node appears in the list of flattened edge endpoints
    counts = Counter(edges_flat)
    streets_per_node = {node: counts[node] for node in nodes}
    utils.log("Counted undirected street segments incident to each node")
    return streets_per_node


def remove_isolated_nodes(G):
    """
    Remove from a graph all nodes that have no incident edges.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        graph from which to remove isolated nodes

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with all isolated nodes removed
    """
    isolated_nodes = [node for node, degree in G.degree() if degree < 1]
    G.remove_nodes_from(isolated_nodes)
    utils.log(f"Removed {len(isolated_nodes)} isolated nodes")
    return G


def _is_duplicate_edge(data, data_other):
    """
    Check if two edge data dicts are the same based on OSM ID and geometry.

    Parameters
    ----------
    data : dict
        the first edge's data
    data_other : dict
        the second edge's data

    Returns
    -------
    is_dupe : bool
    """
    is_dupe = False

    # if either edge's OSM ID contains multiple values (due to simplification), we want
    # to compare as sets so they are order-invariant, otherwise uv does not match vu
    osmid = set(data["osmid"]) if isinstance(data["osmid"], list) else data["osmid"]
    osmid_other = (
        set(data_other["osmid"]) if isinstance(data_other["osmid"], list) else data_other["osmid"]
    )

    if osmid == osmid_other:
        # if they contain the same OSM ID or set of OSM IDs (due to simplification)
        if ("geometry" in data) and ("geometry" in data_other):
            # if both edges have a geometry attribute
            if _is_same_geometry(data["geometry"], data_other["geometry"]):
                # if their edge geometries have the same coordinates
                is_dupe = True
        elif ("geometry" in data) and ("geometry" in data_other):
            # if neither edge has a geometry attribute
            is_dupe = True
        else:
            # if one edge has geometry attribute but the other doesn't, keep it
            pass

    return is_dupe


def _is_same_geometry(ls1, ls2):
    """
    Check if LineString geometries in two edges are the same.

    Check both normal and reversed order of constituent points.

    Parameters
    ----------
    ls1 : shapely.geometry.LineString
        the first edge's geometry
    ls2 : shapely.geometry.LineString
        the second edge's geometry

    Returns
    -------
    bool
    """
    # extract geometries from each edge data dict
    geom1 = [list(coords) for coords in ls1.xy]
    geom2 = [list(coords) for coords in ls2.xy]

    # reverse the first edge's list of x's and y's to look for a match in
    # either order
    geom1_r = [list(reversed(list(coords))) for coords in ls1.xy]

    # if the edge's geometry matches its reverse's geometry in either order,
    # return True
    return geom1 == geom2 or geom1_r == geom2


def _update_edge_keys(G):
    """
    Update keys of edges that share u, v with other edge but differ in geometry.

    For example, two one-way streets from u to v that bow away from each other
    as separate streets, rather than opposite direction edges of a single
    street.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    # identify all the edges that are duplicates based on a sorted combination
    # of their origin, destination, and key. that is, edge uv will match edge vu
    # as a duplicate, but only if they have the same key
    edges = graph_to_gdfs(G, nodes=False, fill_edge_geometry=False)
    edges["uvk"] = edges.apply(
        lambda row: "_".join(sorted([str(row["u"]), str(row["v"])]) + [str(row["key"])]), axis=1
    )
    edges["dupe"] = edges["uvk"].duplicated(keep=False)
    dupes = edges[edges["dupe"]].dropna(subset=["geometry"])

    different_streets = []
    groups = dupes[["geometry", "uvk", "u", "v", "key", "dupe"]].groupby("uvk")

    # for each set of duplicate edges
    for _, group in groups:

        # if there are more than 2 edges here, make sure to compare all
        if len(group["geometry"]) > 2:
            li = group["geometry"].tolist()
            li.append(li[0])
            geom_pairs = list(zip(li[:-1], li[1:]))
        # otherwise, just compare the first edge to the second edge
        else:
            geom_pairs = [(group["geometry"].iloc[0], group["geometry"].iloc[1])]

        # for each pair of edges to compare
        for geom1, geom2 in geom_pairs:
            # if they don't have the same geometry, flag them as different streets
            if not _is_same_geometry(geom1, geom2):
                # add edge uvk, but not edge vuk, otherwise we'll iterate both their keys
                # and they'll still duplicate each other at the end of this process
                different_streets.append(
                    (group["u"].iloc[0], group["v"].iloc[0], group["key"].iloc[0])
                )

    # for each unique different street, iterate its key + 1 so it's unique
    for u, v, k in set(different_streets):
        # filter out key if it appears in data dict as we'll pass it explicitly
        attributes = {k: v for k, v in G[u][v][k].items() if k != "key"}
        G.add_edge(u, v, key=k + 1, **attributes)
        G.remove_edge(u, v, key=k)

    return G


def get_digraph(G, weight="length"):
    """
    Convert MultiDiGraph to DiGraph.

    Chooses between parallel edges by minimizing `weight` attribute value.
    Note: see also `get_undirected` to convert MultiDiGraph to MultiGraph.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    weight : string
        attribute value to minimize when choosing between parallel edges

    Returns
    -------
    networkx.DiGraph
    """
    # make a copy to not edit the original graph object the caller passed in
    G = G.copy()
    to_remove = []

    # identify all the parallel edges in the MultiDiGraph
    parallels = ((u, v) for u, v, k in G.edges(keys=True) if k > 0)

    # remove the parallel edge with greater "weight" attribute value
    for u, v in set(parallels):
        k, d = max(G.get_edge_data(u, v).items(), key=lambda x: x[1][weight])
        to_remove.append((u, v, k))

    G.remove_edges_from(to_remove)
    return nx.DiGraph(G)


def get_undirected(G):
    """
    Convert MultiDiGraph to MultiGraph.

    Maintains parallel edges only if their geometries differ. Note: see also
    `get_digraph` to convert MultiDiGraph to DiGraph.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph

    Returns
    -------
    networkx.MultiGraph
    """
    # make a copy to not edit the original graph object the caller passed in
    G = G.copy()

    # set from/to nodes before making graph undirected
    for u, v, k, data in G.edges(keys=True, data=True):
        G.edges[u, v, k]["from"] = u
        G.edges[u, v, k]["to"] = v

        # add geometry if it doesn't already exist, to retain parallel
        # edges' distinct geometries
        if "geometry" not in data:
            point_u = Point((G.nodes[u]["x"], G.nodes[u]["y"]))
            point_v = Point((G.nodes[v]["x"], G.nodes[v]["y"]))
            data["geometry"] = LineString([point_u, point_v])

    # update edge keys so we don't retain only one edge of sets of parallel edges
    # when we convert from a multidigraph to a multigraph
    G = _update_edge_keys(G)

    # now convert multidigraph to a multigraph, retaining all edges in both
    # directions for now, as well as all graph attributes
    H = nx.MultiGraph()
    H.add_nodes_from(G.nodes(data=True))
    H.add_edges_from(G.edges(keys=True, data=True))
    H.graph = G.graph

    # the previous operation added all directed edges from G as undirected
    # edges in H. this means we have duplicate edges for every bi-directional
    # street. so, look through the edges and remove any duplicates
    duplicate_edges = []
    for u, v, key, data in H.edges(keys=True, data=True):

        # if we haven't already flagged this edge as a duplicate
        if not (u, v, key) in duplicate_edges:

            # look at every other edge between u and v, one at a time
            for key_other in H[u][v]:

                # don't compare this edge to itself
                if not key_other == key:

                    # compare the first edge's data to the second's to see if
                    # they are duplicates
                    data_other = H.edges[u, v, key_other]
                    if _is_duplicate_edge(data, data_other):

                        # if they match up, flag the duplicate for removal
                        duplicate_edges.append((u, v, key_other))

    H.remove_edges_from(duplicate_edges)
    utils.log("Converted MultiDiGraph to undirected MultiGraph")

    return H


def add_edge_lengths(G, precision=3):
    """
    Add `length` (meters) attribute to each edge.

    Calculated via great-circle distance between each edge's incident nodes,
    so ensure graph is in unprojected coordinates.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    precision : int
        decimal precision to round lengths

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with edge length attributes
    """
    # first load all the edges' origin and destination coordinates as a
    # dataframe indexed by u, v, key
    try:
        coords = np.array(
            [
                (u, v, k, G.nodes[u]["y"], G.nodes[u]["x"], G.nodes[v]["y"], G.nodes[v]["x"])
                for u, v, k in G.edges(keys=True)
            ]
        )
    except KeyError:  # pragma: no cover
        missing_nodes = {
            str(i)
            for u, v, _ in G.edges(keys=True)
            if not (G.nodes[u] or G.nodes[u])
            for i in (u, v)
            if not G.nodes[i]
        }
        missing_str = ", ".join(missing_nodes)
        raise KeyError(f"Edge(s) missing nodes {missing_str} possibly due to clipping issue")

    df_coords = pd.DataFrame(coords, columns=["u", "v", "k", "u_y", "u_x", "v_y", "v_x"])
    df_coords[["u", "v", "k"]] = df_coords[["u", "v", "k"]].astype(np.int64)
    df_coords = df_coords.set_index(["u", "v", "k"])

    # then calculate the great circle distance with the vectorized function
    gc_distances = distance.great_circle_vec(
        lat1=df_coords["u_y"], lng1=df_coords["u_x"], lat2=df_coords["v_y"], lng2=df_coords["v_x"]
    )

    # fill nulls with zeros and round
    gc_distances = gc_distances.fillna(value=0).round(precision)
    nx.set_edge_attributes(G, name="length", values=gc_distances.to_dict())

    utils.log("Added edge lengths to graph")
    return G
