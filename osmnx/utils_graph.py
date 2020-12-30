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
    Convert a MultiDiGraph to node and/or edge GeoDataFrames.

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

        gdf_nodes.index.rename("osmid", inplace=True)
        utils.log("Created nodes GeoDataFrame from graph")

    if edges:

        if not G.edges:
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

        # add u, v, key attributes as index
        gdf_edges["u"] = u
        gdf_edges["v"] = v
        gdf_edges["key"] = k
        gdf_edges.set_index(["u", "v", "key"], inplace=True)

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

    This function is the inverse of `graph_to_gdfs` and is designed to work in
    conjunction with it. However, you can convert arbitrary node and edge
    GeoDataFrames as long as gdf_nodes is uniquely indexed by `osmid` and
    gdf_edges is uniquely multi-indexed by `u`, `v`, `key` (following normal
    MultiDiGraph structure). This allows you to load any node/edge shapefiles
    or GeoPackage layers as GeoDataFrames then convert them to a MultiDiGraph
    for graph analysis.

    Parameters
    ----------
    gdf_nodes : geopandas.GeoDataFrame
        GeoDataFrame of graph nodes uniquely indexed by osmid
    gdf_edges : geopandas.GeoDataFrame
        GeoDataFrame of graph edges uniquely multi-indexed by u, v, key
    graph_attrs : dict
        the new G.graph attribute dict. if None, use crs from gdf_edges as the
        only graph-level attribute (gdf_edges must have crs attribute set)

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    if graph_attrs is None:
        graph_attrs = {"crs": gdf_edges.crs}
    G = nx.MultiDiGraph(**graph_attrs)

    # add edges and their attributes to graph, but filter out null attribute
    # values so that edges only get attributes with non-null values
    attr_names = gdf_edges.columns.to_list()
    for (u, v, k), attr_vals in zip(gdf_edges.index, gdf_edges.values):
        data_all = zip(attr_names, attr_vals)
        data = {name: val for name, val in data_all if isinstance(val, list) or pd.notnull(val)}
        G.add_edge(u, v, key=k, **data)

    # add nodes' attributes to graph
    for col in gdf_nodes.columns:
        nx.set_node_attributes(G, name=col, values=gdf_nodes[col].dropna())

    utils.log("Created graph from node/edge GeoDataFrames")
    return G


def add_edge_lengths(G, precision=3):
    """
    Add `length` attribute (in meters) to each edge.

    Calculated via great-circle distance between each edge's incident nodes,
    so ensure graph is in unprojected coordinates. Graph should be
    unsimplified to get accurate distances. Note: this function is run by all
    the `graph.graph_from_x` functions automatically to add `length`
    attributes to all edges.

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
    # extract the edges' endpoint nodes' coordinates
    try:
        coords = (
            (u, v, k, G.nodes[u]["y"], G.nodes[u]["x"], G.nodes[v]["y"], G.nodes[v]["x"])
            for u, v, k in G.edges
        )
    except KeyError:  # pragma: no cover
        raise KeyError("Some edges missing nodes, possibly due to input data clipping issue")

    # turn the coordinates into a DataFrame indexed by u, v, k
    cols = ["u", "v", "k", "u_y", "u_x", "v_y", "v_x"]
    df = pd.DataFrame(coords, columns=cols).set_index(["u", "v", "k"])

    # calculate great circle distances, fill nulls with zeros, then round
    dists = distance.great_circle_vec(df["u_y"], df["u_x"], df["v_y"], df["v_x"])
    dists = dists.fillna(value=0).round(precision)
    nx.set_edge_attributes(G, name="length", values=dists)

    utils.log("Added edge lengths to graph")
    return G


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

    utils.log("Counted undirected street segments incident to each node")
    return streets_per_node


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
    # make a copy to not mutate original graph object caller passed in
    G = G.copy()

    # get the set of all isolated nodes, then remove them
    isolated_nodes = {node for node, degree in G.degree() if degree < 1}
    G.remove_nodes_from(isolated_nodes)
    utils.log(f"Removed {len(isolated_nodes)} isolated nodes")
    return G


def get_largest_component(G, strongly=False):
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
        utils.log(f"Got largest {kind} connected component ({len(G)} of {n} total nodes)")

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
    # make a copy to not mutate original graph object caller passed in
    G = G.copy()
    to_remove = []

    # identify all the parallel edges in the MultiDiGraph
    parallels = ((u, v) for u, v, k in G.edges(keys=True) if k > 0)

    # remove the parallel edge with greater "weight" attribute value
    for u, v in set(parallels):
        k, _ = max(G.get_edge_data(u, v).items(), key=lambda x: x[1][weight])
        to_remove.append((u, v, k))

    G.remove_edges_from(to_remove)
    utils.log("Converted MultiDiGraph to DiGraph")

    return nx.DiGraph(G)


def get_undirected(G):
    """
    Convert MultiDiGraph to undirected MultiGraph.

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
    H = nx.MultiGraph(**G.graph)
    H.add_nodes_from(G.nodes(data=True))
    H.add_edges_from(G.edges(keys=True, data=True))

    # the previous operation added all directed edges from G as undirected
    # edges in H. we now have duplicate edges for every bidirectional parallel
    # edge or self-loop. so, look through the edges and remove any duplicates.
    duplicate_edges = set()
    for u1, v1, key1, data1 in H.edges(keys=True, data=True):

        # if we haven't already flagged this edge as a duplicate
        if (u1, v1, key1) not in duplicate_edges:

            # look at every other edge between u and v, one at a time
            for key2 in H[u1][v1]:

                # don't compare this edge to itself
                if key1 != key2:

                    # compare the first edge's data to the second's
                    # if they match up, flag the duplicate for removal
                    data2 = H.edges[u1, v1, key2]
                    if _is_duplicate_edge(data1, data2):
                        duplicate_edges.add((u1, v1, key2))

    H.remove_edges_from(duplicate_edges)
    utils.log("Converted MultiDiGraph to undirected MultiGraph")

    return H


def _is_duplicate_edge(data1, data2):
    """
    Check if two graph edge data dicts have the same osmid and geometry.

    Parameters
    ----------
    data1: dict
        the first edge's data
    data2 : dict
        the second edge's data

    Returns
    -------
    is_dupe : bool
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


def _is_same_geometry(ls1, ls2):
    """
    Determine if two LineString geometries are the same (in either direction).

    Check both the normal and reversed orders of their constituent points.

    Parameters
    ----------
    ls1 : shapely.geometry.LineString
        the first LineString geometry
    ls2 : shapely.geometry.LineString
        the second LineString geometry

    Returns
    -------
    bool
    """
    # extract coordinates from each LineString geometry
    geom1 = [tuple(coords) for coords in ls1.xy]
    geom2 = [tuple(coords) for coords in ls2.xy]

    # reverse the first LineString's coordinates' direction
    geom1_r = [tuple(reversed(coords)) for coords in ls1.xy]

    # if first geometry matches second in either direction, return True
    return geom1 == geom2 or geom1_r == geom2


def _update_edge_keys(G):
    """
    Increment key of one edge of parallel edges that differ in geometry.

    For example, two streets from u to v that bow away from each other as
    separate streets, rather than opposite direction edges of a single street.
    Increment one of these edge's keys so that they do not match across u, v,
    k or v, u, k so we can add both to an undirected MultiGraph.

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
    edges["uvk"] = ["_".join(sorted([str(u), str(v)]) + [str(k)]) for u, v, k in edges.index]
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
