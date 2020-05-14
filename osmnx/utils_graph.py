################################################################################
# Module: utils_graph.py
# Description: Nnetwork utility functions.
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
import time
from collections import Counter
from itertools import chain
from shapely.geometry import LineString
from shapely.geometry import Point

from . import settings
from .utils import log

# scipy and sklearn are optional dependencies for faster nearest node search
try:
    from scipy.spatial import cKDTree
except ImportError as e:
    cKDTree = None
try:
    from sklearn.neighbors import BallTree
except ImportError as e:
    BallTree = None



def induce_subgraph(G, node_subset):
    """
    Induce a subgraph of G.

    Parameters
    ----------
    G : networkx multidigraph
    node_subset : list-like
        the subset of nodes to induce a subgraph of G

    Returns
    -------
    G2 : networkx multidigraph
        the subgraph of G induced by node_subset
    """

    node_subset = set(node_subset)

    # copy nodes into new graph
    G2 = G.__class__()
    G2.add_nodes_from((n, G.nodes[n]) for n in node_subset)

    # copy edges to new graph, including parallel edges
    if G2.is_multigraph:
        G2.add_edges_from((n, nbr, key, d)
            for n, nbrs in G.adj.items() if n in node_subset
            for nbr, keydict in nbrs.items() if nbr in node_subset
            for key, d in keydict.items())
    else:
        G2.add_edges_from((n, nbr, d)
            for n, nbrs in G.adj.items() if n in node_subset
            for nbr, d in nbrs.items() if nbr in node_subset)

    # update graph attribute dict, and return graph
    G2.graph.update(G.graph)
    return G2


def get_largest_component(G, strongly=False):
    """
    Return a subgraph of the largest weakly or strongly connected component
    from a directed graph.

    Parameters
    ----------
    G : networkx multidigraph
    strongly : bool
        if True, return the largest strongly instead of weakly connected
        component

    Returns
    -------
    G : networkx multidigraph
        the largest connected component subgraph from the original graph
    """

    start_time = time.time()
    original_len = len(list(G.nodes()))

    if strongly:
        # if the graph is not connected retain only the largest strongly connected component
        if not nx.is_strongly_connected(G):

            # get all the strongly connected components in graph then identify the largest
            sccs = nx.strongly_connected_components(G)
            largest_scc = max(sccs, key=len)
            G = induce_subgraph(G, largest_scc)

            msg = ('Graph was not connected, retained only the largest strongly '
                   'connected component ({:,} of {:,} total nodes) in {:.2f} seconds')
            log(msg.format(len(list(G.nodes())), original_len, time.time()-start_time))
    else:
        # if the graph is not connected retain only the largest weakly connected component
        if not nx.is_weakly_connected(G):

            # get all the weakly connected components in graph then identify the largest
            wccs = nx.weakly_connected_components(G)
            largest_wcc = max(wccs, key=len)
            G = induce_subgraph(G, largest_wcc)

            msg = ('Graph was not connected, retained only the largest weakly '
                   'connected component ({:,} of {:,} total nodes) in {:.2f} seconds')
            log(msg.format(len(list(G.nodes())), original_len, time.time()-start_time))

    return G



def get_route_edge_attributes(G, route, attribute=None, minimize_key='length', retrieve_default=None):
    """
    Get a list of attribute values for each edge in a path.

    Parameters
    ----------
    G : networkx multidigraph
    route : list
        list of nodes in the path
    attribute : string
        the name of the attribute to get the value of for each edge.
        If not specified, the complete data dict is returned for each edge.
    minimize_key : string
        if there are parallel edges between two nodes, select the one with the
        lowest value of minimize_key
    retrieve_default : Callable[Tuple[Any, Any], Any]
        Function called with the edge nodes as parameters to retrieve a default value, if the edge does not
        contain the given attribute. Per default, a `KeyError` is raised
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
    Count how many street segments emanate from each node (i.e., intersections and dead-ends) in this graph.

    If nodes is passed, then only count the nodes in the graph with those IDs.

    Parameters
    ----------
    G : networkx multidigraph
    nodes : iterable
        the set of node IDs to get counts for

    Returns
    ----------
    streets_per_node : dict
        counts of how many streets emanate from each node with keys=node id and values=count
    """

    start_time = time.time()

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
    non_self_loop_edges = [e for e in all_edges if not e[0]==e[1]]

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
    edges_flat = list(chain.from_iterable(edges))

    # count how often each node appears in the list of flattened edge endpoints
    counts = Counter(edges_flat)
    streets_per_node = {node:counts[node] for node in nodes}
    msg = ('Got the counts of undirected street segments incident to each node '
           '(before removing peripheral edges) in {:,.2f} seconds')
    log(msg.format(time.time()-start_time))
    return streets_per_node





def graph_to_gdfs(G, nodes=True, edges=True, node_geometry=True, fill_edge_geometry=True):
    """
    Convert a graph into node and/or edge GeoDataFrames

    Parameters
    ----------
    G : networkx multidigraph
    nodes : bool
        if True, convert graph nodes to a GeoDataFrame and return it
    edges : bool
        if True, convert graph edges to a GeoDataFrame and return it
    node_geometry : bool
        if True, create a geometry column from node x and y data
    fill_edge_geometry : bool
        if True, fill in missing edge geometry fields using origin and
        destination nodes

    Returns
    -------
    GeoDataFrame or tuple
        gdf_nodes or gdf_edges or both as a tuple
    """

    if not (nodes or edges):
        raise ValueError('You must request nodes or edges, or both.')

    to_return = []

    if nodes:

        start_time = time.time()

        nodes, data = zip(*G.nodes(data=True))
        gdf_nodes = gpd.GeoDataFrame(list(data), index=nodes)
        if node_geometry:
            gdf_nodes['geometry'] = gdf_nodes.apply(lambda row: Point(row['x'], row['y']), axis=1)
            gdf_nodes.set_geometry('geometry', inplace=True)
        gdf_nodes.crs = G.graph['crs']
        gdf_nodes.gdf_name = '{}_nodes'.format(G.graph['name'])

        to_return.append(gdf_nodes)
        log('Created GeoDataFrame "{}" from graph in {:,.2f} seconds'.format(gdf_nodes.gdf_name, time.time()-start_time))

    if edges:

        start_time = time.time()

        # create a list to hold our edges, then loop through each edge in the
        # graph
        edges = []
        for u, v, key, data in G.edges(keys=True, data=True):

            # for each edge, add key and all attributes in data dict to the
            # edge_details
            edge_details = {'u':u, 'v':v, 'key':key}
            for attr_key in data:
                edge_details[attr_key] = data[attr_key]

            # if edge doesn't already have a geometry attribute, create one now
            # if fill_edge_geometry==True
            if 'geometry' not in data:
                if fill_edge_geometry:
                    point_u = Point((G.nodes[u]['x'], G.nodes[u]['y']))
                    point_v = Point((G.nodes[v]['x'], G.nodes[v]['y']))
                    edge_details['geometry'] = LineString([point_u, point_v])
                else:
                    edge_details['geometry'] = np.nan

            edges.append(edge_details)

        # create a GeoDataFrame from the list of edges and set the CRS
        gdf_edges = gpd.GeoDataFrame(edges)
        gdf_edges.crs = G.graph['crs']
        gdf_edges.gdf_name = '{}_edges'.format(G.graph['name'])

        to_return.append(gdf_edges)
        log('Created GeoDataFrame "{}" from graph in {:,.2f} seconds'.format(gdf_edges.gdf_name, time.time()-start_time))

    if len(to_return) > 1:
        return tuple(to_return)
    else:
        return to_return[0]


def gdfs_to_graph(gdf_nodes, gdf_edges):
    """
    Convert node and edge GeoDataFrames into a graph

    Parameters
    ----------
    gdf_nodes : GeoDataFrame
    gdf_edges : GeoDataFrame

    Returns
    -------
    networkx multidigraph
    """

    G = nx.MultiDiGraph()
    G.graph['crs'] = gdf_nodes.crs
    G.graph['name'] = gdf_nodes.gdf_name.rstrip('_nodes')

    # add the nodes and their attributes to the graph
    G.add_nodes_from(gdf_nodes.index)
    attributes = gdf_nodes.to_dict()
    for attribute_name in gdf_nodes.columns:
        # only add this attribute to nodes which have a non-null value for it
        attribute_values = {k:v for k, v in attributes[attribute_name].items() if pd.notnull(v)}
        nx.set_node_attributes(G, name=attribute_name, values=attribute_values)

    # add the edges and attributes that are not u, v, key (as they're added
    # separately) or null
    for _, row in gdf_edges.iterrows():
        attrs = {}
        for label, value in row.iteritems():
            if (label not in ['u', 'v', 'key']) and (isinstance(value, list) or pd.notnull(value)):
                attrs[label] = value
        G.add_edge(row['u'], row['v'], key=row['key'], **attrs)

    return G
