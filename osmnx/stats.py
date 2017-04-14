###################################################################################################
# Module: stats.py
# Description: Calculate graph-theoretic topological and metric measures for a network
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
###################################################################################################

from __future__ import division
from itertools import chain
from collections import Counter
import time
import warnings
import networkx as nx
import numpy as np
import pandas as pd

from .utils import log, get_largest_component, great_circle_vec


def basic_stats(G, area=None):
    """
    Calculate basic descriptive stats and metrics for a graph.

    Parameters
    ----------
    G : networkx multidigraph
    area : numeric
        the area covered by the street network, in square meters (typically land area); if none, will skip all density-based metrics

    Returns
    -------
    stats : dict
        dictionary of network measures containing the following elements:

          - n = number of nodes in the graph
          - m = number of edges in the graph
          - k_avg = average node degree of the graph
          - count_intersections = number of intersections in graph, that is, nodes with >1 street emanating from them
          - streets_per_node_avg = how many streets (edges in the undirected representation of the graph) emanate from each node (ie, intersection or dead-end) on average (mean)
          - streets_per_node_counts = dict, with keys of number of streets emanating from the node, and values of number of nodes with this count
          - streets_per_node_proportion = dict, same as previous, but as a proportion of the total, rather than counts
          - edge_length_total = sum of all edge lengths in the graph, in meters
          - edge_length_avg = mean edge length in the graph, in meters
          - street_length_total = sum of all edges in the undirected representation of the graph
          - street_length_avg = mean edge length in the undirected representation of the graph, in meters
          - street_segments_count = number of edges in the undirected representation of the graph
          - node_density_km = n divided by area in square kilometers
          - intersection_density_km = count_intersections divided by area in square kilometers
          - edge_density_km = edge_length_total divided by area in square kilometers
          - street_density_km = street_length_total divided by area in square kilometers
          - circuity_avg = edge_length_total divided by the sum of the great circle distances between the nodes of each edge
          - self_loop_proportion = proportion of edges that have a single node as its two endpoints (ie, the edge links nodes u and v, and u==v)

    """

    sq_m_in_sq_km = 1e6 #there are 1 million sq meters in 1 sq km
    G_undirected = None

    # calculate the number of nodes, n, and the number of edges, m, in the graph
    n = len(list(G.nodes()))
    m = len(list(G.edges()))

    # calculate the average degree of the graph
    k_avg = 2 * m / n

    if 'streets_per_node' in G.graph:
        # get the degrees saved as a graph attribute (from an undirected representation of the graph)
        # this is not the degree of the nodes in the directed graph, but rather represents the number of streets (unidirected edges) emanating from each node. see count_streets_per_node function.
        streets_per_node = G.graph['streets_per_node']
    else:
        # count how many street segments emanate from each node in this graph
        streets_per_node = count_streets_per_node(G)

    # count number of intersections in graph, as nodes with >1 street emanating from them
    node_ids = set(G.nodes())
    count_intersections = len([True for node, count in streets_per_node.items() if (count > 1) and (node in node_ids)])

    # calculate the average number of streets (unidirected edges) incident to each node
    streets_per_node_avg = sum(streets_per_node.values()) / n

    # create a dict where key = number of streets (unidirected edges) incident to each node, and value = how many nodes are of this number in the graph
    streets_per_node_counts = {num:list(streets_per_node.values()).count(num) for num in range(max(streets_per_node.values()) + 1)}

    # degree proportions: dict where key = each degree and value = what proportion of nodes are of this degree in the graph
    streets_per_node_proportion = {num:count/n for num, count in streets_per_node_counts.items()}

    # calculate the total and average edge lengths
    edge_length_total = sum([d['length'] for u, v, d in G.edges(data=True)])
    edge_length_avg = edge_length_total / m

    # calculate the total and average street segment lengths (so, edges without double-counting two-way streets)
    if G_undirected is None:
        G_undirected = G.to_undirected(reciprocal=False)
    street_length_total = sum([d['length'] for u, v, d in G_undirected.edges(data=True)])
    street_segments_count = len(list(G_undirected.edges(keys=True)))
    street_length_avg = street_length_total / street_segments_count

    # we can calculate density metrics only if area is not null
    if not area is None:
        area_km = area / sq_m_in_sq_km

        # calculate node density as nodes per sq km
        node_density_km = n / area_km

        # calculate intersection density as nodes with >1 street emanating from them, per sq km
        intersection_density_km = count_intersections / area_km

        # calculate edge density as linear meters per sq km
        edge_density_km = edge_length_total / area_km

        # calculate street density as linear meters per sq km
        street_density_km = street_length_total / area_km
    else:
        # if area is None, then we cannot calculate density
        node_density_km = None
        intersection_density_km = None
        edge_density_km = None
        street_density_km = None

    # average circuity: sum of edge lengths divided by sum of great circle distance between edge endpoints
    # first load all the edges origin and destination coordinates as a dataframe, then calculate the great circle distance with the vectorized function
    coords = np.array([[G.node[u]['y'], G.node[u]['x'], G.node[v]['y'], G.node[v]['x']] for u, v, k in G.edges(keys=True)])
    df_coords = pd.DataFrame(coords, columns=['u_y', 'u_x', 'v_y', 'v_x'])

    # ignore warnings during this calculation because numpy warns it cannot calculate arccos for self-loops since u==v
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        gc_distances = great_circle_vec(lat1=df_coords['u_y'],
                                        lng1=df_coords['u_x'],
                                        lat2=df_coords['v_y'],
                                        lng2=df_coords['v_x'])
    gc_distances = gc_distances.fillna(value=0)
    try:
        circuity_avg = edge_length_total / gc_distances.sum()
    except ZeroDivisionError:
        circuity_avg = np.nan

    # percent of edges that are self-loops, ie both endpoints are the same node
    self_loops = [True for u, v, k in G.edges(keys=True) if u == v]
    self_loops_count = len(self_loops)
    self_loop_proportion = self_loops_count / m

    # assemble the results
    stats = {'n':n,
             'm':m,
             'k_avg':k_avg,
             'count_intersections':count_intersections,
             'streets_per_node_avg':streets_per_node_avg,
             'streets_per_node_counts':streets_per_node_counts,
             'streets_per_node_proportion':streets_per_node_proportion,
             'edge_length_total':edge_length_total,
             'edge_length_avg':edge_length_avg,
             'street_length_total':street_length_total,
             'street_length_avg':street_length_avg,
             'street_segments_count':street_segments_count,
             'node_density_km':node_density_km,
             'intersection_density_km':intersection_density_km,
             'edge_density_km':edge_density_km,
             'street_density_km':street_density_km,
             'circuity_avg':circuity_avg,
             'self_loop_proportion':self_loop_proportion}

    # return the results
    return stats


def extended_stats(G, connectivity=False, anc=False, ecc=False, bc=False, cc=False):
    """
    Calculate extended topological stats and metrics for a graph.

    Global topological analysis of large complex networks is extremely
    time consuming and may exhaust computer memory. Consider using function
    arguments to not run metrics that require computation of
    a full matrix of paths if they will not be needed.

    Parameters
    ----------
    G : networkx multidigraph
    connectivity : bool
        if True, calculate node and edge connectivity
    anc : bool, if True
        calculate average node connectivity
    ecc : bool
        if True, calculate shortest paths, eccentricity, and topological metrics that use eccentricity
    bc : bool
        if True, calculate node betweenness centrality
    cc : bool
        if True, calculate node closeness centrality

    Returns
    -------
    stats : dict
        dictionary of network measures containing the following elements (some only calculated/returned optionally, based on passed parameters):

          - avg_neighbor_degree
          - avg_neighbor_degree_avg
          - avg_weighted_neighbor_degree
          - avg_weighted_neighbor_degree_avg
          - degree_centrality
          - degree_centrality_avg
          - clustering_coefficient
          - clustering_coefficient_avg
          - clustering_coefficient_weighted
          - clustering_coefficient_weighted_avg
          - pagerank
          - pagerank_max_node
          - pagerank_max
          - pagerank_min_node
          - pagerank_min
          - node_connectivity
          - node_connectivity_avg
          - edge_connectivity
          - eccentricity
          - diameter
          - radius
          - center
          - periphery
          - closeness_centrality
          - closeness_centrality_avg
          - betweenness_centrality
          - betweenness_centrality_avg

    """

    stats = {}
    full_start_time = time.time()

    # create a DiGraph from the MultiDiGraph, for those metrics that require it
    G_dir = nx.DiGraph(G)

    # create an undirected Graph from the MultiDiGraph, for those metrics that require it
    G_undir = nx.Graph(G)

    # get the largest strongly connected component, for those metrics that require strongly connected graphs
    G_strong = get_largest_component(G, strongly=True)

    # average degree of the neighborhood of each node, and average for the graph
    avg_neighbor_degree = nx.average_neighbor_degree(G)
    stats['avg_neighbor_degree'] = avg_neighbor_degree
    stats['avg_neighbor_degree_avg'] = sum(avg_neighbor_degree.values())/len(avg_neighbor_degree)

    # average weighted degree of the neighborhood of each node, and average for the graph
    avg_weighted_neighbor_degree = nx.average_neighbor_degree(G, weight='length')
    stats['avg_weighted_neighbor_degree'] = avg_weighted_neighbor_degree
    stats['avg_weighted_neighbor_degree_avg'] = sum(avg_weighted_neighbor_degree.values())/len(avg_weighted_neighbor_degree)

    # degree centrality for a node is the fraction of nodes it is connected to
    degree_centrality = nx.degree_centrality(G)
    stats['degree_centrality'] = degree_centrality
    stats['degree_centrality_avg'] = sum(degree_centrality.values())/len(degree_centrality)

    # calculate clustering coefficient for the nodes
    stats['clustering_coefficient'] = nx.clustering(G_undir)

    # average clustering coefficient for the graph
    stats['clustering_coefficient_avg'] = nx.average_clustering(G_undir)

    # calculate weighted clustering coefficient for the nodes
    stats['clustering_coefficient_weighted'] = nx.clustering(G_undir, weight='length')

    # average clustering coefficient (weighted) for the graph
    stats['clustering_coefficient_weighted_avg'] = nx.average_clustering(G_undir, weight='length')

    # pagerank: a ranking of the nodes in the graph based on the structure of the incoming links
    pagerank = nx.pagerank(G_dir, weight='length')
    stats['pagerank'] = pagerank

    # node with the highest page rank, and its value
    pagerank_max_node = max(pagerank, key=lambda x: pagerank[x])
    stats['pagerank_max_node'] = pagerank_max_node
    stats['pagerank_max'] = pagerank[pagerank_max_node]

    # node with the lowest page rank, and its value
    pagerank_min_node = min(pagerank, key=lambda x: pagerank[x])
    stats['pagerank_min_node'] = pagerank_min_node
    stats['pagerank_min'] = pagerank[pagerank_min_node]

    # if True, calculate node and edge connectivity
    if connectivity:
        start_time = time.time()

        # node connectivity is the minimum number of nodes that must be removed to disconnect G or render it trivial
        stats['node_connectivity'] = nx.node_connectivity(G_strong)

        # edge connectivity is equal to the minimum number of edges that must be removed to disconnect G or render it trivial
        stats['edge_connectivity'] = nx.edge_connectivity(G_strong)
        log('Calculated node and edge connectivity in {:,.2f} seconds'.format(time.time() - start_time))

    # if True, calculate average node connectivity
    if anc:
        # mean number of internally node-disjoint paths between each pair of nodes in G
        # i.e., the expected number of nodes that must be removed to disconnect a randomly selected pair of non-adjacent nodes
        start_time = time.time()
        stats['node_connectivity_avg'] = nx.average_node_connectivity(G)
        log('Calculated average node connectivity in {:,.2f} seconds'.format(time.time() - start_time))

    # if True, calculate shortest paths, eccentricity, and topological metrics that use eccentricity
    if ecc:
        # precompute shortest paths between all nodes for eccentricity-based stats
        start_time = time.time()
        sp = {source:dict(nx.single_source_dijkstra_path_length(G_strong, source, weight='length')) for source in G_strong.nodes()}

        log('Calculated shortest path lengths in {:,.2f} seconds'.format(time.time() - start_time))

        # eccentricity of a node v is the maximum distance from v to all other nodes in G
        eccentricity = nx.eccentricity(G_strong, sp=sp)
        stats['eccentricity'] = eccentricity

        # diameter is the maximum eccentricity
        diameter = nx.diameter(G_strong, e=eccentricity)
        stats['diameter'] = diameter

        # radius is the minimum eccentricity
        radius = nx.radius(G_strong, e=eccentricity)
        stats['radius'] = radius

        # center is the set of nodes with eccentricity equal to radius
        center = nx.center(G_strong, e=eccentricity)
        stats['center'] = center

        # periphery is the set of nodes with eccentricity equal to the diameter
        periphery = nx.periphery(G_strong, e=eccentricity)
        stats['periphery'] = periphery

    # if True, calculate node closeness centrality
    if cc:
        # closeness centrality of a node is the reciprocal of the sum of the shortest path distances from u to all other nodes
        start_time = time.time()
        closeness_centrality = nx.closeness_centrality(G, distance='length')
        stats['closeness_centrality'] = closeness_centrality
        stats['closeness_centrality_avg'] = sum(closeness_centrality.values())/len(closeness_centrality)
        log('Calculated closeness centrality in {:,.2f} seconds'.format(time.time() - start_time))

    # if True, calculate node betweenness centrality
    if bc:
        # betweenness centrality of a node is the sum of the fraction of all-pairs shortest paths that pass through node
        start_time = time.time()
        betweenness_centrality = nx.betweenness_centrality(G, weight='length')
        stats['betweenness_centrality'] = betweenness_centrality
        stats['betweenness_centrality_avg'] = sum(betweenness_centrality.values())/len(betweenness_centrality)
        log('Calculated betweenness centrality in {:,.2f} seconds'.format(time.time() - start_time))

    log('Calculated extended stats in {:,.2f} seconds'.format(time.time()-full_start_time))
    return stats




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

    # to calculate the counts, get undirected representation of the graph. for each node, get the list of the set of unique u,v,key edges, including parallel edges but excluding self-loop parallel edges
    # (this is necessary because bi-directional self-loops will appear twice in the undirected graph as you have u,v,key0 and u,v,key1 where u==v when you convert from MultiDiGraph to MultiGraph - BUT,
    # one-way self-loops will appear only once. to get consistent accurate counts of physical streets, ignoring directionality, we need the list of the set of unique edges...).
    # then, count how many times the node appears in the u,v tuples in the list. this is the count of how many street segments emanate from this node.
    # finally create a dict of node id:count
    G_undir = G.to_undirected(reciprocal=False)
    all_edges = G_undir.edges(keys=False)
    if nodes is None:
        nodes = G_undir.nodes()

    # get all unique edges - this throws away any parallel edges (including those in self-loops)
    all_unique_edges = set(all_edges)

    # get all edges (including parallel edges) that are not self-loops
    non_self_loop_edges = [e for e in all_edges if not e[0]==e[1]]

    # get a single copy of each self-loop edge (ie, if it's bi-directional, we ignore the parallel edge going the reverse direction and keep only one copy)
    set_non_self_loop_edges = set(non_self_loop_edges)
    self_loop_edges = [e for e in all_unique_edges if e not in set_non_self_loop_edges]

    # final list contains all unique edges, including each parallel edge, unless the parallel edge is a self-loop
    # in which case it doesn't double-count the self-loop
    edges = non_self_loop_edges + self_loop_edges

    # flatten the list of (u,v) tuples
    edges_flat = list(chain.from_iterable(edges))

    # count how often each node appears in the list of flattened edge endpoints
    counts = Counter(edges_flat)
    streets_per_node = {node:counts[node] for node in nodes}
    log('Got the counts of undirected street segments incident to each node (before removing peripheral edges) in {:,.2f} seconds'.format(time.time()-start_time))
    return streets_per_node
