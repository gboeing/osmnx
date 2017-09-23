################################################################################
# Module: stats.py
# Description: Calculate graph-theoretic topological and metric measures for a
#              network
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

from __future__ import division
import time
import networkx as nx
import numpy as np
import pandas as pd

from .simplify import clean_intersections
from .utils import log
from .utils import get_largest_component
from .utils import great_circle_vec
from .utils import count_streets_per_node
from .utils import euclidean_dist_vec


def basic_stats(G, area=None, clean_intersects=False, tolerance=15,
                circuity_dist='gc'):
    """
    Calculate basic descriptive metric and topological stats for a graph.

    For an unprojected lat-lng graph, tolerance and graph units should be in
    degrees, and circuity_dist should be 'gc'. For a projected graph, tolerance
    and graph units should be in meters (or similar) and circuity_dist should be
    'euclidean'.

    Parameters
    ----------
    G : networkx multidigraph
    area : numeric
        the area covered by the street network, in square meters (typically land
        area); if none, will skip all density-based metrics
    clean_intersects : bool
        if True, calculate clean intersections count (and density, if area is
        provided)
    tolerance : numeric
        tolerance value passed along if clean_intersects=True, see
        clean_intersections() function documentation for details and usage
    circuity_dist : str
        'gc' or 'euclidean', how to calculate straight-line distances for
        circuity measurement; use former for lat-lng networks and latter for
        projected networks

    Returns
    -------
    stats : dict
        dictionary of network measures containing the following elements (some
        keys may not be present, based on the arguments passed into the function):

          - n = number of nodes in the graph
          - m = number of edges in the graph
          - k_avg = average node degree of the graph
          - intersection_count = number of intersections in graph, that is,
                nodes with >1 street emanating from them
          - streets_per_node_avg = how many streets (edges in the undirected
                representation of the graph) emanate from each node (ie,
                intersection or dead-end) on average (mean)
          - streets_per_node_counts = dict, with keys of number of streets
                emanating from the node, and values of number of nodes with this
                count
          - streets_per_node_proportion = dict, same as previous, but as a
                proportion of the total, rather than counts
          - edge_length_total = sum of all edge lengths in the graph, in meters
          - edge_length_avg = mean edge length in the graph, in meters
          - street_length_total = sum of all edges in the undirected
                representation of the graph
          - street_length_avg = mean edge length in the undirected
                representation of the graph, in meters
          - street_segments_count = number of edges in the undirected
                representation of the graph
          - node_density_km = n divided by area in square kilometers
          - intersection_density_km = intersection_count divided by area in
                square kilometers
          - edge_density_km = edge_length_total divided by area in square
                kilometers
          - street_density_km = street_length_total divided by area in square
                kilometers
          - circuity_avg = edge_length_total divided by the sum of the great
                circle distances between the nodes of each edge
          - self_loop_proportion = proportion of edges that have a single node
                as its two endpoints (ie, the edge links nodes u and v, and u==v)
          - clean_intersection_count = number of intersections in street
                network, merging complex ones into single points
          - clean_intersection_density_km = clean_intersection_count divided by
                area in square kilometers
    """

    sq_m_in_sq_km = 1e6 #there are 1 million sq meters in 1 sq km
    G_undirected = None

    # calculate the number of nodes, n, and the number of edges, m, in the graph
    n = len(list(G.nodes()))
    m = len(list(G.edges()))

    # calculate the average degree of the graph
    k_avg = 2 * m / n

    if 'streets_per_node' in G.graph:
        # get the degrees saved as a graph attribute (from an undirected
        # representation of the graph). this is not the degree of the nodes in
        # the directed graph, but rather represents the number of streets
        # (unidirected edges) emanating from each node. see
        # count_streets_per_node function.
        streets_per_node = G.graph['streets_per_node']
    else:
        # count how many street segments emanate from each node in this graph
        streets_per_node = count_streets_per_node(G)

    # count number of intersections in graph, as nodes with >1 street emanating
    # from them
    node_ids = set(G.nodes())
    intersection_count = len([True for node, count in streets_per_node.items() if (count > 1) and (node in node_ids)])

    # calculate the average number of streets (unidirected edges) incident to
    # each node
    streets_per_node_avg = sum(streets_per_node.values()) / n

    # create a dict where key = number of streets (unidirected edges) incident
    # to each node, and value = how many nodes are of this number in the graph
    streets_per_node_counts = {num:list(streets_per_node.values()).count(num) for num in range(max(streets_per_node.values()) + 1)}

    # degree proportions: dict where key = each degree and value = what
    # proportion of nodes are of this degree in the graph
    streets_per_node_proportion = {num:count/n for num, count in streets_per_node_counts.items()}

    # calculate the total and average edge lengths
    edge_length_total = sum([d['length'] for u, v, d in G.edges(data=True)])
    edge_length_avg = edge_length_total / m

    # calculate the total and average street segment lengths (so, edges without
    # double-counting two-way streets)
    if G_undirected is None:
        G_undirected = G.to_undirected(reciprocal=False)
    street_length_total = sum([d['length'] for u, v, d in G_undirected.edges(data=True)])
    street_segments_count = len(list(G_undirected.edges(keys=True)))
    street_length_avg = street_length_total / street_segments_count

    # calculate clean intersection counts
    if clean_intersects:
        clean_intersection_points = clean_intersections(G, tolerance=tolerance, dead_ends=False )
        clean_intersection_count = len(clean_intersection_points)
    else:
        clean_intersection_count = None

    # we can calculate density metrics only if area is not null
    if area is not None:
        area_km = area / sq_m_in_sq_km

        # calculate node density as nodes per sq km
        node_density_km = n / area_km

        # calculate intersection density as nodes with >1 street emanating from
        # them, per sq km
        intersection_density_km = intersection_count / area_km

        # calculate edge density as linear meters per sq km
        edge_density_km = edge_length_total / area_km

        # calculate street density as linear meters per sq km
        street_density_km = street_length_total / area_km

        if clean_intersects:
            clean_intersection_density_km = clean_intersection_count / area_km
        else:
            clean_intersection_density_km = None
    else:
        # if area is None, then we cannot calculate density
        node_density_km = None
        intersection_density_km = None
        edge_density_km = None
        street_density_km = None
        clean_intersection_density_km = None

    # average circuity: sum of edge lengths divided by sum of straight-line
    # distance between edge endpoints. first load all the edges origin and
    # destination coordinates as a dataframe, then calculate the straight-line
    # distance
    coords = np.array([[G.nodes[u]['y'], G.nodes[u]['x'], G.nodes[v]['y'], G.nodes[v]['x']] for u, v, k in G.edges(keys=True)])
    df_coords = pd.DataFrame(coords, columns=['u_y', 'u_x', 'v_y', 'v_x'])
    if circuity_dist == 'gc':
        gc_distances = great_circle_vec(lat1=df_coords['u_y'],
                                        lng1=df_coords['u_x'],
                                        lat2=df_coords['v_y'],
                                        lng2=df_coords['v_x'])
    elif circuity_dist == 'euclidean':
        gc_distances = euclidean_dist_vec(y1=df_coords['u_y'],
                                          x1=df_coords['u_x'],
                                          y2=df_coords['v_y'],
                                          x2=df_coords['v_x'])
    else:
        raise ValueError('circuity_dist must be "gc" or "euclidean"')

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
             'intersection_count':intersection_count,
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
             'self_loop_proportion':self_loop_proportion,
             'clean_intersection_count':clean_intersection_count,
             'clean_intersection_density_km':clean_intersection_density_km}

    # return the results
    return stats


def extended_stats(G, connectivity=False, anc=False, ecc=False, bc=False, cc=False):
    """
    Calculate extended topological stats and metrics for a graph.

    Many of these algorithms have an inherently high time complexity. Global
    topological analysis of large complex networks is extremely time consuming
    and may exhaust computer memory. Consider using function arguments to not
    run metrics that require computation of a full matrix of paths if they
    will not be needed.

    Parameters
    ----------
    G : networkx multidigraph
    connectivity : bool
        if True, calculate node and edge connectivity
    anc : bool
        if True, calculate average node connectivity
    ecc : bool
        if True, calculate shortest paths, eccentricity, and topological metrics
        that use eccentricity
    bc : bool
        if True, calculate node betweenness centrality
    cc : bool
        if True, calculate node closeness centrality

    Returns
    -------
    stats : dict
        dictionary of network measures containing the following elements (some
        only calculated/returned optionally, based on passed parameters):

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

    # create an undirected Graph from the MultiDiGraph, for those metrics that
    # require it
    G_undir = nx.Graph(G)

    # get the largest strongly connected component, for those metrics that
    # require strongly connected graphs
    G_strong = get_largest_component(G, strongly=True)

    # average degree of the neighborhood of each node, and average for the graph
    avg_neighbor_degree = nx.average_neighbor_degree(G)
    stats['avg_neighbor_degree'] = avg_neighbor_degree
    stats['avg_neighbor_degree_avg'] = sum(avg_neighbor_degree.values())/len(avg_neighbor_degree)

    # average weighted degree of the neighborhood of each node, and average for
    # the graph
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

    # pagerank: a ranking of the nodes in the graph based on the structure of
    # the incoming links
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

        # node connectivity is the minimum number of nodes that must be removed
        # to disconnect G or render it trivial
        stats['node_connectivity'] = nx.node_connectivity(G_strong)

        # edge connectivity is equal to the minimum number of edges that must be
        # removed to disconnect G or render it trivial
        stats['edge_connectivity'] = nx.edge_connectivity(G_strong)
        log('Calculated node and edge connectivity in {:,.2f} seconds'.format(time.time() - start_time))

    # if True, calculate average node connectivity
    if anc:
        # mean number of internally node-disjoint paths between each pair of
        # nodes in G, i.e., the expected number of nodes that must be removed to
        # disconnect a randomly selected pair of non-adjacent nodes
        start_time = time.time()
        stats['node_connectivity_avg'] = nx.average_node_connectivity(G)
        log('Calculated average node connectivity in {:,.2f} seconds'.format(time.time() - start_time))

    # if True, calculate shortest paths, eccentricity, and topological metrics
    # that use eccentricity
    if ecc:
        # precompute shortest paths between all nodes for eccentricity-based
        # stats
        start_time = time.time()
        sp = {source:dict(nx.single_source_dijkstra_path_length(G_strong, source, weight='length')) for source in G_strong.nodes()}

        log('Calculated shortest path lengths in {:,.2f} seconds'.format(time.time() - start_time))

        # eccentricity of a node v is the maximum distance from v to all other
        # nodes in G
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
        # closeness centrality of a node is the reciprocal of the sum of the
        # shortest path distances from u to all other nodes
        start_time = time.time()
        closeness_centrality = nx.closeness_centrality(G, distance='length')
        stats['closeness_centrality'] = closeness_centrality
        stats['closeness_centrality_avg'] = sum(closeness_centrality.values())/len(closeness_centrality)
        log('Calculated closeness centrality in {:,.2f} seconds'.format(time.time() - start_time))

    # if True, calculate node betweenness centrality
    if bc:
        # betweenness centrality of a node is the sum of the fraction of
        # all-pairs shortest paths that pass through node
        start_time = time.time()
        betweenness_centrality = nx.betweenness_centrality(G, weight='length')
        stats['betweenness_centrality'] = betweenness_centrality
        stats['betweenness_centrality_avg'] = sum(betweenness_centrality.values())/len(betweenness_centrality)
        log('Calculated betweenness centrality in {:,.2f} seconds'.format(time.time() - start_time))

    log('Calculated extended stats in {:,.2f} seconds'.format(time.time()-full_start_time))
    return stats
