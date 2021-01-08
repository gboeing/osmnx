"""Calculate geometric and topological network measures."""

import logging as lg

import networkx as nx
import numpy as np
import pandas as pd

from . import distance
from . import simplification
from . import utils
from . import utils_graph


def basic_stats(G, area=None, clean_intersects=False, tolerance=15, circuity_dist="gc"):
    """
    Calculate basic descriptive geometric and topological stats for a graph.

    For an unprojected lat-lng graph, tolerance and graph units should be in
    degrees, and circuity_dist should be 'gc'. For a projected graph,
    tolerance and graph units should be in meters (or similar) and
    circuity_dist should be 'euclidean'.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    area : numeric
        the land area of this study site, in square meters. must be greater
        than 0. if None, will skip all density-based metrics.
    clean_intersects : bool
        if True, calculate consolidated intersections count (and density, if
        area is provided) via consolidate_intersections function
    tolerance : numeric
        tolerance value passed along if clean_intersects=True, see
        consolidate_intersections function documentation for details and usage
    circuity_dist : string
        'gc' or 'euclidean', how to calculate straight-line distances for
        circuity measurement; use former for lat-lng networks and latter for
        projected networks

    Returns
    -------
    stats : dict
        network measures containing the following elements (some keys may not
        be present, based on the arguments passed into the function):

          - n = number of nodes in the graph
          - m = number of edges in the graph
          - k_avg = average node degree of the graph
          - intersection_count = number of intersections in graph, that is,
                nodes with >1 physical street connected to them
          - streets_per_node_avg = how many physical streets (edges in the
                undirected representation of the graph) connect to each node
                (ie, intersection or dead-end) on average (mean)
          - streets_per_node_counts = dict with keys of number of physical
                streets connecting to a node, and values of number of nodes
                with this count
          - streets_per_node_proportion = dict, same as previous, but as a
                proportion of the total, rather than counts
          - edge_length_total = sum of all edge lengths in graph, in meters
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
                as its endpoints (ie, the edge links nodes u and v, and u==v)
          - clean_intersection_count = number of intersections in street
                network, merging complex ones into single points
          - clean_intersection_density_km = clean_intersection_count divided
                by area in square kilometers
    """
    sq_m_in_sq_km = 1_000_000  # there are 1 million sq meters in 1 sq km
    Gu = None
    node_ids = set(G.nodes)

    # calculate number of nodes, n, and number of edges, m, in the graph
    n = len(G)
    m = len(G.edges)

    # calculate the average degree of the graph
    k_avg = 2 * m / n

    # get number of physical streets (undirected edges) connected to each node
    spn = nx.get_node_attributes(G, "street_count")
    if set(spn) != node_ids:
        utils.log("Graph nodes changed since `street_count`s were calculated", level=lg.WARN)

    # count number of intersections in graph, as nodes with >1 street
    intersect_count = sum(count > 1 and node in node_ids for node, count in spn.items())

    # calculate streets-per-node average: the average number of streets
    # (unidirected edges) incident to each node
    spna = sum(spn.values()) / n

    # calculate streets-per-node counts
    # create a dict where key = number of streets (unidirected edges) incident
    # to each node, and value = how many nodes are of this number in the graph
    spnc = {num: list(spn.values()).count(num) for num in range(max(spn.values()) + 1)}

    # calculate streets-per-node proportion
    # degree proportions: dict where key = each degree and value = what
    # proportion of nodes are of this degree in the graph
    spnp = {num: count / n for num, count in spnc.items()}

    # calculate the total and average edge lengths
    edge_length_total = sum(d["length"] for u, v, d in G.edges(data=True))
    edge_length_avg = edge_length_total / m

    # calculate total and average street segment lengths (so, edges without
    # double-counting two-way streets)
    if Gu is None:
        Gu = utils_graph.get_undirected(G)
    street_length_total = sum(d["length"] for u, v, d in Gu.edges(data=True))
    street_segments_count = len(Gu.edges(keys=True))
    street_length_avg = street_length_total / street_segments_count

    # calculate clean intersection counts
    if clean_intersects:
        points = simplification.consolidate_intersections(G, tolerance, False, False)
        clean_intersect_count = len(points)
    else:
        clean_intersect_count = None

    # we can calculate density metrics only if area is not null
    if area is not None:
        area_km = area / sq_m_in_sq_km

        # calculate node density as nodes per sq km
        node_density_km = n / area_km

        # calculate intersection density as nodes with >1 physical street
        # connected to them, per sq km
        intersect_density_km = intersect_count / area_km

        # calculate edge density as linear meters per sq km
        edge_density_km = edge_length_total / area_km

        # calculate street density as linear meters per sq km
        street_density_km = street_length_total / area_km

        if clean_intersects:
            clean_intersect_density_km = clean_intersect_count / area_km
        else:
            clean_intersect_density_km = None
    else:
        # if area is None, then we cannot calculate density
        node_density_km = None
        intersect_density_km = None
        edge_density_km = None
        street_density_km = None
        clean_intersect_density_km = None

    # average circuity: sum of edge lengths divided by sum of straight-line
    # distance between edge endpoints. first load all the edges origin and
    # destination coordinates as a dataframe, then calculate the straight-line
    # distance
    coords = (
        (G.nodes[u]["y"], G.nodes[u]["x"], G.nodes[v]["y"], G.nodes[v]["x"]) for u, v, k in G.edges
    )
    df_coords = pd.DataFrame(coords, columns=["u_y", "u_x", "v_y", "v_x"])
    if circuity_dist == "gc":
        gc_distances = distance.great_circle_vec(
            lat1=df_coords["u_y"],
            lng1=df_coords["u_x"],
            lat2=df_coords["v_y"],
            lng2=df_coords["v_x"],
        )
    elif circuity_dist == "euclidean":
        gc_distances = distance.euclidean_dist_vec(
            y1=df_coords["u_y"], x1=df_coords["u_x"], y2=df_coords["v_y"], x2=df_coords["v_x"]
        )
    else:
        raise ValueError('circuity_dist argument must be "gc" or "euclidean"')

    gc_distances = gc_distances.fillna(value=0)
    try:
        circuity_avg = edge_length_total / gc_distances.sum()
    except ZeroDivisionError:
        circuity_avg = np.nan

    # percent of edges that are self-loops, ie both endpoints are same node
    self_loop_proportion = sum(u == v for u, v, k in G.edges) / m

    # assemble the results
    stats = {
        "n": n,
        "m": m,
        "k_avg": k_avg,
        "intersection_count": intersect_count,
        "streets_per_node_avg": spna,
        "streets_per_node_counts": spnc,
        "streets_per_node_proportion": spnp,
        "edge_length_total": edge_length_total,
        "edge_length_avg": edge_length_avg,
        "street_length_total": street_length_total,
        "street_length_avg": street_length_avg,
        "street_segments_count": street_segments_count,
        "node_density_km": node_density_km,
        "intersection_density_km": intersect_density_km,
        "edge_density_km": edge_density_km,
        "street_density_km": street_density_km,
        "circuity_avg": circuity_avg,
        "self_loop_proportion": self_loop_proportion,
        "clean_intersection_count": clean_intersect_count,
        "clean_intersection_density_km": clean_intersect_density_km,
    }

    # return the results
    return stats


def extended_stats(G, connectivity=False, anc=False, ecc=False, bc=False, cc=False):
    """
    Calculate extended topological measures for a graph.

    Many of these algorithms have an inherently high time complexity. Global
    topological analysis of large complex networks is extremely time consuming
    and may exhaust computer memory. Consider using function arguments to not
    run metrics that require computation of a full matrix of paths if they
    will not be needed.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    connectivity : bool
        if True, calculate node and edge connectivity
    anc : bool
        if True, calculate average node connectivity
    ecc : bool
        if True, calculate shortest paths, eccentricity, and topological
        metrics that use eccentricity
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
    stats = dict()

    # create DiGraph from the MultiDiGraph, for those metrics that need it
    D = utils_graph.get_digraph(G, weight="length")

    # create undirected Graph from the DiGraph, for those metrics that need it
    Gu = nx.Graph(D)

    # get largest strongly connected component, for those metrics that require
    # strongly connected graphs
    Gs = utils_graph.get_largest_component(G, strongly=True)

    # average degree of the neighborhood of each node, and average for graph
    avg_neighbor_degree = nx.average_neighbor_degree(G)
    stats["avg_neighbor_degree"] = avg_neighbor_degree
    stats["avg_neighbor_degree_avg"] = sum(avg_neighbor_degree.values()) / len(avg_neighbor_degree)

    # avg weighted degree of neighborhood of each node, and average for graph
    avg_wtd_nbr_deg = nx.average_neighbor_degree(G, weight="length")
    stats["avg_weighted_neighbor_degree"] = avg_wtd_nbr_deg
    stats["avg_weighted_neighbor_degree_avg"] = sum(avg_wtd_nbr_deg.values()) / len(avg_wtd_nbr_deg)

    # degree centrality for a node is the fraction of nodes it is connected to
    degree_centrality = nx.degree_centrality(G)
    stats["degree_centrality"] = degree_centrality
    stats["degree_centrality_avg"] = sum(degree_centrality.values()) / len(degree_centrality)

    # calculate clustering coefficient for the nodes
    stats["clustering_coefficient"] = nx.clustering(Gu)

    # average clustering coefficient for the graph
    stats["clustering_coefficient_avg"] = nx.average_clustering(Gu)

    # calculate weighted clustering coefficient for the nodes
    stats["clustering_coefficient_weighted"] = nx.clustering(Gu, weight="length")

    # average clustering coefficient (weighted) for the graph
    stats["clustering_coefficient_weighted_avg"] = nx.average_clustering(Gu, weight="length")

    # pagerank: a ranking of the nodes in the graph based on the structure of
    # the incoming links
    pagerank = nx.pagerank(D, weight="length")
    stats["pagerank"] = pagerank

    # node with the highest page rank, and its value
    pagerank_max_node = max(pagerank, key=lambda x: pagerank[x])
    stats["pagerank_max_node"] = pagerank_max_node
    stats["pagerank_max"] = pagerank[pagerank_max_node]

    # node with the lowest page rank, and its value
    pagerank_min_node = min(pagerank, key=lambda x: pagerank[x])
    stats["pagerank_min_node"] = pagerank_min_node
    stats["pagerank_min"] = pagerank[pagerank_min_node]

    # if True, calculate node and edge connectivity
    if connectivity:

        # node connectivity is minimum number of nodes that must be removed
        # to disconnect G or render it trivial
        stats["node_connectivity"] = nx.node_connectivity(Gs)

        # edge connectivity is equal to minimum number of edges that must be
        # removed to disconnect G or render it trivial
        stats["edge_connectivity"] = nx.edge_connectivity(Gs)
        utils.log("Calculated node and edge connectivity")

    # if True, calculate average node connectivity
    if anc:
        # mean number of internally node-disjoint paths between each pair of
        # nodes in G, i.e., expected number of nodes that must be removed to
        # disconnect a randomly selected pair of non-adjacent nodes
        stats["node_connectivity_avg"] = nx.average_node_connectivity(G)
        utils.log("Calculated average node connectivity")

    # if True, calculate shortest paths, eccentricity, and topological metrics
    # that use eccentricity
    if ecc:
        # precompute shortest paths between all nodes for eccentricity-based
        # stats
        length_func = nx.single_source_dijkstra_path_length
        sp = {source: dict(length_func(Gs, source, weight="length")) for source in Gs.nodes}

        utils.log("Calculated shortest path lengths")

        # eccentricity of a node v is the maximum distance from v to all other
        # nodes in G
        eccentricity = nx.eccentricity(Gs, sp=sp)
        stats["eccentricity"] = eccentricity

        # diameter is the maximum eccentricity
        diameter = nx.diameter(Gs, e=eccentricity)
        stats["diameter"] = diameter

        # radius is the minimum eccentricity
        radius = nx.radius(Gs, e=eccentricity)
        stats["radius"] = radius

        # center is the set of nodes with eccentricity equal to radius
        center = nx.center(Gs, e=eccentricity)
        stats["center"] = center

        # periphery is the set of nodes with eccentricity equal to diameter
        periphery = nx.periphery(Gs, e=eccentricity)
        stats["periphery"] = periphery

    # if True, calculate node closeness centrality
    if cc:
        # closeness centrality of a node is the reciprocal of the sum of the
        # shortest path distances from u to all other nodes
        close_cent = nx.closeness_centrality(G, distance="length")
        stats["closeness_centrality"] = close_cent
        stats["closeness_centrality_avg"] = sum(close_cent.values()) / len(close_cent)
        utils.log("Calculated closeness centrality")

    # if True, calculate node betweenness centrality
    if bc:
        # betweenness centrality of a node is the sum of the fraction of
        # all-pairs shortest paths that pass through node. nx2.4+
        # implementation cannot run on Multi(Di)Graphs, so use DiGraph
        btwn_cent = nx.betweenness_centrality(D, weight="length")
        stats["betweenness_centrality"] = btwn_cent
        stats["betweenness_centrality_avg"] = sum(btwn_cent.values()) / len(btwn_cent)
        utils.log("Calculated betweenness centrality")

    utils.log("Calculated extended stats")
    return stats
