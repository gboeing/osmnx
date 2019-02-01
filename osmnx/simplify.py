################################################################################
# Module: simplify.py
# Description: Simplify and correct network topology
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

import time
import collections
import logging as lg
import numpy   as np
import pandas  as pd
import networkx as nx
import geopandas as gpd
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import LineString

from .save_load import graph_to_gdfs
from .utils import log
from .utils import count_streets_per_node

# scipy and sklearn are optional dependencies for faster nearest node search
try:
    from scipy.spatial import cKDTree
except ImportError as e:
    cKDTree = None

def is_endpoint(G, node, strict=True):
    """
    Return True if the node is a "real" endpoint of an edge in the network, \
    otherwise False. OSM data includes lots of nodes that exist only as points \
    to help streets bend around curves. An end point is a node that either: \
    1) is its own neighbor, ie, it self-loops. \
    2) or, has no incoming edges or no outgoing edges, ie, all its incident \
        edges point inward or all its incident edges point outward. \
    3) or, it does not have exactly two neighbors and degree of 2 or 4. \
    4) or, if strict mode is false, if its edges have different OSM IDs. \

    Parameters
    ----------
    G : networkx multidigraph

    node : int
        the node to examine
    strict : bool
        if False, allow nodes to be end points even if they fail all other rules \
        but have edges with different OSM IDs

    Returns
    -------
    bool

    """
    neighbors = set(list(G.predecessors(node)) + list(G.successors(node)))
    n = len(neighbors)
    d = G.degree(node)

    if node in neighbors:
        # if the node appears in its list of neighbors, it self-loops. this is
        # always an endpoint.
        return True

    # if node has no incoming edges or no outgoing edges, it must be an endpoint
    elif G.out_degree(node)==0 or G.in_degree(node)==0:
        return True

    elif not (n==2 and (d==2 or d==4)):
        # else, if it does NOT have 2 neighbors AND either 2 or 4 directed
        # edges, it is an endpoint. either it has 1 or 3+ neighbors, in which
        # case it is a dead-end or an intersection of multiple streets or it has
        # 2 neighbors but 3 degree (indicating a change from oneway to twoway)
        # or more than 4 degree (indicating a parallel edge) and thus is an
        # endpoint
        return True

    elif not strict:
        # non-strict mode
        osmids = []

        # add all the edge OSM IDs for incoming edges
        for u in G.predecessors(node):
            for key in G[u][node]:
                osmids.append(G.edges[u, node, key]['osmid'])

        # add all the edge OSM IDs for outgoing edges
        for v in G.successors(node):
            for key in G[node][v]:
                osmids.append(G.edges[node, v, key]['osmid'])

        # if there is more than 1 OSM ID in the list of edge OSM IDs then it is
        # an endpoint, if not, it isn't
        return len(set(osmids)) > 1

    else:
        # if none of the preceding rules returned true, then it is not an endpoint
        return False


def build_path(G, node, endpoints, path):
    """
    Recursively build a path of nodes until you hit an endpoint node.

    Parameters
    ----------
    G : networkx multidigraph
    node : int
        the current node to start from
    endpoints : set
        the set of all nodes in the graph that are endpoints
    path : list
        the list of nodes in order in the path so far

    Returns
    -------
    paths_to_simplify : list
    """
    # for each successor in the passed-in node
    for successor in G.successors(node):
        if successor not in path:
            # if this successor is already in the path, ignore it, otherwise add
            # it to the path
            path.append(successor)
            if successor not in endpoints:
                # if this successor is not an endpoint, recursively call
                # build_path until you find an endpoint
                path = build_path(G, successor, endpoints, path)
            else:
                # if this successor is an endpoint, we've completed the path,
                # so return it
                return path

    if (path[-1] not in endpoints) and (path[0] in G.successors(path[-1])):
        # if the end of the path is not actually an endpoint and the path's
        # first node is a successor of the path's final node, then this is
        # actually a self loop, so add path's first node to end of path to
        # close it
        path.append(path[0])

    return path


def get_paths_to_simplify(G, strict=True):
    """
    Create a list of all the paths to be simplified between endpoint nodes.

    The path is ordered from the first endpoint, through the interstitial nodes,
    to the second endpoint. If your street network is in a rural area with many
    interstitial nodes between true edge endpoints, you may want to increase
    your system's recursion limit to avoid recursion errors.

    Parameters
    ----------
    G : networkx multidigraph
    strict : bool
        if False, allow nodes to be end points even if they fail all other rules
        but have edges with different OSM IDs

    Returns
    -------
    paths_to_simplify : list
    """

    # first identify all the nodes that are endpoints
    start_time = time.time()
    endpoints = set([node for node in G.nodes() if is_endpoint(G, node, strict=strict)])
    log('Identified {:,} edge endpoints in {:,.2f} seconds'.format(len(endpoints), time.time()-start_time))

    start_time = time.time()
    paths_to_simplify = []

    # for each endpoint node, look at each of its successor nodes
    for node in endpoints:
        for successor in G.successors(node):
            if successor not in endpoints:
                # if the successor is not an endpoint, build a path from the
                # endpoint node to the next endpoint node
                try:
                    path = build_path(G, successor, endpoints, path=[node, successor])
                    paths_to_simplify.append(path)
                except RuntimeError:
                    log('Recursion error: exceeded max depth, moving on to next endpoint successor', level=lg.WARNING)
                    # recursion errors occur if some connected component is a
                    # self-contained ring in which all nodes are not end points.
                    # could also occur in extremely long street segments (eg, in
                    # rural areas) with too many nodes between true endpoints.
                    # handle it by just ignoring that component and letting its
                    # topology remain intact (this should be a rare occurrence)
                    # RuntimeError is what Python <3.5 will throw, Py3.5+ throws
                    # RecursionError but it is a subtype of RuntimeError so it
                    # still gets handled

    log('Constructed all paths to simplify in {:,.2f} seconds'.format(time.time()-start_time))
    return paths_to_simplify


def is_simplified(G):
    """
    Determine if a graph has already had its topology simplified.

    If any of its edges have a geometry attribute, we know that it has
    previously been simplified.

    Parameters
    ----------
    G : networkx multidigraph

    Returns
    -------
    bool
    """

    return 'simplified' in G.graph and G.graph['simplified']


def simplify_graph(G, strict=True):
    """
    Simplify a graph's topology by removing all nodes that are not intersections
    or dead-ends.

    Create an edge directly between the end points that encapsulate them,
    but retain the geometry of the original edges, saved as attribute in new
    edge.

    Parameters
    ----------
    G : networkx multidigraph
    strict : bool
        if False, allow nodes to be end points even if they fail all other rules
        but have edges with different OSM IDs

    Returns
    -------
    networkx multidigraph
    """

    if is_simplified(G):
        raise Exception('This graph has already been simplified, cannot simplify it again.')

    log('Begin topologically simplifying the graph...')
    G = G.copy()
    initial_node_count = len(list(G.nodes()))
    initial_edge_count = len(list(G.edges()))
    all_nodes_to_remove = []
    all_edges_to_add = []

    # construct a list of all the paths that need to be simplified
    paths = get_paths_to_simplify(G, strict=strict)

    start_time = time.time()
    for path in paths:

        # add the interstitial edges we're removing to a list so we can retain
        # their spatial geometry
        edge_attributes = {}
        for u, v in zip(path[:-1], path[1:]):

            # there shouldn't be multiple edges between interstitial nodes
            if not G.number_of_edges(u, v) == 1:
                log('Multiple edges between "{}" and "{}" found when simplifying'.format(u, v), level=lg.WARNING)

            # the only element in this list as long as above check is True
            # (MultiGraphs use keys (the 0 here), indexed with ints from 0 and
            # up)
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
            if len(set(edge_attributes[key])) == 1 and not key == 'length':
                # if there's only 1 unique value in this attribute list,
                # consolidate it to the single value (the zero-th)
                edge_attributes[key] = edge_attributes[key][0]
            elif not key == 'length':
                # otherwise, if there are multiple values, keep one of each value
                edge_attributes[key] = list(set(edge_attributes[key]))

        # construct the geometry and sum the lengths of the segments
        edge_attributes['geometry'] = LineString([Point((G.nodes[node]['x'], G.nodes[node]['y'])) for node in path])
        edge_attributes['length'] = sum(edge_attributes['length'])

        # add the nodes and edges to their lists for processing at the end
        all_nodes_to_remove.extend(path[1:-1])
        all_edges_to_add.append({'origin':path[0],
                                 'destination':path[-1],
                                 'attr_dict':edge_attributes})

    # for each edge to add in the list we assembled, create a new edge between
    # the origin and destination
    for edge in all_edges_to_add:
        G.add_edge(edge['origin'], edge['destination'], **edge['attr_dict'])

    # finally remove all the interstitial nodes between the new edges
    G.remove_nodes_from(set(all_nodes_to_remove))

    G.graph['simplified'] = True

    msg = 'Simplified graph (from {:,} to {:,} nodes and from {:,} to {:,} edges) in {:,.2f} seconds'
    log(msg.format(initial_node_count, len(list(G.nodes())), initial_edge_count, len(list(G.edges())), time.time()-start_time))
    return G


def clean_intersections(G, tolerance=15, dead_ends=False):
    """
    Clean-up intersections comprising clusters of nodes by merging them and
    returning their centroids.

    Divided roads are represented by separate centerline edges. The intersection
    of two divided roads thus creates 4 nodes, representing where each edge
    intersects a perpendicular edge. These 4 nodes represent a single
    intersection in the real world. This function cleans them up by buffering
    their points to an arbitrary distance, merging overlapping buffers, and
    taking their centroid. For best results, the tolerance argument should be
    adjusted to approximately match street design standards in the specific
    street network.

    Parameters
    ----------
    G : networkx multidigraph
    tolerance : float
        nodes within this distance (in graph's geometry's units) will be
        dissolved into a single intersection
    dead_ends : bool
        if False, discard dead-end nodes to return only street-intersection
        points

    Returns
    ----------
    G__ : networkx multidigraph
        a new spatial graph, topologically equal to G (approximately), but with
        clusters of nodes merged into a single centroid node.
    """

    start_time = time.time()

    if 'proj' not in G.graph['crs'] and G.graph['proj'] != "utm":
        log("The provided graph might not have been projected to utm.",
            level = lg.WARNING)

    # Let's make a copy of G
    G = G.copy()

    # if dead_ends is False, discard dead-end nodes to only work with edge
    # intersections
    if not dead_ends:
        if 'streets_per_node' in G.graph:
            streets_per_node = G.graph['streets_per_node']
        else:
            streets_per_node = count_streets_per_node(G)

        dead_end_nodes = [node for node, count in streets_per_node.items() if count <= 1]
        G.remove_nodes_from(dead_end_nodes)

    # create a GeoDataFrame of nodes, buffer to passed-in distance, merge
    # overlaps
    gdf_nodes = graph_to_gdfs(G, edges=False)
    buffered_nodes = gdf_nodes.buffer(tolerance).unary_union
    if isinstance(buffered_nodes, Polygon):
        # if only a single node results, make it iterable so we can turn it into
        # a GeoSeries
        buffered_nodes = [buffered_nodes]

    # get the centroids of the merged intersection polygons
    unified_intersections = gpd.GeoSeries(list(buffered_nodes))
    intersection_centroids = unified_intersections.centroid

    # To make things simpler, every edge should have a geometry
    # (this avoids KeyError later when choosing between several
    #  geometries or when the centroid inherits a geometry directly)
    for u, v, data in G.edges(keys=False, data=True):
        if 'geometry' not in data:
            # if it doesn't have a geometry attribute, the edge is a straight
            # line from node to node
            x1 = G.nodes[u]['x']
            y1 = G.nodes[u]['y']
            x2 = G.nodes[v]['x']
            y2 = G.nodes[v]['y']
            data['geometry'] = LineString([(x1, y1), (x2, y2)])

    # First, we need to find the nearest centroid to each node
    if cKDTree:
        # the search is faster if we use a cKDTree
        # build a k-d tree for euclidean nearest node search
        centroids = pd.DataFrame({'x' : intersection_centroids.apply(lambda pt: pt.x),
                                  'y' : intersection_centroids.apply(lambda pt: pt.y)})

        tree = cKDTree(data=centroids, compact_nodes=True, balanced_tree=True)

        X = [ x for node, x in G.nodes(data = 'x') ]
        Y = [ y for node, y in G.nodes(data = 'y') ]

        # query the tree for nearest centroid to each point
        points = np.array([X, Y]).T
        dist, idx = tree.query(points, k=1)
        nearest_centroid = centroids.iloc[idx].index

        # We build the mapping: osmid -> centroid
        node_centroids = { osmid : nearest_centroid[i]
                           for i, osmid in enumerate(G.nodes(data = False))}
    else:
        # else a lot slower if scipy is not available
        nodes_centroids = {
            osmid : np.argmin([node['geometry'].distance(centroid) \
                               for centroid in intersection_centroids]) \
            for osmid, node in gdf_nodes.iterrows()
        }

    # We also build the reverse mapping: centroid -> osmids
    centroid_elements = {i : [] for i in range(len(intersection_centroids))}
    for i, osmid in enumerate(G.nodes(data = False)):
        centroid = nearest_centroid[i]
        centroid_elements[centroid].append(osmid)

    # We're ready to start building the new graph
    G__ = nx.DiGraph(name = G.graph['name'],
                     crs = G.graph['crs'])

    # Add the nodes first
    # The centroids are given new ids = 0 .. total_centroids-1
    for centroid, elements in centroid_elements.items():
        x = intersection_centroids[centroid].x
        y = intersection_centroids[centroid].y

        osmid = elements[0] if len(elements) == 1 else elements
        G__.add_node(centroid, osmid = osmid, x = x, y = y)

    for centroid, elements in centroid_elements.items():
            # get the outgoing edges for all elements in the centroid
            out_edges = []
            for node in elements:
                out_edges += list(G.out_edges(node))

            # select the neighbor nodes of each element in the centroid
            out_nodes = map(lambda edge: edge[1], out_edges)

            # and get the corresponding centroid of each element
            neighbors = list(map(lambda osmid: node_centroids[osmid],
                                 out_nodes))

            # remove duplicate neighbor centroids
            neighbors = set(neighbors)
            # and self-loops
            if centroid in neighbors:
                neighbors.remove(centroid)

            # helper for fixing geometries
            centroid_point = Point(
                intersection_centroids[centroid].x,
                intersection_centroids[centroid].y)

            # add an outgoing edge to each resulting neighbor
            for neighbor in neighbors:
                # retrieve the elements of this neighbor
                neighbor_elements = centroid_elements[neighbor]
                # select the edges (c1,c2) where c1 is an element of this centroid,
                # and c2 is an element of this neighbor's centroid
                out_edges_neighbor = set(filter(lambda x: x[1] in neighbor_elements, out_edges))

                if len(out_edges_neighbor) > 1:
                    # If there are multiple edges that fit this criteria, then we need
                    # to combine their attributes, except for geometry and length
                    # for which we need to pick a single value

                    edges_attrs = [ G.edges[edge[0], edge[1], 0] for edge in out_edges_neighbor ]

                    # Not every edge has the same set of keys, so we compute the union of the individual
                    # sets, and also count their occurrence for initialisation purposes
                    all_keys = []
                    for edge_attr in edges_attrs:
                        all_keys += edge_attr.keys()

                    key_counter = collections.Counter(all_keys)

                    # If if the key counter is greater than 1, we initialise the new
                    # dictionary with empty lists, so that we can later call append()
                    attr = {key : [] if count > 1 else "" for key, count in key_counter.items()}

                    for edge_attr in edges_attrs:
                        for key in edge_attr.keys():
                            val = edge_attr[key]

                            # This key only occurs once, so we can set its value
                            if key_counter[key] == 1:
                                attr[key] = val
                            # Otherwise, append
                            else:
                                attr[key].append(val)

                    # We pick the geometry with largest length
                    # but, there may be a better way to make this decision
                    whichmax = np.argmax(attr['length'])
                    attr['length'] = attr['length'][whichmax]
                    attr['geometry'] = attr['geometry'][whichmax]
                    length = attr['length']

                # If there is only one edge to this centroid, then we can just
                # reuse the existing attributes
                else:
                    edge = list(out_edges_neighbor)[0]
                    attr = G.edges[edge[0], edge[1], 0]
                    length = attr['length']

                # fix edge geometries - make sure they intersect their centroids,
                # otherwise we may have gaps on the resulting graph when plotted
                neighbor_point = Point(
                    intersection_centroids[neighbor].x,
                    intersection_centroids[neighbor].y)

                # make sure the resulting geometry intersects this centroid
                coords = attr['geometry'].coords[:]
                if not attr['geometry'].intersects(centroid_point):
                    coords = centroid_point.coords[:] + coords

                # make sure the resulting geometry intersects the neighbor's centroid
                if not attr['geometry'].intersects(neighbor_point):
                    coords = coords + neighbor_point.coords[:]

                # update geometry
                attr['geometry'] = LineString(coords)
                # update length based on geometry
                attr['length'] = attr['geometry'].length

                # Is the new length larger than expected? theoretically this
                # shouldn't happen, but there are a few instances where it does.
                # For now, we just log these cases
                if attr['length'] > (length + 2 * tolerance):
                    log("New edge ({},{}) has length larger than expected: {:,.1f} > {:,.1f} + {}."\
                            .format(centroid, neighbor, attr['length'], length, 2*tolerance),
                        level = lg.WARNING)

                G__.add_edge(centroid, neighbor, **attr)

    msg = 'Cleaned intersections of graph (from {:,} to {:,} nodes and from {:,} to {:,} edges) in {:,.2f} seconds'
    log(msg.format(len(list(G.nodes())), len(list(G__.nodes())),
                   len(list(G.edges())), len(list(G__.edges())), time.time()-start_time))

    return nx.MultiDiGraph(G__)
