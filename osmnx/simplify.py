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
from .utils import get_nearest_nodes
from . import settings

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


def clean_intersections(G, tolerance=15, dead_ends=False, method='kdtree'):
    """
    Clean-up intersections comprising clusters of nodes by merging them and
    returning their centroids.

    Divided roads are represented by separate centerline edges. The intersection
    of two divided roads thus creates 4 nodes, representing where each edge
    intersects a perpendicular edge. These 4 nodes represent a single
    intersection in the real world. This function cleans them up by buffering
    their points to an arbitrary distance, merging overlapping buffers, and
    taking their centroid. Then, it constructs a new graph with the centroids
    as nodes. Edges between centroids are calculated such that the new graph
    remains topologically (approximately) equal to G.

    Cleaning occurs in the graph's current units, but the use of unprojected
    units is not recommended. For best results, the tolerance argument should be
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
    method : str {None, 'kdtree', 'balltree'}
        Which method to use for finding nearest node to each point.
        If None, we manually find each node one at a time using
        osmnx.utils.get_nearest_node and haversine. If 'kdtree' we use
        scipy.spatial.cKDTree for very fast euclidean search. If
        'balltree', we use sklearn.neighbors.BallTree for fast
        haversine search.

    Returns
    ----------
    G2 : networkx multidigraph
        the new cleaned graph
    """

    start_time = time.time()

    if G.graph['crs'] == settings.default_crs:
        log(("The graph seems to be using unprojected units which is not recommended."),
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

    log('Finding centroids (checkpoint 0-1) took {:,.2f} seconds'\
         .format(time.time()-start_time))

    start_time2 = time.time()

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

    # Let's first create the Graph first
    G__ = nx.MultiDiGraph(name = G.graph['name'],
                          crs  = G.graph['crs'])

    # And add the nodes (without the osmid attributes just yet)
    # The centroids are given new ids = 0 .. total_centroids-1
    for i in range(len(intersection_centroids)):
        G__.add_node(i,
                     x = intersection_centroids[i].x,
                     y = intersection_centroids[i].y)

    # We can now run ox.get_nearest_nodes to calculate the closest centroid
    # to each node
    osmids = [ osmid for osmid in G.nodes()           ]
    X      = [ x     for _, x  in G.nodes(data = 'x') ]
    Y      = [ y     for _, y  in G.nodes(data = 'y') ]

    # Mapping: i -> centroid
    nearest_centroid = get_nearest_nodes(G__, X, Y, method = method)

    # For ease of use, we build the mappings:
    # osmid -> centroid
    nearest_centroid = { osmids[i] : centroid for i, centroid in enumerate(nearest_centroid) }

    # centroid -> osmids
    centroid_elements = {i : [] for i in range(len(intersection_centroids))}
    for osmid, centroid in nearest_centroid.items():
        centroid_elements[centroid].append(osmid)

    # Add osmids as attribute
    nx.set_node_attributes(G__, centroid_elements, 'osmids')
    # For compability, set attribute the 'osmid' with a single value
    nx.set_node_attributes(G__, { centroid : elements[0] for centroid, elements in centroid_elements.items() } , 'osmid')

    log('Getting nearest nodes and computing the mapping and reverse mapping of '
        'centroid->elements (checkpoint 1-2) took {:,.2f} seconds'\
        .format(time.time()-start_time2))

    start_time2 = time.time()

    # --------
    # Independently of the cases below, in the main loop, we need to
    # fix edge geometries and make sure they intersect their centroids,
    # otherwise we may have gaps on the resulting graph when plotted.
    # --------
    # We define this behavior as an inner function here because it will be
    # extensively used in the main loop, thus avoiding clutter in the main loop.
    # --------
    def fix_geometry(centroid_p, neighbor_p, geom):
        # make sure the resulting geometry intersects the centroid
        coords = geom.coords[:]
        if not geom.intersects(centroid_point):
            coords = centroid_point.coords[:] + coords

        # make sure the resulting geometry intersects the neighbor's centroid
        if not geom.intersects(neighbor_point):
            coords = coords + neighbor_point.coords[:]

        # return updated geometry
        return LineString(coords)

    # Now we're ready to start adding the edges to the cleaned graph
    for centroid, elements in centroid_elements.items():
            # get the outgoing edges for all elements in the centroid
            out_edges = []
            for node in elements:
                out_edges += list(G.out_edges(node))

            # select the neighbor nodes of each element in the centroid
            out_nodes = map(lambda edge: edge[1], out_edges)

            # and get the corresponding centroid of each element
            neighbors = list(map(lambda osmid: nearest_centroid[osmid],
                                 out_nodes))

            # remove duplicate neighbor centroids
            neighbors = set(neighbors)
            # and self-loops
            if centroid in neighbors:
                neighbors.remove(centroid)

            # helper for fixing geometries
            # initialised outside the for loop even thoug it's used later on
            # to avoid multiple redundant calls
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

                # We're gonna need this to fix any problems with resulting geometries
                neighbor_point = Point(
                    intersection_centroids[neighbor].x,
                    intersection_centroids[neighbor].y)

                # --------
                # Cases 1 and 2:
                # --------
                # This centroid and my neighbor centroid have only one element each,
                # OR
                # There are multiple elements in one or both of this centroid and
                # its neighbor, but only a single edge exists between them.
                #
                # This means that we can keep any parallel edges between them and
                # reuse any existing attributes. Case 2 may still need corrections
                # for the resulting geometry (see below, after case 3)
                if (len(elements) == 1 and len(neighbor_elements) == 1) or \
                   (len(out_edges_neighbor) == 1):
                   # There is only one edge in either case
                   u, v = list(out_edges_neighbor)[0]
                   # Iterate through existing keys
                   for key in list(G[u][v].keys()):
                       attr = G.edges[u, v, key]
                       # Case 1 does not need to fix_geometry but case 2 might.
                       # fix_geometry does not affect an already correct geometry.
                       attr['geometry'] = fix_geometry(centroid_point,
                                                       neighbor_point,
                                                       attr['geometry'])
                       # update length based on geometry
                       attr['length'] = attr['geometry'].length
                       # Time to add the edge
                       G__.add_edge(centroid, neighbor, key = key, **attr)

                # --------
                # Case 3:
                # --------
                # There are multiple elements in one or both of this centroid and
                # its neighbor, and multiple edges between them.
                # The result is a single edge in the new graph, whose attributes
                # are combined among the existing edges, except for geometry
                # and length, whose value must be a single value (and not a list)
                else:
                    edges_attrs = [ G.edges[edge[0], edge[1], 0] for edge in out_edges_neighbor ]

                    # Not every edge has the same set of keys, so we compute the union of the individual
                    # sets, and also count their occurrence for initialisation purposes
                    all_keys = []
                    for edge_attr in edges_attrs:
                        all_keys += edge_attr.keys()

                    key_counter = collections.Counter(all_keys)

                    # If if the key counter is greater than 1, we initialise the new
                    # dictionary with empty lists, so that we can later call append()
                    attr = {key : [] if count > 1 else '' for key, count in key_counter.items()}

                    for edge_attr in edges_attrs:
                        for key in edge_attr.keys():
                            val = edge_attr[key]

                            # This key only occurs once, so we can set its value
                            if key_counter[key] == 1:
                                attr[key] = val
                            # Otherwise, append
                            elif isinstance(val, list):
                                attr[key].extend(val)
                            else:
                                attr[key].append(val)

                    # We pick the geometry with largest length
                    # but, there may be a better way to make this decision
                    whichmax = np.argmax(attr['length'])
                    geom = attr['geometry'][whichmax]
                    attr['geometry'] = fix_geometry(centroid_point, neighbor_point, geom)

                    # update length based on geometry
                    attr['length'] = attr['geometry'].length

                    # Time to add the edge
                    G__.add_edge(centroid, neighbor, **attr)

    log('Adding edges to the new graph (checkpoint 2-3) took {:,.2f} seconds'\
        .format(time.time()-start_time2))

    msg = 'Cleaned intersections of graph (from {:,} to {:,} nodes and from {:,} to {:,} edges) in {:,.2f} seconds'
    log(msg.format(len(list(G.nodes())), len(list(G__.nodes())),
                   len(list(G.edges())), len(list(G__.edges())), time.time()-start_time))

    return G__
