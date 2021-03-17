"""Calculate distances and shortest paths and find nearest node/edge(s) to point(s)."""

import itertools
import warnings

import networkx as nx
import numpy as np
import pandas as pd
from rtree.index import Index as RTreeIndex
from scipy.spatial import cKDTree
from shapely.geometry import Point

from . import projection
from . import utils_geo
from . import utils_graph

# scikit-learn is optional dependency for unprojected nearest-neighbor search
try:
    from sklearn.neighbors import BallTree
except ImportError:
    BallTree = None

EARTH_RADIUS_M = 6_371_009


def great_circle_vec(lat1, lng1, lat2, lng2, earth_radius=EARTH_RADIUS_M):
    """
    Calculate great-circle distances between pairs of points.

    Vectorized function to calculate the great-circle distance between two
    points' coordinates or between arrays of points' coordinates using the
    haversine formula. Expects coordinates in decimal degrees.

    Parameters
    ----------
    lat1 : float or numpy.array of float
        first point's latitude coordinate
    lng1 : float or numpy.array of float
        first point's longitude coordinate
    lat2 : float or numpy.array of float
        second point's latitude coordinate
    lng2 : float or numpy.array of float
        second point's longitude coordinate
    earth_radius : int or float
        radius of earth in units in which distance will be returned
        (default is meters)

    Returns
    -------
    dist : float or numpy.array of float
        distance from each (lat1, lng1) to each (lat2, lng2) in units of
        earth_radius
    """
    y1 = np.deg2rad(lat1)
    y2 = np.deg2rad(lat2)
    dy = y2 - y1

    x1 = np.deg2rad(lng1)
    x2 = np.deg2rad(lng2)
    dx = x2 - x1

    h = np.sin(dy / 2) ** 2 + np.cos(y1) * np.cos(y2) * np.sin(dx / 2) ** 2
    h = np.minimum(1, h)  # protect against floating point errors
    arc = 2 * np.arcsin(np.sqrt(h))

    # return distance in units of earth_radius
    return arc * earth_radius


def euclidean_dist_vec(y1, x1, y2, x2):
    """
    Calculate Euclidean distances between pairs of points.

    Vectorized function to calculate the Euclidean distance between two
    points' coordinates or between arrays of points' coordinates. For accurate
    results, use projected coordinates rather than decimal degrees.

    Parameters
    ----------
    y1 : float or numpy.array of float
        first point's y coordinate
    x1 : float or numpy.array of float
        first point's x coordinate
    y2 : float or numpy.array of float
        second point's y coordinate
    x2 : float or numpy.array of float
        second point's x coordinate

    Returns
    -------
    dist : float or numpy.array of float
        distance from each (x1, y1) to each (x2, y2) in coordinates' units
    """
    # pythagorean theorem
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5


def nearest_nodes(G, X, Y, return_dist=False):
    """
    Find the nearest node(s) to some point(s).

    If the graph is projected, this uses a k-d tree for euclidean nearest
    neighbor search (the fastest method). If it is unprojected, this uses a
    ball tree for haversine nearest neighbor search, which requires that
    scikit-learn is installed as an optional dependency.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        graph in which to find nearest nodes
    X : list
        points' x or longitude coordinates, in same CRS/units as graph and
        containing no nulls
    Y : list
        points' y or latitude coordinates, in same CRS/units as graph and
        containing no nulls
    return_dist : bool
        optionally also return distance between points and nearest nodes

    Returns
    -------
    nn or (nn, dist): numpy.array or tuple of numpy.array
        nearest node IDs or optionally a tuple of arrays where `dist` contains
        distances between the points and their nearest nodes
    """
    if pd.Series(X).isna().any() or pd.Series(Y).isna().any():
        raise ValueError("`X` and `Y` cannot contain nulls")
    nodes = utils_graph.graph_to_gdfs(G, edges=False, node_geometry=False)[["x", "y"]]

    if projection.is_projected(G.graph["crs"]):
        # if projected, use k-d tree for euclidean nearest-neighbor search
        dist, pos = cKDTree(nodes).query(np.array([X, Y]).T, k=1)
        nn = nodes.index[pos].values

    else:
        # if unprojected, use ball tree for haversine nearest-neighbor search
        if BallTree is None:
            raise ImportError("scikit-learn must be installed to search an unprojected graph")
        # haversine requires lat, lng coords in radians
        nodes_rad = np.deg2rad(nodes[["y", "x"]])
        points_rad = np.deg2rad(np.array([Y, X]).T)
        dist, pos = BallTree(nodes_rad, metric="haversine").query(points_rad, k=1)
        dist = dist[:, 0] * EARTH_RADIUS_M  # convert radians -> meters
        nn = nodes.index[pos[:, 0]].values

    if return_dist:
        return nn, dist
    else:
        return nn


def nearest_edges(G, X, Y, interpolate=None, return_dist=False):
    """
    Find the nearest edge(s) to some point(s).

    If `interpolate` is None, search for the nearest edge to each point, one
    at a time, using an r-tree and minimizing the euclidean distances from the
    point to the possible matches. For accuracy, use a projected graph and
    points. This method is precise and also fastest if searching for few
    points relative to the graph's size.

    For a faster method if searching for many points relative to the graph's
    size, use the `interpolate` argument to interpolate points along the edges
    and index them. If the graph is projected, this uses a k-d tree for
    euclidean nearest neighbor search. If graph is unprojected, this uses a
    ball tree for haversine nearest neighbor search, which requires that
    scikit-learn is installed as an optional dependency.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        graph in which to find nearest edges
    X : list
        points' x or longitude coordinates, in same CRS/units as graph and
        containing no nulls
    Y : list
        points' y or latitude coordinates, in same CRS/units as graph and
        containing no nulls
    interpolate : float
        spacing distance between interpolated points, in same units as graph.
        smaller values generate more points.
    return_dist : bool
        optionally also return distance between points and nearest edges

    Returns
    -------
    ne or (ne, dist): numpy.array or tuple of numpy.array
        nearest edges as [u, v, key] or optionally a tuple of arrays where
        `dist` contains distances between the points and their nearest edges
    """
    if (pd.isnull(X) | pd.isnull(Y)).any():
        raise ValueError("`X` and `Y` cannot contain nulls")

    geoms = utils_graph.graph_to_gdfs(G, nodes=False)["geometry"]

    # if no interpolation distance was provided
    if interpolate is None:

        # build the r-tree spatial index by position for subsequent iloc
        rtree = RTreeIndex()
        for pos, bounds in enumerate(geoms.bounds.values):
            rtree.insert(pos, bounds)

        # use r-tree to find possible nearest neighbors, one point at a time,
        # then minimize euclidean distance from point to the possible matches
        ne_dist = list()
        for xy in zip(X, Y):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                dists = geoms.iloc[list(rtree.nearest(xy))].distance(Point(xy))
            ne_dist.append((dists.idxmin(), dists.min()))
        ne, dist = zip(*ne_dist)

    # otherwise, if interpolation distance was provided
    else:

        # interpolate points along edges to index with a k-d tree or ball tree
        uvk_xy = list()
        for uvk, geom in zip(geoms.index, geoms.values):
            uvk_xy.extend((uvk, xy) for xy in utils_geo.interpolate_points(geom, interpolate))
        labels, xy = zip(*uvk_xy)
        vertices = pd.DataFrame(xy, index=labels, columns=["x", "y"])

        if projection.is_projected(G.graph["crs"]):
            # if projected, use k-d tree for euclidean nearest-neighbor search
            dist, pos = cKDTree(vertices).query(np.array([X, Y]).T, k=1)
            ne = vertices.index[pos]

        else:
            # if unprojected, use ball tree for haversine nearest-neighbor search
            if BallTree is None:
                raise ImportError("scikit-learn must be installed to search an unprojected graph")
            # haversine requires lat, lng coords in radians
            vertices_rad = np.deg2rad(vertices[["y", "x"]])
            points_rad = np.deg2rad(np.array([Y, X]).T)
            dist, pos = BallTree(vertices_rad, metric="haversine").query(points_rad, k=1)
            dist = dist[:, 0] * EARTH_RADIUS_M  # convert radians -> meters
            ne = vertices.index[pos[:, 0]]

    if return_dist:
        return np.array(ne), np.array(dist)
    else:
        return np.array(ne)


def get_nearest_node(G, point, method="haversine", return_dist=False):
    """
    Find the nearest node to a point.

    Return the graph node nearest to some (lat, lng) or (y, x) point and
    optionally the distance between the node and the point. This function can
    use either the haversine formula or Euclidean distance.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    point : tuple
        The (lat, lng) or (y, x) point for which we will find the nearest node
        in the graph
    method : string {'haversine', 'euclidean'}
        Which method to use for calculating distances to find nearest node.
        If 'haversine', graph nodes' coordinates must be in units of decimal
        degrees. If 'euclidean', graph nodes' coordinates must be projected.
    return_dist : bool
        Optionally also return the distance (in meters if haversine, or graph
        node coordinate units if euclidean) between the point and the nearest
        node

    Returns
    -------
    int or tuple of (int, float)
        Nearest node ID or optionally a tuple of (node ID, dist), where dist
        is the distance (in meters if haversine, or graph node coordinate
        units if euclidean) between the point and nearest node
    """
    nn, dist = nearest_nodes(G, X=[point[1]], Y=[point[0]], return_dist=True)
    if return_dist:
        return nn[0], dist[0]
    else:
        return nn[0]


def get_nearest_edge(G, point, return_geom=False, return_dist=False):
    """
    Find the nearest edge to a point by minimum Euclidean distance.

    For best results, both G and point should be projected.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    point : tuple
        the (lat, lng) or (y, x) point for which we will find the nearest edge
        in the graph
    return_geom : bool
        Optionally return the geometry of the nearest edge
    return_dist : bool
        Optionally return the distance in graph's coordinates' units between
        the point and the nearest edge

    Returns
    -------
    tuple
        Graph edge unique identifier as a tuple of (u, v, key).
        Or a tuple of (u, v, key, geom) if return_geom is True.
        Or a tuple of (u, v, key, dist) if return_dist is True.
        Or a tuple of (u, v, key, geom, dist) if return_geom and return_dist are True.
    """
    ne, dist = nearest_edges(G, X=[point[1]], Y=[point[0]], return_dist=True)
    u, v, key = ne[0]
    geom = utils_graph.graph_to_gdfs(G, nodes=False).loc[(u, v, key), "geometry"]
    if return_dist and return_geom:
        return u, v, key, geom, dist[0]
    elif return_dist:
        return u, v, key, dist[0]
    elif return_geom:
        return u, v, key, geom
    else:
        return u, v, key


def get_nearest_nodes(G, X, Y, method=None, return_dist=False):
    """
    Find the nearest node to each point in a list of points.

    Pass in points as separate lists of X and Y coordinates. The 'kdtree'
    method is by far the fastest with large data sets, but only finds
    approximate nearest nodes if working in unprojected coordinates like
    lat-lng (it precisely finds the nearest node if working in projected
    coordinates). The 'balltree' method is second fastest with large data
    sets but it is precise if working in unprojected coordinates like lat-lng.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    X : list-like
        the longitudes or x coordinates for which we will find the nearest
        node in the graph
    Y : list-like
        the latitudes or y coordinates for which we will find the nearest node
        in the graph
    method : string {None, 'kdtree', 'balltree'}
        Which method to use for finding the nearest node to each point. If
        None, we manually find each node one at a time using
        utils.get_nearest_node and haversine. If 'kdtree' we use
        scipy.spatial.cKDTree for very fast euclidean search. If
        'balltree', we use sklearn.neighbors.BallTree for fast
        haversine search.
    return_dist : bool
        Optionally also return the distance (in meters if haversine, or graph
        node coordinate units if euclidean) between the point and the nearest
        node

    Returns
    -------
    nn or (dist, nn): np.array or tuple of (np.array<float>, np.array<int>)
        array of node IDs representing the node nearest to each point in the
        passed-in list of points, or optionally tuple of arrays (dist,nn) where dist
        is an array of the distances (in meters if haversine, or graph node coordinate
        units if euclidean) between the point and nearest node
    """
    return nearest_nodes(G, X=X, Y=Y, return_dist=return_dist)


def get_nearest_edges(G, X, Y, method=None, dist=0.0001):
    """
    Find the nearest edge to each point in a list of points.

    Pass in points as separate lists of X and Y coordinates. The 'kdtree'
    method is by far the fastest with large data sets, but only finds
    approximate nearest edges if working in unprojected coordinates like
    lat-lng (it precisely finds the nearest edge if working in projected
    coordinates). The 'balltree' method is second fastest with large data
    sets, but it is precise if working in unprojected coordinates like
    lat-lng. As a rule of thumb, if you have a small graph just use
    method=None. If you have a large graph with lat-lng coordinates, use
    method='balltree'. If you have a large graph with projected coordinates,
    use method='kdtree'. Note that if you are working in units of lat-lng,
    the X vector corresponds to longitude and the Y vector corresponds
    to latitude. The method creates equally distanced points along the edges
    of the network. Then, these points are used in a kdTree or BallTree search
    to identify which is nearest. Note that this method will not give exact
    perpendicular point along the edge, but the smaller the *dist* parameter,
    the closer (but slower) the solution will be.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    X : list-like
        the longitudes or x coordinates for which we will find the nearest
        edge in the graph. For projected graphs use the projected coordinates,
        usually in meters.
    Y : list-like
        the latitudes or y coordinates for which we will find the nearest edge
        in the graph. For projected graphs use the projected coordinates,
        usually in meters.
    method : string {None, 'kdtree', 'balltree'}
        Which method to use for finding nearest edge to each point. If None,
        we manually find each edge one at a time using get_nearest_edge. If
        'kdtree' we use scipy.spatial.cKDTree for very fast euclidean search.
        Recommended for projected graphs. If 'balltree', we use
        sklearn.neighbors.BallTree for fast haversine search. Recommended for
        unprojected graphs.
    dist : float
        spacing length along edges. Units are the same as the graph's
        geometries. The smaller the value, the more points are created.

    Returns
    -------
    ne : np.array
        array of edge IDs representing the edge nearest to each point in the
        passed-in list of points. Edge IDs are represented by u, v, key where
        u and v the node IDs of the nodes the edge links.
    """
    return nearest_edges(G, X, Y)


def shortest_path(G, orig, dest, weight="length"):
    """
    Get shortest path from origin node to destination node.

    See also `k_shortest_paths` to get multiple shortest paths.

    This function is a convenience wrapper around networkx.shortest_path. For
    more functionality or different algorithms, use networkx directly.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    orig : int
        origin node ID
    dest : int
        destination node ID
    weight : string
        edge attribute to minimize when solving shortest path. default is edge
        length in meters.

    Returns
    -------
    path : list
        list of node IDs consituting the shortest path
    """
    path = nx.shortest_path(G, orig, dest, weight=weight)
    return path


def k_shortest_paths(G, orig, dest, k, weight="length"):
    """
    Get `k` shortest paths from origin node to destination node.

    See also `shortest_path` to get just the one shortest path.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    orig : int
        origin node ID
    dest : int
        destination node ID
    k : int
        number of shortest paths to get
    weight : string
        edge attribute to minimize when solving shortest paths. default is
        edge length in meters.

    Returns
    -------
    paths : generator
        a generator of `k` shortest paths ordered by total weight. each path
        is a list of node IDs.
    """
    paths_gen = nx.shortest_simple_paths(utils_graph.get_digraph(G, weight), orig, dest, weight)
    for path in itertools.islice(paths_gen, 0, k):
        yield path
