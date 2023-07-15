"""Calculate distances and shortest paths and find nearest node/edge(s) to point(s)."""

import itertools
import multiprocessing as mp
from warnings import warn

import networkx as nx
import numpy as np
import pandas as pd
from shapely.geometry import Point
from shapely.strtree import STRtree

from . import projection
from . import utils
from . import utils_geo
from . import utils_graph

# scipy is optional dependency for projected nearest-neighbor search
try:
    from scipy.spatial import cKDTree
except ImportError:  # pragma: no cover
    cKDTree = None

# scikit-learn is optional dependency for unprojected nearest-neighbor search
try:
    from sklearn.neighbors import BallTree
except ImportError:  # pragma: no cover
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
    earth_radius : float
        earth's radius in units in which distance will be returned (default is
        meters)

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


def add_edge_lengths(G, precision=None, edges=None):
    """
    Add `length` attribute (in meters) to each edge.

    Vectorized function to calculate great-circle distance between each edge's
    incident nodes. Ensure graph is in unprojected coordinates, and
    unsimplified to get accurate distances.

    Note: this function is run by all the `graph.graph_from_x` functions
    automatically to add `length` attributes to all edges. It calculates edge
    lengths as the great-circle distance from node `u` to node `v`. When
    OSMnx automatically runs this function upon graph creation, it does it
    before simplifying the graph: thus it calculates the straight-line lengths
    of edge segments that are themselves all straight. Only after
    simplification do edges take on a (potentially) curvilinear geometry. If
    you wish to calculate edge lengths later, you are calculating
    straight-line distances which necessarily ignore the curvilinear geometry.
    You only want to run this function on a graph with all straight edges
    (such as is the case with an unsimplified graph).

    Parameters
    ----------
    G : networkx.MultiDiGraph
        unprojected, unsimplified input graph
    precision : int
        deprecated, do not use
    edges : tuple
        tuple of (u, v, k) tuples representing subset of edges to add length
        attributes to. if None, add lengths to all edges.

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with edge length attributes
    """
    if precision is None:
        precision = 3
    else:
        warn(
            "the `precision` parameter is deprecated and will be removed in a future release",
            stacklevel=2,
        )

    if edges is None:
        uvk = tuple(G.edges)
    else:
        uvk = edges

    # extract edge IDs and corresponding coordinates from their nodes
    x = G.nodes(data="x")
    y = G.nodes(data="y")
    try:
        # two-dimensional array of coordinates: y0, x0, y1, x1
        c = np.array([(y[u], x[u], y[v], x[v]) for u, v, k in uvk])
        # ensure all coordinates can be converted to float and are non-null
        assert not np.isnan(c.astype(float)).any()
    except (AssertionError, KeyError) as e:  # pragma: no cover
        msg = "some edges missing nodes, possibly due to input data clipping issue"
        raise ValueError(msg) from e

    # calculate great circle distances, round, and fill nulls with zeros
    dists = great_circle_vec(c[:, 0], c[:, 1], c[:, 2], c[:, 3]).round(precision)
    dists[np.isnan(dists)] = 0
    nx.set_edge_attributes(G, values=dict(zip(uvk, dists)), name="length")

    utils.log("Added length attributes to graph edges")
    return G


def nearest_nodes(G, X, Y, return_dist=False):
    """
    Find the nearest node to a point or to each of several points.

    If `X` and `Y` are single coordinate values, this will return the nearest
    node to that point. If `X` and `Y` are lists of coordinate values, this
    will return the nearest node to each point.

    If the graph is projected, this uses a k-d tree for euclidean nearest
    neighbor search, which requires that scipy is installed as an optional
    dependency. If it is unprojected, this uses a ball tree for haversine
    nearest neighbor search, which requires that scikit-learn is installed as
    an optional dependency.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        graph in which to find nearest nodes
    X : float or list
        points' x (longitude) coordinates, in same CRS/units as graph and
        containing no nulls
    Y : float or list
        points' y (latitude) coordinates, in same CRS/units as graph and
        containing no nulls
    return_dist : bool
        optionally also return distance between points and nearest nodes

    Returns
    -------
    nn or (nn, dist) : int/list or tuple
        nearest node IDs or optionally a tuple where `dist` contains distances
        between the points and their nearest nodes
    """
    is_scalar = False
    if not (hasattr(X, "__iter__") and hasattr(Y, "__iter__")):
        # make coordinates arrays if user passed non-iterable values
        is_scalar = True
        X = np.array([X])
        Y = np.array([Y])

    if np.isnan(X).any() or np.isnan(Y).any():  # pragma: no cover
        msg = "`X` and `Y` cannot contain nulls"
        raise ValueError(msg)
    nodes = utils_graph.graph_to_gdfs(G, edges=False, node_geometry=False)[["x", "y"]]

    if projection.is_projected(G.graph["crs"]):
        # if projected, use k-d tree for euclidean nearest-neighbor search
        if cKDTree is None:  # pragma: no cover
            msg = "scipy must be installed to search a projected graph"
            raise ImportError(msg)
        dist, pos = cKDTree(nodes).query(np.array([X, Y]).T, k=1)
        nn = nodes.index[pos]

    else:
        # if unprojected, use ball tree for haversine nearest-neighbor search
        if BallTree is None:  # pragma: no cover
            msg = "scikit-learn must be installed to search an unprojected graph"
            raise ImportError(msg)
        # haversine requires lat, lng coords in radians
        nodes_rad = np.deg2rad(nodes[["y", "x"]])
        points_rad = np.deg2rad(np.array([Y, X]).T)
        dist, pos = BallTree(nodes_rad, metric="haversine").query(points_rad, k=1)
        dist = dist[:, 0] * EARTH_RADIUS_M  # convert radians -> meters
        nn = nodes.index[pos[:, 0]]

    # convert results to correct types for return
    nn = nn.tolist()
    dist = dist.tolist()
    if is_scalar:
        nn = nn[0]
        dist = dist[0]

    if return_dist:
        return nn, dist
    else:
        return nn


def nearest_edges(G, X, Y, interpolate=None, return_dist=False):
    """
    Find the nearest edge to a point or to each of several points.

    If `X` and `Y` are single coordinate values, this will return the nearest
    edge to that point. If `X` and `Y` are lists of coordinate values, this
    will return the nearest edge to each point. This function uses an R-tree
    spatial index and minimizes the euclidean distance from each point to the
    possible matches. For accurate results, use a projected graph and points.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        graph in which to find nearest edges
    X : float or list
        points' x (longitude) coordinates, in same CRS/units as graph and
        containing no nulls
    Y : float or list
        points' y (latitude) coordinates, in same CRS/units as graph and
        containing no nulls
    interpolate : float
        deprecated, do not use
    return_dist : bool
        optionally also return distance between points and nearest edges

    Returns
    -------
    ne or (ne, dist) : tuple or list
        nearest edges as (u, v, key) or optionally a tuple where `dist`
        contains distances between the points and their nearest edges
    """
    is_scalar = False
    if not (hasattr(X, "__iter__") and hasattr(Y, "__iter__")):
        # make coordinates arrays if user passed non-iterable values
        is_scalar = True
        X = np.array([X])
        Y = np.array([Y])

    if np.isnan(X).any() or np.isnan(Y).any():  # pragma: no cover
        msg = "`X` and `Y` cannot contain nulls"
        raise ValueError(msg)
    geoms = utils_graph.graph_to_gdfs(G, nodes=False)["geometry"]

    # if no interpolation distance was provided
    if interpolate is None:
        # build an r-tree spatial index by position for subsequent iloc
        rtree = STRtree(geoms)

        # use the r-tree to find each point's nearest neighbor and distance
        points = [Point(xy) for xy in zip(X, Y)]
        pos, dist = rtree.query_nearest(points, all_matches=False, return_distance=True)

        # if user passed X/Y lists, the 2nd subarray contains geom indices
        if len(pos.shape) > 1:
            pos = pos[1]
        ne = geoms.iloc[pos].index

    # otherwise, if interpolation distance was provided
    else:
        warn(
            "The `interpolate` parameter has been deprecated and will be removed in a future release",
            stacklevel=2,
        )

        # interpolate points along edges to index with k-d tree or ball tree
        uvk_xy = []
        for uvk, geom in zip(geoms.index, geoms.values):
            uvk_xy.extend((uvk, xy) for xy in utils_geo.interpolate_points(geom, interpolate))
        labels, xy = zip(*uvk_xy)
        vertices = pd.DataFrame(xy, index=labels, columns=["x", "y"])

        if projection.is_projected(G.graph["crs"]):
            # if projected, use k-d tree for euclidean nearest-neighbor search
            if cKDTree is None:  # pragma: no cover
                msg = "scipy must be installed to search a projected graph"
                raise ImportError(msg)
            dist, pos = cKDTree(vertices).query(np.array([X, Y]).T, k=1)
            ne = vertices.index[pos]

        else:
            # if unprojected, use ball tree for haversine nearest-neighbor search
            if BallTree is None:  # pragma: no cover
                msg = "scikit-learn must be installed to search an unprojected graph"
                raise ImportError(msg)
            # haversine requires lat, lng coords in radians
            vertices_rad = np.deg2rad(vertices[["y", "x"]])
            points_rad = np.deg2rad(np.array([Y, X]).T)
            dist, pos = BallTree(vertices_rad, metric="haversine").query(points_rad, k=1)
            dist = dist[:, 0] * EARTH_RADIUS_M  # convert radians -> meters
            ne = vertices.index[pos[:, 0]]

    # convert results to correct types for return
    ne = list(ne)
    dist = list(dist)
    if is_scalar:
        ne = ne[0]
        dist = dist[0]

    if return_dist:
        return ne, dist
    else:
        return ne


def _single_shortest_path(G, orig, dest, weight):
    """
    Solve the shortest path from an origin node to a destination node.

    This function is a convenience wrapper around networkx.shortest_path, with
    exception handling for unsolvable paths. It uses Dijkstra's algorithm.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    orig : int
        origin node ID
    dest : int
        destination node ID
    weight : string
        edge attribute to minimize when solving shortest path

    Returns
    -------
    path : list
        list of node IDs constituting the shortest path
    """
    try:
        return nx.shortest_path(G, orig, dest, weight=weight, method="dijkstra")
    except nx.exception.NetworkXNoPath:  # pragma: no cover
        utils.log(f"Cannot solve path from {orig} to {dest}")
        return None


def shortest_path(G, orig, dest, weight="length", cpus=1):
    """
    Solve shortest path from origin node(s) to destination node(s).

    Uses Dijkstra's algorithm. If `orig` and `dest` are single node IDs, this
    will return a list of the nodes constituting the shortest path between
    them. If `orig` and `dest` are lists of node IDs, this will return a list
    of lists of the nodes constituting the shortest path between each
    origin-destination pair. If a path cannot be solved, this will return None
    for that path. You can parallelize solving multiple paths with the `cpus`
    parameter, but be careful to not exceed your available RAM.

    See also `k_shortest_paths` to solve multiple shortest paths between a
    single origin and destination. For additional functionality or different
    solver algorithms, use NetworkX directly.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    orig : int or list
        origin node ID, or a list of origin node IDs
    dest : int or list
        destination node ID, or a list of destination node IDs
    weight : string
        edge attribute to minimize when solving shortest path
    cpus : int
        how many CPU cores to use; if None, use all available

    Returns
    -------
    path : list
        list of node IDs constituting the shortest path, or, if orig and dest
        are lists, then a list of path lists
    """
    _verify_edge_attribute(G, weight)

    # if neither orig nor dest is iterable, just return the shortest path
    if not (hasattr(orig, "__iter__") or hasattr(dest, "__iter__")):
        return _single_shortest_path(G, orig, dest, weight)

    # if both orig and dest are iterables, ensure they have same lengths
    elif hasattr(orig, "__iter__") and hasattr(dest, "__iter__"):
        if len(orig) != len(dest):  # pragma: no cover
            msg = "orig and dest must contain same number of elements"
            raise ValueError(msg)

        if cpus is None:
            cpus = mp.cpu_count()
        cpus = min(cpus, mp.cpu_count())
        utils.log(f"Solving {len(orig)} paths with {cpus} CPUs...")

        # if single-threading, calculate each shortest path one at a time
        if cpus == 1:
            paths = [_single_shortest_path(G, o, d, weight) for o, d in zip(orig, dest)]

        # if multi-threading, calculate shortest paths in parallel
        else:
            args = ((G, o, d, weight) for o, d in zip(orig, dest))
            pool = mp.Pool(cpus)
            sma = pool.starmap_async(_single_shortest_path, args)
            paths = sma.get()
            pool.close()
            pool.join()

        return paths

    # if only one of orig or dest is iterable and the other is not
    else:  # pragma: no cover
        msg = "orig and dest must either both be iterable or neither must be iterable"
        raise ValueError(msg)


def k_shortest_paths(G, orig, dest, k, weight="length"):
    """
    Solve `k` shortest paths from an origin node to a destination node.

    Uses Yen's algorithm. See also `shortest_path` to solve just the one
    shortest path.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    orig : int
        origin node ID
    dest : int
        destination node ID
    k : int
        number of shortest paths to solve
    weight : string
        edge attribute to minimize when solving shortest paths. default is
        edge length in meters.

    Yields
    ------
    path : list
        a generator of `k` shortest paths ordered by total weight. each path
        is a list of node IDs.
    """
    _verify_edge_attribute(G, weight)
    paths_gen = nx.shortest_simple_paths(utils_graph.get_digraph(G, weight), orig, dest, weight)
    yield from itertools.islice(paths_gen, 0, k)


def _verify_edge_attribute(G, attr):
    """
    Verify attribute values are numeric and non-null across graph edges.

    Raises a `ValueError` if attribute contains non-numeric values and raises
    a warning if attribute is missing or null on any edges.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    attr : string
        edge attribute to verify

    Returns
    -------
    None
    """
    try:
        values = np.array(tuple(G.edges(data=attr)))[:, 2]
        values_float = values.astype(float)
        if np.isnan(values_float).any():
            warn(f"The attribute {attr!r} is missing or null on some edges.", stacklevel=2)
    except ValueError as e:
        msg = f"The edge attribute {attr!r} contains non-numeric values."
        raise ValueError(msg) from e
