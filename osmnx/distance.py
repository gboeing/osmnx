"""Calculate distances and shortest paths and find nearest node/edge(s) to point(s)."""

import itertools
import multiprocessing as mp
import warnings

import networkx as nx
import numpy as np
import pandas as pd
from rtree.index import Index as RTreeIndex
from shapely.geometry import Point

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


def add_edge_lengths(G, precision=3):
    """
    Add `length` attribute (in meters) to each edge.

    Vectorized function to calculate great-circle distance between each edge's
    incident nodes. Ensure graph is in unprojected coordinates, and
    unsimplified to get accurate distances. Note: this function is run by all
    the `graph.graph_from_x` functions automatically to add `length`
    attributes to all edges.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        unprojected, unsimplified input graph
    precision : int
        decimal precision to round lengths

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with edge length attributes
    """
    # extract edge IDs and corresponding coordinates from their nodes
    uvk = tuple(G.edges)
    x = G.nodes(data="x")
    y = G.nodes(data="y")
    try:
        # two-dimensional array of coordinates: y0, x0, y1, x1
        c = np.array([(y[u], x[u], y[v], x[v]) for u, v, k in uvk])
    except KeyError:  # pragma: no cover
        raise KeyError("some edges missing nodes, possibly due to input data clipping issue")

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
        raise ValueError("`X` and `Y` cannot contain nulls")
    nodes = utils_graph.graph_to_gdfs(G, edges=False, node_geometry=False)[["x", "y"]]

    if projection.is_projected(G.graph["crs"]):
        # if projected, use k-d tree for euclidean nearest-neighbor search
        if cKDTree is None:  # pragma: no cover
            raise ImportError("scipy must be installed to search a projected graph")
        dist, pos = cKDTree(nodes).query(np.array([X, Y]).T, k=1)
        nn = nodes.index[pos]

    else:
        # if unprojected, use ball tree for haversine nearest-neighbor search
        if BallTree is None:  # pragma: no cover
            raise ImportError("scikit-learn must be installed to search an unprojected graph")
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
    will return the nearest edge to each point.

    If `interpolate` is None, search for the nearest edge to each point, one
    at a time, using an r-tree and minimizing the euclidean distances from the
    point to the possible matches. For accuracy, use a projected graph and
    points. This method is precise and also fastest if searching for few
    points relative to the graph's size.

    For a faster method if searching for many points relative to the graph's
    size, use the `interpolate` argument to interpolate points along the edges
    and index them. If the graph is projected, this uses a k-d tree for
    euclidean nearest neighbor search, which requires that scipy is installed
    as an optional dependency. If graph is unprojected, this uses a ball tree
    for haversine nearest neighbor search, which requires that scikit-learn is
    installed as an optional dependency.

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
        spacing distance between interpolated points, in same units as graph.
        smaller values generate more points.
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
            dists = geoms.iloc[list(rtree.nearest(xy))].distance(Point(xy))
            ne_dist.append((dists.idxmin(), dists.min()))
        ne, dist = zip(*ne_dist)

    # otherwise, if interpolation distance was provided
    else:

        # interpolate points along edges to index with k-d tree or ball tree
        uvk_xy = list()
        for uvk, geom in zip(geoms.index, geoms.values):
            uvk_xy.extend((uvk, xy) for xy in utils_geo.interpolate_points(geom, interpolate))
        labels, xy = zip(*uvk_xy)
        vertices = pd.DataFrame(xy, index=labels, columns=["x", "y"])

        if projection.is_projected(G.graph["crs"]):
            # if projected, use k-d tree for euclidean nearest-neighbor search
            if cKDTree is None:  # pragma: no cover
                raise ImportError("scipy must be installed to search a projected graph")
            dist, pos = cKDTree(vertices).query(np.array([X, Y]).T, k=1)
            ne = vertices.index[pos]

        else:
            # if unprojected, use ball tree for haversine nearest-neighbor search
            if BallTree is None:  # pragma: no cover
                raise ImportError("scikit-learn must be installed to search an unprojected graph")
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


def get_nearest_node(G, point, method=None, return_dist=False):
    """
    Do not use, deprecated.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        deprecated, do not use
    point : tuple
        deprecated, do not use
    method : string
        deprecated, do not use
    return_dist : bool
        deprecated, do not use

    Returns
    -------
    int or tuple
    """
    msg = (
        "The `get_nearest_node` function has been deprecated and will be removed in a "
        "future release. Use the more efficient `distance.nearest_nodes` instead."
    )
    warnings.warn(msg)
    nn, dist = nearest_nodes(G, X=[point[1]], Y=[point[0]], return_dist=True)
    if return_dist:
        return nn[0], dist[0]
    else:
        return nn[0]


def get_nearest_edge(G, point, return_geom=False, return_dist=False):
    """
    Do not use, deprecated.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        deprecated, do not use
    point : tuple
        deprecated, do not use
    return_geom : bool
        deprecated, do not use
    return_dist : bool
        deprecated, do not use

    Returns
    -------
    tuple
    """
    msg = (
        "The `get_nearest_edge` function has been deprecated and will be removed in a "
        "future release. Use the more efficient `distance.nearest_edges` instead."
    )
    warnings.warn(msg)
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
    Do not use, deprecated.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        deprecated, do not use
    X : list
        deprecated, do not use
    Y : list
        deprecated, do not use
    method : string
        deprecated, do not use
    return_dist : bool
        deprecated, do not use

    Returns
    -------
    numpy.array or tuple of numpy.array
    """
    msg = (
        "The `get_nearest_nodes` function has been deprecated and will be removed in a "
        "future release. Use the more efficient `distance.nearest_nodes` instead."
    )
    warnings.warn(msg)
    return nearest_nodes(G, X=X, Y=Y, return_dist=return_dist)


def get_nearest_edges(G, X, Y, method=None, dist=None):
    """
    Do not use, deprecated.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        deprecated, do not use
    X : list-like
        deprecated, do not use
    Y : list-like
        deprecated, do not use
    method : string
        deprecated, do not use
    dist : float
        deprecated, do not use

    Returns
    -------
    numpy.array
    """
    msg = (
        "The `get_nearest_edges` function has been deprecated and will be removed in a "
        "future release. Use the more efficient `distance.nearest_edges` instead."
    )
    warnings.warn(msg)
    return nearest_edges(G, X, Y, dist)


def _single_shortest_path(G, orig, dest, weight):
    """
    Solve the shortest path from an origin node to a destination node.

    This function is a convenience wrapper around networkx.shortest_path, with
    exception handling for unsolvable paths.

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
        return nx.shortest_path(G, orig, dest, weight=weight)
    except nx.exception.NetworkXNoPath:  # pragma: no cover
        utils.log(f"Cannot solve path from {orig} to {dest}")
        return None


def shortest_path(G, orig, dest, weight="length", cpus=1):
    """
    Solve shortest path from origin node(s) to destination node(s).

    If `orig` and `dest` are single node IDs, this will return a list of the
    nodes constituting the shortest path between them.  If `orig` and `dest`
    are lists of node IDs, this will return a list of lists of the nodes
    constituting the shortest path between each origin-destination pair. If a
    path cannot be solved, this will return None for that path. You can
    parallelize solving multiple paths with the `cpus` parameter, but be
    careful to not exceed your available RAM.

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
    if not (hasattr(orig, "__iter__") or hasattr(dest, "__iter__")):
        # if neither orig nor dest is iterable, just return the shortest path
        return _single_shortest_path(G, orig, dest, weight)

    elif hasattr(orig, "__iter__") and hasattr(dest, "__iter__"):
        # if both orig and dest are iterables ensure they have same lengths
        if len(orig) != len(dest):  # pragma: no cover
            raise ValueError("orig and dest must contain same number of elements")

        if cpus is None:
            cpus = mp.cpu_count()
        cpus = min(cpus, mp.cpu_count())
        utils.log(f"Solving {len(orig)} paths with {cpus} CPUs...")

        if cpus == 1:
            # if single-threading, calculate each shortest path one at a time
            paths = [_single_shortest_path(G, o, d, weight) for o, d in zip(orig, dest)]
        else:
            # if multi-threading, calculate shortest paths in parallel
            args = ((G, o, d, weight) for o, d in zip(orig, dest))
            pool = mp.Pool(cpus)
            sma = pool.starmap_async(_single_shortest_path, args)
            paths = sma.get()
            pool.close()
            pool.join()

        return paths

    else:
        # if only one of orig or dest is iterable and the other is not
        raise ValueError("orig and dest must either both be iterable or neither must be iterable")


def k_shortest_paths(G, orig, dest, k, weight="length"):
    """
    Solve `k` shortest paths from an origin node to a destination node.

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
        number of shortest paths to solve
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
