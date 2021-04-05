"""Calculate distances and shortest paths and find nearest node/edge(s) to point(s)."""

import itertools
import warnings

import networkx as nx
import numpy as np
import pandas as pd
from rtree.index import Index as RTreeIndex
from shapely.geometry import Point

from . import distance
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
        raise KeyError("some edges missing nodes, possibly due to input data clipping issue")

    # turn the coordinates into a DataFrame indexed by u, v, k
    cols = ["u", "v", "k", "u_y", "u_x", "v_y", "v_x"]
    df = pd.DataFrame(coords, columns=cols).set_index(["u", "v", "k"])

    # calculate great circle distances, fill nulls with zeros, then round
    dists = distance.great_circle_vec(df["u_y"], df["u_x"], df["v_y"], df["v_x"])
    dists = dists.fillna(value=0).round(precision)
    nx.set_edge_attributes(G, name="length", values=dists)

    utils.log("Added edge lengths to graph")
    return G


def nearest_nodes(G, X, Y, return_dist=False):
    """
    Find the nearest nodes to some points.

    If the graph is projected, this uses a k-d tree for euclidean nearest
    neighbor search, which requires that scipy is installed as an optional
    dependency. If it is unprojected, this uses a ball tree for haversine
    nearest neighbor search, which requires that scikit-learn is installed as
    an optional dependency.

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
    if pd.Series(X).isna().any() or pd.Series(Y).isna().any():  # pragma: no cover
        raise ValueError("`X` and `Y` cannot contain nulls")
    nodes = utils_graph.graph_to_gdfs(G, edges=False, node_geometry=False)[["x", "y"]]

    if projection.is_projected(G.graph["crs"]):
        # if projected, use k-d tree for euclidean nearest-neighbor search
        if cKDTree is None:  # pragma: no cover
            raise ImportError("scipy must be installed to search a projected graph")
        dist, pos = cKDTree(nodes).query(np.array([X, Y]).T, k=1)
        nn = nodes.index[pos].values

    else:
        # if unprojected, use ball tree for haversine nearest-neighbor search
        if BallTree is None:  # pragma: no cover
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
    euclidean nearest neighbor search, which requires that scipy is installed
    as an optional dependency. If graph is unprojected, this uses a ball tree
    for haversine nearest neighbor search, which requires that scikit-learn is
    installed as an optional dependency.

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
    if (pd.isnull(X) | pd.isnull(Y)).any():  # pragma: no cover
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

    if return_dist:
        return np.array(ne), np.array(dist)
    else:
        return np.array(ne)


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
