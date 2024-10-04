"""Calculate distances and find nearest graph node/edge(s) to point(s)."""

from __future__ import annotations

import logging as lg
from collections.abc import Iterable
from typing import Literal
from typing import overload

import networkx as nx
import numpy as np
import numpy.typing as npt
from shapely import Point
from shapely.strtree import STRtree

from . import convert
from . import projection
from . import utils

# scipy is optional dependency for projected nearest-neighbor search
try:
    from scipy.spatial import cKDTree
except ImportError:  # pragma: no cover
    cKDTree = None  # noqa: N816

# scikit-learn is optional dependency for unprojected nearest-neighbor search
try:
    from sklearn.neighbors import BallTree
except ImportError:  # pragma: no cover
    BallTree = None

EARTH_RADIUS_M = 6_371_009


# if coords are all floats, return float
@overload
def great_circle(lat1: float, lon1: float, lat2: float, lon2: float) -> float: ...


# if coords are all floats (and optional arg is provided), return float
@overload
def great_circle(
    lat1: float,
    lon1: float,
    lat2: float,
    lon2: float,
    earth_radius: float,
) -> float: ...


# if coords are all arrays, return array
@overload
def great_circle(
    lat1: npt.NDArray[np.float64],
    lon1: npt.NDArray[np.float64],
    lat2: npt.NDArray[np.float64],
    lon2: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]: ...


# if coords are all arrays (and optional arg is provided), return array
@overload
def great_circle(
    lat1: npt.NDArray[np.float64],
    lon1: npt.NDArray[np.float64],
    lat2: npt.NDArray[np.float64],
    lon2: npt.NDArray[np.float64],
    earth_radius: float,
) -> npt.NDArray[np.float64]: ...


def great_circle(
    lat1: float | npt.NDArray[np.float64],
    lon1: float | npt.NDArray[np.float64],
    lat2: float | npt.NDArray[np.float64],
    lon2: float | npt.NDArray[np.float64],
    earth_radius: float = EARTH_RADIUS_M,
) -> float | npt.NDArray[np.float64]:
    """
    Calculate great-circle distances between pairs of points.

    Vectorized function to calculate the great-circle distance between two
    points' coordinates or between arrays of points' coordinates using the
    haversine formula. Expects coordinates in decimal degrees.

    Parameters
    ----------
    lat1
        First point's latitude coordinate(s).
    lon1
        First point's longitude coordinate(s).
    lat2
        Second point's latitude coordinate(s).
    lon2
        Second point's longitude coordinate(s).
    earth_radius
        Earth's radius in units in which distance will be returned (default
        represents meters).

    Returns
    -------
    dist
        Distance from each `(lat1, lon1)` point to each `(lat2, lon2)` point
        in units of `earth_radius`.
    """
    y1 = np.deg2rad(lat1)
    y2 = np.deg2rad(lat2)
    delta_y = y2 - y1

    x1 = np.deg2rad(lon1)
    x2 = np.deg2rad(lon2)
    delta_x = x2 - x1

    h = np.sin(delta_y / 2) ** 2 + np.cos(y1) * np.cos(y2) * np.sin(delta_x / 2) ** 2
    h = np.minimum(1, h)  # protect against floating point errors
    arc = 2 * np.arcsin(np.sqrt(h))

    # return distance in units of earth_radius
    dist: float | npt.NDArray[np.float64] = arc * earth_radius
    return dist


# if coords are all floats, return float
@overload
def euclidean(y1: float, x1: float, y2: float, x2: float) -> float: ...


# if coords are all arrays, return array
@overload
def euclidean(
    y1: npt.NDArray[np.float64],
    x1: npt.NDArray[np.float64],
    y2: npt.NDArray[np.float64],
    x2: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]: ...


def euclidean(
    y1: float | npt.NDArray[np.float64],
    x1: float | npt.NDArray[np.float64],
    y2: float | npt.NDArray[np.float64],
    x2: float | npt.NDArray[np.float64],
) -> float | npt.NDArray[np.float64]:
    """
    Calculate Euclidean distances between pairs of points.

    Vectorized function to calculate the Euclidean distance between two
    points' coordinates or between arrays of points' coordinates. For accurate
    results, use projected coordinates rather than decimal degrees.

    Parameters
    ----------
    y1
        First point's y coordinate(s).
    x1
        First point's x coordinate(s).
    y2
        Second point's y coordinate(s).
    x2
        Second point's x coordinate(s).

    Returns
    -------
    dist
        Distance from each `(x1, y1)` point to each `(x2, y2)` point in same
        units as the points' coordinates.
    """
    # pythagorean theorem
    dist: float | npt.NDArray[np.float64] = ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5
    return dist


def add_edge_lengths(
    G: nx.MultiDiGraph,
    *,
    edges: Iterable[tuple[int, int, int]] | None = None,
) -> nx.MultiDiGraph:
    """
    Calculate and add `length` attribute (in meters) to each edge.

    Vectorized function to calculate great-circle distance between each edge's
    incident nodes. Ensure graph is unprojected and unsimplified to calculate
    accurate distances.

    Note: this function is run by all the `graph.graph_from_x` functions
    automatically to add `length` attributes to all edges. It calculates edge
    lengths as the great-circle distance from node `u` to node `v`. When
    OSMnx automatically runs this function upon graph creation, it does it
    before simplifying the graph: thus it calculates the straight-line lengths
    of edge segments that are themselves all straight. Only after
    simplification do edges take on (potentially) curvilinear geometry. If you
    wish to calculate edge lengths later, note that you will be calculating
    straight-line distances which necessarily ignore the curvilinear geometry.
    Thus you only want to run this function on a graph with all straight edges
    (such as is the case with an unsimplified graph).

    Parameters
    ----------
    G
        Unprojected and unsimplified input graph.
    edges
        The subset of edges to add `length` attributes to, as `(u, v, k)`
        tuples. If None, add lengths to all edges.

    Returns
    -------
    G
        Graph with `length` attributes on the edges.
    """
    uvk = G.edges if edges is None else edges

    # extract edge IDs and corresponding coordinates from their nodes
    x = G.nodes(data="x")
    y = G.nodes(data="y")
    msg = "Some edges missing nodes, possibly due to input data clipping issue."
    try:
        # two-dimensional array of coordinates: y0, x0, y1, x1
        c = np.array([(y[u], x[u], y[v], x[v]) for u, v, k in uvk])
    except KeyError as e:  # pragma: no cover
        raise ValueError(msg) from e
    else:
        # ensure all coordinates can be converted to float and are non-null
        if np.isnan(c.astype(float)).any():
            raise ValueError(msg)

    # calculate great circle distances, round, and fill nulls with zeros
    dists = great_circle(c[:, 0], c[:, 1], c[:, 2], c[:, 3])
    dists[np.isnan(dists)] = 0
    nx.set_edge_attributes(G, values=dict(zip(uvk, dists)), name="length")

    msg = "Added length attributes to graph edges"
    utils.log(msg, level=lg.INFO)
    return G


# if X and Y are floats and return_dist is not provided (defaults False)
@overload
def nearest_nodes(G: nx.MultiDiGraph, X: float, Y: float) -> int: ...


# if X and Y are floats and return_dist is provided/False
@overload
def nearest_nodes(
    G: nx.MultiDiGraph,
    X: float,
    Y: float,
    *,
    return_dist: Literal[False],
) -> int: ...


# if X and Y are floats and return_dist is provided/True
@overload
def nearest_nodes(
    G: nx.MultiDiGraph,
    X: float,
    Y: float,
    *,
    return_dist: Literal[True],
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]: ...


# if X and Y are iterable and return_dist is not provided (defaults False)
@overload
def nearest_nodes(
    G: nx.MultiDiGraph,
    X: Iterable[float],
    Y: Iterable[float],
) -> npt.NDArray[np.int64]: ...


# if X and Y are iterable and return_dist is provided/False
@overload
def nearest_nodes(
    G: nx.MultiDiGraph,
    X: Iterable[float],
    Y: Iterable[float],
    *,
    return_dist: Literal[False],
) -> npt.NDArray[np.int64]: ...


# if X and Y are iterable and return_dist is provided/True
@overload
def nearest_nodes(
    G: nx.MultiDiGraph,
    X: Iterable[float],
    Y: Iterable[float],
    *,
    return_dist: Literal[True],
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]: ...


def nearest_nodes(
    G: nx.MultiDiGraph,
    X: float | Iterable[float],
    Y: float | Iterable[float],
    *,
    return_dist: bool = False,
) -> (
    int
    | npt.NDArray[np.int64]
    | tuple[int, float]
    | tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]
):
    """
    Find the nearest node to a point or to each of several points.

    If `X` and `Y` are single coordinate values, this function will return the
    nearest node to that point. If `X` and `Y` are iterables of coordinate
    values, it will return the nearest node to each point.

    This function is vectorized: if you have many points to search for, pass
    them in one call as numpy arrays (avoid using loops) to maximize runtime
    speed. If the graph is projected, it uses a k-d tree for Euclidean nearest
    neighbor search, which requires that scipy is installed as an optional
    dependency. If the graph is unprojected, it uses a ball tree for haversine
    nearest neighbor search, which requires that scikit-learn is installed as
    an optional dependency.

    Parameters
    ----------
    G
        Graph in which to find nearest nodes.
    X
        The points' x (longitude) coordinates, in same CRS/units as graph and
        containing no nulls.
    Y
        The points' y (latitude) coordinates, in same CRS/units as graph and
        containing no nulls.
    return_dist
        If True, optionally also return the distance(s) between point(s) and
        nearest node(s).

    Returns
    -------
    nn or (nn, dist)
        Nearest node ID(s) or optionally a tuple of ID(s) and distance(s)
        between each point and its nearest node.
    """
    # make coordinates arrays whether user passed iterable values or not
    if not (isinstance(X, Iterable) and isinstance(Y, Iterable)):
        is_scalar = True
        X_arr = np.array([X])
        Y_arr = np.array([Y])
    else:
        is_scalar = False
        X_arr = np.array(X)
        Y_arr = np.array(Y)

    if np.isnan(X_arr).any() or np.isnan(Y_arr).any():  # pragma: no cover
        msg = "`X` and `Y` cannot contain nulls."
        raise ValueError(msg)

    nodes = convert.graph_to_gdfs(G, edges=False, node_geometry=False)[["x", "y"]]
    nn_array: npt.NDArray[np.int64]
    dist_array: npt.NDArray[np.float64]

    if projection.is_projected(G.graph["crs"]):
        # if projected, use k-d tree for euclidean nearest-neighbor search
        if cKDTree is None:  # pragma: no cover
            msg = "scipy must be installed as an optional dependency to search a projected graph."
            raise ImportError(msg)
        dist_array, pos = cKDTree(nodes).query(np.array([X_arr, Y_arr]).T, k=1)
        nn_array = nodes.index[pos].to_numpy()

    else:
        # if unprojected, use ball tree for haversine nearest-neighbor search
        if BallTree is None:  # pragma: no cover
            msg = "scikit-learn must be installed as an optional dependency to search an unprojected graph."
            raise ImportError(msg)
        # haversine requires lat, lon coords in radians
        nodes_rad = np.deg2rad(nodes[["y", "x"]])
        points_rad = np.deg2rad(np.array([Y_arr, X_arr]).T)
        dist_array, pos = BallTree(nodes_rad, metric="haversine").query(points_rad, k=1)
        dist_array = dist_array[:, 0] * EARTH_RADIUS_M  # convert radians -> meters
        nn_array = nodes.index[pos[:, 0]].to_numpy()

    # convert results to correct types for return
    if is_scalar:
        nn = int(nn_array[0])
        dist = float(dist_array[0])
        if return_dist:
            return nn, dist
        # otherwise
        return nn

    # otherwise
    if return_dist:
        return nn_array, dist_array
    # otherwise
    return nn_array


# if X and Y are floats and return_dist is not provided (defaults False)
@overload
def nearest_edges(G: nx.MultiDiGraph, X: float, Y: float) -> tuple[int, int, int]: ...


# if X and Y are floats and return_dist is provided/False
@overload
def nearest_edges(
    G: nx.MultiDiGraph,
    X: float,
    Y: float,
    *,
    return_dist: Literal[False],
) -> tuple[int, int, int]: ...


# if X and Y are floats and return_dist is provided/True
@overload
def nearest_edges(
    G: nx.MultiDiGraph,
    X: float,
    Y: float,
    *,
    return_dist: Literal[True],
) -> tuple[tuple[int, int, int], float]: ...


# if X and Y are iterable and return_dist is not provided (defaults False)
@overload
def nearest_edges(
    G: nx.MultiDiGraph,
    X: Iterable[float],
    Y: Iterable[float],
) -> npt.NDArray[np.object_]: ...


# if X and Y are iterable and return_dist is provided/False
@overload
def nearest_edges(
    G: nx.MultiDiGraph,
    X: Iterable[float],
    Y: Iterable[float],
    *,
    return_dist: Literal[False],
) -> npt.NDArray[np.object_]: ...


# if X and Y are iterable and return_dist is provided/True
@overload
def nearest_edges(
    G: nx.MultiDiGraph,
    X: Iterable[float],
    Y: Iterable[float],
    *,
    return_dist: Literal[True],
) -> tuple[npt.NDArray[np.object_], npt.NDArray[np.float64]]: ...


def nearest_edges(
    G: nx.MultiDiGraph,
    X: float | Iterable[float],
    Y: float | Iterable[float],
    *,
    return_dist: bool = False,
) -> (
    tuple[int, int, int]
    | npt.NDArray[np.object_]
    | tuple[tuple[int, int, int], float]
    | tuple[npt.NDArray[np.object_], npt.NDArray[np.float64]]
):
    """
    Find the nearest edge to a point or to each of several points.

    If `X` and `Y` are single coordinate values, this function will return the
    nearest edge to that point. If `X` and `Y` are iterables of coordinate
    values, it will return the nearest edge to each point.

    This function is vectorized: if you have many points to search for, pass
    them in one call as numpy arrays (avoid using loops) to maximize runtime
    speed. It uses an R-tree spatial index and minimizes the Euclidean
    distance from each point to the possible matches. For accurate results,
    use a projected graph and projected points.

    Parameters
    ----------
    G
        Graph in which to find nearest edges.
    X
        The points' x (longitude) coordinates, in same CRS/units as graph and
        containing no nulls.
    Y
        The points' y (latitude) coordinates, in same CRS/units as graph and
        containing no nulls.
    return_dist
        If True, optionally also return the distance(s) between point(s) and
        nearest edge(s), in same units as graph and points.

    Returns
    -------
    ne or (ne, dist)
        Nearest edge ID(s) as `(u, v, k)` tuples, or optionally a tuple of
        ID(s) and distance(s) between each point and its nearest edge.
    """
    # make coordinates arrays whether user passed iterable values or not
    if not (isinstance(X, Iterable) and isinstance(Y, Iterable)):
        is_scalar = True
        X_arr = np.array([X])
        Y_arr = np.array([Y])
    else:
        is_scalar = False
        X_arr = np.array(X)
        Y_arr = np.array(Y)

    if np.isnan(X_arr).any() or np.isnan(Y_arr).any():  # pragma: no cover
        msg = "`X` and `Y` cannot contain nulls."
        raise ValueError(msg)
    geoms = convert.graph_to_gdfs(G, nodes=False)["geometry"]
    ne_array: npt.NDArray[np.object_]  # array of tuple[int, int, int]
    dist_array: npt.NDArray[np.float64]

    # build an r-tree spatial index by position for subsequent iloc
    rtree = STRtree(geoms)

    # use the r-tree to find each point's nearest neighbor and distance
    points = [Point(xy) for xy in zip(X_arr, Y_arr)]
    pos, dist_array = rtree.query_nearest(points, all_matches=False, return_distance=True)

    # if user passed X/Y lists, the 2nd subarray contains geom indices
    if len(pos.shape) > 1:
        pos = pos[1]
    ne_array = geoms.iloc[pos].index.to_numpy()

    # convert results to correct types for return
    if is_scalar:
        ne: tuple[int, int, int] = ne_array[0]
        dist = float(dist_array[0])
        if return_dist:
            return ne, dist
        # otherwise
        return ne

    # otherwise
    if return_dist:
        return ne_array, dist_array
    # otherwise
    return ne_array
