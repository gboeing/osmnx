"""Calculate weighted shortest paths between graph nodes."""

import itertools
import multiprocessing as mp
from warnings import warn

import networkx as nx
import numpy as np

from . import utils
from . import utils_graph


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

    # if only 1 of orig or dest is iterable and the other is not, raise error
    if not (hasattr(orig, "__iter__") and hasattr(dest, "__iter__")):
        msg = "orig and dest must either both be iterable or neither must be iterable"
        raise ValueError(msg)

    # if both orig and dest are iterable, ensure they have same lengths
    if len(orig) != len(dest):  # pragma: no cover
        msg = "orig and dest must be of equal length"
        raise ValueError(msg)

    # determine how many cpu cores to use
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
        with mp.get_context("spawn").Pool(cpus) as pool:
            paths = pool.starmap_async(_single_shortest_path, args).get()

    return paths


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
