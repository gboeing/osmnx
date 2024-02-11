"""Calculate weighted shortest paths between graph nodes."""

from __future__ import annotations

import itertools
import logging as lg
import multiprocessing as mp
from collections.abc import Iterable
from collections.abc import Iterator
from typing import overload
from warnings import warn

import networkx as nx
import numpy as np

from . import utils
from . import utils_graph


# orig/dest int, weight present, cpus present
@overload  # pragma: no cover
def shortest_path(
    G: nx.MultiDiGraph,
    orig: int,
    dest: int,
    weight: str,
    cpus: int | None,
) -> list[int] | None:
    ...


# orig/dest int, weight missing, cpus present
@overload  # pragma: no cover
def shortest_path(
    G: nx.MultiDiGraph,
    orig: int,
    dest: int,
    *,
    cpus: int | None,
) -> list[int] | None:
    ...


# orig/dest int, weight present, cpus missing
@overload  # pragma: no cover
def shortest_path(
    G: nx.MultiDiGraph,
    orig: int,
    dest: int,
    weight: str,
) -> list[int] | None:
    ...


# orig/dest int, weight missing, cpus missing
@overload  # pragma: no cover
def shortest_path(
    G: nx.MultiDiGraph,
    orig: int,
    dest: int,
) -> list[int] | None:
    ...


# orig/dest Iterable, weight present, cpus present
@overload  # pragma: no cover
def shortest_path(
    G: nx.MultiDiGraph,
    orig: Iterable[int],
    dest: Iterable[int],
    weight: str,
    cpus: int | None,
) -> list[list[int] | None]:
    ...


# orig/dest Iterable, weight missing, cpus present
@overload  # pragma: no cover
def shortest_path(
    G: nx.MultiDiGraph,
    orig: Iterable[int],
    dest: Iterable[int],
    *,
    cpus: int | None,
) -> list[list[int] | None]:
    ...


# orig/dest Iterable, weight present, cpus missing
@overload  # pragma: no cover
def shortest_path(
    G: nx.MultiDiGraph,
    orig: Iterable[int],
    dest: Iterable[int],
    weight: str,
) -> list[list[int] | None]:
    ...


# orig/dest Iterable, weight missing, cpus missing
@overload  # pragma: no cover
def shortest_path(
    G: nx.MultiDiGraph,
    orig: Iterable[int],
    dest: Iterable[int],
) -> list[list[int] | None]:
    ...


def shortest_path(
    G: nx.MultiDiGraph,
    orig: int | Iterable[int],
    dest: int | Iterable[int],
    weight: str = "length",
    cpus: int | None = 1,
) -> list[int] | None | list[list[int] | None]:
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
    G
        Input graph,
    orig
        Origin node ID(s).
    dest
        Destination node ID(s).
    weight
        Edge attribute to minimize when solving shortest path.
    cpus
        How many CPU cores to use. If None, use all available.

    Returns
    -------
    path
        The node IDs constituting the shortest path, or, if `orig` and `dest`
        are both iterable, then a list of such paths.
    """
    _verify_edge_attribute(G, weight)

    # if neither orig nor dest is iterable, just return the shortest path
    if not (isinstance(orig, Iterable) or isinstance(dest, Iterable)):
        return _single_shortest_path(G, orig, dest, weight)

    # if only 1 of orig or dest is iterable and the other is not, raise error
    if not (isinstance(orig, Iterable) and isinstance(dest, Iterable)):
        msg = "orig and dest must either both be iterable or neither must be iterable"
        raise TypeError(msg)

    # if both orig and dest are iterable, make them lists (so we're guaranteed
    # to be able to get their sizes) then ensure they have same lengths
    orig = list(orig)
    dest = list(dest)
    if len(orig) != len(dest):  # pragma: no cover
        msg = "orig and dest must be of equal length"
        raise ValueError(msg)

    # determine how many cpu cores to use
    if cpus is None:
        cpus = mp.cpu_count()
    cpus = min(cpus, mp.cpu_count())

    msg = f"Solving {len(orig)} paths with {cpus} CPUs..."
    utils.log(msg, level=lg.INFO)

    # if single-threading, calculate each shortest path one at a time
    if cpus == 1:
        paths = [_single_shortest_path(G, o, d, weight) for o, d in zip(orig, dest)]

    # if multi-threading, calculate shortest paths in parallel
    else:
        args = ((G, o, d, weight) for o, d in zip(orig, dest))
        with mp.get_context("spawn").Pool(cpus) as pool:
            paths = pool.starmap_async(_single_shortest_path, args).get()

    return paths


def k_shortest_paths(
    G: nx.MultiDiGraph,
    orig: int,
    dest: int,
    k: int,
    weight: str = "length",
) -> Iterator[list[int]]:
    """
    Solve `k` shortest paths from an origin node to a destination node.

    Uses Yen's algorithm. See also `shortest_path` to solve just the one
    shortest path.

    Parameters
    ----------
    G
        Input graph.
    orig
        Origin node ID.
    dest
        Destination node ID.
    k
        Number of shortest paths to solve.
    weight
        Edge attribute to minimize when solving shortest paths.

    Yields
    ------
    path
        The node IDs constituting the next-shortest path.
    """
    _verify_edge_attribute(G, weight)
    paths_gen = nx.shortest_simple_paths(utils_graph.get_digraph(G, weight), orig, dest, weight)
    yield from itertools.islice(paths_gen, 0, k)


def _single_shortest_path(
    G: nx.MultiDiGraph,
    orig: int,
    dest: int,
    weight: str,
) -> list[int] | None:
    """
    Solve the shortest path from an origin node to a destination node.

    This function uses Dijkstra's algorithm. It is a convenience wrapper
    around `networkx.shortest_path`, with exception handling for unsolvable
    paths. If the path is unsolvable, it returns None.

    Parameters
    ----------
    G
        Input graph.
    orig
        Origin node ID.
    dest
        Destination node ID.
    weight
        Edge attribute to minimize when solving shortest path.

    Returns
    -------
    path
        The node IDs constituting the shortest path.
    """
    try:
        return list(nx.shortest_path(G, orig, dest, weight=weight, method="dijkstra"))
    except nx.exception.NetworkXNoPath:  # pragma: no cover
        msg = f"Cannot solve path from {orig} to {dest}"
        utils.log(msg, level=lg.WARNING)
        return None


def _verify_edge_attribute(G: nx.MultiDiGraph, attr: str) -> None:
    """
    Verify attribute values are numeric and non-null across graph edges.

    Raises a ValueError if this attribute contains non-numeric values, and
    issues a UserWarning if this attribute is missing or null on any edges.

    Parameters
    ----------
    G
        Input graph.
    attr
        Name of the edge attribute to verify.

    Returns
    -------
    None
    """
    try:
        values_float = (np.array(tuple(G.edges(data=attr)))[:, 2]).astype(float)
        if np.isnan(values_float).any():
            msg = f"The attribute {attr!r} is missing or null on some edges."
            warn(msg, category=UserWarning, stacklevel=2)
    except ValueError as e:
        msg = f"The edge attribute {attr!r} contains non-numeric values."
        raise ValueError(msg) from e
