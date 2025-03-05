"""Calculate graph edge bearings and orientation entropy."""

from __future__ import annotations

from typing import overload
from warnings import warn

import networkx as nx
import numpy as np
import numpy.typing as npt

from . import projection

# scipy is an optional dependency for entropy calculation
try:
    import scipy
except ImportError:  # pragma: no cover
    scipy = None


# if coords are all floats, return float
@overload
def calculate_bearing(
    lat1: float,
    lon1: float,
    lat2: float,
    lon2: float,
) -> float: ...


# if coords are all arrays, return array
@overload
def calculate_bearing(
    lat1: npt.NDArray[np.float64],
    lon1: npt.NDArray[np.float64],
    lat2: npt.NDArray[np.float64],
    lon2: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]: ...


def calculate_bearing(
    lat1: float | npt.NDArray[np.float64],
    lon1: float | npt.NDArray[np.float64],
    lat2: float | npt.NDArray[np.float64],
    lon2: float | npt.NDArray[np.float64],
) -> float | npt.NDArray[np.float64]:
    """
    Calculate the compass bearing(s) between pairs of lat-lon points.

    Vectorized function to calculate initial bearings between two points'
    coordinates or between arrays of points' coordinates. Expects coordinates
    in decimal degrees. The bearing represents the clockwise angle in degrees
    between north and the geodesic line from `(lat1, lon1)` to `(lat2, lon2)`.

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

    Returns
    -------
    bearing
        The bearing(s) in decimal degrees.
    """
    # get the latitudes and the difference in longitudes, all in radians
    lat1 = np.deg2rad(lat1)
    lat2 = np.deg2rad(lat2)
    delta_lon = np.deg2rad(lon2 - lon1)

    # calculate initial bearing from -180 degrees to +180 degrees
    y = np.sin(delta_lon) * np.cos(lat2)
    x = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(delta_lon)
    initial_bearing = np.rad2deg(np.arctan2(y, x))

    # normalize to 0-360 degrees to get compass bearing
    bearing: float | npt.NDArray[np.float64] = initial_bearing % 360
    return bearing


def add_edge_bearings(G: nx.MultiDiGraph) -> nx.MultiDiGraph:
    """
    Calculate and add compass `bearing` attributes to all graph edges.

    Vectorized function to calculate (initial) bearing from origin node to
    destination node for each edge in a directed, unprojected graph then add
    these bearings as new `bearing` edge attributes. Bearing represents angle
    in degrees (clockwise) between north and the geodesic line from the origin
    node to the destination node. Ignores self-loop edges as their bearings
    are undefined.

    Parameters
    ----------
    G
        Unprojected graph.

    Returns
    -------
    G
        Graph with `bearing` attributes on the edges.
    """
    if projection.is_projected(G.graph["crs"]):  # pragma: no cover
        msg = "Graph must be unprojected to add edge bearings."
        raise ValueError(msg)

    # extract edge IDs and corresponding coordinates from their nodes
    uvk = [(u, v, k) for u, v, k in G.edges if u != v]
    x = G.nodes(data="x")
    y = G.nodes(data="y")
    coords = np.array([(y[u], x[u], y[v], x[v]) for u, v, k in uvk])

    # calculate bearings then set as edge attributes
    bearings = calculate_bearing(coords[:, 0], coords[:, 1], coords[:, 2], coords[:, 3])
    values = zip(uvk, bearings)
    nx.set_edge_attributes(G, dict(values), name="bearing")

    return G


def orientation_entropy(
    G: nx.MultiGraph | nx.MultiDiGraph,
    *,
    num_bins: int = 36,
    min_length: float = 0,
    weight: str | None = None,
) -> float:
    """
    Calculate graph's orientation entropy.

    Orientation entropy is the Shannon entropy of the graphs' edges' bearings
    across evenly spaced bins. Ignores self-loop edges as their bearings are
    undefined. If `G` is a MultiGraph, all edge bearings will be bidirectional
    (ie, two reciprocal bearings per undirected edge). If `G` is a
    MultiDiGraph, all edge bearings will be directional (ie, one bearing per
    directed edge).

    For more info see: Boeing, G. 2019. "Urban Spatial Order: Street Network
    Orientation, Configuration, and Entropy." Applied Network Science, 4 (1),
    67. https://doi.org/10.1007/s41109-019-0189-1

    Parameters
    ----------
    G
        Unprojected graph with `bearing` attributes on each edge.
    num_bins
        Number of bins. For example, if `num_bins=36` is provided, then each
        bin will represent 10 degrees around the compass.
    min_length
        Ignore edges with "length" attributes less than `min_length`. Useful
        to ignore the noise of many very short edges.
    weight
        If None, apply equal weight for each bearing. Otherwise, weight edges'
        bearings by this (non-null) edge attribute. For example, if "length"
        is provided, each edge's bearing observation will be weighted by its
        "length" attribute value.

    Returns
    -------
    entropy
        The orientation entropy of `G`.
    """
    # check if we were able to import scipy
    if scipy is None:  # pragma: no cover
        msg = "scipy must be installed as an optional dependency to calculate entropy."
        raise ImportError(msg)
    bin_counts, _ = _bearings_distribution(G, num_bins, min_length, weight)
    entropy: float = scipy.stats.entropy(bin_counts)
    return entropy


def _extract_edge_bearings(
    G: nx.MultiGraph | nx.MultiDiGraph,
    min_length: float,
    weight: str | None,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """
    Extract graph's edge bearings.

    Ignores self-loop edges as their bearings are undefined. If `G` is a
    MultiGraph, all edge bearings will be bidirectional (ie, two reciprocal
    bearings per undirected edge). If `G` is a MultiDiGraph, all edge bearings
    will be directional (ie, one bearing per directed edge). For example, if
    an undirected edge has a bearing of 90 degrees then we will record
    bearings of both 90 degrees and 270 degrees for this edge.

    Parameters
    ----------
    G
        Unprojected graph with `bearing` attributes on each edge.
    min_length
        Ignore edges with `length` attributes less than `min_length`. Useful
        to ignore the noise of many very short edges.
    weight
        If None, apply equal weight for each bearing. Otherwise, weight edges'
        bearings by this (non-null) edge attribute. For example, if "length"
        is provided, each edge's bearing observation will be weighted by its
        "length" attribute value.

    Returns
    -------
    bearings, weights
        The edge bearings of `G` and their corresponding weights.
    """
    if projection.is_projected(G.graph["crs"]):  # pragma: no cover
        msg = "Graph must be unprojected to analyze edge bearings."
        raise ValueError(msg)
    bearings = []
    weights = []
    for u, v, data in G.edges(data=True):
        # ignore self-loops and any edges below min_length
        if u != v and data["length"] >= min_length:
            bearings.append(data["bearing"])
            weights.append(data[weight] if weight is not None else 1.0)

    # drop any nulls
    bearings_array = np.array(bearings)
    weights_array = np.array(weights)
    keep_idx = ~np.isnan(bearings_array)
    bearings_array = bearings_array[keep_idx]
    weights_array = weights_array[keep_idx]
    if nx.is_directed(G):
        msg = (
            "`G` is a MultiDiGraph, so edge bearings will be directional (one per "
            "edge). If you want bidirectional edge bearings (two reciprocal bearings "
            "per edge), pass a MultiGraph instead. Use `convert.to_undirected`."
        )
        warn(msg, category=UserWarning, stacklevel=2)
        return bearings_array, weights_array
    # for undirected graphs, add reverse bearings
    bearings_array = np.concatenate([bearings_array, (bearings_array - 180) % 360])
    weights_array = np.concatenate([weights_array, weights_array])
    return bearings_array, weights_array


def _bearings_distribution(
    G: nx.MultiGraph | nx.MultiDiGraph,
    num_bins: int,
    min_length: float,
    weight: str | None,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """
    Compute distribution of bearings across evenly spaced bins.

    Prevents bin-edge effects around common values like 0 degrees and 90
    degrees by initially creating twice as many bins as desired, then merging
    them in pairs. For example, if `num_bins=36` is provided, then each bin
    will represent 10 degrees around the compass, with the first bin
    representing 355 degrees to 5 degrees.

    Parameters
    ----------
    G
        Unprojected graph with `bearing` attributes on each edge.
    num_bins
        Number of bins for the bearing histogram.
    min_length
        Ignore edges with `length` attributes less than `min_length`. Useful
        to ignore the noise of many very short edges.
    weight
        If None, apply equal weight for each bearing. Otherwise, weight edges'
        bearings by this (non-null) edge attribute. For example, if "length"
        is provided, each edge's bearing observation will be weighted by its
        "length" attribute value.

    Returns
    -------
    bin_counts, bin_centers
        Counts of bearings per bin and the bins' centers in degrees. Both
        arrays are of length `num_bins`.
    """
    # Split bins in half to prevent bin-edge effects around common values.
    # Bins will be merged in pairs after the histogram is computed. The last
    # bin edge is the same as the first (i.e., 0 degrees = 360 degrees).
    num_split_bins = num_bins * 2
    split_bin_edges = np.linspace(0, 360, num_split_bins + 1)

    bearings, weights = _extract_edge_bearings(G, min_length, weight)
    split_bin_counts, split_bin_edges = np.histogram(
        bearings,
        bins=split_bin_edges,
        weights=weights,
    )

    # Move last bin to front, so eg 0.01 degrees and 359.99 degrees will be
    # binned together. Then combine counts from pairs of split bins.
    split_bin_counts = np.roll(split_bin_counts, 1)
    bin_counts = split_bin_counts[::2] + split_bin_counts[1::2]

    # Every other edge of the split bins is the center of a merged bin.
    bin_centers = split_bin_edges[range(0, num_split_bins - 1, 2)]
    return bin_counts, bin_centers
