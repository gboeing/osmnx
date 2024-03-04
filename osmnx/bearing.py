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
    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)
    delta_lon = np.radians(lon2 - lon1)

    # calculate initial bearing from -180 degrees to +180 degrees
    y = np.sin(delta_lon) * np.cos(lat2)
    x = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(delta_lon)
    initial_bearing = np.degrees(np.arctan2(y, x))

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
    G: nx.MultiGraph,
    *,
    num_bins: int = 36,
    min_length: float = 0,
    weight: str | None = None,
) -> float:
    """
    Calculate graph's orientation entropy.

    Orientation entropy is the Shannon entropy of the graphs' edges'
    bearings across evenly spaced bins. Ignores self-loop edges
    as their bearings are undefined.

    For MultiGraph input, calculates entropy of bidirectional bearings.
    For MultiDiGraph input, calculates entropy of directional bearings.

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
        If not None, weight edges' bearings by this (non-null) edge attribute.
        For example, if "length" is provided, this will return 1 bearing
        observation per meter per street.

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
    G: nx.MultiGraph,
    min_length: float,
    weight: str | None,
) -> npt.NDArray[np.float64]:
    """
    Extract graph's edge bearings.

    A MultiGraph input receives bidirectional bearings.
    For example, if an undirected edge has a bearing of 90 degrees then we will record
    bearings of both 90 degrees and 270 degrees for this edge.
    For MultiDiGraph input, record only one bearing per edge.

    Parameters
    ----------
    G
        Unprojected graph with `bearing` attributes on each edge.
    min_length
        Ignore edges with `length` attributes less than `min_length`. Useful
        to ignore the noise of many very short edges.
    weight
        If not None, weight edges' bearings by this (non-null) edge attribute.
        For example, if "length" is provided, this will return 1 bearing
        observation per meter per street (which could result in a very large
        `bearings` array).

    Returns
    -------
    bearings
        The edge bearings of `Gu`.
    """
    if projection.is_projected(G.graph["crs"]):  # pragma: no cover
        msg = "Graph must be unprojected to analyze edge bearings."
        raise ValueError(msg)
    bearings = []
    for u, v, data in G.edges(data=True):
        # ignore self-loops and any edges below min_length
        if u != v and data["length"] >= min_length:
            if weight:
                # weight edges' bearings by some edge attribute value
                bearings.extend([data["bearing"]] * int(data[weight]))
            else:
                # don't weight bearings, just take one value per edge
                bearings.append(data["bearing"])

    # drop any nulls
    bearings_array = np.array(bearings)
    bearings_array = bearings_array[~np.isnan(bearings_array)]
    if nx.is_directed(G):
        # https://github.com/gboeing/osmnx/issues/1137
        msg = (
            "Extracting directional bearings (one bearing per edge) due to MultiDiGraph input. "
            "To extract bidirectional bearings (two bearings per edge, including the reverse bearing), "
            "supply an undirected graph instead via `osmnx.get_undirected(G)`."
        )
        warn(msg, category=UserWarning, stacklevel=2)
        return bearings_array
    # for undirected graphs, add reverse bearings and return
    bearings_array_r = (bearings_array - 180) % 360
    return np.concatenate([bearings_array, bearings_array_r])


def _bearings_distribution(
    G: nx.MultiGraph,
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
        If not None, weight edges' bearings by this (non-null) edge attribute.
        For example, if "length" is provided, this will return 1 bearing
        observation per meter per street (which could result in a very large
        `bearings` array).

    Returns
    -------
    bin_counts, bin_edges
        Counts of bearings per bin and the bins edges.
    """
    n = num_bins * 2
    bins = np.arange(n + 1) * 360 / n

    bearings = _extract_edge_bearings(G, min_length, weight)
    count, bin_edges = np.histogram(bearings, bins=bins)

    # move last bin to front, so eg 0.01 degrees and 359.99 degrees will be
    # binned together
    count = np.roll(count, 1)
    bin_counts = count[::2] + count[1::2]

    # because we merged the bins, their edges are now only every other one
    bin_edges = bin_edges[range(0, len(bin_edges), 2)]
    return bin_counts, bin_edges
