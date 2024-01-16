"""Calculate graph edge bearings."""
from __future__ import annotations

from typing import overload
from warnings import warn

import networkx as nx
import numpy as np
from numpy.typing import NDArray

from . import plot
from . import projection

# scipy is an optional dependency for entropy calculation
try:
    import scipy
except ImportError:  # pragma: no cover
    scipy = None


# if coords are all floats, return float
@overload  # pragma: no cover
def calculate_bearing(
    lat1: float,
    lon1: float,
    lat2: float,
    lon2: float,
) -> float:
    ...


# if coords are all arrays, return array
@overload  # pragma: no cover
def calculate_bearing(
    lat1: NDArray[np.float64],
    lon1: NDArray[np.float64],
    lat2: NDArray[np.float64],
    lon2: NDArray[np.float64],
) -> NDArray[np.float64]:
    ...


def calculate_bearing(
    lat1: float | NDArray[np.float64],
    lon1: float | NDArray[np.float64],
    lat2: float | NDArray[np.float64],
    lon2: float | NDArray[np.float64],
) -> float | NDArray[np.float64]:
    """
    Calculate the compass bearing(s) between pairs of lat-lon points.

    Vectorized function to calculate initial bearings between two points'
    coordinates or between arrays of points' coordinates. Expects coordinates
    in decimal degrees. Bearing represents the clockwise angle in degrees
    between north and the geodesic line from (lat1, lon1) to (lat2, lon2).

    Parameters
    ----------
    lat1 : float or numpy.array of float
        first point's latitude coordinate
    lon1 : float or numpy.array of float
        first point's longitude coordinate
    lat2 : float or numpy.array of float
        second point's latitude coordinate
    lon2 : float or numpy.array of float
        second point's longitude coordinate

    Returns
    -------
    bearing : float or numpy.array of float
        the bearing(s) in decimal degrees
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
    bearing: float | NDArray[np.float64] = initial_bearing % 360
    return bearing


def add_edge_bearings(G: nx.MultiDiGraph, precision: int | None = None) -> nx.MultiDiGraph:
    """
    Add compass `bearing` attributes to all graph edges.

    Vectorized function to calculate (initial) bearing from origin node to
    destination node for each edge in a directed, unprojected graph then add
    these bearings as new edge attributes. Bearing represents angle in degrees
    (clockwise) between north and the geodesic line from the origin node to
    the destination node. Ignores self-loop edges as their bearings are
    undefined.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        unprojected graph
    precision : int
        deprecated, do not use

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with edge bearing attributes
    """
    if precision is None:
        precision = 1
    else:
        warn(
            "the `precision` parameter is deprecated and will be removed in a future release",
            stacklevel=2,
        )

    if projection.is_projected(G.graph["crs"]):  # pragma: no cover
        msg = "graph must be unprojected to add edge bearings"
        raise ValueError(msg)

    # extract edge IDs and corresponding coordinates from their nodes
    uvk = [(u, v, k) for u, v, k in G.edges if u != v]
    x = G.nodes(data="x")
    y = G.nodes(data="y")
    coords = np.array([(y[u], x[u], y[v], x[v]) for u, v, k in uvk])

    # calculate bearings then set as edge attributes
    bearings = calculate_bearing(coords[:, 0], coords[:, 1], coords[:, 2], coords[:, 3])
    values = zip(uvk, bearings.round(precision))
    nx.set_edge_attributes(G, dict(values), name="bearing")

    return G


def orientation_entropy(
    Gu: nx.MultiGraph, num_bins: int = 36, min_length: float = 0, weight: str | None = None
) -> float:
    """
    Calculate undirected graph's orientation entropy.

    Orientation entropy is the entropy of its edges' bidirectional bearings
    across evenly spaced bins. Ignores self-loop edges as their bearings are
    undefined.

    For more info see: Boeing, G. 2019. "Urban Spatial Order: Street Network
    Orientation, Configuration, and Entropy." Applied Network Science, 4 (1),
    67. https://doi.org/10.1007/s41109-019-0189-1

    Parameters
    ----------
    Gu : networkx.MultiGraph
        undirected, unprojected graph with `bearing` attributes on each edge
    num_bins : int
        number of bins; for example, if `num_bins=36` is provided, then each
        bin will represent 10 degrees around the compass
    min_length : float
        ignore edges with `length` attributes less than `min_length`; useful
        to ignore the noise of many very short edges
    weight : string
        if not None, weight edges' bearings by this (non-null) edge attribute.
        for example, if "length" is provided, this will return 1 bearing
        observation per meter per street, which could result in a very large
        `bearings` array.

    Returns
    -------
    entropy : float
        the graph's orientation entropy
    """
    # check if we were able to import scipy
    if scipy is None:  # pragma: no cover
        msg = "scipy must be installed to calculate entropy"
        raise ImportError(msg)
    bin_counts, _ = _bearings_distribution(Gu, num_bins, min_length, weight)
    entropy: float = scipy.stats.entropy(bin_counts)
    return entropy


def _extract_edge_bearings(
    Gu: nx.MultiGraph, min_length: float = 0, weight: str | None = None
) -> NDArray[np.float64]:
    """
    Extract undirected graph's bidirectional edge bearings.

    For example, if an edge has a bearing of 90 degrees then we will record
    bearings of both 90 degrees and 270 degrees for this edge.

    Parameters
    ----------
    Gu : networkx.MultiGraph
        undirected, unprojected graph with `bearing` attributes on each edge
    min_length : float
        ignore edges with `length` attributes less than `min_length`; useful
        to ignore the noise of many very short edges
    weight : string
        if not None, weight edges' bearings by this (non-null) edge attribute.
        for example, if "length" is provided, this will return 1 bearing
        observation per meter per street, which could result in a very large
        `bearings` array.

    Returns
    -------
    bearings : numpy.array of float
        the graph's bidirectional edge bearings
    """
    if nx.is_directed(Gu) or projection.is_projected(Gu.graph["crs"]):  # pragma: no cover
        msg = "graph must be undirected and unprojected to analyze edge bearings"
        raise ValueError(msg)
    bearings = []
    for u, v, data in Gu.edges(data=True):
        # ignore self-loops and any edges below min_length
        if u != v and data["length"] >= min_length:
            if weight:
                # weight edges' bearings by some edge attribute value
                bearings.extend([data["bearing"]] * int(data[weight]))
            else:
                # don't weight bearings, just take one value per edge
                bearings.append(data["bearing"])

    # drop any nulls, calculate reverse bearings, concatenate and return
    bearings_array = np.array(bearings)
    bearings_array = bearings_array[~np.isnan(bearings_array)]
    bearings_array_r = (bearings_array - 180) % 360
    return np.concatenate([bearings_array, bearings_array_r])


def _bearings_distribution(
    Gu: nx.MultiGraph, num_bins: int, min_length: float = 0, weight: str | None = None
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Compute distribution of bearings across evenly spaced bins.

    Prevents bin-edge effects around common values like 0 degrees and 90
    degrees by initially creating twice as many bins as desired, then merging
    them in pairs. For example, if `num_bins=36` is provided, then each bin
    will represent 10 degrees around the compass, with the first bin
    representing 355 degrees to 5 degrees.

    Parameters
    ----------
    Gu : networkx.MultiGraph
        undirected, unprojected graph with `bearing` attributes on each edge
    num_bins : int
        number of bins for the bearings histogram
    min_length : float
        ignore edges with `length` attributes less than `min_length`; useful
        to ignore the noise of many very short edges
    weight : string
        if not None, weight edges' bearings by this (non-null) edge attribute.
        for example, if "length" is provided, this will return 1 bearing
        observation per meter per street, which could result in a very large
        `bearings` array.

    Returns
    -------
    bin_counts, bin_edges : tuple of numpy.array of float
        counts of bearings per bin and the bins edges
    """
    n = num_bins * 2
    bins = np.arange(n + 1) * 360 / n

    bearings = _extract_edge_bearings(Gu, min_length, weight)
    count, bin_edges = np.histogram(bearings, bins=bins)

    # move last bin to front, so eg 0.01 degrees and 359.99 degrees will be
    # binned together
    count = np.roll(count, 1)
    bin_counts = count[::2] + count[1::2]

    # because we merged the bins, their edges are now only every other one
    bin_edges = bin_edges[range(0, len(bin_edges), 2)]
    return bin_counts, bin_edges


def plot_orientation(  # type: ignore[no-untyped-def]
    Gu,
    num_bins=36,
    min_length=0,
    weight=None,
    ax=None,
    figsize=(5, 5),
    area=True,
    color="#003366",
    edgecolor="k",
    linewidth=0.5,
    alpha=0.7,
    title=None,
    title_y=1.05,
    title_font=None,
    xtick_font=None,
):
    """
    Do not use: deprecated.

    The plot_orientation function moved to the plot module. Calling it via the
    bearing module will raise an error in a future release.

    Parameters
    ----------
    Gu : networkx.MultiGraph
        deprecated, do not use
    num_bins : int
        deprecated, do not use
    min_length : float
        deprecated, do not use
    weight : string
        deprecated, do not use
    ax : matplotlib.axes.PolarAxesSubplot
        deprecated, do not use
    figsize : tuple
        deprecated, do not use
    area : bool
        deprecated, do not use
    color : string
        deprecated, do not use
    edgecolor : string
        deprecated, do not use
    linewidth : float
        deprecated, do not use
    alpha : float
        deprecated, do not use
    title : string
        deprecated, do not use
    title_y : float
        deprecated, do not use
    title_font : dict
        deprecated, do not use
    xtick_font : dict
        deprecated, do not use

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axis
    """
    warn(
        "The `plot_orientation` function moved to the `plot` module. Calling it "
        "via the `bearing` module will cause an exception in a future release.",
        stacklevel=2,
    )
    return plot.plot_orientation(
        Gu,
        num_bins=num_bins,
        min_length=min_length,
        weight=weight,
        ax=ax,
        figsize=figsize,
        area=area,
        color=color,
        edgecolor=edgecolor,
        linewidth=linewidth,
        alpha=alpha,
        title=title,
        title_y=title_y,
        title_font=title_font,
        xtick_font=xtick_font,
    )
