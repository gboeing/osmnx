"""Calculate graph edge bearings."""

import warnings

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from . import projection

# scipy is an optional dependency for entropy calculation
try:
    import scipy
except ImportError:  # pragma: no cover
    scipy = None


def calculate_bearing(lat1, lng1, lat2, lng2):
    """
    Calculate the compass bearing(s) between pairs of lat-lng points.

    Vectorized function to calculate (initial) bearings between two points'
    coordinates or between arrays of points' coordinates. Expects coordinates
    in decimal degrees. Bearing represents angle in degrees (clockwise)
    between north and the geodesic line from point 1 to point 2.

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

    Returns
    -------
    bearing : float or numpy.array of float
        the bearing(s) in decimal degrees
    """
    # get the latitudes and the difference in longitudes, in radians
    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)
    d_lng = np.radians(lng2 - lng1)

    # calculate initial bearing from -180 degrees to +180 degrees
    y = np.sin(d_lng) * np.cos(lat2)
    x = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(d_lng)
    initial_bearing = np.degrees(np.arctan2(y, x))

    # normalize to 0-360 degrees to get compass bearing
    return initial_bearing % 360


def get_bearing(origin_point, destination_point):
    """
    Do not use, deprecated.

    Parameters
    ----------
    origin_point : tuple
        deprecated, do not use
    destination_point : tuple
        deprecated, do not use

    Returns
    -------
    bearing : float
        deprecated, do not use
    """
    msg = (
        "The `get_bearing` function has been deprecated and will be removed in a "
        "future release. Use the vectorized `bearing.calculate_bearing` instead."
    )
    warnings.warn(msg)
    return calculate_bearing(
        origin_point[0], origin_point[1], destination_point[0], destination_point[1]
    )


def add_edge_bearings(G, precision=1):
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
        decimal precision to round bearing

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with edge bearing attributes
    """
    if projection.is_projected(G.graph["crs"]):  # pragma: no cover
        raise ValueError("graph must be unprojected to add edge bearings")

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


def orientation_entropy(Gu, num_bins=36, min_length=0, weight=None):
    """
    Calculate undirected graph's orientation entropy.

    Orientation entropy is the entropy of its edges' bidirectional bearings
    across evenly spaced bins. Ignores self-loop edges as their bearings are
    undefined.

    Parameters
    ----------
    Gu : networkx.MultiGraph
        undirected, unprojected graph with `bearing` attributes on each edge
    num_bins : int
        number of bins; for example, if `num_bins=36` is provided, then each
        bin will represent 10° around the compass
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
        raise ImportError("scipy must be installed to calculate entropy")
    bin_counts, _ = _bearings_distribution(Gu, num_bins, min_length, weight)
    return scipy.stats.entropy(bin_counts)


def _extract_edge_bearings(Gu, min_length=0, weight=None):
    """
    Extract undirected graph's bidirectional edge bearings.

    For example, if an edge has a bearing of 90° then we will record bearings
    of both 90° and 270° for this edge.

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
    bearings : numpy.array
        the graph's bidirectional edge bearings
    """
    if nx.is_directed(Gu) or projection.is_projected(Gu.graph["crs"]):  # pragma: no cover
        raise ValueError("graph must be undirected and unprojected to analyze edge bearings")
    bearings = list()
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
    bearings = np.array(bearings)
    bearings = bearings[~np.isnan(bearings)]
    bearings_r = (bearings - 180) % 360
    return np.concatenate([bearings, bearings_r])


def _bearings_distribution(Gu, num_bins, min_length=0, weight=None):
    """
    Compute distribution of bearings across evenly spaced bins.

    Prevents bin-edge effects around common values like 0° and 90° by
    initially creating twice as many bins as desired, then merging them in
    pairs. For example, if `num_bins=36` is provided, then each bin will
    represent 10° around the compass, with the first bin representing 355°-5°.

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
    bin_counts, bin_edges : tuple of numpy.array
        counts of bearings per bin and the bins edges
    """
    n = num_bins * 2
    bins = np.arange(n + 1) * 360 / n

    bearings = _extract_edge_bearings(Gu, min_length, weight)
    count, bin_edges = np.histogram(bearings, bins=bins)

    # move last bin to front, so eg 0.01° and 359.99° will be binned together
    count = np.roll(count, 1)
    bin_counts = count[::2] + count[1::2]

    # because we merged the bins, their edges are now only every other one
    bin_edges = bin_edges[range(0, len(bin_edges), 2)]
    return bin_counts, bin_edges


def plot_orientation(
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
    Plot a polar histogram of a spatial network's bidirectional edge bearings.

    Ignores self-loop edges as their bearings are undefined.

    For more info see: Boeing, G. 2019. "Urban Spatial Order: Street Network
    Orientation, Configuration, and Entropy." Applied Network Science, 4 (1),
    67. https://doi.org/10.1007/s41109-019-0189-1

    Parameters
    ----------
    Gu : networkx.MultiGraph
        undirected, unprojected graph with `bearing` attributes on each edge
    num_bins : int
        number of bins; for example, if `num_bins=36` is provided, then each
        bin will represent 10° around the compass
    min_length : float
        ignore edges with `length` attributes less than `min_length`
    weight : string
        if not None, weight edges' bearings by this (non-null) edge attribute
    ax : matplotlib.axes.PolarAxesSubplot
        if not None, plot on this preexisting axis; must have projection=polar
    figsize : tuple
        if ax is None, create new figure with size (width, height)
    area : bool
        if True, set bar length so area is proportional to frequency,
        otherwise set bar length so height is proportional to frequency
    color : string
        color of histogram bars
    edgecolor : string
        color of histogram bar edges
    linewidth : float
        width of histogram bar edges
    alpha : float
        opacity of histogram bars
    title : string
        title for plot
    title_y : float
        y position to place title
    title_font : dict
        the title's fontdict to pass to matplotlib
    xtick_font : dict
        the xtick labels' fontdict to pass to matplotlib

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axis
    """
    if title_font is None:
        title_font = {"family": "DejaVu Sans", "size": 24, "weight": "bold"}
    if xtick_font is None:
        xtick_font = {
            "family": "DejaVu Sans",
            "size": 10,
            "weight": "bold",
            "alpha": 1.0,
            "zorder": 3,
        }

    # get the bearings' distribution's bin counts and edges
    bin_counts, bin_edges = _bearings_distribution(Gu, num_bins, min_length, weight)

    # positions: where to center each bar. ignore the last bin edge, because
    # it's the same as the first (i.e., 0° = 360°)
    positions = np.radians(bin_edges[:-1])

    # width: make bars fill the circumference without gaps or overlaps
    width = 2 * np.pi / num_bins

    # radius: how long to make each bar
    bin_frequency = bin_counts / bin_counts.sum()
    if area:
        # set bar length so area is proportional to frequency
        radius = np.sqrt(bin_frequency)
    else:
        # set bar length so height is proportional to frequency
        radius = bin_frequency

    # create ax (if necessary) then set N at top and go clockwise
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={"projection": "polar"})
    else:
        fig = ax.figure
    ax.set_theta_zero_location("N")
    ax.set_theta_direction("clockwise")
    ax.set_ylim(top=radius.max())

    # configure the y-ticks and remove their labels
    ax.set_yticks(np.linspace(0, radius.max(), 5))
    ax.set_yticklabels(labels="")

    # configure the x-ticks and their labels
    xticklabels = ["N", "", "E", "", "S", "", "W", ""]
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(labels=xticklabels, fontdict=xtick_font)
    ax.tick_params(axis="x", which="major", pad=-2)

    # draw the bars
    ax.bar(
        positions,
        height=radius,
        width=width,
        align="center",
        bottom=0,
        zorder=2,
        color=color,
        edgecolor=edgecolor,
        linewidth=linewidth,
        alpha=alpha,
    )

    if title:
        ax.set_title(title, y=title_y, fontdict=title_font)
    fig.tight_layout()
    return fig, ax
