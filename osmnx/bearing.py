"""Calculate graph edge bearings."""

import math

import networkx as nx
import numpy as np


def get_bearing(origin_point, destination_point):
    """
    Calculate the bearing between two lat-lng points.

    Each argument tuple should represent (lat, lng) as decimal degrees.
    Bearing represents angle in degrees (clockwise) between north and the
    direction from the origin point to the destination point.

    Parameters
    ----------
    origin_point : tuple
        (lat, lng)
    destination_point : tuple
        (lat, lng)

    Returns
    -------
    bearing : float
        the compass bearing in decimal degrees from the origin point to the
        destination point
    """
    if not (isinstance(origin_point, tuple) and isinstance(destination_point, tuple)):
        raise TypeError("origin_point and destination_point must be (lat, lng) tuples")

    # get latitudes and the difference in longitude, as radians
    lat1 = math.radians(origin_point[0])
    lat2 = math.radians(destination_point[0])
    diff_lng = math.radians(destination_point[1] - origin_point[1])

    # calculate initial bearing from -180 degrees to +180 degrees
    x = math.sin(diff_lng) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - (math.sin(lat1) * math.cos(lat2) * math.cos(diff_lng))
    initial_bearing = math.atan2(x, y)

    # normalize initial bearing to 0-360 degrees to get compass bearing
    initial_bearing = math.degrees(initial_bearing)
    bearing = initial_bearing % 360

    return bearing


def add_edge_bearings(G, precision=1):
    """
    Add `bearing` attributes to all graph edges.

    Calculate the compass bearing from origin node to destination node for
    each edge in the directed graph then add each bearing as a new edge
    attribute. Bearing represents angle in degrees (clockwise) between north
    and the direction from the origin node to the destination node.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    precision : int
        decimal precision to round bearing

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with edge bearing attributes
    """
    for u, v, data in G.edges(keys=False, data=True):

        if u == v:
            # a self-loop has an undefined compass bearing
            data["bearing"] = np.nan

        else:
            # calculate bearing from edge's origin to its destination
            origin_point = (G.nodes[u]["y"], G.nodes[u]["x"])
            destination_point = (G.nodes[v]["y"], G.nodes[v]["x"])
            bearing = get_bearing(origin_point, destination_point)
            data["bearing"] = round(bearing, precision)

    return G


def _extract_edge_bearings(Gu, min_length=0, weight=None):
    """
    Extract undirected graph's bidirectional edge bearings.

    For example, if an edge has a bearing of 90° then we will record bearings
    of both 90° and 270° for this edge.

    Parameters
    ----------
    Gu : networkx.MultiGraph
        undirected input graph
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
    if nx.is_directed(Gu):
        raise ValueError("`Gu` must be undirected")
    bearings = list()
    for _, _, data in add_edge_bearings(Gu).edges(data=True):
        if data["length"] >= min_length:
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


def _bearings_distribution(bearings, num_bins):
    """
    Compute distribution of bearings across evenly spaced bins.

    Prevents bin-edge effects around common values like 0° and 90° by
    initially creating twice as many bins as desired, then merging them in
    pairs. For example, if `num_bins=36` is provided, then each bin will
    represent 10° around the compass, with the first bin representing 355°-5°.

    Parameters
    ----------
    bearings : numpy.array
        the graph's bidirectional edge bearings
    num_bins : int
        number of bins for the bearings histogram

    Returns
    -------
    bin_counts : numpy.array
        counts of bearings per bin
    """
    n = num_bins * 2
    bin_edges = np.arange(n + 1) * 360 / n
    count, _ = np.histogram(bearings, bins=bin_edges)

    # move last bin to front, so eg 0.01° and 359.99° will be binned together
    count = np.roll(count, 1)
    return count[::2] + count[1::2]
