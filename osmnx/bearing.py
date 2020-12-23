"""Calculate graph edge bearings."""

import math

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
    bearing = (initial_bearing + 360) % 360

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
