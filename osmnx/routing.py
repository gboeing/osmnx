"""Calculate weighted shortest paths between graph nodes."""

import itertools
import multiprocessing as mp
import re
from warnings import warn

import networkx as nx
import numpy as np
import pandas as pd

from . import convert
from . import utils


def route_to_gdf(G, route, weight="length"):
    """
    Return a GeoDataFrame of the edges in a path, in order.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    route : list
        list of node IDs constituting the path
    weight : string
        if there are parallel edges between two nodes, choose lowest weight

    Returns
    -------
    gdf_edges : geopandas.GeoDataFrame
        GeoDataFrame of the edges
    """
    pairs = zip(route[:-1], route[1:])
    uvk = ((u, v, min(G[u][v].items(), key=lambda i: i[1][weight])[0]) for u, v in pairs)
    return convert.graph_to_gdfs(G.subgraph(route), nodes=False).loc[uvk]


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
    paths_gen = nx.shortest_simple_paths(convert.to_digraph(G, weight), orig, dest, weight)
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


def add_edge_speeds(G, hwy_speeds=None, fallback=None, precision=None, agg=np.mean):
    """
    Add edge speeds (km per hour) to graph as new `speed_kph` edge attributes.

    By default, this imputes free-flow travel speeds for all edges via the
    mean `maxspeed` value of the edges of each highway type. For highway types
    in the graph that have no `maxspeed` value on any edge, it assigns the
    mean of all `maxspeed` values in graph.

    This default mean-imputation can obviously be imprecise, and the user can
    override it by passing in `hwy_speeds` and/or `fallback` arguments that
    correspond to local speed limit standards. The user can also specify a
    different aggregation function (such as the median) to impute missing
    values from the observed values.

    If edge `maxspeed` attribute has "mph" in it, value will automatically be
    converted from miles per hour to km per hour. Any other speed units should
    be manually converted to km per hour prior to running this function,
    otherwise there could be unexpected results. If "mph" does not appear in
    the edge's maxspeed attribute string, then function assumes kph, per OSM
    guidelines: https://wiki.openstreetmap.org/wiki/Map_Features/Units

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    hwy_speeds : dict
        dict keys = OSM highway types and values = typical speeds (km per
        hour) to assign to edges of that highway type for any edges missing
        speed data. Any edges with highway type not in `hwy_speeds` will be
        assigned the mean preexisting speed value of all edges of that highway
        type.
    fallback : numeric
        default speed value (km per hour) to assign to edges whose highway
        type did not appear in `hwy_speeds` and had no preexisting speed
        values on any edge
    precision : int
        deprecated, do not use
    agg : function
        aggregation function to impute missing values from observed values.
        the default is numpy.mean, but you might also consider for example
        numpy.median, numpy.nanmedian, or your own custom function

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with speed_kph attributes on all edges
    """
    if precision is None:
        precision = 1
    else:
        warn(
            "The `precision` parameter is deprecated and will be removed in the v2.0.0 release. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123",
            FutureWarning,
            stacklevel=2,
        )

    if fallback is None:
        fallback = np.nan

    edges = convert.graph_to_gdfs(G, nodes=False, fill_edge_geometry=False)

    # collapse any highway lists (can happen during graph simplification)
    # into string values simply by keeping just the first element of the list
    edges["highway"] = edges["highway"].map(lambda x: x[0] if isinstance(x, list) else x)

    if "maxspeed" in edges.columns:
        # collapse any maxspeed lists (can happen during graph simplification)
        # into a single value
        edges["maxspeed"] = edges["maxspeed"].apply(_collapse_multiple_maxspeed_values, agg=agg)

        # create speed_kph by cleaning maxspeed strings and converting mph to
        # kph if necessary
        edges["speed_kph"] = edges["maxspeed"].astype(str).map(_clean_maxspeed).astype(float)
    else:
        # if no edges in graph had a maxspeed attribute
        edges["speed_kph"] = None

    # if user provided hwy_speeds, use them as default values, otherwise
    # initialize an empty series to populate with values
    hwy_speed_avg = pd.Series(dtype=float) if hwy_speeds is None else pd.Series(hwy_speeds).dropna()

    # for each highway type that caller did not provide in hwy_speeds, impute
    # speed of type by taking the mean of the preexisting speed values of that
    # highway type
    for hwy, group in edges.groupby("highway"):
        if hwy not in hwy_speed_avg:
            hwy_speed_avg.loc[hwy] = agg(group["speed_kph"])

    # if any highway types had no preexisting speed values, impute their speed
    # with fallback value provided by caller. if fallback=np.nan, impute speed
    # as the mean speed of all highway types that did have preexisting values
    hwy_speed_avg = hwy_speed_avg.fillna(fallback).fillna(agg(hwy_speed_avg))

    # for each edge missing speed data, assign it the imputed value for its
    # highway type
    speed_kph = (
        edges[["highway", "speed_kph"]].set_index("highway").iloc[:, 0].fillna(hwy_speed_avg)
    )

    # all speeds will be null if edges had no preexisting maxspeed data and
    # caller did not pass in hwy_speeds or fallback arguments
    if pd.isna(speed_kph).all():
        msg = (
            "this graph's edges have no preexisting `maxspeed` attribute "
            "values so you must pass `hwy_speeds` or `fallback` arguments."
        )
        raise ValueError(msg)

    # add speed kph attribute to graph edges
    edges["speed_kph"] = speed_kph.round(precision).to_numpy()
    nx.set_edge_attributes(G, values=edges["speed_kph"], name="speed_kph")

    return G


def add_edge_travel_times(G, precision=None):
    """
    Add edge travel time (seconds) to graph as new `travel_time` edge attributes.

    Calculates free-flow travel time along each edge, based on `length` and
    `speed_kph` attributes. Note: run `add_edge_speeds` first to generate the
    `speed_kph` attribute. All edges must have `length` and `speed_kph`
    attributes and all their values must be non-null.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    precision : int
        deprecated, do not use

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with travel_time attributes on all edges
    """
    if precision is None:
        precision = 1
    else:
        warn(
            "The `precision` parameter is deprecated and will be removed in the v2.0.0 release. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123",
            FutureWarning,
            stacklevel=2,
        )

    edges = convert.graph_to_gdfs(G, nodes=False)

    # verify edge length and speed_kph attributes exist
    if not ("length" in edges.columns and "speed_kph" in edges.columns):  # pragma: no cover
        msg = "all edges must have `length` and `speed_kph` attributes."
        raise KeyError(msg)

    # verify edge length and speed_kph attributes contain no nulls
    if pd.isna(edges["length"]).any() or pd.isna(edges["speed_kph"]).any():  # pragma: no cover
        msg = "edge `length` and `speed_kph` values must be non-null."
        raise ValueError(msg)

    # convert distance meters to km, and speed km per hour to km per second
    distance_km = edges["length"] / 1000
    speed_km_sec = edges["speed_kph"] / (60 * 60)

    # calculate edge travel time in seconds
    travel_time = distance_km / speed_km_sec

    # add travel time attribute to graph edges
    edges["travel_time"] = travel_time.round(precision).to_numpy()
    nx.set_edge_attributes(G, values=edges["travel_time"], name="travel_time")

    return G


def _clean_maxspeed(maxspeed, agg=np.mean, convert_mph=True):
    """
    Clean a maxspeed string and convert mph to kph if necessary.

    If present, splits maxspeed on "|" (which denotes that the value contains
    different speeds per lane) then aggregates the resulting values. Invalid
    inputs return None. See https://wiki.openstreetmap.org/wiki/Key:maxspeed
    for details on values and formats.

    Parameters
    ----------
    maxspeed : string
        a valid OpenStreetMap way maxspeed value
    agg : function
        aggregation function if maxspeed contains multiple values (default
        is numpy.mean)
    convert_mph : bool
        if True, convert miles per hour to km per hour

    Returns
    -------
    clean_value : string
    """
    MILES_TO_KM = 1.60934
    # regex adapted from OSM wiki
    pattern = "^([0-9][\\.,0-9]+?)(?:[ ]?(?:km/h|kmh|kph|mph|knots))?$"
    values = re.split(r"\|", maxspeed)  # creates a list even if it's a single value
    try:
        clean_values = []
        for value in values:
            match = re.match(pattern, value)
            clean_value = float(match.group(1).replace(",", "."))
            if convert_mph and "mph" in maxspeed.lower():
                clean_value = clean_value * MILES_TO_KM
            clean_values.append(clean_value)
        return agg(clean_values)

    except (ValueError, AttributeError):
        # if invalid input, return None
        return None


def _collapse_multiple_maxspeed_values(value, agg):
    """
    Collapse a list of maxspeed values to a single value.

    Parameters
    ----------
    value : list or string
        an OSM way maxspeed value, or a list of them
    agg : function
        the aggregation function to reduce the list to a single value

    Returns
    -------
    agg_value : int
        an integer representation of the aggregated value in the list,
        converted to kph if original value was in mph.
    """
    # if this isn't a list, just return it right back to the caller
    if not isinstance(value, list):
        return value

    # otherwise, if it is a list, process it
    try:
        # clean each value in list and convert to kph if it is mph then
        # return a single aggregated value
        values = [_clean_maxspeed(x) for x in value]
        return int(agg(pd.Series(values).dropna()))
    except ValueError:
        return None
