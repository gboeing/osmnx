"""Calculate graph edge speeds and travel times."""

import re
from warnings import warn

import networkx as nx
import numpy as np
import pandas as pd

from . import utils_graph


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
            "the `precision` parameter is deprecated and will be removed in a future release",
            stacklevel=2,
        )

    if fallback is None:
        fallback = np.nan

    edges = utils_graph.graph_to_gdfs(G, nodes=False, fill_edge_geometry=False)

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
    if hwy_speeds is None:
        hwy_speed_avg = pd.Series(dtype=float)
    else:
        hwy_speed_avg = pd.Series(hwy_speeds).dropna()

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
    if pd.isnull(speed_kph).all():
        msg = (
            "this graph's edges have no preexisting `maxspeed` attribute "
            "values so you must pass `hwy_speeds` or `fallback` arguments."
        )
        raise ValueError(msg)

    # add speed kph attribute to graph edges
    edges["speed_kph"] = speed_kph.round(precision).values
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
            "the `precision` parameter is deprecated and will be removed in a future release",
            stacklevel=2,
        )

    edges = utils_graph.graph_to_gdfs(G, nodes=False)

    # verify edge length and speed_kph attributes exist
    if not ("length" in edges.columns and "speed_kph" in edges.columns):  # pragma: no cover
        msg = "all edges must have `length` and `speed_kph` attributes."
        raise KeyError(msg)

    # verify edge length and speed_kph attributes contain no nulls
    if pd.isnull(edges["length"]).any() or pd.isnull(edges["speed_kph"]).any():  # pragma: no cover
        msg = "edge `length` and `speed_kph` values must be non-null."
        raise ValueError(msg)

    # convert distance meters to km, and speed km per hour to km per second
    distance_km = edges["length"] / 1000
    speed_km_sec = edges["speed_kph"] / (60 * 60)

    # calculate edge travel time in seconds
    travel_time = distance_km / speed_km_sec

    # add travel time attribute to graph edges
    edges["travel_time"] = travel_time.round(precision).values
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

    else:
        try:
            # clean each value in list and convert to kph if it is mph then
            # return a single aggregated value
            values = [_clean_maxspeed(x) for x in value]
            return int(agg(pd.Series(values).dropna()))
        except ValueError:
            return None
