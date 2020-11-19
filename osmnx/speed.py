"""Calculate graph edge speeds and travel times."""

import re

import networkx as nx
import numpy as np
import pandas as pd

from . import utils_graph


def add_edge_speeds(G, hwy_speeds=None, fallback=None, precision=1):
    """
    Add edge speeds (km per hour) to graph as new `speed_kph` edge attributes.

    Imputes free-flow travel speeds for all edges based on mean `maxspeed`
    value of edges, per highway type. For highway types in graph that have no
    `maxspeed` value on any edge, function assigns the mean of all `maxspeed`
    values in graph.

    This mean-imputation can obviously be imprecise, and the caller can
    override it by passing in `hwy_speeds` and/or `fallback` arguments that
    correspond to local speed limit standards.

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
        decimal precision to round speed_kph

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with speed_kph attributes on all edges
    """
    if fallback is None:
        fallback = np.nan

    edges = utils_graph.graph_to_gdfs(G, nodes=False, fill_edge_geometry=False)

    # collapse any highway lists (can happen during graph simplification)
    # into string values simply by keeping just the first element of the list
    edges["highway"] = edges["highway"].map(lambda x: x[0] if isinstance(x, list) else x)

    if "maxspeed" in edges.columns:
        # collapse any maxspeed lists (can happen during graph simplification)
        # into a single value
        edges["maxspeed"] = edges["maxspeed"].map(_collapse_multiple_maxspeed_values)

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
            hwy_speed_avg.loc[hwy] = group["speed_kph"].mean()

    # if any highway types had no preexisting speed values, impute their speed
    # with fallback value provided by caller. if fallback=np.nan, impute speed
    # as the mean speed of all highway types that did have preexisting values
    hwy_speed_avg = hwy_speed_avg.fillna(fallback).fillna(hwy_speed_avg.mean())

    # for each edge missing speed data, assign it the imputed value for its
    # highway type
    speed_kph = (
        edges[["highway", "speed_kph"]].set_index("highway").iloc[:, 0].fillna(hwy_speed_avg)
    )

    # all speeds will be null if edges had no preexisting maxspeed data and
    # caller did not pass in hwy_speeds or fallback arguments
    if pd.isnull(speed_kph).all():
        raise ValueError(
            (
                "this graph's edges have no preexisting `maxspeed` "
                "attribute values so you must pass `hwy_speeds` or "
                "`fallback` arguments."
            )
        )

    # add speed kph attribute to graph edges
    edges["speed_kph"] = speed_kph.round(precision).values
    nx.set_edge_attributes(G, values=edges["speed_kph"], name="speed_kph")

    return G


def add_edge_travel_times(G, precision=1):
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
        decimal precision to round travel_time

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with travel_time attributes on all edges
    """
    edges = utils_graph.graph_to_gdfs(G, nodes=False)

    # verify edge length and speed_kph attributes exist and contain no nulls
    if not ("length" in edges.columns and "speed_kph" in edges.columns):
        raise KeyError("all edges must have `length` and `speed_kph` attributes.")
    else:
        if pd.isnull(edges["length"]).any() or pd.isnull(edges["speed_kph"]).any():
            raise ValueError("edge `length` and `speed_kph` values must be non-null.")

    # convert distance km to meters, and speed km per hour to km per second
    distance_km = edges["length"] / 1000
    speed_km_sec = edges["speed_kph"] / (60 * 60)

    # calculate edge travel time in seconds
    travel_time = distance_km / speed_km_sec

    # add travel time attribute to graph edges
    edges["travel_time"] = travel_time.round(precision).values
    nx.set_edge_attributes(G, values=edges["travel_time"], name="travel_time")

    return G


def _clean_maxspeed(value, convert_mph=True):
    """
    Clean a maxspeed string and convert mph to kph if necessary.

    Parameters
    ----------
    value : string
        an OSM way maxspeed value
    convert_mph : bool
        if True, convert mph to kph

    Returns
    -------
    value_clean : string
    """
    MPH_TO_KPH = 1.60934
    pattern = re.compile(r"[^\d\.,;]")

    try:
        # strip out everything but numbers, periods, commas, semicolons
        value_clean = float(re.sub(pattern, "", value).replace(",", "."))
        if convert_mph and "mph" in value.lower():
            value_clean = value_clean * MPH_TO_KPH
        return value_clean

    except ValueError:
        return None


def _collapse_multiple_maxspeed_values(value):
    """
    Collapse a list of maxspeed values into its mean value.

    Parameters
    ----------
    value : list or string
        an OSM way maxspeed value, or a list of them

    Returns
    -------
    mean_value : int
        an integer representation of the mean value in the list, converted
        to kph if original value was in mph.
    """
    # if this isn't a list, just return it right back to the caller
    if not isinstance(value, list):
        return value

    else:
        try:
            # clean each value in list and convert to kph if it is mph then
            # return mean value
            values = [_clean_maxspeed(x) for x in value]
            mean_value = int(pd.Series(values).dropna().mean())
            return mean_value
        except ValueError:
            return None
