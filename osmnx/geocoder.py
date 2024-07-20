"""
Geocode place names or addresses or retrieve OSM elements by place name or ID.

This module uses the Nominatim API's "search" and "lookup" endpoints. For more
details see https://wiki.openstreetmap.org/wiki/Elements and
https://nominatim.org/.
"""

from __future__ import annotations

import logging as lg
from collections import OrderedDict
from typing import Any

import geopandas as gpd
import pandas as pd

from . import _nominatim
from . import settings
from . import utils
from ._errors import InsufficientResponseError


def geocode(query: str) -> tuple[float, float]:
    """
    Geocode place names or addresses to `(lat, lon)` with the Nominatim API.

    This geocodes the query via the Nominatim "search" endpoint.

    Parameters
    ----------
    query
        The query string to geocode.

    Returns
    -------
    point
        The `(lat, lon)` coordinates returned by the geocoder.
    """
    # define the parameters
    params: OrderedDict[str, int | str] = OrderedDict()
    params["format"] = "json"
    params["limit"] = 1
    params["dedupe"] = 0  # prevent deduping to get precise number of results
    params["q"] = query
    response_json = _nominatim._nominatim_request(params=params)

    # if results were returned, parse lat and lon out of the result
    if response_json and "lat" in response_json[0] and "lon" in response_json[0]:
        lat = float(response_json[0]["lat"])
        lon = float(response_json[0]["lon"])
        point = (lat, lon)

        msg = f"Geocoded {query!r} to {point}"
        utils.log(msg, level=lg.INFO)
        return point

    # otherwise we got no results back
    msg = f"Nominatim could not geocode query {query!r}."
    raise InsufficientResponseError(msg)


def geocode_to_gdf(
    query: str | dict[str, str] | list[str | dict[str, str]],
    *,
    which_result: int | None | list[int | None] = None,
    by_osmid: bool = False,
) -> gpd.GeoDataFrame:
    """
    Retrieve OSM elements by place name or OSM ID with the Nominatim API.

    If searching by place name, the `query` argument can be a string or
    structured dict, or a list of such strings/dicts to send to the geocoder.
    This uses the Nominatim "search" endpoint to geocode the place name to the
    best-matching OSM element, then returns that element and its attribute
    data.

    You can instead query by OSM ID by passing `by_osmid=True`. This uses the
    Nominatim "lookup" endpoint to retrieve the OSM element with that ID. In
    this case, the function treats the `query` argument as an OSM ID (or list
    of OSM IDs), which must be prepended with their types: node (N), way (W),
    or relation (R) in accordance with the Nominatim API format. For example,
    `query=["R2192363", "N240109189", "W427818536"]`.

    If `query` is a list, then `which_result` must be either an int or a list
    with the same length as `query`. The queries you provide must be
    resolvable to elements in the Nominatim database. The resulting
    GeoDataFrame's geometry column contains place boundaries if they exist.

    Parameters
    ----------
    query
        The query string(s) or structured dict(s) to geocode.
    which_result
        Which search result to return. If None, auto-select the first
        (Multi)Polygon or raise an error if OSM doesn't return one. To get
        the top match regardless of geometry type, set `which_result=1`.
        Ignored if `by_osmid=True`.
    by_osmid
        If True, treat query as an OSM ID lookup rather than text search.

    Returns
    -------
    gdf
        GeoDataFrame with one row for each query result.
    """
    if isinstance(query, list):
        # if query is a list of queries but which_result is int/None, then
        # turn which_result into a list with same length as query list
        q_list = query
        wr_list = which_result if isinstance(which_result, list) else [which_result] * len(query)
    else:
        # if query is not already a list, turn it into one
        # if which_result was a list, take 0th element, otherwise make it list
        q_list = [query]
        wr_list = [which_result[0]] if isinstance(which_result, list) else [which_result]

    # ensure same length
    if len(q_list) != len(wr_list):  # pragma: no cover
        msg = "`which_result` length must equal `query` length."
        raise ValueError(msg)

    # geocode each query, concat as GeoDataFrame rows, then set the CRS
    results = (_geocode_query_to_gdf(q, wr, by_osmid) for q, wr in zip(q_list, wr_list))
    gdf = pd.concat(results, ignore_index=True).set_crs(settings.default_crs)

    msg = f"Created GeoDataFrame with {len(gdf)} rows from {len(q_list)} queries"
    utils.log(msg, level=lg.INFO)
    return gdf


def _geocode_query_to_gdf(
    query: str | dict[str, str],
    which_result: int | None,
    by_osmid: bool,  # noqa: FBT001
) -> gpd.GeoDataFrame:
    """
    Geocode a single place query to a GeoDataFrame.

    Parameters
    ----------
    query
        Query string or structured dict to geocode.
    which_result
        Which search result to return. If None, auto-select the first
        (Multi)Polygon or raise an error if OSM doesn't return one. To get
        the top match regardless of geometry type, set `which_result=1`.
        Ignored if `by_osmid=True`.
    by_osmid
        If True, treat query as an OSM ID lookup rather than text search.

    Returns
    -------
    gdf
        GeoDataFrame with one row containing the geocoding result.
    """
    limit = 50 if which_result is None else which_result
    results = _nominatim._download_nominatim_element(query, by_osmid=by_osmid, limit=limit)

    # choose the right result from the JSON response
    if len(results) == 0:
        # if no results were returned, raise error
        msg = f"Nominatim geocoder returned 0 results for query {query!r}."
        raise InsufficientResponseError(msg)

    if by_osmid:
        # if searching by OSM ID, always take the first (ie, only) result
        result = results[0]

    elif which_result is None:
        # else, if which_result=None, auto-select the first (Multi)Polygon
        try:
            result = _get_first_polygon(results)
        except TypeError as e:
            msg = f"Nominatim did not geocode query {query!r} to a geometry of type (Multi)Polygon."
            raise TypeError(msg) from e

    elif len(results) >= which_result:
        # else, if we got at least which_result results, choose that one
        result = results[which_result - 1]

    else:  # pragma: no cover
        # else, we got fewer results than which_result, raise error
        msg = f"Nominatim returned {len(results)} result(s) but `which_result={which_result}`."
        raise InsufficientResponseError(msg)

    # if we got a non (Multi)Polygon geometry type (like a point), log warning
    geom_type = result["geojson"]["type"]
    if geom_type not in {"Polygon", "MultiPolygon"}:
        msg = f"Nominatim geocoder returned a {geom_type} as the geometry for query {query!r}"
        utils.log(msg, level=lg.WARNING)

    # build the GeoJSON feature from the chosen result
    bottom, top, left, right = result["boundingbox"]
    feature = {
        "type": "Feature",
        "geometry": result["geojson"],
        "properties": {
            "bbox_west": left,
            "bbox_south": bottom,
            "bbox_east": right,
            "bbox_north": top,
        },
    }

    # add the other attributes we retrieved
    for attr in result:
        if attr not in {"address", "boundingbox", "geojson", "icon", "licence"}:
            feature["properties"][attr] = result[attr]

    # create and return the GeoDataFrame
    gdf = gpd.GeoDataFrame.from_features([feature])
    cols = ["lat", "lon", "bbox_north", "bbox_south", "bbox_east", "bbox_west"]
    gdf[cols] = gdf[cols].astype(float)
    return gdf


def _get_first_polygon(results: list[dict[str, Any]]) -> dict[str, Any]:
    """
    Choose first result of geometry type (Multi)Polygon from list of results.

    Parameters
    ----------
    results
        Results from the Nominatim API.

    Returns
    -------
    result
        The chosen result.
    """
    polygon_types = {"Polygon", "MultiPolygon"}

    for result in results:
        if "geojson" in result and result["geojson"]["type"] in polygon_types:
            return result

    # if we never found a polygon, raise an error
    raise TypeError
