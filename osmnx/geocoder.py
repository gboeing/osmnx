"""
Geocode place names or addresses or retrieve OSM elements by place name or ID.

This module uses the Nominatim API's "search" and "lookup" endpoints. For more
details see https://wiki.openstreetmap.org/wiki/Elements and
https://nominatim.org/.
"""

import logging as lg
from collections import OrderedDict
from warnings import warn

import geopandas as gpd
import pandas as pd

from . import _nominatim
from . import projection
from . import settings
from . import utils
from ._errors import InsufficientResponseError


def geocode(query):
    """
    Geocode place names or addresses to (lat, lon) with the Nominatim API.

    This geocodes the query via the Nominatim "search" endpoint.

    Parameters
    ----------
    query : string
        the query string to geocode

    Returns
    -------
    point : tuple
        the (lat, lon) coordinates returned by the geocoder
    """
    # define the parameters
    params = OrderedDict()
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
        utils.log(f"Geocoded {query!r} to {point}")
        return point

    # otherwise we got no results back
    msg = f"Nominatim could not geocode query {query!r}"
    raise InsufficientResponseError(msg)


def geocode_to_gdf(query, which_result=None, by_osmid=False, buffer_dist=None):
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

    If `query` is a list, then `which_result` must be either a single value or
    a list with the same length as `query`. The queries you provide must be
    resolvable to elements in the Nominatim database. The resulting
    GeoDataFrame's geometry column contains place boundaries if they exist.

    Parameters
    ----------
    query : string or dict or list of strings/dicts
        query string(s) or structured dict(s) to geocode
    which_result : int
        which search result to return. if None, auto-select the first
        (Multi)Polygon or raise an error if OSM doesn't return one. to get
        the top match regardless of geometry type, set which_result=1.
        ignored if by_osmid=True.
    by_osmid : bool
        if True, treat query as an OSM ID lookup rather than text search
    buffer_dist : float
        deprecated, do not use

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        a GeoDataFrame with one row for each query
    """
    if buffer_dist is not None:
        warn(
            "The buffer_dist argument has been deprecated and will be removed "
            "in the v2.0.0 release. Buffer your results directly, if desired. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123",
            FutureWarning,
            stacklevel=2,
        )

    if not isinstance(query, (str, dict, list)):  # pragma: no cover
        msg = "query must be a string or dict or list"
        raise TypeError(msg)

    # if caller passed a list of queries but a scalar which_result value, then
    # turn which_result into a list with same length as query list
    if isinstance(query, list) and (isinstance(which_result, int) or which_result is None):
        which_result = [which_result] * len(query)

    # turn query and which_result into lists if they're not already
    if not isinstance(query, list):
        query = [query]
    if not isinstance(which_result, list):
        which_result = [which_result]

    # ensure same length
    if len(query) != len(which_result):  # pragma: no cover
        msg = "which_result length must equal query length"
        raise ValueError(msg)

    # ensure query type of each item
    for q in query:
        if not isinstance(q, (str, dict)):  # pragma: no cover
            msg = "each query must be a dict or a string"
            raise TypeError(msg)

    # geocode each query and add to GeoDataFrame as a new row
    gdf = gpd.GeoDataFrame()
    for q, wr in zip(query, which_result):
        gdf = pd.concat([gdf, _geocode_query_to_gdf(q, wr, by_osmid)])

    # reset GeoDataFrame index and set its CRS
    gdf = gdf.reset_index(drop=True)
    gdf = gdf.set_crs(settings.default_crs)

    # if buffer_dist was passed in, project the geometry to UTM, buffer it in
    # meters, then project it back to lat-lon
    if buffer_dist is not None and len(gdf) > 0:
        gdf_utm = projection.project_gdf(gdf)
        gdf_utm["geometry"] = gdf_utm["geometry"].buffer(buffer_dist)
        gdf = projection.project_gdf(gdf_utm, to_latlong=True)
        utils.log(f"Buffered GeoDataFrame to {buffer_dist} meters")

    utils.log(f"Created GeoDataFrame with {len(gdf)} rows from {len(query)} queries")
    return gdf


def _geocode_query_to_gdf(query, which_result, by_osmid):
    """
    Geocode a single place query to a GeoDataFrame.

    Parameters
    ----------
    query : string or dict
        query string or structured dict to geocode
    which_result : int
        which geocoding result to use. if None, auto-select the first
        (Multi)Polygon or raise an error if OSM doesn't return one. to get
        the top match regardless of geometry type, set which_result=1.
        ignored if by_osmid=True.
    by_osmid : bool
        if True, handle query as an OSM ID for lookup rather than text search

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        a GeoDataFrame with one row containing the result of geocoding
    """
    limit = 50 if which_result is None else which_result

    results = _nominatim._download_nominatim_element(query, by_osmid=by_osmid, limit=limit)

    # choose the right result from the JSON response
    if not results:
        # if no results were returned, raise error
        msg = f"Nominatim geocoder returned 0 results for query {query!r}"
        raise InsufficientResponseError(msg)

    if by_osmid:
        # if searching by OSM ID, always take the first (ie, only) result
        result = results[0]

    elif which_result is None:
        # else, if which_result=None, auto-select the first (Multi)Polygon
        result = _get_first_polygon(results, query)

    elif len(results) >= which_result:
        # else, if we got at least which_result results, choose that one
        result = results[which_result - 1]

    else:  # pragma: no cover
        # else, we got fewer results than which_result, raise error
        msg = f"Nominatim returned {len(results)} result(s) but which_result={which_result}"
        raise InsufficientResponseError(msg)

    # if we got a non (Multi)Polygon geometry type (like a point), log warning
    geom_type = result["geojson"]["type"]
    if geom_type not in {"Polygon", "MultiPolygon"}:
        msg = f"Nominatim geocoder returned a {geom_type} as the geometry for query {query!r}"
        utils.log(msg, level=lg.WARNING)

    # build the GeoJSON feature from the chosen result
    south, north, west, east = result["boundingbox"]
    feature = {
        "type": "Feature",
        "geometry": result["geojson"],
        "properties": {
            "bbox_north": north,
            "bbox_south": south,
            "bbox_east": east,
            "bbox_west": west,
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


def _get_first_polygon(results, query):
    """
    Choose first result of geometry type (Multi)Polygon from list of results.

    Parameters
    ----------
    results : list
        list of results from _downloader._osm_place_download
    query : str
        the query string or structured dict that was geocoded

    Returns
    -------
    result : dict
        the chosen result
    """
    polygon_types = {"Polygon", "MultiPolygon"}
    for result in results:
        if "geojson" in result and result["geojson"]["type"] in polygon_types:
            return result

    # if we never found a polygon, throw an error
    msg = f"Nominatim could not geocode query {query!r} to a geometry of type (Multi)Polygon"
    raise TypeError(msg)
