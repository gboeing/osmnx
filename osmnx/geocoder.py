"""Geocode queries and create GeoDataFrames of place boundaries."""

import logging as lg
from collections import OrderedDict

import geopandas as gpd

from . import downloader
from . import projection
from . import settings
from . import utils


def geocode(query):
    """
    Geocode a query string to (lat, lng) with the Nominatim geocoder.

    Parameters
    ----------
    query : string
        the query string to geocode

    Returns
    -------
    point : tuple
        the (lat, lng) coordinates returned by the geocoder
    """
    # define the parameters
    params = OrderedDict()
    params["format"] = "json"
    params["limit"] = 1
    params["dedupe"] = 0  # prevent deduping to get precise number of results
    params["q"] = query
    response_json = downloader.nominatim_request(params=params)

    # if results were returned, parse lat and lng out of the result
    if response_json and "lat" in response_json[0] and "lon" in response_json[0]:
        lat = float(response_json[0]["lat"])
        lng = float(response_json[0]["lon"])
        point = (lat, lng)
        utils.log(f'Geocoded "{query}" to {point}')
        return point
    else:
        raise Exception(f'Nominatim geocoder returned no results for query "{query}"')


def geocode_to_gdf(query, which_result=None, by_osmid=False, buffer_dist=None):
    """
    Retrieve place(s) by name or ID from the Nominatim API as a GeoDataFrame.

    You can query by place name or OSM ID. If querying by place name, the
    query argument can be a string or structured dict, or a list of such
    strings/dicts to send to geocoder. You can instead query by OSM ID by
    setting `by_osmid=True`. In this case, geocode_to_gdf treats the query
    argument as an OSM ID (or list of OSM IDs) for Nominatim lookup rather
    than text search. OSM IDs must be prepended with their types: node (N),
    way (W), or relation (R), in accordance with the Nominatim format. For
    example, `query=["R2192363", "N240109189", "W427818536"]`.

    If query argument is a list, then which_result should be either a single
    value or a list with the same length as query. The queries you provide
    must be resolvable to places in the Nominatim database. The resulting
    GeoDataFrame's geometry column contains place boundaries if they exist in
    OpenStreetMap.

    Parameters
    ----------
    query : string or dict or list
        query string(s) or structured dict(s) to geocode
    which_result : int
        which geocoding result to use. if None, auto-select the first
        (Multi)Polygon or raise an error if OSM doesn't return one. to get
        the top match regardless of geometry type, set which_result=1
    by_osmid : bool
        if True, handle query as an OSM ID for lookup rather than text search
    buffer_dist : float
        distance to buffer around the place geometry, in meters

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        a GeoDataFrame with one row for each query
    """
    if not isinstance(query, (str, dict, list)):
        raise ValueError("query must be a string or dict or list")

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
    if len(query) != len(which_result):
        raise ValueError("which_result length must equal query length")

    # ensure query type of each item
    for q in query:
        if not isinstance(q, (str, dict)):
            raise ValueError("each query must be a dict or a string")

    # geocode each query and add to GeoDataFrame as a new row
    gdf = gpd.GeoDataFrame()
    for q, wr in zip(query, which_result):
        gdf = gdf.append(_geocode_query_to_gdf(q, wr, by_osmid))

    # reset GeoDataFrame index and set its CRS
    gdf = gdf.reset_index(drop=True)
    gdf.crs = settings.default_crs

    # if buffer_dist was passed in, project the geometry to UTM, buffer it in
    # meters, then project it back to lat-lng
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
        the top match regardless of geometry type, set which_result=1
    by_osmid : bool
        if True, handle query as an OSM ID for lookup rather than text search

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        a GeoDataFrame with one row containing the result of geocoding
    """
    if which_result is None:
        limit = 50
    else:
        limit = which_result

    results = downloader._osm_place_download(query, by_osmid=by_osmid, limit=limit)

    # choose the right result from the JSON response
    if not results:
        # if no results were returned, raise error
        raise ValueError(f'OSM returned no results for query "{query}"')

    elif by_osmid:
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
        msg = f'OSM returned fewer than `which_result={which_result}` results for query "{query}"'
        raise ValueError(msg)

    # if we got a non (Multi)Polygon geometry type (like a point), log warning
    geom_type = result["geojson"]["type"]
    if geom_type not in {"Polygon", "MultiPolygon"}:
        msg = f'OSM returned a {geom_type} as the geometry for query "{query}"'
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
        list of results from downloader._osm_place_download
    query : str
        the query string or structured dict that was geocoded

    Returns
    -------
    result : dict
        the chosen result
    """
    polygon_types = {"Polygon", "MultiPolygon"}
    for result in results:
        if result["geojson"]["type"] in polygon_types:
            return result

    # if we never found a polygon, throw an error
    raise ValueError(f'OSM did not return any polygonal geometries for query "{query}"')
