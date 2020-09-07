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
    params["dedupe"] = 0  # prevent OSM deduping results so we get precisely 'limit' # of results
    params["q"] = query
    response_json = downloader.nominatim_request(params=params)

    # if results were returned, parse lat and lng out of the result
    if len(response_json) > 0 and "lat" in response_json[0] and "lon" in response_json[0]:
        lat = float(response_json[0]["lat"])
        lng = float(response_json[0]["lon"])
        point = (lat, lng)
        utils.log(f'Geocoded "{query}" to {point}')
        return point
    else:
        raise Exception(f'Nominatim geocoder returned no results for query "{query}"')


def geocode_to_gdf(query, which_result=None, buffer_dist=None):
    """
    Geocode a query or queries to a GeoDataFrame with the Nominatim geocoder.

    Geometry column contains place boundaries if they exist in OpenStreetMap.
    Query can be a string or dict, or a list of strings/dicts to send to the
    geocoder. If query is a list, then which_result should be either a single
    value or a list of the same length as query.

    Parameters
    ----------
    query : string or dict or list
        query string(s) or structured dict(s) to geocode
    which_result : int
        which geocoding result to use. if None, auto-select the first
        multi/polygon or raise an error if OSM doesn't return one.
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
        gdf = gdf.append(_geocode_query_to_gdf(q, wr))

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


def _geocode_query_to_gdf(query, which_result):
    """
    Geocode a single place query to a GeoDataFrame.

    Parameters
    ----------
    query : string or dict
        query string or structured dict to geocode
    which_result : int
        which geocoding result to use. if None, auto-select the first
        multi/polygon or raise an error if OSM doesn't return one.

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        a GeoDataFrame with one row containing the result of geocoding
    """
    if which_result is None:
        limit = 50
    else:
        limit = which_result

    results = downloader._osm_polygon_download(query, limit=limit)

    # choose the right result from the JSON response
    if len(results) == 0:
        # if no results were returned, raise error
        raise ValueError(f'OSM returned no results for query "{query}"')

    elif which_result is None:
        # else, if which_result=None, auto-select the first multi/polygon
        result = _get_first_polygon(results, query)

    elif len(results) >= which_result:
        # else, if we got at least which_result results, choose that one
        result = results[which_result - 1]

    else:
        # else, we got fewer results than which_result, raise error
        msg = f'OSM returned fewer than `which_result={which_result}` results for query "{query}"'
        raise ValueError(msg)

    # build the geojson feature from the chosen result
    south, north, west, east = [float(x) for x in result["boundingbox"]]
    place = result["display_name"]
    geometry = result["geojson"]
    features = [
        {
            "type": "Feature",
            "geometry": geometry,
            "properties": {
                "place_name": place,
                "bbox_north": north,
                "bbox_south": south,
                "bbox_east": east,
                "bbox_west": west,
            },
        }
    ]

    # if we got a non multi/polygon geometry type (like a point), log warning
    if geometry["type"] not in {"Polygon", "MultiPolygon"}:
        msg = f'OSM returned a {geometry["type"]} as the geometry for query "{query}"'
        utils.log(msg, level=lg.WARNING)

    # create and return the GeoDataFrame
    return gpd.GeoDataFrame.from_features(features)


def _get_first_polygon(results, query):
    """
    Choose first result with geometry type multi/polygon from list of results.

    Parameters
    ----------
    results : list
        list of results from downloader._osm_polygon_download
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
