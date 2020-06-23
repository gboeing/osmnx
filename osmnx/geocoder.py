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


def geocode_to_gdf(query, which_result=1, buffer_dist=None):
    """
    Geocode a query or queries to a GeoDataFrame with the Nominatim geocoder.

    Geometry column contains place boundaries if they exist in OpenStreetMap.
    Query can be a string or dict, or a list of strings/dicts to send to the
    geocoder. If query is a list, then which_result should be a list of the
    same length.

    Parameters
    ----------
    query : string or dict or list
        query string or structured dict to geocode/download
    which_result : int or list
        max number of results to return and which to process upon receipt; if
        passing a list then it must be same length as query list
    buffer_dist : float
        distance to buffer around the place geometry, in meters

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    if not isinstance(query, (str, dict, list)):
        raise ValueError("query must be a string or dict or list")

    # if caller passed a list of queries but a scalar which_result value, then
    # turn which_result into a list with same length as query list
    if isinstance(query, list) and isinstance(which_result, int):
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


def _geocode_query_to_gdf(query, which_result=1):
    """
    Geocode a single place query to a GeoDataFrame.

    Parameters
    ----------
    query : string or dict
        query string or structured dict to geocode/download
    which_result : int
        max number of results to return and which to process upon receipt

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    data = downloader._osm_polygon_download(query, limit=which_result)

    if len(data) >= which_result:
        # extract data elements from the JSON response
        result = data[which_result - 1]
        bbox_south, bbox_north, bbox_west, bbox_east = [float(x) for x in result["boundingbox"]]
        geometry = result["geojson"]
        place = result["display_name"]
        features = [
            {
                "type": "Feature",
                "geometry": geometry,
                "properties": {
                    "place_name": place,
                    "bbox_north": bbox_north,
                    "bbox_south": bbox_south,
                    "bbox_east": bbox_east,
                    "bbox_west": bbox_west,
                },
            }
        ]

        # if we got an unexpected geometry type (like a point), log a warning
        if geometry["type"] not in {"Polygon", "MultiPolygon"}:
            utils.log(f'OSM returned a {geometry["type"]} as the geometry', level=lg.WARNING)

        # create the GeoDataFrame
        return gpd.GeoDataFrame.from_features(features)
    else:
        # if no data returned (or fewer results than which_result)
        msg = f'OSM returned no results (or fewer than which_result) for query "{query}"'
        utils.log(msg, level=lg.WARNING)
        return gpd.GeoDataFrame()
