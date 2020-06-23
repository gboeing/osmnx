"""Create GeoDataFrames of place boundaries."""

import logging as lg
import warnings

import geopandas as gpd

from . import downloader
from . import projection
from . import settings
from . import utils


def gdf_from_place(query, which_result=1, buffer_dist=None):
    """
    Use `geocoder.geocode_to_gdf()` instead (deprecated).

    Parameters
    ----------
    query : string or dict
        query string or structured query dict to geocode/download
    which_result : int
        max number of results to return and which to process upon receipt
    buffer_dist : float
        distance to buffer around the place geometry, in meters

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    msg = (
        "The `boundaries` module has been deprecated and will be removed "
        "in a future relase. Use the `geocoder` module's `geocode_to_gdf` "
        "function instead."
    )
    warnings.warn(msg)

    # ensure query type
    if not isinstance(query, (str, dict)):
        raise ValueError("query must be a dict or a string")

    # get the data from OSM
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

        # create the GeoDataFrame, name it, and set its original CRS to default_crs
        gdf = gpd.GeoDataFrame.from_features(features)
        gdf.crs = settings.default_crs

        # if buffer_dist was passed in, project the geometry to UTM, buffer it
        # in meters, then project it back to lat-lng
        if buffer_dist is not None:
            gdf_utm = projection.project_gdf(gdf)
            gdf_utm["geometry"] = gdf_utm["geometry"].buffer(buffer_dist)
            gdf = projection.project_gdf(gdf_utm, to_latlong=True)
            utils.log(f"Buffered GeoDataFrame to {buffer_dist} meters")

        # return the gdf
        utils.log(f'Created GeoDataFrame with {len(gdf)} row for query "{query}"')
        return gdf
    else:
        # if no data returned (or fewer results than which_result)
        utils.log(
            f'OSM returned no results (or fewer than which_result) for query "{query}"',
            level=lg.WARNING,
        )
        return gpd.GeoDataFrame()


def gdf_from_places(queries, which_results=None, buffer_dist=None):
    """
    Use `geocoder.geocode_to_gdf()` instead (deprecated).

    Parameters
    ----------
    queries : list
        list of query strings or structured query dicts to geocode/download,
        one at a time
    which_results : list
        if not None, a list of max number of results to return and which to
        process upon receipt, for each query in queries
    buffer_dist : float
        distance to buffer around the place geometry, in meters

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    msg = (
        "The `boundaries` module has been deprecated and will be removed "
        "in a future relase. Use the `geocoder` module's `geocode_to_gdf` "
        "function instead."
    )
    warnings.warn(msg)

    # create an empty GeoDataFrame then append each result as a new row,
    # checking for the presence of which_results
    gdf = gpd.GeoDataFrame()
    if which_results is not None:

        if len(queries) != len(which_results):
            raise ValueError("which_results length must equal queries length")

        for query, which_result in zip(queries, which_results):
            gdf_tmp = gdf_from_place(query, buffer_dist=buffer_dist, which_result=which_result)
            gdf = gdf.append(gdf_tmp)
    else:
        for query in queries:
            gdf = gdf.append(gdf_from_place(query, buffer_dist=buffer_dist))

    # reset the index
    gdf = gdf.reset_index(drop=True)

    # set the original CRS of the GeoDataFrame to default_crs, and return it
    gdf.crs = settings.default_crs
    utils.log(f"Finished creating GeoDataFrame with {len(gdf)} rows from {len(queries)} queries")
    return gdf
