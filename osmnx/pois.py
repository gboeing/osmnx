"""
Download points of interests (POIs) from OpenStreetMap.

Deprecated: use the new `geometries` module instead.
"""

import warnings

from . import geometries


def pois_from_point(point, tags, dist=1000):
    """
    Get point of interests (POIs) within some distance N, S, E, W of a point.

    Deprecated: use geometries module instead.

    Parameters
    ----------
    point : tuple
        a (lat, lng) point
    tags : dict
        Dict of tags used for finding POIs from the selected area. Results
        returned are the union, not intersection of each individual tag.
        Each result matches at least one tag given. The dict keys should be
        OSM tags, (e.g., `amenity`, `landuse`, `highway`, etc) and the dict
        values should be either `True` to retrieve all items with the given
        tag, or a string to get a single tag-value combination, or a list of
        strings to get multiple values for the given tag. For example,
        `tags = {'amenity':True, 'landuse':['retail','commercial'],
        'highway':'bus_stop'}` would return all amenities, landuse=retail,
        landuse=commercial, and highway=bus_stop.
    dist : numeric
        distance in meters

    Returns
    -------
    gdf : geopandas.GeoDataFrame

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """
    msg = (
        "The `pois` module has been deprecated and will be removed in a "
        "future release. Use the `geometries` module's `geometries_from_point` "
        "function instead."
    )
    warnings.warn(msg)
    return geometries.geometries_from_point(point, tags, dist)


def pois_from_address(address, tags, dist=1000):
    """
    Get point of interests (POIs) within some distance N, S, E, W of address.

    Deprecated: use geometries module instead.

    Parameters
    ----------
    address : string
        the address to geocode to a lat-lng point
    tags : dict
        Dict of tags used for finding POIs from the selected area. Results
        returned are the union, not intersection of each individual tag.
        Each result matches at least one tag given. The dict keys should be
        OSM tags, (e.g., `amenity`, `landuse`, `highway`, etc) and the dict
        values should be either `True` to retrieve all items with the given
        tag, or a string to get a single tag-value combination, or a list of
        strings to get multiple values for the given tag. For example,
        `tags = {'amenity':True, 'landuse':['retail','commercial'],
        'highway':'bus_stop'}` would return all amenities, landuse=retail,
        landuse=commercial, and highway=bus_stop.
    dist : numeric
        distance in meters

    Returns
    -------
    gdf : geopandas.GeoDataFrame

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """
    msg = (
        "The `pois` module has been deprecated and will be removed in a "
        "future release. Use the `geometries` module's `geometries_from_address` "
        "function instead."
    )
    warnings.warn(msg)
    return geometries.geometries_from_address(address, tags, dist)


def pois_from_place(place, tags, which_result=None):
    """
    Get points of interest (POIs) within the boundaries of some place.

    Deprecated: use geometries module instead.

    Parameters
    ----------
    place : string
        the query to geocode to get place boundary polygon
    tags : dict
        Dict of tags used for finding POIs from the selected area. Results
        returned are the union, not intersection of each individual tag.
        Each result matches at least one tag given. The dict keys should be
        OSM tags, (e.g., `amenity`, `landuse`, `highway`, etc) and the dict
        values should be either `True` to retrieve all items with the given
        tag, or a string to get a single tag-value combination, or a list of
        strings to get multiple values for the given tag. For example,
        `tags = {'amenity':True, 'landuse':['retail','commercial'],
        'highway':'bus_stop'}` would return all amenities, landuse=retail,
        landuse=commercial, and highway=bus_stop.
    which_result : int
        which geocoding result to use. if None, auto-select the first
        multi/polygon or raise an error if OSM doesn't return one.

    Returns
    -------
    gdf : geopandas.GeoDataFrame

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """
    msg = (
        "The `pois` module has been deprecated and will be removed in a "
        "future release. Use the `geometries` module's `geometries_from_place` "
        "function instead."
    )
    warnings.warn(msg)
    return geometries.geometries_from_place(place, tags, which_result=which_result)


def pois_from_polygon(polygon, tags):
    """
    Get point of interests (POIs) within some polygon.

    Deprecated: use geometries module instead.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon
        geographic boundaries to fetch POIs within
    tags : dict
        Dict of tags used for finding POIs from the selected area. Results
        returned are the union, not intersection of each individual tag.
        Each result matches at least one tag given. The dict keys should be
        OSM tags, (e.g., `amenity`, `landuse`, `highway`, etc) and the dict
        values should be either `True` to retrieve all items with the given
        tag, or a string to get a single tag-value combination, or a list of
        strings to get multiple values for the given tag. For example,
        `tags = {'amenity':True, 'landuse':['retail','commercial'],
        'highway':'bus_stop'}` would return all amenities, landuse=retail,
        landuse=commercial, and highway=bus_stop.

    Returns
    -------
    gdf : geopandas.GeoDataFrame

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """
    msg = (
        "The `pois` module has been deprecated and will be removed in a "
        "future release. Use the `geometries` module's `geometries_from_polygon` "
        "function instead."
    )
    warnings.warn(msg)
    return geometries.geometries_from_polygon(polygon, tags)
