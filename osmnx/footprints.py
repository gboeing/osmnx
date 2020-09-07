"""
Download footprints from OpenStreetMap.

Deprecated: use the new `geometries` module instead.
"""

import warnings

from . import geometries


def footprints_from_point(point, dist=1000, footprint_type="building", retain_invalid=False):
    """
    Get footprints within some distance N, S, E, W of a lat-lng point.

    Deprecated: use geometries module instead.

    Parameters
    ----------
    point : tuple
        a lat-lng point
    dist : numeric
        distance in meters
    footprint_type : string
        type of footprint to be downloaded. OSM tag key e.g. 'building',
        'landuse', 'place', etc.
    retain_invalid : bool
        deprecated, is ignored

    Returns
    -------
    geopandas.GeoDataFrame

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """
    msg = (
        "The `footprints` module has been deprecated and will be removed in a "
        "future release. Instead, use the `geometries` module's "
        "`geometries_from_point` function, passing `tags={'building':True}`."
    )
    warnings.warn(msg)

    tags = {footprint_type: True}
    return geometries.geometries_from_point(point, tags, dist)


def footprints_from_address(address, dist=1000, footprint_type="building", retain_invalid=False):
    """
    Get footprints within some distance N, S, E, W of an address.

    Deprecated: use geometries module instead.

    Parameters
    ----------
    address : string
        the address to geocode to a lat-lng point
    dist : numeric
        distance in meters
    footprint_type : string
        type of footprint to be downloaded. OSM tag key e.g. 'building',
        'landuse', 'place', etc.
    retain_invalid : bool
        deprecated, is ignored

    Returns
    -------
    geopandas.GeoDataFrame

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """
    msg = (
        "The `footprints` module has been deprecated and will be removed in a "
        "future release. Instead, use the `geometries` module's "
        "`geometries_from_address` function, passing `tags={'building':True}`."
    )
    warnings.warn(msg)

    tags = {footprint_type: True}
    return geometries.geometries_from_address(address, tags, dist)


def footprints_from_place(
    place, footprint_type="building", retain_invalid=False, which_result=None
):
    """
    Get footprints within the boundaries of some place.

    Deprecated: use geometries module instead.

    The query must be geocodable and OSM must have polygon boundaries for the
    geocode result. If OSM does not have a polygon for this place, you can
    instead get its footprints using the footprints_from_address function,
    which geocodes the place name to a point and gets the footprints within
    some distance of that point.

    Parameters
    ----------
    place : string
        the query to geocode to get place boundary polygon
    footprint_type : string
        type of footprint to be downloaded. OSM tag key e.g. 'building',
        'landuse', 'place', etc.
    retain_invalid : bool
        deprecated, is ignored
    which_result : int
        which geocoding result to use. if None, auto-select the first
        multi/polygon or raise an error if OSM doesn't return one.

    Returns
    -------
    geopandas.GeoDataFrame

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """
    msg = (
        "The `footprints` module has been deprecated and will be removed in a "
        "future release. Instead, use the `geometries` module's "
        "`geometries_from_place` function, passing `tags={'building':True}`."
    )
    warnings.warn(msg)

    tags = {footprint_type: True}
    return geometries.geometries_from_place(place, tags, which_result=which_result)


def footprints_from_polygon(polygon, footprint_type="building", retain_invalid=False):
    """
    Get footprints within some polygon.

    Deprecated: use geometries module instead.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        the shape to get data within. coordinates should be in units of
        latitude-longitude degrees.
    footprint_type : string
        type of footprint to be downloaded. OSM tag key e.g. 'building',
        'landuse', 'place', etc.
    retain_invalid : bool
        deprecated, is ignored

    Returns
    -------
    geopandas.GeoDataFrame

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """
    msg = (
        "The `footprints` module has been deprecated and will be removed in a "
        "future release. Instead, use the `geometries` module's "
        "`geometries_from_polygon` function, passing `tags={'building':True}`."
    )
    warnings.warn(msg)

    tags = {footprint_type: True}
    return geometries.geometries_from_polygon(polygon, tags)
