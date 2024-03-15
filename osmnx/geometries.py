"""
Do not use: deprecated.

The `geometries` module has been renamed the `features` module. The
`geometries` module is deprecated and will be removed in the v2.0.0 release.
"""

from warnings import warn

from . import features

DEP_MSG = (
    "The `geometries` module and `geometries_from_X` functions have been "
    "renamed the `features` module and `features_from_X` functions. Use these "
    "instead. The `geometries` module and function names are deprecated and "
    "will be removed in the v2.0.0 release. "
    "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
)


def geometries_from_bbox(north, south, east, west, tags):
    """
    Do not use: deprecated.

    The `geometries` module and `geometries_from_X` functions have been
    renamed the `features` module and `features_from_X` functions. Use these
    instead. The `geometries` module and functions are deprecated and will be
    removed in the v2.0.0 release.

    Parameters
    ----------
    north : float
        Do not use: deprecated.
    south : float
        Do not use: deprecated.
    east : float
        Do not use: deprecated.
    west : float
        Do not use: deprecated.
    tags : dict
        Do not use: deprecated.

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    warn(DEP_MSG, FutureWarning, stacklevel=2)
    return features.features_from_bbox(north, south, east, west, tags=tags)


def geometries_from_point(center_point, tags, dist=1000):
    """
    Do not use: deprecated.

    The `geometries` module and `geometries_from_X` functions have been
    renamed the `features` module and `features_from_X` functions. Use these
    instead. The `geometries` module and functions are deprecated and will be
    removed in the v2.0.0 release.

    Parameters
    ----------
    center_point : tuple
        Do not use: deprecated.
    tags : dict
        Do not use: deprecated.
    dist : numeric
        Do not use: deprecated.

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    warn(DEP_MSG, FutureWarning, stacklevel=2)
    return features.features_from_point(center_point, tags, dist)


def geometries_from_address(address, tags, dist=1000):
    """
    Do not use: deprecated.

    The `geometries` module and `geometries_from_X` functions have been
    renamed the `features` module and `features_from_X` functions. Use these
    instead. The `geometries` module and functions are deprecated and will be
    removed in the v2.0.0 release.

    Parameters
    ----------
    address : string
        Do not use: deprecated.
    tags : dict
        Do not use: deprecated.
    dist : numeric
        Do not use: deprecated.

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    warn(DEP_MSG, FutureWarning, stacklevel=2)
    return features.features_from_address(address, tags, dist)


def geometries_from_place(query, tags, which_result=None, buffer_dist=None):
    """
    Do not use: deprecated.

    The `geometries` module and `geometries_from_X` functions have been
    renamed the `features` module and `features_from_X` functions. Use these
    instead. The `geometries` module and functions are deprecated and will be
    removed in the v2.0.0 release.

    Parameters
    ----------
    query : string or dict or list
        Do not use: deprecated.
    tags : dict
        Do not use: deprecated.
    which_result : int
        Do not use: deprecated.
    buffer_dist : float
        Do not use: deprecated.

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    warn(DEP_MSG, FutureWarning, stacklevel=2)
    return features.features_from_place(query, tags, which_result, buffer_dist)


def geometries_from_polygon(polygon, tags):
    """
    Do not use: deprecated.

    The `geometries` module and `geometries_from_X` functions have been
    renamed the `features` module and `features_from_X` functions. Use these
    instead. The `geometries` module and functions are deprecated and will be
    removed in the v2.0.0 release.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        Do not use: deprecated.
    tags : dict
        Do not use: deprecated.

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    warn(DEP_MSG, FutureWarning, stacklevel=2)
    return features.features_from_polygon(polygon, tags)


def geometries_from_xml(filepath, polygon=None, tags=None):
    """
    Do not use: deprecated.

    The `geometries` module and `geometries_from_X` functions have been
    renamed the `features` module and `features_from_X` functions. Use these
    instead. The `geometries` module and functions are deprecated and will be
    removed in the v2.0.0 release.

    Parameters
    ----------
    filepath : string or pathlib.Path
        Do not use: deprecated.
    polygon : shapely.geometry.Polygon
        Do not use: deprecated.
    tags : dict
        Do not use: deprecated.

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    warn(DEP_MSG, FutureWarning, stacklevel=2)
    return features.features_from_xml(filepath, polygon, tags)
