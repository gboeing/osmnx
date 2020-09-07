"""Geospatial utility functions."""

import math

import numpy as np
from shapely.geometry import LineString
from shapely.geometry import MultiLineString
from shapely.geometry import MultiPoint
from shapely.geometry import MultiPolygon
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.ops import split

from . import projection
from . import settings
from . import utils


def redistribute_vertices(geom, dist):
    """
    Redistribute the vertices on a projected LineString or MultiLineString.

    The distance argument is only approximate since the total distance of the
    linestring may not be a multiple of the preferred distance. This function
    works on only (Multi)LineString geometry types.

    Parameters
    ----------
    geom : shapely.geometry.LineString or shapely.geometry.MultiLineString
        a Shapely geometry (should be projected)
    dist : float
        spacing length along edges. Units are same as the geom: degrees for
        unprojected geometries and meters for projected geometries. The
        smaller the dist value, the more points are created.

    Returns
    -------
    list or shapely.geometry.MultiLineString
        the redistributed vertices as a list if geom is a LineString or
        MultiLineString if geom is a MultiLineString
    """
    if geom.geom_type == "LineString":
        num_vert = int(round(geom.length / dist))
        if num_vert == 0:
            num_vert = 1
        return [geom.interpolate(float(n) / num_vert, normalized=True) for n in range(num_vert + 1)]
    elif geom.geom_type == "MultiLineString":
        parts = [redistribute_vertices(part, dist) for part in geom]
        return type(geom)([p for p in parts if not p])
    else:
        raise ValueError(f"unhandled geometry {geom.geom_type}")


def _round_polygon_coords(p, precision):
    """
    Round the coordinates of a shapely Polygon to some decimal precision.

    Parameters
    ----------
    p : shapely.geometry.Polygon
        the polygon to round the coordinates of
    precision : int
        decimal precision to round coordinates to

    Returns
    -------
    new_poly : shapely.geometry.Polygon
        the polygon with rounded coordinates
    """
    # round the coordinates of the Polygon exterior
    new_exterior = [[round(x, precision) for x in c] for c in p.exterior.coords]

    # round the coordinates of the (possibly multiple, possibly none) Polygon interior(s)
    new_interiors = []
    for interior in p.interiors:
        new_interiors.append([[round(x, precision) for x in c] for c in interior.coords])

    # construct a new Polygon with the rounded coordinates
    # buffer by zero to clean self-touching or self-crossing polygons
    new_poly = Polygon(shell=new_exterior, holes=new_interiors).buffer(0)
    return new_poly


def _round_multipolygon_coords(mp, precision):
    """
    Round the coordinates of a shapely MultiPolygon to some decimal precision.

    Parameters
    ----------
    mp : shapely.geometry.MultiPolygon
        the MultiPolygon to round the coordinates of
    precision : int
        decimal precision to round coordinates to

    Returns
    -------
    shapely.geometry.MultiPolygon
    """
    return MultiPolygon([_round_polygon_coords(p, precision) for p in mp])


def _round_point_coords(pt, precision):
    """
    Round the coordinates of a shapely Point to some decimal precision.

    Parameters
    ----------
    pt : shapely.geometry.Point
        the Point to round the coordinates of
    precision : int
        decimal precision to round coordinates to

    Returns
    -------
    shapely.geometry.Point
    """
    return Point([round(x, precision) for x in pt.coords[0]])


def _round_multipoint_coords(mpt, precision):
    """
    Round the coordinates of a shapely MultiPoint to some decimal precision.

    Parameters
    ----------
    mpt : shapely.geometry.MultiPoint
        the MultiPoint to round the coordinates of
    precision : int
        decimal precision to round coordinates to

    Returns
    -------
    shapely.geometry.MultiPoint
    """
    return MultiPoint([_round_point_coords(pt, precision) for pt in mpt])


def _round_linestring_coords(ls, precision):
    """
    Round the coordinates of a shapely LineString to some decimal precision.

    Parameters
    ----------
    ls : shapely.geometry.LineString
        the LineString to round the coordinates of
    precision : int
        decimal precision to round coordinates to

    Returns
    -------
    shapely.geometry.LineString
    """
    return LineString([[round(x, precision) for x in c] for c in ls.coords])


def _round_multilinestring_coords(mls, precision):
    """
    Round the coordinates of a shapely MultiLineString to some decimal precision.

    Parameters
    ----------
    mls : shapely.geometry.MultiLineString
        the MultiLineString to round the coordinates of
    precision : int
        decimal precision to round coordinates to

    Returns
    -------
    shapely.geometry.MultiLineString
    """
    return MultiLineString([_round_linestring_coords(ls, precision) for ls in mls])


def round_geometry_coords(shape, precision):
    """
    Round the coordinates of a shapely geometry to some decimal precision.

    Parameters
    ----------
    shape : shapely.geometry.geometry
        the geometry to round the coordinates of; must be one of {Point,
        MultiPoint, LineString, MultiLineString, Polygon, MultiPolygon}
    precision : int
        decimal precision to round coordinates to

    Returns
    -------
    shapely.geometry.geometry
    """
    if isinstance(shape, Point):
        return _round_point_coords(shape, precision)

    elif isinstance(shape, MultiPoint):
        return _round_multipoint_coords(shape, precision)

    elif isinstance(shape, LineString):
        return _round_linestring_coords(shape, precision)

    elif isinstance(shape, MultiLineString):
        return _round_multilinestring_coords(shape, precision)

    elif isinstance(shape, Polygon):
        return _round_polygon_coords(shape, precision)

    elif isinstance(shape, MultiPolygon):
        return _round_multipolygon_coords(shape, precision)

    else:
        raise TypeError(f"cannot round coordinates of unhandled geometry type: {type(shape)}")


def _consolidate_subdivide_geometry(geometry, max_query_area_size=None):
    """
    Consolidate and subdivide some geometry.

    Consolidate a geometry into a convex hull, then subdivide it into smaller
    sub-polygons if its area exceeds max size (in geometry's units). Configure
    the max size via max_query_area_size in the settings module.

    Parameters
    ----------
    geometry : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        the geometry to consolidate and subdivide
    max_query_area_size : int
        maximum area for any part of the geometry in meters: any polygon
        bigger than this will get divided up for multiple queries to API
        (default 50km x 50km). if None, use settings.max_query_area_size

    Returns
    -------
    geometry : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
    """
    if max_query_area_size is None:
        max_query_area_size = settings.max_query_area_size

    # let the linear length of the quadrats (with which to subdivide the
    # geometry) be the square root of max area size
    quadrat_width = math.sqrt(max_query_area_size)

    if not isinstance(geometry, (Polygon, MultiPolygon)):
        raise TypeError("Geometry must be a shapely Polygon or MultiPolygon")

    # if geometry is a MultiPolygon OR a single Polygon whose area exceeds the
    # max size, get the convex hull around the geometry
    if isinstance(geometry, MultiPolygon) or (
        isinstance(geometry, Polygon) and geometry.area > max_query_area_size
    ):
        geometry = geometry.convex_hull

    # if geometry area exceeds max size, subdivide it into smaller sub-polygons
    if geometry.area > max_query_area_size:
        geometry = _quadrat_cut_geometry(geometry, quadrat_width=quadrat_width)

    if isinstance(geometry, Polygon):
        geometry = MultiPolygon([geometry])

    return geometry


def _get_polygons_coordinates(geometry):
    """
    Extract exterior coordinates from polygon(s) to pass to OSM.

    Ignore the interior ("holes") coordinates.

    Parameters
    ----------
    geometry : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        the geometry to extract exterior coordinates from

    Returns
    -------
    polygon_coord_strs : list
    """
    # extract the exterior coordinates of the geometry to pass to the API later
    polygons_coords = []
    if isinstance(geometry, Polygon):
        x, y = geometry.exterior.xy
        polygons_coords.append(list(zip(x, y)))
    elif isinstance(geometry, MultiPolygon):
        for polygon in geometry:
            x, y = polygon.exterior.xy
            polygons_coords.append(list(zip(x, y)))
    else:
        raise TypeError("Geometry must be a shapely Polygon or MultiPolygon")

    # convert the exterior coordinates of the polygon(s) to the string format
    # the API expects
    polygon_coord_strs = []
    for coords in polygons_coords:
        s = ""
        separator = " "
        for coord in list(coords):
            # round floating point lats and longs to 6 decimal places (ie, ~100 mm),
            # so we can hash and cache strings consistently
            s = f"{s}{separator}{coord[1]:.6f}{separator}{coord[0]:.6f}"
        polygon_coord_strs.append(s.strip(separator))

    return polygon_coord_strs


def _quadrat_cut_geometry(geometry, quadrat_width, min_num=3):
    """
    Split a Polygon or MultiPolygon up into sub-polygons of a specified size.

    Parameters
    ----------
    geometry : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        the geometry to split up into smaller sub-polygons
    quadrat_width : numeric
        the linear width of the quadrats with which to cut up the geometry (in
        the units the geometry is in)
    min_num : int
        the minimum number of linear quadrat lines (e.g., min_num=3 would
        produce a quadrat grid of 4 squares)

    Returns
    -------
    geometry : shapely.geometry.MultiPolygon
    """
    # create n evenly spaced points between the min and max x and y bounds
    west, south, east, north = geometry.bounds
    x_num = math.ceil((east - west) / quadrat_width) + 1
    y_num = math.ceil((north - south) / quadrat_width) + 1
    x_points = np.linspace(west, east, num=max(x_num, min_num))
    y_points = np.linspace(south, north, num=max(y_num, min_num))

    # create a quadrat grid of lines at each of the evenly spaced points
    vertical_lines = [LineString([(x, y_points[0]), (x, y_points[-1])]) for x in x_points]
    horizont_lines = [LineString([(x_points[0], y), (x_points[-1], y)]) for y in y_points]
    lines = vertical_lines + horizont_lines

    # recursively split the geometry by each quadrat line
    for line in lines:
        geometry = MultiPolygon(split(geometry, line))

    return geometry


def _intersect_index_quadrats(geometries, polygon, quadrat_width=0.05, min_num=3):
    """
    Identify geometries that intersect a (multi)polygon.

    Use an r-tree spatial index and cut the polygon up into smaller
    sub-polygons for r-tree acceleration.

    Parameters
    ----------
    geometries : geopandas.GeoSeries
        the geometries to intersect with the polygon
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        the polygon to intersect with the geometries
    quadrat_width : numeric
        the linear length (in units the polygon is in) of the quadrats with
        which to cut up the polygon (default = 0.05 degrees, approx 4km at
        NYC's latitude)
    min_num : int
        the minimum number of linear quadrat lines (e.g., min_num=3 would
        produce a quadrat grid of 4 squares)

    Returns
    -------
    geoms_in_poly : set
        index labels of geometries that intersected polygon
    """
    # create an r-tree spatial index for the geometries
    sindex = geometries.sindex
    utils.log(f"Created r-tree spatial index for {len(geometries)} geometries")

    # cut the polygon into chunks for spatial index intersecting
    multipoly = _quadrat_cut_geometry(polygon, quadrat_width=quadrat_width, min_num=min_num)
    geoms_in_poly = set()

    # loop through each chunk of the polygon to find intersecting geometries
    for poly in multipoly:
        # first find approximate matches with spatial index, then precise
        # matches from those approximate ones
        poly = poly.buffer(0)
        if poly.is_valid and poly.area > 0:
            possible_matches_iloc = sindex.intersection(poly.bounds)
            possible_matches = geometries.iloc[list(possible_matches_iloc)]
            precise_matches = possible_matches[possible_matches.intersects(poly)]
            geoms_in_poly.update(precise_matches.index)

    utils.log(f"Identified {len(geoms_in_poly)} geometries inside polygon")
    return geoms_in_poly


def bbox_from_point(point, dist=1000, project_utm=False, return_crs=False):
    """
    Create a bounding box from a (lat, lng) center point.

    Create a bounding box some distance in each direction (north, south, east,
    and west) from the center point and optionally project it.

    Parameters
    ----------
    point : tuple
        the (lat, lng) center point to create the bounding box around
    dist : int
        bounding box distance in meters from the center point
    project_utm : bool
        if True, return bounding box as UTM-projected coordinates
    return_crs : bool
        if True, and project_utm=True, return the projected CRS too

    Returns
    -------
    tuple
        (north, south, east, west) or (north, south, east, west, crs_proj)
    """
    earth_radius = 6371009  # meters
    lat, lng = point

    delta_lat = (dist / earth_radius) * (180 / math.pi)
    delta_lng = (dist / earth_radius) * (180 / math.pi) / math.cos(lat * math.pi / 180)
    north = lat + delta_lat
    south = lat - delta_lat
    east = lng + delta_lng
    west = lng - delta_lng

    if project_utm:
        bbox_poly = bbox_to_poly(north, south, east, west)
        bbox_proj, crs_proj = projection.project_geometry(bbox_poly)
        west, south, east, north = bbox_proj.bounds

    utils.log(f"Created bbox {dist} m from {point}: {north},{south},{east},{west}")

    if project_utm and return_crs:
        return north, south, east, west, crs_proj
    else:
        return north, south, east, west


def bbox_to_poly(north, south, east, west):
    """
    Convert bounding box coordinates to shapely Polygon.

    Parameters
    ----------
    north : float
        northern coordinate
    south : float
        southern coordinate
    east : float
        eastern coordinate
    west : float
        western coordinate

    Returns
    -------
    shapely.geometry.Polygon
    """
    return Polygon([(west, south), (east, south), (east, north), (west, north)])
