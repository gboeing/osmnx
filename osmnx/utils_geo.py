"""Geospatial utility functions."""

from warnings import warn

import networkx as nx
import numpy as np
from shapely.geometry import LineString
from shapely.geometry import MultiLineString
from shapely.geometry import MultiPoint
from shapely.geometry import MultiPolygon
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.ops import split

from . import convert
from . import projection
from . import settings
from . import utils


def sample_points(G, n):
    """
    Randomly sample points constrained to a spatial graph.

    This generates a graph-constrained uniform random sample of points. Unlike
    typical spatially uniform random sampling, this method accounts for the
    graph's geometry. And unlike equal-length edge segmenting, this method
    guarantees uniform randomness.

    Parameters
    ----------
    G : networkx.MultiGraph
        graph from which to sample points. should be undirected (to avoid
        oversampling bidirectional edges) and projected (for accurate point
        interpolation)
    n : int
        how many points to sample

    Returns
    -------
    points : geopandas.GeoSeries
        the sampled points, multi-indexed by (u, v, key) of the edge from
        which each point was drawn
    """
    if nx.is_directed(G):  # pragma: no cover
        warn("graph should be undirected to avoid oversampling bidirectional edges", stacklevel=2)
    gdf_edges = convert.graph_to_gdfs(G, nodes=False)[["geometry", "length"]]
    weights = gdf_edges["length"] / gdf_edges["length"].sum()
    idx = np.random.default_rng().choice(gdf_edges.index, size=n, p=weights)
    lines = gdf_edges.loc[idx, "geometry"]
    return lines.interpolate(np.random.default_rng().random(n), normalized=True)


def interpolate_points(geom, dist):
    """
    Interpolate evenly spaced points along a LineString.

    The spacing is approximate because the LineString's length may not be
    evenly divisible by it.

    Parameters
    ----------
    geom : shapely.geometry.LineString
        a LineString geometry
    dist : float
        spacing distance between interpolated points, in same units as `geom`.
        smaller values accordingly generate more points.

    Yields
    ------
    points : generator
        tuples of (x, y) floats of the interpolated points' coordinates
    """
    if isinstance(geom, LineString):
        num_vert = max(round(geom.length / dist), 1)
        for n in range(num_vert + 1):
            point = geom.interpolate(n / num_vert, normalized=True)
            yield point.x, point.y
    else:  # pragma: no cover
        msg = f"unhandled geometry type {geom.geom_type}"
        raise TypeError(msg)


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
    shapely.geometry.Polygon
    """
    # round coords of Polygon exterior
    shell = [[round(x, precision) for x in c] for c in p.exterior.coords]

    # round coords of (possibly multiple, possibly none) Polygon interior(s)
    holes = [[[round(x, precision) for x in c] for c in i.coords] for i in p.interiors]

    # construct new Polygon with rounded coordinates and buffer by zero to
    # clean self-touching or self-crossing polygons
    return Polygon(shell=shell, holes=holes).buffer(0)


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
    return MultiPolygon([_round_polygon_coords(p, precision) for p in mp.geoms])


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
    return MultiPoint([_round_point_coords(pt, precision) for pt in mpt.geoms])


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
    return MultiLineString([_round_linestring_coords(ls, precision) for ls in mls.geoms])


def round_geometry_coords(geom, precision):
    """
    Do not use: deprecated.

    Parameters
    ----------
    geom : shapely.geometry.geometry {Point, MultiPoint, LineString, MultiLineString, Polygon, MultiPolygon}
        deprecated, do not use
    precision : int
        deprecated, do not use

    Returns
    -------
    shapely.geometry.geometry
    """
    warn(
        "The `round_geometry_coords` function is deprecated and will be "
        "removed in the v2.0.0 release. "
        "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123",
        FutureWarning,
        stacklevel=2,
    )

    if isinstance(geom, Point):
        return _round_point_coords(geom, precision)

    if isinstance(geom, MultiPoint):
        return _round_multipoint_coords(geom, precision)

    if isinstance(geom, LineString):
        return _round_linestring_coords(geom, precision)

    if isinstance(geom, MultiLineString):
        return _round_multilinestring_coords(geom, precision)

    if isinstance(geom, Polygon):
        return _round_polygon_coords(geom, precision)

    if isinstance(geom, MultiPolygon):
        return _round_multipolygon_coords(geom, precision)

    # otherwise
    msg = f"cannot round coordinates of unhandled geometry type: {type(geom)}"
    raise TypeError(msg)


def _consolidate_subdivide_geometry(geometry):
    """
    Consolidate and subdivide some geometry.

    Consolidate a geometry into a convex hull, then subdivide it into smaller
    sub-polygons if its area exceeds max size (in geometry's units). Configure
    the max size via max_query_area_size in the settings module.

    When the geometry has a very large area relative to its vertex count,
    the resulting MultiPolygon's boundary may differ somewhat from the input,
    due to the way long straight lines are projected. You can interpolate
    additional vertices along your input geometry's exterior to mitigate this.

    Parameters
    ----------
    geometry : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        the projected (in meter units) geometry to consolidate and subdivide

    Returns
    -------
    geometry : shapely.geometry.MultiPolygon
    """
    if not isinstance(geometry, (Polygon, MultiPolygon)):  # pragma: no cover
        msg = "Geometry must be a shapely Polygon or MultiPolygon"
        raise TypeError(msg)

    # if geometry is either 1) a Polygon whose area exceeds the max size, or
    # 2) a MultiPolygon, then get the convex hull around the geometry
    mqas = settings.max_query_area_size
    if isinstance(geometry, MultiPolygon) or (
        isinstance(geometry, Polygon) and geometry.area > mqas
    ):
        geometry = geometry.convex_hull

    # warn user if they passed a geometry with area much larger than max size
    ratio = int(geometry.area / mqas)
    warning_threshold = 10
    if ratio > warning_threshold:
        msg = (
            f"This area is {ratio:,} times your configured Overpass max query "
            "area size. It will automatically be divided up into multiple "
            "sub-queries accordingly. This may take a long time."
        )
        warn(msg, stacklevel=2)

    # if geometry area exceeds max size, subdivide it into smaller subpolygons
    # that are no greater than settings.max_query_area_size in size
    if geometry.area > mqas:
        geometry = _quadrat_cut_geometry(geometry, quadrat_width=np.sqrt(mqas))

    if isinstance(geometry, Polygon):
        geometry = MultiPolygon([geometry])

    return geometry


def _quadrat_cut_geometry(geometry, quadrat_width):
    """
    Split a Polygon or MultiPolygon up into sub-polygons of a specified size.

    Parameters
    ----------
    geometry : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        the geometry to split up into smaller sub-polygons
    quadrat_width : float
        width (in geometry's units) of quadrat squares with which to split up
        the geometry

    Returns
    -------
    geometry : shapely.geometry.MultiPolygon
    """
    # min number of dividing lines (3 produces a grid of 4 quadrat squares)
    min_num = 3

    # create n evenly spaced points between the min and max x and y bounds
    west, south, east, north = geometry.bounds
    x_num = int(np.ceil((east - west) / quadrat_width) + 1)
    y_num = int(np.ceil((north - south) / quadrat_width) + 1)
    x_points = np.linspace(west, east, num=max(x_num, min_num))
    y_points = np.linspace(south, north, num=max(y_num, min_num))

    # create a quadrat grid of lines at each of the evenly spaced points
    vertical_lines = [LineString([(x, y_points[0]), (x, y_points[-1])]) for x in x_points]
    horizont_lines = [LineString([(x_points[0], y), (x_points[-1], y)]) for y in y_points]
    lines = vertical_lines + horizont_lines

    # recursively split the geometry by each quadrat line
    geometries = [geometry]
    for line in lines:
        # split polygon by line if they intersect, otherwise just keep it
        split_geoms = [split(g, line).geoms if g.intersects(line) else [g] for g in geometries]
        # now flatten the list and process these split geoms on the next line in the list of lines
        geometries = [g for g_list in split_geoms for g in g_list]

    return MultiPolygon(geometries)


def _intersect_index_quadrats(geometries, polygon):
    """
    Identify geometries that intersect a (Multi)Polygon.

    Uses an r-tree spatial index and cuts polygon up into smaller sub-polygons
    for r-tree acceleration. Ensure that geometries and polygon are in the
    same coordinate reference system.

    Parameters
    ----------
    geometries : geopandas.GeoSeries
        the geometries to intersect with the polygon
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        the polygon to intersect with the geometries

    Returns
    -------
    geoms_in_poly : set
        set of the index labels of the geometries that intersected the polygon
    """
    # create an r-tree spatial index for the geometries
    rtree = geometries.sindex
    utils.log(f"Built r-tree spatial index for {len(geometries):,} geometries")

    # cut polygon into chunks for faster spatial index intersecting. specify a
    # sensible quadrat_width to balance performance (eg, 0.1 degrees is approx
    # 8 km at NYC's latitude) with either projected or unprojected coordinates
    quadrat_width = max(0.1, np.sqrt(polygon.area) / 10)
    multipoly = _quadrat_cut_geometry(polygon, quadrat_width)
    utils.log(f"Accelerating r-tree with {len(multipoly.geoms)} quadrats")

    # loop through each chunk of the polygon to find intersecting geometries
    # first find approximate matches with spatial index, then precise matches
    # from those approximate ones
    geoms_in_poly = set()
    for poly in multipoly.geoms:
        poly_buff = poly.buffer(0)
        if poly_buff.is_valid and poly_buff.area > 0:
            possible_matches_iloc = rtree.intersection(poly_buff.bounds)
            possible_matches = geometries.iloc[list(possible_matches_iloc)]
            precise_matches = possible_matches[possible_matches.intersects(poly_buff)]
            geoms_in_poly.update(precise_matches.index)

    utils.log(f"Identified {len(geoms_in_poly):,} geometries inside polygon")
    return geoms_in_poly


def bbox_from_point(point, dist=1000, project_utm=False, return_crs=False):
    """
    Create a bounding box around a (lat, lon) point.

    Create a bounding box some distance (in meters) in each direction (north,
    south, east, and west) from the center point and optionally project it.

    Parameters
    ----------
    point : tuple
        the (lat, lon) center point to create the bounding box around
    dist : int
        bounding box distance in meters from the center point
    project_utm : bool
        if True, return bounding box as UTM-projected coordinates
    return_crs : bool
        if True, and project_utm=True, return the projected CRS too

    Returns
    -------
    bbox or bbox, crs: tuple or tuple, crs
        (north, south, east, west) or ((north, south, east, west), crs)
    """
    EARTH_RADIUS_M = 6_371_009  # meters
    lat, lon = point

    delta_lat = (dist / EARTH_RADIUS_M) * (180 / np.pi)
    delta_lon = (dist / EARTH_RADIUS_M) * (180 / np.pi) / np.cos(lat * np.pi / 180)
    north = lat + delta_lat
    south = lat - delta_lat
    east = lon + delta_lon
    west = lon - delta_lon

    if project_utm:
        bbox_poly = bbox_to_poly(bbox=(north, south, east, west))
        bbox_proj, crs_proj = projection.project_geometry(bbox_poly)
        west, south, east, north = bbox_proj.bounds

    utils.log(f"Created bbox {dist} m from {point}: {north},{south},{east},{west}")

    if project_utm and return_crs:
        return (north, south, east, west), crs_proj

    # otherwise
    return north, south, east, west


def bbox_to_poly(north=None, south=None, east=None, west=None, bbox=None):
    """
    Convert bounding box coordinates to shapely Polygon.

    Parameters
    ----------
    north : float
        deprecated, do not use
    south : float
        deprecated, do not use
    east : float
        deprecated, do not use
    west : float
        deprecated, do not use
    bbox : tuple of floats
        bounding box as (north, south, east, west)

    Returns
    -------
    shapely.geometry.Polygon
    """
    if not (north is None and south is None and east is None and west is None):
        msg = (
            "The `north`, `south`, `east`, and `west` parameters are deprecated and "
            "will be removed in the v2.0.0 release. Use the `bbox` parameter instead. "
            "Note that the expected order of coordinates in `bbox` will change in the "
            "v2.0.0 release to `(left, bottom, right, top)`. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
        )
        warn(msg, FutureWarning, stacklevel=2)
    else:
        north, south, east, west = bbox
    return Polygon([(west, south), (east, south), (east, north), (west, north)])
