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

from . import projection
from . import settings
from . import utils
from . import utils_graph


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
        graph to sample points from; should be undirected (to not oversample
        bidirectional edges) and projected (for accurate point interpolation)
    n : int
        how many points to sample

    Returns
    -------
    points : geopandas.GeoSeries
        the sampled points, multi-indexed by (u, v, key) of the edge from
        which each point was drawn
    """
    if nx.is_directed(G):  # pragma: no cover
        warn("graph should be undirected to not oversample bidirectional edges", stacklevel=2)
    gdf_edges = utils_graph.graph_to_gdfs(G, nodes=False)[["geometry", "length"]]
    weights = gdf_edges["length"] / gdf_edges["length"].sum()
    idx = np.random.choice(gdf_edges.index, size=n, p=weights)
    lines = gdf_edges.loc[idx, "geometry"]
    return lines.interpolate(np.random.rand(n), normalized=True)


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
        smaller values generate more points.

    Yields
    ------
    point : tuple of floats
        a generator of (x, y) tuples of the interpolated points' coordinates
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
        "the `round_geometry_coords` function is deprecated and will be removed in a future release",
        stacklevel=2,
    )

    if isinstance(geom, Point):
        return _round_point_coords(geom, precision)

    elif isinstance(geom, MultiPoint):
        return _round_multipoint_coords(geom, precision)

    elif isinstance(geom, LineString):
        return _round_linestring_coords(geom, precision)

    elif isinstance(geom, MultiLineString):
        return _round_multilinestring_coords(geom, precision)

    elif isinstance(geom, Polygon):
        return _round_polygon_coords(geom, precision)

    elif isinstance(geom, MultiPolygon):
        return _round_multipolygon_coords(geom, precision)

    else:  # pragma: no cover
        msg = f"cannot round coordinates of unhandled geometry type: {type(geom)}"
        raise TypeError(msg)


def _consolidate_subdivide_geometry(geometry, max_query_area_size=None):
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
        the geometry to consolidate and subdivide
    max_query_area_size : int
        maximum area for any part of the geometry in meters: any polygon
        bigger than this will get divided up for multiple queries to API
        (default 50km x 50km). if None, use settings.max_query_area_size

    Returns
    -------
    geometry : shapely.geometry.MultiPolygon
    """
    if max_query_area_size is None:
        max_query_area_size = settings.max_query_area_size

    # let the linear length of the quadrats (with which to subdivide the
    # geometry) be the square root of max area size
    quadrat_width = np.sqrt(max_query_area_size)

    if not isinstance(geometry, (Polygon, MultiPolygon)):  # pragma: no cover
        msg = "Geometry must be a shapely Polygon or MultiPolygon"
        raise TypeError(msg)

    # if geometry is either 1) a Polygon whose area exceeds the max size, or
    # 2) a MultiPolygon, then get the convex hull around the geometry
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
    if not isinstance(geometry, MultiPolygon):  # pragma: no cover
        msg = "Geometry must be a shapely MultiPolygon"
        raise TypeError(msg)

    # extract geometry's exterior coords
    polygons_coords = []
    for polygon in geometry.geoms:
        x, y = polygon.exterior.xy
        polygons_coords.append(list(zip(x, y)))

    # convert exterior coords to the string format the API expects
    polygon_coord_strs = []
    for coords in polygons_coords:
        s = ""
        separator = " "
        for coord in list(coords):
            # round floating point lats and longs to 6 decimals (ie, ~100 mm)
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


def _intersect_index_quadrats(geometries, polygon, quadrat_width=0.05, min_num=3):
    """
    Identify geometries that intersect a (multi)polygon.

    Uses an r-tree spatial index and cuts polygon up into smaller sub-polygons
    for r-tree acceleration. Ensure that geometries and polygon are in the
    same coordinate reference system.

    Parameters
    ----------
    geometries : geopandas.GeoSeries
        the geometries to intersect with the polygon
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        the polygon to intersect with the geometries
    quadrat_width : numeric
        linear length (in polygon's units) of quadrat lines with which to cut
        up the polygon (default = 0.05 degrees, approx 4km at NYC's latitude)
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
    utils.log(f"Created r-tree spatial index for {len(geometries):,} geometries")

    # cut the polygon into chunks for spatial index intersecting
    multipoly = _quadrat_cut_geometry(polygon, quadrat_width=quadrat_width, min_num=min_num)
    geoms_in_poly = set()

    # loop through each chunk of the polygon to find intersecting geometries
    for poly in multipoly.geoms:
        # first find approximate matches with spatial index, then precise
        # matches from those approximate ones
        poly_buff = poly.buffer(0)
        if poly_buff.is_valid and poly_buff.area > 0:
            possible_matches_iloc = sindex.intersection(poly_buff.bounds)
            possible_matches = geometries.iloc[list(possible_matches_iloc)]
            precise_matches = possible_matches[possible_matches.intersects(poly_buff)]
            geoms_in_poly.update(precise_matches.index)

    utils.log(f"Identified {len(geoms_in_poly):,} geometries inside polygon")
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
    earth_radius = 6_371_009  # meters
    lat, lng = point

    delta_lat = (dist / earth_radius) * (180 / np.pi)
    delta_lng = (dist / earth_radius) * (180 / np.pi) / np.cos(lat * np.pi / 180)
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
