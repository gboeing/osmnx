"""Geospatial utility functions."""

from __future__ import annotations

import logging as lg
from typing import TYPE_CHECKING
from typing import Any
from typing import Literal
from typing import overload
from warnings import warn

import networkx as nx
import numpy as np
from shapely import Geometry
from shapely import LineString
from shapely import MultiPolygon
from shapely import Polygon
from shapely.ops import split

from . import convert
from . import projection
from . import settings
from . import utils

if TYPE_CHECKING:
    from collections.abc import Iterator

    import geopandas as gpd


def buffer_geometry(geom: Geometry, dist: float) -> Geometry:
    """
    Buffer an unprojected Shapely geometry by some distance in meters.

    Parameters
    ----------
    geom
        The geometry to be buffered. Coordinates should be in unprojected
        latitude-longitude degrees (EPSG:4326).
    dist
        The buffer distance in meters.

    Returns
    -------
    geometry_buff
        The (also unprojected) buffered geometry.
    """
    geom_proj, crs_proj = projection.project_geometry(geom)
    geom_buff, _ = projection.project_geometry(
        geometry=geom_proj.buffer(dist),
        crs=crs_proj,
        to_latlong=True,
    )
    return geom_buff


def sample_points(G: nx.MultiGraph, n: int) -> gpd.GeoSeries:
    """
    Randomly sample points constrained to a spatial graph.

    This generates a graph-constrained uniform random sample of points. Unlike
    typical spatially uniform random sampling, this method accounts for the
    graph's geometry. And unlike equal-length edge segmenting, this method
    guarantees uniform randomness.

    Parameters
    ----------
    G
        Graph from which to sample points. Should be undirected (to avoid
        oversampling bidirectional edges) and projected (for accurate point
        interpolation).
    n
        How many points to sample.

    Returns
    -------
    points
        The sampled points, multi-indexed by `(u, v, key)` of the edge from
        which each point was sampled.
    """
    if nx.is_directed(G):  # pragma: no cover
        msg = "`G` should be undirected to avoid oversampling bidirectional edges."
        warn(msg, category=UserWarning, stacklevel=2)
    gdf_edges = convert.graph_to_gdfs(G, nodes=False)[["geometry", "length"]]
    weights = gdf_edges["length"] / gdf_edges["length"].sum()
    idx = np.random.default_rng().choice(gdf_edges.index, size=n, p=weights)
    lines = gdf_edges.loc[idx, "geometry"]
    return lines.interpolate(np.random.default_rng().random(n), normalized=True)


def interpolate_points(geom: LineString, dist: float) -> Iterator[tuple[float, float]]:
    """
    Interpolate evenly spaced points along a LineString.

    The spacing is approximate because the LineString's length may not be
    evenly divisible by it.

    Parameters
    ----------
    geom
        A LineString geometry.
    dist
        Spacing distance between interpolated points, in same units as `geom`.
        Smaller values accordingly generate more points.

    Yields
    ------
    point
        Interpolated point's `(x, y)` coordinates.
    """
    if isinstance(geom, LineString):
        num_vert = max(round(geom.length / dist), 1)
        for n in range(num_vert + 1):
            point = geom.interpolate(n / num_vert, normalized=True)
            yield point.x, point.y
    else:  # pragma: no cover
        msg = "`geom` must be a LineString."
        raise TypeError(msg)


def _consolidate_subdivide_geometry(geom: Polygon | MultiPolygon) -> MultiPolygon:
    """
    Consolidate and subdivide some (projected) geometry.

    Consolidate a geometry into a convex hull, then subdivide it into smaller
    sub-polygons if its area exceeds max size (in geometry's units). Configure
    the max size via the `settings` module's `max_query_area_size`. Geometries
    with areas much larger than `max_query_area_size` may take a long time to
    process.

    When the geometry has a very large area relative to its vertex count,
    the resulting MultiPolygon's boundary may differ somewhat from the input,
    due to the way long straight lines are projected. You can interpolate
    additional vertices along your input geometry's exterior to mitigate this
    if necessary.

    Parameters
    ----------
    geom
        The projected (in meter units) geometry to consolidate and subdivide.

    Returns
    -------
    geom
    """
    if not isinstance(geom, (Polygon, MultiPolygon)):  # pragma: no cover
        msg = "Geometry must be a shapely Polygon or MultiPolygon."
        raise TypeError(msg)

    # if geometry is either 1) a Polygon whose area exceeds the max size, or
    # 2) a MultiPolygon, then get the convex hull around the geometry
    mqas = settings.max_query_area_size
    if isinstance(geom, MultiPolygon) or (isinstance(geom, Polygon) and geom.area > mqas):
        geom = geom.convex_hull

    # warn user if they passed a geometry with area much larger than max size
    ratio = int(geom.area / mqas)
    warning_threshold = 10
    if ratio > warning_threshold:
        msg = (
            f"This area is {ratio:,} times your configured Overpass max query "
            "area size. It will automatically be divided up into multiple "
            "sub-queries accordingly. This may take a long time."
        )
        warn(msg, category=UserWarning, stacklevel=2)

    # if geometry area exceeds max size, subdivide it into smaller subpolygons
    # that are no greater than settings.max_query_area_size in size
    if geom.area > mqas:
        geom = _quadrat_cut_geometry(geom, quadrat_width=np.sqrt(mqas))

    if isinstance(geom, Polygon):
        geom = MultiPolygon([geom])

    return geom


def _quadrat_cut_geometry(geom: Polygon | MultiPolygon, quadrat_width: float) -> MultiPolygon:
    """
    Split a Polygon or MultiPolygon up into sub-polygons of a specified size.

    Parameters
    ----------
    geom
        The geometry to split up into smaller sub-polygons.
    quadrat_width
        Width (in geometry's units) of quadrat squares with which to split up
        the geometry.

    Returns
    -------
    geom
    """
    # min number of dividing lines (3 produces a grid of 4 quadrat squares)
    min_num = 3

    # create n evenly spaced points between the min and max x and y bounds
    left, bottom, right, top = geom.bounds
    x_num = int(np.ceil((right - left) / quadrat_width) + 1)
    y_num = int(np.ceil((top - bottom) / quadrat_width) + 1)
    x_points = np.linspace(left, right, num=max(x_num, min_num))
    y_points = np.linspace(bottom, top, num=max(y_num, min_num))

    # create a quadrat grid of lines at each of the evenly spaced points
    vertical_lines = [LineString([(x, y_points[0]), (x, y_points[-1])]) for x in x_points]
    horizont_lines = [LineString([(x_points[0], y), (x_points[-1], y)]) for y in y_points]
    lines = vertical_lines + horizont_lines

    # recursively split the geometry by each quadrat line
    geoms = [geom]
    for line in lines:
        # split polygon by line if they intersect, otherwise just keep it
        split_geoms = [split(g, line).geoms if g.intersects(line) else [g] for g in geoms]
        # now flatten the list and process these split geoms on the next line in the list of lines
        geoms = [g for g_list in split_geoms for g in g_list]

    return MultiPolygon(geoms)


def _intersect_index_quadrats(
    geoms: gpd.GeoSeries,
    polygon: Polygon | MultiPolygon,
) -> set[Any]:
    """
    Identify geometries that intersect a (Multi)Polygon.

    Uses an r-tree spatial index and cuts polygon up into smaller sub-polygons
    for r-tree acceleration. Ensure that geometries and polygon are in the
    same coordinate reference system.

    Parameters
    ----------
    geoms
        The geometries to intersect with the polygon.
    polygon
        The polygon to intersect with the geometries.

    Returns
    -------
    geoms_in_poly
        The index labels of the geometries that intersected the polygon.
    """
    # create an r-tree spatial index for the geometries
    rtree = geoms.sindex
    msg = f"Built r-tree spatial index for {len(geoms):,} geometries"
    utils.log(msg, level=lg.INFO)

    # cut polygon into chunks for faster spatial index intersecting. specify a
    # sensible quadrat_width to balance performance (eg, 0.1 degrees is approx
    # 8 km at NYC's latitude) with either projected or unprojected coordinates
    quadrat_width = max(0.1, np.sqrt(polygon.area) / 10)
    multipoly = _quadrat_cut_geometry(polygon, quadrat_width)
    msg = f"Accelerating r-tree with {len(multipoly.geoms)} quadrats"
    utils.log(msg, level=lg.INFO)

    # loop through each chunk of the polygon to find intersecting geometries
    # first find approximate matches with spatial index, then precise matches
    # from those approximate ones
    geoms_in_poly = set()
    for poly in multipoly.geoms:
        poly_buff = poly.buffer(0)
        if poly_buff.is_valid and poly_buff.area > 0:
            possible_matches_iloc = rtree.intersection(poly_buff.bounds)
            possible_matches = geoms.iloc[list(possible_matches_iloc)]
            precise_matches = possible_matches[possible_matches.intersects(poly_buff)]
            geoms_in_poly.update(precise_matches.index)

    msg = f"Identified {len(geoms_in_poly):,} geometries inside polygon"
    utils.log(msg, level=lg.INFO)
    return geoms_in_poly


# dist present, project_utm missing/False, return_crs missing/False
@overload
def bbox_from_point(
    point: tuple[float, float],
    dist: float,
) -> tuple[float, float, float, float]: ...


# dist present, project_utm missing/False, return_crs present/True
@overload
def bbox_from_point(
    point: tuple[float, float],
    dist: float,
    *,
    return_crs: Literal[True],
) -> tuple[float, float, float, float]: ...


# dist present, project_utm missing/False, return_crs present/False
@overload
def bbox_from_point(
    point: tuple[float, float],
    dist: float,
    *,
    return_crs: Literal[False],
) -> tuple[float, float, float, float]: ...


# dist present, project_utm present/True, return_crs missing/False
@overload
def bbox_from_point(
    point: tuple[float, float],
    dist: float,
    *,
    project_utm: Literal[True],
) -> tuple[float, float, float, float]: ...


# dist present, project_utm present/True, return_crs present/True
@overload
def bbox_from_point(
    point: tuple[float, float],
    dist: float,
    *,
    project_utm: Literal[True],
    return_crs: Literal[True],
) -> tuple[tuple[float, float, float, float], Any]: ...


# dist present, project_utm present/True, return_crs present/False
@overload
def bbox_from_point(
    point: tuple[float, float],
    dist: float,
    *,
    project_utm: Literal[True],
    return_crs: Literal[False],
) -> tuple[float, float, float, float]: ...


# dist present, project_utm present/False, return_crs missing/False
@overload
def bbox_from_point(
    point: tuple[float, float],
    dist: float,
    *,
    project_utm: Literal[False],
) -> tuple[float, float, float, float]: ...


# dist present, project_utm present/False, return_crs present/True
@overload
def bbox_from_point(
    point: tuple[float, float],
    dist: float,
    *,
    project_utm: Literal[False],
    return_crs: Literal[True],
) -> tuple[float, float, float, float]: ...


# dist present, project_utm present/False, return_crs present/False
@overload
def bbox_from_point(
    point: tuple[float, float],
    dist: float,
    *,
    project_utm: Literal[False],
    return_crs: Literal[False],
) -> tuple[float, float, float, float]: ...


def bbox_from_point(
    point: tuple[float, float],
    dist: float,
    *,
    project_utm: bool = False,
    return_crs: bool = False,
) -> tuple[float, float, float, float] | tuple[tuple[float, float, float, float], Any]:
    """
    Create a bounding box around a (lat, lon) point.

    Create a bounding box some distance (in meters) in each direction (top,
    bottom, right, and left) from the center point and optionally project it.

    Parameters
    ----------
    point
        The `(lat, lon)` center point to create the bounding box around.
    dist
        Bounding box distance in meters from the center point.
    project_utm
        If True, return bounding box as UTM-projected coordinates.
    return_crs
        If True, and `project_utm` is True, then return the projected CRS too.

    Returns
    -------
    bbox or bbox, crs
        `(left, bottom, right, top)` or `((left, bottom, right, top), crs)`.
    """
    EARTH_RADIUS_M = 6_371_009  # meters
    lat, lon = point

    delta_lat = (dist / EARTH_RADIUS_M) * (180 / np.pi)
    delta_lon = (dist / EARTH_RADIUS_M) * (180 / np.pi) / np.cos(lat * np.pi / 180)
    top = lat + delta_lat
    bottom = lat - delta_lat
    right = lon + delta_lon
    left = lon - delta_lon
    bbox = left, bottom, right, top

    if project_utm:
        bbox_poly = bbox_to_poly(bbox=bbox)
        bbox_proj, crs_proj = projection.project_geometry(bbox_poly)
        bbox = bbox_proj.bounds

    msg = f"Created bbox {dist} meters from {point}: {bbox}"
    utils.log(msg, level=lg.INFO)

    if project_utm and return_crs:
        return bbox, crs_proj

    # otherwise
    return bbox


def bbox_to_poly(bbox: tuple[float, float, float, float]) -> Polygon:
    """
    Convert bounding box coordinates to Shapely Polygon.

    Parameters
    ----------
    bbox
        Bounding box as `(left, bottom, right, top)`.

    Returns
    -------
    polygon
    """
    left, bottom, right, top = bbox
    return Polygon([(left, bottom), (right, bottom), (right, top), (left, top)])
