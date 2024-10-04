"""
Download and create GeoDataFrames from OpenStreetMap geospatial features.

Retrieve points of interest, building footprints, transit lines/stops, or any
other map features from OSM, including their geometries and attribute data,
then construct a GeoDataFrame of them. You can use this module to query for
nodes, ways, and relations (the latter of type "multipolygon" or "boundary"
only) by passing a dictionary of desired OSM tags.

For more details, see https://wiki.openstreetmap.org/wiki/Map_features and
https://wiki.openstreetmap.org/wiki/Elements

Refer to the Getting Started guide for usage limitations.
"""

from __future__ import annotations

import logging as lg
from typing import TYPE_CHECKING
from typing import Any

import geopandas as gpd
import pandas as pd
from shapely import LineString
from shapely import MultiLineString
from shapely import MultiPolygon
from shapely import Point
from shapely import Polygon
from shapely import prepare
from shapely.errors import GEOSException
from shapely.ops import linemerge
from shapely.ops import polygonize
from shapely.ops import unary_union

from . import _osm_xml
from . import _overpass
from . import geocoder
from . import settings
from . import utils
from . import utils_geo
from ._errors import CacheOnlyInterruptError
from ._errors import InsufficientResponseError

if TYPE_CHECKING:
    from collections.abc import Iterable
    from pathlib import Path

# define what types of OSM relations we currently handle
_RELATION_TYPES = {"boundary", "multipolygon"}

# OSM tags to determine if closed ways should be polygons, based on JSON from
# https://wiki.openstreetmap.org/wiki/Overpass_turbo/Polygon_Features
_POLYGON_FEATURES: dict[str, dict[str, str | set[str]]] = {
    "aeroway": {"polygon": "blocklist", "values": {"taxiway"}},
    "amenity": {"polygon": "all"},
    "area": {"polygon": "all"},
    "area:highway": {"polygon": "all"},
    "barrier": {
        "polygon": "passlist",
        "values": {"city_wall", "ditch", "hedge", "retaining_wall", "spikes"},
    },
    "boundary": {"polygon": "all"},
    "building": {"polygon": "all"},
    "building:part": {"polygon": "all"},
    "craft": {"polygon": "all"},
    "golf": {"polygon": "all"},
    "highway": {"polygon": "passlist", "values": {"elevator", "escape", "rest_area", "services"}},
    "historic": {"polygon": "all"},
    "indoor": {"polygon": "all"},
    "landuse": {"polygon": "all"},
    "leisure": {"polygon": "all"},
    "man_made": {"polygon": "blocklist", "values": {"cutline", "embankment", "pipeline"}},
    "military": {"polygon": "all"},
    "natural": {
        "polygon": "blocklist",
        "values": {"arete", "cliff", "coastline", "ridge", "tree_row"},
    },
    "office": {"polygon": "all"},
    "place": {"polygon": "all"},
    "power": {"polygon": "passlist", "values": {"generator", "plant", "substation", "transformer"}},
    "public_transport": {"polygon": "all"},
    "railway": {
        "polygon": "passlist",
        "values": {"platform", "roundhouse", "station", "turntable"},
    },
    "ruins": {"polygon": "all"},
    "shop": {"polygon": "all"},
    "tourism": {"polygon": "all"},
    "waterway": {"polygon": "passlist", "values": {"boatyard", "dam", "dock", "riverbank"}},
}


def features_from_bbox(
    bbox: tuple[float, float, float, float],
    tags: dict[str, bool | str | list[str]],
) -> gpd.GeoDataFrame:
    """
    Download OSM features within a lat-lon bounding box.

    You can use the `settings` module to retrieve a snapshot of historical OSM
    data as of a certain date, or to configure the Overpass server timeout,
    memory allocation, and other custom settings. This function searches for
    features using tags. For more details, see:
    https://wiki.openstreetmap.org/wiki/Map_features

    Parameters
    ----------
    bbox
        Bounding box as `(left, bottom, right, top)`. Coordinates should be in
        unprojected latitude-longitude degrees (EPSG:4326).
    tags
        Tags for finding elements in the selected area. Results are the union,
        not intersection of the tags and each result matches at least one tag.
        The keys are OSM tags (e.g. `building`, `landuse`, `highway`, etc) and
        the values can be either `True` to retrieve all elements matching the
        tag, or a string to retrieve a single `tag:value` combination, or a
        list of strings to retrieve multiple values for the tag. For example,
        `tags = {'building': True}` would return all buildings in the area.
        Or, `tags = {'amenity':True, 'landuse':['retail','commercial'],
        'highway':'bus_stop'}` would return all amenities, any landuse=retail,
        any landuse=commercial, and any highway=bus_stop.

    Returns
    -------
    gdf
        The features, multi-indexed by element type and OSM ID.
    """
    # convert bbox to polygon then create GeoDataFrame of features within it
    polygon = utils_geo.bbox_to_poly(bbox)
    return features_from_polygon(polygon, tags)


def features_from_point(
    center_point: tuple[float, float],
    tags: dict[str, bool | str | list[str]],
    dist: float,
) -> gpd.GeoDataFrame:
    """
    Download OSM features within some distance of a lat-lon point.

    You can use the `settings` module to retrieve a snapshot of historical OSM
    data as of a certain date, or to configure the Overpass server timeout,
    memory allocation, and other custom settings. This function searches for
    features using tags. For more details, see:
    https://wiki.openstreetmap.org/wiki/Map_features

    Parameters
    ----------
    center_point
        The `(lat, lon)` center point around which to retrieve the features.
        Coordinates should be in unprojected latitude-longitude degrees
        (EPSG:4326).
    tags
        Tags for finding elements in the selected area. Results are the union,
        not intersection of the tags and each result matches at least one tag.
        The keys are OSM tags (e.g. `building`, `landuse`, `highway`, etc) and
        the values can be either `True` to retrieve all elements matching the
        tag, or a string to retrieve a single `tag:value` combination, or a
        list of strings to retrieve multiple values for the tag. For example,
        `tags = {'building': True}` would return all buildings in the area.
        Or, `tags = {'amenity':True, 'landuse':['retail','commercial'],
        'highway':'bus_stop'}` would return all amenities, any landuse=retail,
        any landuse=commercial, and any highway=bus_stop.
    dist
        Distance in meters from `center_point` to create a bounding box to
        query.

    Returns
    -------
    gdf
        The features, multi-indexed by element type and OSM ID.
    """
    # create bbox from point and dist, then create gdf of features within it
    bbox = utils_geo.bbox_from_point(center_point, dist)
    return features_from_bbox(bbox, tags)


def features_from_address(
    address: str,
    tags: dict[str, bool | str | list[str]],
    dist: float,
) -> gpd.GeoDataFrame:
    """
    Download OSM features within some distance of an address.

    You can use the `settings` module to retrieve a snapshot of historical OSM
    data as of a certain date, or to configure the Overpass server timeout,
    memory allocation, and other custom settings. This function searches for
    features using tags. For more details, see:
    https://wiki.openstreetmap.org/wiki/Map_features

    Parameters
    ----------
    address
        The address to geocode and use as the center point around which to
        retrieve the features.
    tags
        Tags for finding elements in the selected area. Results are the union,
        not intersection of the tags and each result matches at least one tag.
        The keys are OSM tags (e.g. `building`, `landuse`, `highway`, etc) and
        the values can be either `True` to retrieve all elements matching the
        tag, or a string to retrieve a single `tag:value` combination, or a
        list of strings to retrieve multiple values for the tag. For example,
        `tags = {'building': True}` would return all buildings in the area.
        Or, `tags = {'amenity':True, 'landuse':['retail','commercial'],
        'highway':'bus_stop'}` would return all amenities, any landuse=retail,
        any landuse=commercial, and any highway=bus_stop.
    dist
        Distance in meters from `address` to create a bounding box to query.

    Returns
    -------
    gdf
        The features, multi-indexed by element type and OSM ID.
    """
    # geocode the address to a point, then create gdf of features around it
    center_point = geocoder.geocode(address)
    return features_from_point(center_point, tags, dist)


def features_from_place(
    query: str | dict[str, str] | list[str | dict[str, str]],
    tags: dict[str, bool | str | list[str]],
    *,
    which_result: int | None | list[int | None] = None,
) -> gpd.GeoDataFrame:
    """
    Download OSM features within the boundaries of some place(s).

    The query must be geocodable and OSM must have polygon boundaries for the
    geocode result. If OSM does not have a polygon for this place, you can
    instead get features within it using the `features_from_address`
    function, which geocodes the place name to a point and gets the features
    within some distance of that point.

    If OSM does have polygon boundaries for this place but you're not finding
    it, try to vary the query string, pass in a structured query dict, or vary
    the `which_result` argument to use a different geocode result. If you know
    the OSM ID of the place, you can retrieve its boundary polygon using the
    `geocode_to_gdf` function, then pass it to the `features_from_polygon`
    function.

    You can use the `settings` module to retrieve a snapshot of historical OSM
    data as of a certain date, or to configure the Overpass server timeout,
    memory allocation, and other custom settings. This function searches for
    features using tags. For more details, see:
    https://wiki.openstreetmap.org/wiki/Map_features

    Parameters
    ----------
    query
        The query or queries to geocode to retrieve place boundary polygon(s).
    tags
        Tags for finding elements in the selected area. Results are the union,
        not intersection of the tags and each result matches at least one tag.
        The keys are OSM tags (e.g. `building`, `landuse`, `highway`, etc) and
        the values can be either `True` to retrieve all elements matching the
        tag, or a string to retrieve a single `tag:value` combination, or a
        list of strings to retrieve multiple values for the tag. For example,
        `tags = {'building': True}` would return all buildings in the area.
        Or, `tags = {'amenity':True, 'landuse':['retail','commercial'],
        'highway':'bus_stop'}` would return all amenities, any landuse=retail,
        any landuse=commercial, and any highway=bus_stop.
    which_result
        Which search result to return. If None, auto-select the first
        (Multi)Polygon or raise an error if OSM doesn't return one.

    Returns
    -------
    gdf
        The features, multi-indexed by element type and OSM ID.
    """
    # extract the geometry from the GeoDataFrame to use in query
    polygon = geocoder.geocode_to_gdf(query, which_result=which_result).union_all()
    msg = "Constructed place geometry polygon(s) to query Overpass"
    utils.log(msg, level=lg.INFO)

    # create GeoDataFrame using this polygon(s) geometry
    return features_from_polygon(polygon, tags)


def features_from_polygon(
    polygon: Polygon | MultiPolygon,
    tags: dict[str, bool | str | list[str]],
) -> gpd.GeoDataFrame:
    """
    Download OSM features within the boundaries of a (Multi)Polygon.

    You can use the `settings` module to retrieve a snapshot of historical OSM
    data as of a certain date, or to configure the Overpass server timeout,
    memory allocation, and other custom settings. This function searches for
    features using tags. For more details, see:
    https://wiki.openstreetmap.org/wiki/Map_features

    Parameters
    ----------
    polygon
        The geometry within which to retrieve features. Coordinates should be
        in unprojected latitude-longitude degrees (EPSG:4326).
    tags
        Tags for finding elements in the selected area. Results are the union,
        not intersection of the tags and each result matches at least one tag.
        The keys are OSM tags (e.g. `building`, `landuse`, `highway`, etc) and
        the values can be either `True` to retrieve all elements matching the
        tag, or a string to retrieve a single `tag:value` combination, or a
        list of strings to retrieve multiple values for the tag. For example,
        `tags = {'building': True}` would return all buildings in the area.
        Or, `tags = {'amenity':True, 'landuse':['retail','commercial'],
        'highway':'bus_stop'}` would return all amenities, any landuse=retail,
        any landuse=commercial, and any highway=bus_stop.

    Returns
    -------
    gdf
        The features, multi-indexed by element type and OSM ID.
    """
    # verify that the geometry is valid and is a Polygon/MultiPolygon
    if not polygon.is_valid:
        msg = "The geometry of `polygon` is invalid."
        raise ValueError(msg)

    if not isinstance(polygon, (Polygon, MultiPolygon)):
        msg = (
            "Boundaries must be a Polygon or MultiPolygon. If you requested "
            "`features_from_place`, ensure your query geocodes to a Polygon "
            "or MultiPolygon. See the documentation for details."
        )
        raise TypeError(msg)

    # retrieve the data from Overpass then turn it into a GeoDataFrame
    response_jsons = _overpass._download_overpass_features(polygon, tags)
    return _create_gdf(response_jsons, polygon, tags)


def features_from_xml(
    filepath: str | Path,
    *,
    polygon: Polygon | MultiPolygon | None = None,
    tags: dict[str, bool | str | list[str]] | None = None,
    encoding: str = "utf-8",
) -> gpd.GeoDataFrame:
    """
    Create a GeoDataFrame of OSM features from data in an OSM XML file.

    Because this function creates a GeoDataFrame of features from an OSM XML
    file that has already been downloaded (i.e., no query is made to the
    Overpass API), the `polygon` and `tags` arguments are optional. If they
    are None, filtering will be skipped.

    Parameters
    ----------
    filepath
        Path to file containing OSM XML data.
    tags
        Query tags to optionally filter the final GeoDataFrame.
    polygon
        Spatial boundaries to optionally filter the final GeoDataFrame.
    encoding
        The OSM XML file's character encoding.

    Returns
    -------
    gdf
        The features, multi-indexed by element type and OSM ID.
    """
    # if tags or polygon is None, create an empty object to skip filtering
    if tags is None:
        tags = {}
    if polygon is None:
        polygon = Polygon()

    # transmogrify OSM XML file to JSON then create GeoDataFrame from it
    response_jsons = [_osm_xml._overpass_json_from_xml(filepath, encoding)]
    gdf = _create_gdf(response_jsons, polygon, tags)

    # drop misc element attrs that might have been added from OSM XML file
    to_drop = set(gdf.columns) & {"changeset", "timestamp", "uid", "user", "version"}
    return gdf.drop(columns=list(to_drop))


def _create_gdf(
    response_jsons: Iterable[dict[str, Any]],
    polygon: Polygon | MultiPolygon,
    tags: dict[str, bool | str | list[str]],
) -> gpd.GeoDataFrame:
    """
    Convert Overpass API JSON responses to a GeoDataFrame of features.

    Parameters
    ----------
    response_jsons
        Iterable of Overpass API JSON responses.
    polygon
        Spatial boundaries to optionally filter the final GeoDataFrame.
    tags
        Query tags to optionally filter the final GeoDataFrame.

    Returns
    -------
    gdf
        GeoDataFrame of features with tags and geometry columns.
    """
    # consume response_jsons generator to download data from server
    elements = []
    response_count = 0
    for response_json in response_jsons:
        response_count += 1
        if not settings.cache_only_mode:
            elements.extend(response_json["elements"])

    msg = f"Retrieved {len(elements):,} elements from API in {response_count} request(s)"
    utils.log(msg, level=lg.INFO)
    if settings.cache_only_mode:
        msg = "Interrupted because `settings.cache_only_mode=True`."
        raise CacheOnlyInterruptError(msg)

    # convert the elements into a GeoDataFrame of features
    gdf = (
        gpd.GeoDataFrame(
            data=_process_features(elements, set(tags.keys())),
            geometry="geometry",
            crs=settings.default_crs,
        )
        .set_index(["element", "id"])
        .sort_index()
    )
    return _filter_features(gdf, polygon, tags)


def _process_features(
    elements: list[dict[str, Any]],
    query_tag_keys: set[str],
) -> list[dict[str, Any]]:
    """
    Convert node/way/relation elements into features with geometries.

    Parameters
    ----------
    elements
        The node/way/relation elements retrieved from the server.
    query_tag_keys
        The keys of the tags used to query for matching features.

    Returns
    -------
    features
    """
    nodes = []  # all nodes, including ones that just compose ways
    feature_nodes = []  # nodes that possibly match our query tags
    node_coords = {}  # hold node lon,lat tuples to create way geoms
    ways = []  # all ways, including ones that just compose relations
    feature_ways = []  # ways that possibly match our query tags
    way_geoms = {}  # hold way geoms to create relation geoms
    relations = []  # all relations

    # sort elements by node, way, and relation. only retain relations that
    # match the relation types we currently handle
    for element in elements:
        et = element["type"]
        if et == "node":
            nodes.append(element)
        elif et == "way":
            ways.append(element)
        elif et == "relation" and element.get("tags", {}).get("type") in _RELATION_TYPES:
            relations.append(element)

    # extract all nodes' coords, then add to features any nodes with tags that
    # match the passed query tags, or with any tags if no query tags passed
    for node in nodes:
        node_coords[node["id"]] = (node["lon"], node["lat"])
        if (len(query_tag_keys) == 0 and len(node.get("tags", {}).keys()) > 0) or (
            len(query_tag_keys & node.get("tags", {}).keys()) > 0
        ):
            node["element"] = node.pop("type")
            node["geometry"] = Point(node.pop("lon"), node.pop("lat"))
            node.update(node.pop("tags"))
            feature_nodes.append(node)

    # build all ways' geometries, then add to features any ways with tags that
    # match the passed query tags, or with any tags if no query tags passed
    for way in ways:
        way["geometry"] = _build_way_geometry(
            way["id"],
            way.pop("nodes"),
            way.get("tags", {}),
            node_coords,
        )
        way_geoms[way["id"]] = way["geometry"]
        if (len(query_tag_keys) == 0 and len(way.get("tags", {}).keys()) > 0) or (
            len(query_tag_keys & way.get("tags", {}).keys()) > 0
        ):
            way["element"] = way.pop("type")
            way.update(way.pop("tags"))
            feature_ways.append(way)

    # process relations and build their geometries
    for relation in relations:
        relation["element"] = "relation"
        relation.update(relation.pop("tags"))
        relation["geometry"] = _build_relation_geometry(relation.pop("members"), way_geoms)

    features = [*feature_nodes, *feature_ways, *relations]
    if len(features) == 0:
        msg = "No matching features. Check query location, tags, and log."
        raise InsufficientResponseError(msg)

    return features


def _build_way_geometry(
    way_id: int,
    way_nodes: list[int],
    way_tags: dict[str, Any],
    node_coords: dict[int, tuple[float, float]],
) -> LineString | Polygon:
    """
    Build a way's geometry from its constituent nodes' coordinates.

    A way can be a LineString (open or closed way) or a Polygon (closed way)
    but multi-geometries and polygons with holes are represented as relations.
    See documentation: https://wiki.openstreetmap.org/wiki/Way#Types_of_ways

    Parameters
    ----------
    way_id
        The way's OSM ID.
    way_nodes
        The way's constituent nodes.
    way_tags
        The way's tags.
    node_coords
        Keyed by OSM node ID with values of `(lat, lon)` coordinate tuples.

    Returns
    -------
    geometry
    """
    # a way is a LineString by default, but if it's a closed way and it's not
    # tagged area=no, check if any of its tags denote it as a polygon instead
    geom_type = LineString
    if way_nodes[0] == way_nodes[-1] and way_tags.get("area") != "no":
        for tag in way_tags.keys() & _POLYGON_FEATURES.keys():
            rule = _POLYGON_FEATURES[tag]["polygon"]
            values = _POLYGON_FEATURES[tag].get("values", set())
            if (
                rule == "all"
                or (rule == "passlist" and way_tags[tag] in values)
                or (rule == "blocklist" and way_tags[tag] not in values)
            ):
                geom_type = Polygon
                break

    # create the way geometry from its constituent nodes' coordinates
    try:
        return geom_type(node_coords[node] for node in way_nodes)
    except (GEOSException, KeyError, ValueError) as e:
        msg = f"Could not build geometry of way {way_id}: {e!r}"
        utils.log(msg, level=lg.WARNING)
        return geom_type()


def _build_relation_geometry(
    members: list[dict[str, Any]],
    way_geoms: dict[int, LineString | Polygon],
) -> Polygon | MultiPolygon:
    """
    Build a relation's geometry from its constituent member ways' geometries.

    OSM represents simple polygons as closed ways (see `_build_way_geometry`),
    but it uses relations to represent multipolygons (with or without holes)
    and polygons with holes. For the former, the relation contains multiple
    members with role "outer". For the latter, the relation contains at least
    one member with role "outer" representing the shell(s), and at least one
    member with role "inner" representing the hole(s). For documentation, see
    https://wiki.openstreetmap.org/wiki/Relation:multipolygon

    Parameters
    ----------
    members
        The members constituting the relation.
    way_geoms
        Keyed by OSM way ID with values of their geometries.

    Returns
    -------
    geometry
    """
    inner_linestrings = []
    outer_linestrings = []
    inner_polygons = []
    outer_polygons = []

    # sort member geometries by member role and geometry type
    for member in members:
        if member["type"] == "way":
            geom = way_geoms[member["ref"]]
            role = member["role"]
            if role == "outer" and geom.geom_type == "LineString":
                outer_linestrings.append(geom)
            elif role == "outer" and geom.geom_type == "Polygon":
                outer_polygons.append(geom)
            elif role == "inner" and geom.geom_type == "LineString":
                inner_linestrings.append(geom)
            elif role == "inner" and geom.geom_type == "Polygon":
                inner_polygons.append(geom)

    # merge/polygonize outer linestring fragments then add to outer polygons
    merged_outer_linestrings = linemerge(outer_linestrings)
    if merged_outer_linestrings.geom_type == "LineString":
        merged_outer_linestrings = MultiLineString([merged_outer_linestrings])
    for merged_outer_linestring in merged_outer_linestrings.geoms:
        outer_polygons += polygonize(merged_outer_linestring)

    # merge/polygonize inner linestring fragments then add to inner polygons
    merged_inner_linestrings = linemerge(inner_linestrings)
    if merged_inner_linestrings.geom_type == "LineString":
        merged_inner_linestrings = MultiLineString([merged_inner_linestrings])
    for merged_inner_linestring in merged_inner_linestrings.geoms:
        inner_polygons += polygonize(merged_inner_linestring)

    # remove holes from polygons, if any, then retun
    return _remove_polygon_holes(outer_polygons, inner_polygons)


def _remove_polygon_holes(
    outer_polygons: list[Polygon],
    inner_polygons: list[Polygon],
) -> Polygon | MultiPolygon:
    """
    Subtract inner holes from outer polygons.

    This allows possible island polygons within a larger polygon's holes.

    Parameters
    ----------
    outer_polygons
        Polygons, including possible islands within a larger polygon's holes.
    inner_polygons
        Inner holes to subtract from the outer polygons that contain them.

    Returns
    -------
    geometry
    """
    if len(inner_polygons) == 0:
        # if there are no holes to remove, geom is the union of outer polygons
        geometry = unary_union(outer_polygons)
    else:
        # otherwise, remove from each outer poly all inner polys it contains
        polygons_with_holes = []
        for outer in outer_polygons:
            prepare(outer)
            holes = [inner for inner in inner_polygons if outer.contains(inner)]
            polygons_with_holes.append(outer.difference(unary_union(holes)))
        geometry = unary_union(polygons_with_holes)

    # ensure returned geometry is a Polygon or MultiPolygon
    if isinstance(geometry, (Polygon, MultiPolygon)):
        return geometry
    return Polygon()


def _filter_features(
    gdf: gpd.GeoDataFrame,
    polygon: Polygon | MultiPolygon,
    tags: dict[str, bool | str | list[str]],
) -> gpd.GeoDataFrame:
    """
    Filter features GeoDataFrame by spatial boundaries and query tags.

    If the `polygon` and `tags` arguments are empty objects, the final
    GeoDataFrame will not be filtered accordingly.

    Parameters
    ----------
    gdf
        Original GeoDataFrame of features.
    polygon
        If not empty, the spatial boundaries to filter the GeoDataFrame.
    tags
        If not empty, the query tags to filter the GeoDataFrame.

    Returns
    -------
    gdf
        Filtered GeoDataFrame of features.
    """
    # remove any null or empty geometries then fix any invalid geometries
    gdf = gdf[~(gdf["geometry"].isna() | gdf["geometry"].is_empty)]
    gdf.loc[:, "geometry"] = gdf["geometry"].make_valid()

    # retain rows with geometries that intersect the polygon
    if polygon.is_empty:
        geom_filter = pd.Series(data=True, index=gdf.index)
    else:
        idx = utils_geo._intersect_index_quadrats(gdf["geometry"], polygon)
        geom_filter = gdf.index.isin(idx)

    # retain rows that have any of their tag filters satisfied
    if len(tags) == 0:
        tags_filter = pd.Series(data=True, index=gdf.index)
    else:
        tags_filter = pd.Series(data=False, index=gdf.index)
        for col in set(gdf.columns) & tags.keys():
            value = tags[col]
            if value is True:
                tags_filter |= gdf[col].notna()
            elif isinstance(value, str):
                tags_filter |= gdf[col] == value
            elif isinstance(value, list):
                tags_filter |= gdf[col].isin(set(value))

    # filter gdf then drop any columns with only nulls left after filtering
    gdf = gdf[geom_filter & tags_filter].dropna(axis="columns", how="all")
    if len(gdf) == 0:  # pragma: no cover
        msg = "No matching features. Check query location, tags, and log."
        raise InsufficientResponseError(msg)

    msg = f"{len(gdf):,} features in the final GeoDataFrame"
    utils.log(msg, level=lg.INFO)
    return gdf
