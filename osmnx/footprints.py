"""Download and plot footprints from OpenStreetMap."""

import geopandas as gpd
from shapely.geometry import LineString
from shapely.geometry import MultiPolygon
from shapely.geometry import Polygon
from shapely.ops import polygonize

from . import downloader
from . import geocoder
from . import projection
from . import settings
from . import utils
from . import utils_geo


def _osm_footprints_download(polygon, footprint_type="building"):
    """
    Download OpenStreetMap footprint data as a list of json responses.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        geographic boundaries to fetch the footprints within
    footprint_type : string
        type of footprint to be downloaded. OSM tag key e.g. 'building', 'landuse', 'place', etc.

    Returns
    -------
    response_jsons : list
        list of response_json dicts
    """
    response_jsons = []
    overpass_settings = downloader._make_overpass_settings()

    # project to utm, divide polygon up into sub-polygons if area exceeds a
    # max size (in meters), project back to lat-lng, then get a list of polygon(s) exterior coordinates
    geometry_proj, crs_proj = projection.project_geometry(polygon)
    gpcs = utils_geo._consolidate_subdivide_geometry(geometry_proj)
    geometry, _ = projection.project_geometry(gpcs, crs=crs_proj, to_latlong=True)
    polygon_coord_strs = utils_geo._get_polygons_coordinates(geometry)
    utils.log(
        f"Requesting footprints within polygon from API in {len(polygon_coord_strs)} request(s)"
    )

    # pass each polygon exterior coordinates in the list to the API, one at
    # a time
    for polygon_coord_str in polygon_coord_strs:
        query_str = (
            f"{overpass_settings};("
            f'way(poly:"{polygon_coord_str}")["{footprint_type}"];(._;>;);'
            f'relation(poly:"{polygon_coord_str}")["{footprint_type}"];(._;>;););out;'
        )
        response_json = downloader.overpass_request(data={"data": query_str})
        response_jsons.append(response_json)
    utils.log(
        f"Got all footprint data within polygon from API in {len(polygon_coord_strs)} request(s)"
    )

    return response_jsons


def _create_footprints_gdf(
    polygon=None, footprint_type="building", retain_invalid=False, responses=None
):
    """
    Get footprint (polygon) data from OSM and convert it into a GeoDataFrame.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        geographic boundaries to fetch the footprints within
    footprint_type : string
        type of footprint to be downloaded. OSM tag key e.g. 'building',
        'landuse', 'place', etc.
    retain_invalid : bool
        if False discard any footprints with an invalid geometry
    responses : list
        list of pre-existing response jsons

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    if polygon is None and responses is None:
        raise ValueError("You must pass a polygon or pre-existing responses.")

    if responses is None:
        responses = _osm_footprints_download(polygon, footprint_type)

    # parse the list of responses into separate dicts of vertices, footprints and relations
    # create a set of ways not directly tagged with footprint_type
    vertices, footprints, relations, untagged_ways = _responses_to_dicts(responses, footprint_type)

    # create simple Shapely geometries (Polygon or LineString) for all of the ways in footprints
    for footprint_key, footprint_val in footprints.items():
        footprint_val["geometry"] = _create_footprint_geometry(
            footprint_key, footprint_val, vertices
        )

    # create a complex Shapely Polygon or MultiPolygon for each relation
    for relation_key, relation_val in relations.items():
        relation_val["geometry"] = _create_relation_geometry(relation_key, relation_val, footprints)

    # merge relations into the footprints dictionary
    footprints.update(relations)

    # delete supporting geometry not directly tagged with footprint_type from the footprints dictionary
    for untagged_way in untagged_ways:
        try:
            del footprints[untagged_way]
        except KeyError:
            utils.log(f"untagged_way {untagged_way} not found in footprints dict")

    # Convert footprints dictionary to a GeoDataFrame
    gdf = gpd.GeoDataFrame.from_dict(footprints, orient="index")
    gdf.crs = settings.default_crs

    # filter the gdf to only include valid Polygons/MultiPolygons if retain_invalid is False
    if not retain_invalid and not gdf.empty:
        filter1 = gdf["geometry"].is_valid
        filter2 = (gdf["geometry"].geom_type == "Polygon") | (
            gdf["geometry"].geom_type == "MultiPolygon"
        )
        filter_combined = filter1 & filter2
        gdf = gdf[filter_combined]

    return gdf


def _responses_to_dicts(responses, footprint_type):
    """
    Parse list of json responses into dicts of vertices, footprints, relations.

    Note: OSM's data model and the Overpass API will return open ways (lines)
    as part of a 'polygon' query. These may be fragments of the inner and
    outer rings of relations or they may be open ways mistakenly tagged with
    'polygon' type tags.

    Ways not directly tagged with the footprint type are added to the
    untagged_ways set for removal from the footprints dictionary at the end of
    the process.

    Some inner ways of relations may be tagged with the footprint type in
    their own right e.g. landuse=meadow as an inner way in a landuse=forest
    relation and need to be kept. These are created here.

    Parameters
    ----------
    responses : list
        list of json responses
    footprint_type : string
        type of footprint downloaded. OSM tag key e.g. 'building', 'landuse',
        'place', etc.

    Returns
    -------
    (vertices, footprints, relations, untagged_footprints) : tuple
        vertices
            dictionary of OSM nodes including their lat, lng coordinates
        footprints
            dictionary of OSM ways including their nodes and tags
        relations
            dictionary of OSM relations including member ids and tags
        untagged_footprints
            set of ids for ways or relations not directly tagged with
            footprint_type
    """
    # create dictionaries to hold vertices, footprints and relations
    vertices = {}
    footprints = {}
    relations = {}

    # create a set to hold the ids of ways not directly tagged as footprint_type
    untagged_footprints = set()

    # loop through each response once adding each element to one of the dicts
    for response in responses:
        for element in response["elements"]:
            # NODES - only keep coordinates
            if "type" in element and element["type"] == "node":
                vertices[element["id"]] = {"lat": element["lat"], "lon": element["lon"]}
            # WAYS - both open and closed
            elif "type" in element and element["type"] == "way":
                footprint = {"nodes": element["nodes"]}
                if "tags" in element:
                    for tag in element["tags"]:
                        footprint[tag] = element["tags"][tag]
                footprints[element["id"]] = footprint
                # add ways not individually tagged with footprint_type to the
                # untagged_footprints set
                if ("tags" not in element) or (footprint_type not in element["tags"]):
                    untagged_footprints.add(element["id"])
            # RELATIONS
            elif "type" in element and element["type"] == "relation":
                relation = {"members": {}}
                for member in element["members"]:
                    if "type" in member and member["type"] == "way":
                        relation["members"].update({member["ref"]: member.get("role")})
                if "tags" in element:
                    for tag in element["tags"]:
                        relation[tag] = element["tags"][tag]
                relations[element["id"]] = relation
                # add relations not individually tagged with footprint_type to
                # the untagged_footprints set
                if ("tags" not in element) or (footprint_type not in element["tags"]):
                    untagged_footprints.add(element["id"])
            else:
                utils.log(f'Element {element["id"]} is not a node, way or relation')

    return vertices, footprints, relations, untagged_footprints


def _create_footprint_geometry(footprint_key, footprint_val, vertices):
    """
    Create geometry for footprint open/closed ways.

    Create Shapely geometry for open or closed ways in the initial footprints
    dictionary. Closed ways are converted directly to Shapely Polygons, open
    ways (fragments that will form the outer and inner rings of relations) are
    converted to LineStrings.

    Parameters
    ----------
    footprint_key : int
        the id of the way/footprint to process
    footprint_val : dict
        the nodes and tags of the footprint
    vertices : dict
        the dictionary of OSM nodes with their coordinates

    Returns
    -------
    Shapely Polygon or LineString
    """
    # CLOSED WAYS
    if footprint_val["nodes"][0] == footprint_val["nodes"][-1]:
        try:
            poly = [
                (vertices[node]["lon"], vertices[node]["lat"]) for node in footprint_val["nodes"]
            ]
            footprint_geometry = Polygon(poly)
        except Exception:
            utils.log(f"Polygon has invalid geometry: {footprint_key}")
    # OPEN WAYS
    else:
        try:
            ls = [(vertices[node]["lon"], vertices[node]["lat"]) for node in footprint_val["nodes"]]
            footprint_geometry = LineString(ls)
        except Exception:
            utils.log(f"LineString has invalid geometry: {footprint_key}")

    return footprint_geometry


def _create_relation_geometry(relation_key, relation_val, footprints):
    """
    Create Shapely geometry for relations: Polygons with holes or MultiPolygons.

    OSM relations are used to define complex polygons - polygons with holes or
    multi-polygons. The polygons' outer and inner rings may be made up of chains
    of LineStrings. https://wiki.openstreetmap.org/wiki/Relation:multipolygon
    requires that multipolygon rings have an outer or inner 'role'.

    OSM's data model allows a polygon type tag e.g. 'building' to be added to
    any OSM element. This can include non-polygon relations e.g. bus routes.
    Relations that do not have at least one closed ring with an outer role
    are filtered out.

    Inner rings that are tagged with the footprint type in their own right e.g.
    landuse=meadow as an inner ring of landuse=forest will have been included in
    the footprints dictionary as part of the original parsing and are not dealt
    with here.

    Parameters
    ----------
    relation_key : int
        the id of the relation to process
    relation_val : dict
        members and tags of the relation
    footprints : dict
        dictionary of all footprints (including open and closed ways)

    Returns
    -------
    Shapely Polygon or MultiPolygon
    """
    # lists to hold member geometries
    outer_polys, outer_lines, inner_polys, inner_lines = _members_geom_lists(
        relation_val, footprints
    )

    # try to polygonize open outer ways and concatenate them to outer_polys
    if len(outer_lines) > 0:
        try:
            result = list(polygonize(outer_lines))
        except Exception:
            utils.log(f"polygonize failed for outer ways in relation: {relation_key}")
        else:
            outer_polys += result

    # try to polygonize open inner ways and concatenate them to inner_polys
    if len(inner_lines) > 0:
        try:
            result = list(polygonize(inner_lines))
        except Exception:
            utils.log(f"polygonize failed for inner ways in relation: {relation_key}")
        else:
            inner_polys += result

    # filter out relations missing both 'outer' and 'inner' polygons or just 'outer'
    multipoly = []
    if len(outer_polys + inner_polys) == 0:
        utils.log(f"Relation {relation_key} missing outer and inner closed ways")
    elif len(outer_polys) == 0:
        utils.log(f"Relation {relation_key} missing outer closed ways")
    # process the others to multipolygons
    else:
        for outer_poly in outer_polys:
            outer_poly = outer_poly.buffer(0)  # fix invalid geometry if present
            temp_poly = outer_poly
            for inner_poly in inner_polys:
                inner_poly = inner_poly.buffer(0)  # fix invalid geometry if present
                if inner_poly.within(outer_poly):
                    temp_poly = temp_poly.difference(inner_poly)
            multipoly.append(temp_poly)

    # return relations with one outer way as Polygons, multiple outer ways
    # as MultiPolygons
    if len(multipoly) == 1:
        return multipoly[0]
    elif len(multipoly) > 1:
        return MultiPolygon(multipoly)
    else:
        utils.log(f"relation {relation_key} could not be converted to a complex footprint")


def _members_geom_lists(relation_val, footprints):
    """
    Add relation members' geoms to lists.

    Parameters
    ----------
    relation_val : dict
        members and tags of the relation
    footprints : dict
        dictionary of all footprints (including open and closed ways)

    Returns
    -------
    tuple of lists
    """
    outer_polys = []
    outer_lines = []
    inner_polys = []
    inner_lines = []

    # add each members geometry to a list according to its role and geometry type
    for member_id, member_role in relation_val["members"].items():
        if member_role == "outer":
            if footprints[member_id]["geometry"].geom_type == "Polygon":
                outer_polys.append(footprints[member_id]["geometry"])
            elif footprints[member_id]["geometry"].geom_type == "LineString":
                outer_lines.append(footprints[member_id]["geometry"])
        elif member_role == "inner":
            if footprints[member_id]["geometry"].geom_type == "Polygon":
                inner_polys.append(footprints[member_id]["geometry"])
            elif footprints[member_id]["geometry"].geom_type == "LineString":
                inner_lines.append(footprints[member_id]["geometry"])

    return outer_polys, outer_lines, inner_polys, inner_lines


def footprints_from_point(point, dist=1000, footprint_type="building", retain_invalid=False):
    """
    Get footprints within some distance N, S, E, W of a lat-lng point.

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
        if False discard any footprints with an invalid geometry

    Returns
    -------
    geopandas.GeoDataFrame

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """
    north, south, east, west = utils_geo.bbox_from_point(point=point, dist=dist)
    polygon = utils_geo.bbox_to_poly(north, south, east, west)
    return footprints_from_polygon(polygon, footprint_type, retain_invalid)


def footprints_from_address(address, dist=1000, footprint_type="building", retain_invalid=False):
    """
    Get footprints within some distance N, S, E, W of an address.

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
        if False discard any footprints with an invalid geometry

    Returns
    -------
    geopandas.GeoDataFrame

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """
    # geocode the address string to a (lat, lng) point
    point = geocoder.geocode(query=address)

    # get footprints within distance of this point
    return footprints_from_point(point, dist, footprint_type, retain_invalid)


def footprints_from_place(place, footprint_type="building", retain_invalid=False, which_result=1):
    """
    Get footprints within the boundaries of some place.

    The query must be geocodable and OSM must have polygon boundaries for the
    geocode result. If OSM does not have a polygon for this place, you can
    instead get its footprints using the footprints_from_address function,
    which geocodes the place name to a point and gets the footprints within
    some distance of that point.

    Parameters
    ----------
    place : string
        the query to geocode to get geojson boundary polygon
    footprint_type : string
        type of footprint to be downloaded. OSM tag key e.g. 'building',
        'landuse', 'place', etc.
    retain_invalid : bool
        if False discard any footprints with an invalid geometry
    which_result : int
        max number of results to return and which to process upon receipt

    Returns
    -------
    geopandas.GeoDataFrame

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """
    city = geocoder.geocode_to_gdf(place, which_result=which_result)
    polygon = city["geometry"].iloc[0]
    return footprints_from_polygon(polygon, footprint_type, retain_invalid)


def footprints_from_polygon(polygon, footprint_type="building", retain_invalid=False):
    """
    Get footprints within some polygon.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        the shape to get data within. coordinates should be in units of
        latitude-longitude degrees.
    footprint_type : string
        type of footprint to be downloaded. OSM tag key e.g. 'building',
        'landuse', 'place', etc.
    retain_invalid : bool
        if False discard any footprints with an invalid geometry

    Returns
    -------
    geopandas.GeoDataFrame

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """
    return _create_footprints_gdf(polygon, footprint_type, retain_invalid)
