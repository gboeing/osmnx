"""Download geometries from OpenStreetMap."""

import geopandas as gpd
from shapely.geometry import MultiPolygon
from shapely.geometry import Point
from shapely.geometry import LineString
from shapely.geometry import Polygon
from shapely.ops import polygonize

from . import downloader
from . import geocoder
from . import settings
from . import utils
from . import utils_geo


def _create_overpass_query(polygon, tags):
    """
    Create an overpass query string based on passed tags.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon
        geographic boundary to fetch geometry within
    tags : dict
        dict of tags used for finding geometry in the selected area

    Returns
    -------
    query : string
    """
    overpass_settings = downloader._make_overpass_settings()

    # make sure every value in dict is bool, str, or list of str
    error_msg = "tags must be a dict with values of bool, str, or list of str"
    if not isinstance(tags, dict):
        raise TypeError(error_msg)

    tags_dict = {}
    for key, value in tags.items():

        if isinstance(value, bool):
            tags_dict[key] = value

        elif isinstance(value, str):
            tags_dict[key] = [value]

        elif isinstance(value, list):
            if not all(isinstance(s, str) for s in value):
                raise TypeError(error_msg)
            tags_dict[key] = value

        else:
            raise TypeError(error_msg)

    # convert the tags dict into a list of {tag:value} dicts
    tags_list = []
    for key, value in tags_dict.items():
        if isinstance(value, bool):
            tags_list.append({key: value})
        else:
            for value_item in value:
                tags_list.append({key: value_item})

    # create query bounding box
    west, south, east, north = polygon.bounds
    bbox = f"({south:.6f},{west:.6f},{north:.6f},{east:.6f})"

    # add node/way/relation query components one at a time
    components = []
    for d in tags_list:
        for key, value in d.items():

            if isinstance(value, bool):
                # if bool (ie, True) just pass the key, no value
                tag_str = f'["{key}"]{bbox};(._;>;);'
            else:
                # otherwise, pass "key"="value"
                tag_str = f'["{key}"="{value}"]{bbox};(._;>;);'

            for kind in ("node", "way", "relation"):
                components.append(f"({kind}{tag_str});")

    # finalize query and return
    components = "".join(components)
    query = f"{overpass_settings};({components});out;"

    return query


def _get_overpass_response(polygon, tags):
    """
    Get geometry from OpenStreetMap based on passed tags.

    Note that if a polygon is passed-in, the query will be limited to its
    bounding box rather than to the shape of the polygon itself.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon
        geographic boundaries to fetch geometry within
    tags : dict
        dict of tags used for finding geometry in the selected area

    Returns
    -------
    responses : dict
        JSON response from Overpass server
    """
    # TODO: add functionality for subdividing search area geometry
    # TODO: add functionality for constraining query to poly rather than its bbox

    # get the geometry
    query = _create_overpass_query(polygon, tags)
    response = downloader.overpass_request(data={"data": query})
    return response


def _parse_node_to_coord(element):
    """
    Parse id and coordinates from a node in the overpass response.

    Some nodes are standalone points of interest, others are vertices in
    other geometry.

    Parameters
    ----------
    element : string
        element type "node" from OSM response

    Returns
    -------
    coord : dict
        dict of node ID and its lat, lng coordinates
    """
    # return the coordinate of a single node element
    coord = {"lat": element["lat"], "lon": element["lon"]}

    return coord


def _parse_node_to_point(element):
    """
    Parse point from OSM node element.

    Parameters
    ----------
    element : JSON string
        element type "node" from overpass response JSON

    Returns
    -------
    point_geometry : dict
        dict of OSM ID, Point geometry, and any tags
    """
    point = {}
    point["osmid"] = element["id"]
    point["element_type"] = "node"

    if "tags" in element:
        for tag in element["tags"]:
            point[tag] = element["tags"][tag]

    geometry = Point(element["lon"], element["lat"])

    point["geometry"] = geometry

    return point


def _parse_way_to_linestring_or_polygon(element, coords):
    """
    Parse linestring or polygon from OSM 'way'

    Parameters
    ----------
    element : JSON string
        element type "way" in OSM response with matching start and end point
    coords : dict
        dict of node IDs and their lat, lng coordinates

    Returns
    -------
    linestring_geometry : dict
        dict of OSM ID, LineString geometry, and any tags or
        None if it cannot
    """
    nodes = element["nodes"]

    linestring_or_polygon = {}
    linestring_or_polygon["osmid"] = element["id"]
    linestring_or_polygon["element_type"] = "way"
    linestring_or_polygon["nodes"] = nodes

    if "tags" in element:
        for tag in element["tags"]:
            linestring_or_polygon[tag] = element["tags"][tag]
    
    # If the way is open, parse it as a LineString
    # If the way is closed, parse it as a Polygon
    if element["nodes"][0] != element["nodes"][-1]:
        geometry = LineString([(coords[node]["lon"], coords[node]["lat"]) for node in nodes])
    elif element["nodes"][0] == element["nodes"][-1]:
        geometry = Polygon([(coords[node]["lon"], coords[node]["lat"]) for node in nodes])

    linestring_or_polygon["geometry"] = geometry
        
    return linestring_or_polygon


def _parse_relation_to_multipolygon(element, linestrings_and_polygons):
    """
    Parse multipolygon from OSM relation (type:MultiPolygon).

    See more information about relations from OSM documentation:
    http://wiki.openstreetmap.org/wiki/Relation

    Parameters
    ----------
    relations : list
        OSM 'relation' items (dictionaries) in a list
    df_osm_ways : geopandas..GeoDataFrame
        OSM 'way' features as a GeoDataFrame that contains all the 'way'
        features that will constitute the multipolygon relations

    Returns
    -------
    df_osm_ways : geopandas.GeoDataFrame
        GeoDataFrame with MultiPolygon representations of the relations and
        the attributes associated with them.
    """
    # Parse member 'way' ids
    member_way_refs = [
        member["ref"] for member in element["members"] if member["type"] == "way"
    ]
    # Extract the ways
    member_ways = [
        linestrings_and_polygons[member_way_ref] for member_way_ref in member_way_refs
    ]
    # Extract the nodes of those ways
    member_nodes = [
        [member_way["nodes"] for member_way in member_ways]
    ]

    multipolygon = {}
    multipolygon["osmid"] = element["id"]
    multipolygon["element_type"] = "relation"
    multipolygon["ways"] = member_way_refs
    multipolygon["nodes"] = member_nodes

    if "tags" in element:
        for tag in element["tags"]:
            multipolygon[tag] = element["tags"][tag]

    try:
        # Assemble MultiPolygon from linestrings and polygons
        geometry = _assemble_multipolygon_geometry(element, linestrings_and_polygons)

        multipolygon["geometry"] = geometry

        return multipolygon

    except Exception as e:
        print(e)
        utils.log(f'Could not parse OSM relation {element["id"]}')


def _assemble_multipolygon_geometry(element, linestrings_and_polygons):
    """

    """
    outer_member_geometries = []
    inner_member_geometries = []

    for member in element["members"]:
        if ("type" in member) and (member["type"] == "way"):
            if ("role" in member) and (member["role"] == "outer"):
                outer_member_geometry = linestrings_and_polygons.get(member["ref"])["geometry"]
                outer_member_geometries.append(outer_member_geometry)
            elif ("role" in member) and (member["role"] == "inner"):
                inner_member_geometry = linestrings_and_polygons.get(member["ref"])["geometry"]
                inner_member_geometries.append(inner_member_geometry)

    outer_polygons = list(polygonize(outer_member_geometries))
    inner_polygons = list(polygonize(inner_member_geometries))

    outer_polygons_with_holes = []

    for outer_polygon in outer_polygons:
        for inner_polygon in inner_polygons:
            if inner_polygon.within(outer_polygon):
                outer_polygon = outer_polygon.difference(inner_polygon)
        outer_polygons_with_holes.append(outer_polygon)

    geometry = MultiPolygon(outer_polygons_with_holes)

    return geometry


def _invalid_multipoly_handler(gdf, relation, way_ids):  # pragma: no cover
    """
    Handle invalid multipolygon geometries.

    For example, when there exists a feature without geometry (geometry==NaN).

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame with Polygon geometries that should be converted into a
        MultiPolygon object.
    relation : dict
        OSM 'relation' dictionary
    way_ids : list
        A list of 'way' ids that should be converted into a MultiPolygon object.

    Returns
    -------
    shapely.geometry.MultiPolygon
    """
    try:
        gdf_clean = gdf.dropna(subset=["geometry"])
        multipoly = MultiPolygon(list(gdf_clean["geometry"]))
        return multipoly

    except Exception:
        utils.log(f'Invalid geometry at relation "{relation["id"]}", way IDs: {way_ids}')
        return None


def _create_gdf(polygon, tags, response=None):
    """
    Create GeoDataFrame from json returned by Overpass API.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon
        geographic boundaries to fetch POIs within
    tags : dict
        dict of tags used for finding POIs from the selected area

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        geometries and their associated tags
    """
    if response is None:
        response = _get_overpass_response(polygon, tags)

    # Dictionaries to hold intermediate geometries
    coords = {}
    points = {}
    linestrings_and_polygons = {}
    multipolygons = {}

    for element in response["elements"]:

        if element["type"] == "node":
            # Parse all nodes to coords
            coord = _parse_node_to_coord(element=element)
            coords[element["id"]] = coord

            # Parse nodes tagged with requested tags to points
            if "tags" in element and (any([key in element["tags"].keys() for key in tags.keys()])):
                point = _parse_node_to_point(element=element)
                points[element["id"]] = point

        elif element["type"] == "way":
            # Parse all ways to linestrings or polygons
            linestring_or_polygon = _parse_way_to_linestring_or_polygon(element=element, coords=coords)
            if linestring_or_polygon:
                linestrings_and_polygons[element["id"]] = linestring_or_polygon

        elif element["type"] == "relation" and element["tags"]["type"] == "multipolygon":
            # Parse all multipolygon relations to multipolygons
            multipolygon = _parse_relation_to_multipolygon(element=element, linestrings_and_polygons=linestrings_and_polygons)
            if multipolygon:
                multipolygons[element["id"]] = multipolygon

    # Create GeoDataFrames
    gdf_points = gpd.GeoDataFrame.from_dict(points, orient="index")
    gdf_lines_and_polygons = gpd.GeoDataFrame.from_dict(linestrings_and_polygons, orient="index")
    gdf_multipolygons = gpd.GeoDataFrame.from_dict(multipolygons, orient="index")

    # Combine GeoDataFrames
    gdf = gdf_points.append(gdf_lines_and_polygons, sort=False)
    gdf = gdf.append(gdf_multipolygons, sort=False)

    # Set default crs
    gdf.crs = settings.default_crs

    # if caller requested pois within a polygon, only retain those that
    # fall within the polygon
    if len(gdf) > 0:
        gdf = gdf.loc[gdf["geometry"].centroid.within(polygon)]

    return gdf


def gdf_from_point(point, tags, dist=1000):
    """
    Get geometry within some distance N, S, E, W of a point.

    Parameters
    ----------
    point : tuple
        a (lat, lng) point
    tags : dict
        Dict of tags used for finding geometry from the selected area. Results
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
    north, south, east, west = utils_geo.bbox_from_point(point=point, dist=dist)
    polygon = utils_geo.bbox_to_poly(north, south, east, west)
    return gdf_from_polygon(polygon, tags)


def gdf_from_address(address, tags, dist=1000):
    """
    Get geometry within some distance N, S, E, W of address.

    Parameters
    ----------
    address : string
        the address to geocode to a lat-lng point
    tags : dict
        Dict of tags used for finding geometry in the selected area. Results
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
    # geocode the address string to a (lat, lng) point
    point = geocoder.geocode(query=address)
    return gdf_from_point(point=point, tags=tags, dist=dist)


def gdf_from_place(place, tags, which_result=1):
    """
    Get geometry within the boundaries of some place.

    Parameters
    ----------
    place : string
        the query to geocode to get boundary polygon
    tags : dict
        Dict of tags used for finding geometry in the selected area. Results
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
        max number of geocoding results to return and which to process

    Returns
    -------
    gdf : geopandas.GeoDataFrame

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """
    city = geocoder.geocode_to_gdf(place, which_result=which_result)
    polygon = city["geometry"].iloc[0]
    return gdf_from_polygon(polygon, tags)


def gdf_from_polygon(polygon, tags):
    """
    Get point of interests (POIs) within some polygon.

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
    return _create_gdf(polygon, tags)