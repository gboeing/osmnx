"""Download geometries from OpenStreetMap."""

import geopandas as gpd
from shapely.geometry import MultiPolygon
from shapely.geometry import Point
from shapely.geometry import LineString
from shapely.geometry import Polygon

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
    try:
        geometry = Point(element["lon"], element["lat"])

        point_geometry = {"osmid": element["id"], "geometry": geometry}

        for tag in element["tags"]:
            point_geometry[tag] = element["tags"][tag]

        # Add element_type
        point_geometry["element_type"] = "node"

    except Exception:
        utils.log(f'OSM node {element["id"]} has invalid point geometry')

    return point_geometry


def _parse_way_to_line_or_polygon(coords, element):
    """
    Parse linestring or polygon from OSM 'way'

    Parameters
    ----------
    coords : dict
        dict of node IDs and their lat, lng coordinates
    element : JSON string
        element type "way" in OSM response with matching start and end point

    Returns
    -------
    linestring_geometry : dict
        dict of OSM ID, LineString geometry, and any tags or
        None if it cannot
    """
    nodes = element["nodes"]
    
    try:
        # If the way is open, parse it as a LineString
        if element["nodes"][0] != element["nodes"][-1]:
            geometry = LineString([(coords[node]["lon"], coords[node]["lat"]) for node in nodes])
        # If the way is closed, parse it as a Polygon
        elif element["nodes"][0] == element["nodes"][-1]:
            geometry = Polygon([(coords[node]["lon"], coords[node]["lat"]) for node in nodes])

        line_or_polygon = {"nodes": nodes, "geometry": geometry, "osmid": element["id"]}

        if "tags" in element:
            for tag in element["tags"]:
                line_or_polygon[tag] = element["tags"][tag]
        
        # Add element_type
        line_or_polygon["element_type"] = "way"

        return line_or_polygon

    except Exception:
        # Don't think this will be triggered
        # By design, Shapely will build invalid Polygons without complaint
        utils.log(f'OSM way {element["id"]} has invalid geometry.')


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


def _parse_osm_relations(relations, df_osm_ways):
    """
    Parse OSM relations (MultiPolygons) from ways and nodes.

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
    gdf_relations = gpd.GeoDataFrame()

    # Iterate over relations and extract the items
    for relation in relations:
        try:
            if relation["tags"]["type"] == "multipolygon":
                # Parse member 'way' ids
                member_way_ids = [
                    member["ref"] for member in relation["members"] if member["type"] == "way"
                ]
                # Extract the ways
                member_ways = df_osm_ways.reindex(member_way_ids)
                # Extract the nodes of those ways
                member_nodes = list(member_ways["nodes"].values)
                try:
                    # Create MultiPolygon from geometries (exclude NaNs)
                    multipoly = MultiPolygon(list(member_ways["geometry"]))
                except Exception:
                    multipoly = _invalid_multipoly_handler(
                        gdf=member_ways, relation=relation, way_ids=member_way_ids
                    )

                if multipoly:
                    # Create GeoDataFrame with the tags and the MultiPolygon and its
                    # 'ways' (ids), and the 'nodes' of those ways
                    geo = gpd.GeoDataFrame(relation["tags"], index=[relation["id"]])
                    # Initialize columns (needed for .loc inserts)
                    geo = geo.assign(
                        geometry=None, ways=None, nodes=None, element_type=None, osmid=None
                    )
                    # Add attributes
                    geo.loc[relation["id"], "geometry"] = multipoly
                    geo.loc[relation["id"], "ways"] = member_way_ids
                    geo.loc[relation["id"], "nodes"] = member_nodes
                    geo.loc[relation["id"], "element_type"] = "relation"
                    geo.loc[relation["id"], "osmid"] = relation["id"]

                    # Append to relation GeoDataFrame
                    gdf_relations = gdf_relations.append(geo, sort=False)
                    # Remove such 'ways' from 'df_osm_ways' that are part of the 'relation'
                    df_osm_ways = df_osm_ways.drop(member_way_ids)
        except Exception:
            utils.log(f'Could not parse OSM relation {relation["id"]}')

    # Merge df_osm_ways and the gdf_relations
    df_osm_ways = df_osm_ways.append(gdf_relations, sort=False)
    return df_osm_ways


def _create_gdf(polygon, tags):
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
    response = _get_overpass_response(polygon, tags)

    # Dict to hold all node coordinates in the response
    coords = {}

    # Points
    point_geometries = {}

    # Lines and Polygons
    line_and_polygon_geometries = {}

    # A list of POI relations
    relations = []

    for element in response["elements"]:

        if element["type"] == "node":
            coord = _parse_node_to_coord(element=element)
            # Add to coords
            coords[element["id"]] = coord

            if "tags" in element:
                point_geometry = _parse_node_to_point(element=element)
                # Add to 'points'
                point_geometries[element["id"]] = point_geometry

        elif element["type"] == "way":

                line_or_polygon = _parse_way_to_line_or_polygon(coords=coords, element=element)
                if line_or_polygon:
                    # Add to 'line_and_polygon_geometries'
                    line_and_polygon_geometries[element["id"]] = line_or_polygon

        elif element["type"] == "relation":
            # Add relation to a relation list (needs to be parsed after
            # all nodes and ways have been parsed)
            relations.append(element)

    # Create GeoDataFrames
    gdf_points = gpd.GeoDataFrame(point_geometries).T
    gdf_points.crs = settings.default_crs

    gdf_lines_and_polygons = gpd.GeoDataFrame(line_and_polygon_geometries).T
    gdf_lines_and_polygons.crs = settings.default_crs

    # Parse relations (MultiPolygons) from 'ways'
    gdf_lines_and_polygons = _parse_osm_relations(relations=relations, df_osm_ways=gdf_lines_and_polygons)

    # Combine GeoDataFrames
    gdf = gdf_points.append(gdf_lines_and_polygons, sort=False)

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