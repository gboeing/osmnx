"""Download geometries from OpenStreetMap."""

import geopandas as gpd
import numpy as np
from shapely.geometry import Point
from shapely.geometry import LineString
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
from shapely.ops import linemerge
from shapely.ops import polygonize
from shapely.geos import TopologicalError

from . import downloader
from . import geocoder
from . import settings
from . import utils
from . import utils_geo
from ._errors import EmptyOverpassResponse
from .polygon_features import polygon_features


def gdf_from_bbox(point, size, tags):
    """
    Get geometry within bounding box around a central point.

    Parameters
    ----------
    point : tuple
        a (lat, lng) center of bounding box
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
    size : tuple
        a (width, height)
        Total width and height of the bounding box in meters

    Returns
    -------
    gdf : geopandas.GeoDataFrame

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """

    lat,   lng               = point
    width, height            = size
    north, south, _   , _    = utils_geo.bbox_from_point((lat,lng), dist=height/2)
    _    , _    , east, west = utils_geo.bbox_from_point((lat,lng), dist=width/2)

    polygon = utils_geo.bbox_to_poly(north, south, east, west)
    return gdf_from_polygon(polygon, tags)


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
    # create a GeoDataFrame with the spatial boundaries of the place(s)
    gdf_place = geocoder.geocode_to_gdf(place, which_result=which_result)

    # extract the geometry from the GeoDataFrame to use in API query
    polygon = gdf_place["geometry"].iloc[0]

    # create geodataframe using this polygon(s) geometry
    GDF = gdf_from_polygon(polygon, tags)

    utils.log(f"gdf_from_place returned geodataframe with {len(GDF)} geometries")
    return GDF


def gdf_from_polygon(polygon, tags):
    """
    Get point of interests (POIs) within some polygon.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon
        geographic boundaries to fetch geometries within
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
    # verify that the geometry is valid and is a shapely Polygon/MultiPolygon
    # before proceeding
    if not polygon.is_valid:
        raise ValueError("The geometry to query within is invalid")
    if not isinstance(polygon, (Polygon, MultiPolygon)):
        raise TypeError(
            "Geometry must be a shapely Polygon or MultiPolygon. If you requested "
            "graph from place name, make sure your query resolves to a Polygon or "
            "MultiPolygon, and not some other geometry, like a Point. See OSMnx "
            "documentation for details."
        )

    # download the geometry data from OSM
    response_jsons = downloader._osm_geometry_download(polygon, tags)

    # create geodataframe from the downloaded data
    GDF = _create_gdf(response_jsons, polygon, tags)

    utils.log(f"gdf_from_polygon returned geodataframe with {len(GDF)} geometries")
    return GDF


def _create_gdf(response_jsons, polygon, tags, saved_response_json=None):
    """
    Requests and parses a JSON response from the Overpass API to a GeoDataFrame.

    Parameters
    ----------
    response_jsons : list
        list of dicts of JSON responses from from the Overpass API
    polygon : shapely.geometry.Polygon
        geographic boundaries to fetch geometries within
    tags : dict
        dict of tags used for finding geometries from the selected area
    saved_response_json : JSON
        allows a single saved response JSON to be passed in for parsing

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        geometries and their associated tags
    """
    # RE/MOVE?
    # If a pre-saved response is passed in
    if saved_response_json is not None:
        # Requests a JSON from the Overpass API
        response_jsons = [saved_response_json]

    utils.log("Creating geodataframe from downloaded OSM data...")

    # make sure we got data back from the server requests
    elements = []
    for response_json in response_jsons:
        elements.extend(response_json["elements"])
    if len(elements) < 1:
        raise EmptyOverpassResponse("There are no data elements in the response JSON")

    # Dictionaries to hold nodes and complete geometries
    coords = {}
    geometries = {}

    # extract geometries from the downloaded osm data
    for response_json in response_jsons:
        # Parses the JSON of OSM nodes, ways and (multipolygon) relations to dictionaries of
        # coordinates, Shapely Points, LineStrings, Polygons and MultiPolygons
        for element in response_json["elements"]:

            # id numbers are only unique within element types
            # create unique id from combination of type and id
            unique_id = f"{element['type']}/{element['id']}"

            if element["type"] == "node":
                # Parse all nodes to coords
                coord = _parse_node_to_coord(element=element)
                coords[element["id"]] = coord
                # Parse nodes with tags to points
                if "tags" in element:
                    point = _parse_node_to_point(element=element)
                    geometries[unique_id] = point

            elif element["type"] == "way":
                # Parse all ways to linestrings or polygons
                linestring_or_polygon = _parse_way_to_linestring_or_polygon(
                    element=element,
                    coords=coords,
                    polygon_features=polygon_features,
                    )
                if linestring_or_polygon:
                    geometries[unique_id] = linestring_or_polygon

            elif element["type"] == "relation" and element.get("tags").get("type") == "multipolygon":
                # Parse all multipolygon relations to multipolygons
                multipolygon = _parse_relation_to_multipolygon(element=element, geometries=geometries)
                if multipolygon:
                    geometries[unique_id] = multipolygon

    # Create GeoDataFrame
    GDF = gpd.GeoDataFrame.from_dict(geometries, orient="index")

    # ensure GDF has a geometry col before assigning crs
    if "geometry" not in GDF.columns:
        # if there is no geometry column, create a null column
        GDF["geometry"] = np.nan
    GDF.set_geometry("geometry")

    # Set default crs
    GDF.crs = settings.default_crs

    # Filter final gdf to requested tags and bounding polygon
    GDF = _filter_final_gdf(gdf=GDF, polygon=polygon, tags=tags)

    return GDF


def _parse_node_to_coord(element):
    """
    Parse coordinates from a node in the overpass response.

    The coords are only used to create LineStrings and Polygons.

    Parameters
    ----------
    element : dict
        element type "node" from overpass response JSON

    Returns
    -------
    coord : dict
        dict of lat and lon coordinates
    """
    # return the coordinate of a single node element
    coord = {"lat": element["lat"], "lon": element["lon"]}

    return coord


def _parse_node_to_point(element):
    """
    Parse point from a tagged node in the overpass response.

    The points are geometries in their own right.

    Parameters
    ----------
    element : dict
        element type "node" from overpass response JSON

    Returns
    -------
    point : dict
        dict of OSM ID, OSM element type, tags and geometry
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


def _parse_way_to_linestring_or_polygon(element, coords, polygon_features):
    """
    Parse open LineString, closed LineString or Polygon from OSM 'way'

    Please see https://wiki.openstreetmap.org/wiki/Overpass_turbo/Polygon_Features
    for more information on which tags should be parsed to polygons

    Parameters
    ----------
    element : dict
        element type "way" from overpass response JSON
    coords : dict
        dict of node IDs and their lat, lon coordinates
    polygon_features : dict
        dict for determining whether closed ways are LineStrings or Polygons

    Returns
    -------
    linestring_or_polygon : dict
        dict of OSM ID, OSM element type, nodes, tags and geometry
    """
    nodes = element["nodes"]

    linestring_or_polygon = {}
    linestring_or_polygon["osmid"] = element["id"]
    linestring_or_polygon["element_type"] = "way"
    linestring_or_polygon["nodes"] = nodes

    if "tags" in element:
        for tag in element["tags"]:
            linestring_or_polygon[tag] = element["tags"][tag]

    # if the OSM element is an open way (i.e. first and last nodes are not the same)
    # the geometry should be a Shapely LineString
    if element["nodes"][0] != element["nodes"][-1]:
        geometry = LineString([(coords[node]["lon"], coords[node]["lat"]) for node in nodes])

    # if the OSM element is a closed way (i.e. first and last nodes are the same)
    # depending upon the tags the geometry could be a Shapely LineString or Polygon
    elif element["nodes"][0] == element["nodes"][-1]:
        # Use the tags to determine whether the way represents a LineString, Polygon
        is_polygon = _closed_way_is_linestring_or_polygon(element, polygon_features)

        if is_polygon:
            geometry = Polygon([(coords[node]["lon"], coords[node]["lat"]) for node in nodes])
        else:
            geometry = LineString([(coords[node]["lon"], coords[node]["lat"]) for node in nodes])

    linestring_or_polygon["geometry"] = geometry
        
    return linestring_or_polygon


def _closed_way_is_linestring_or_polygon(element, polygon_features):
    """
    Determines whether a closed OSM way represents a LineString or Polygon

    Closed OSM ways may represent LineStrings (e.g. a roundabout or hedge
    round a field) or Polygons (e.g. a building footprint or land use area)
    depending on the tags applied to them.

    It is possible for a single closed OSM way to have both LineString and
    Polygon type tags (e.g. both barrier=fence and landuse=agricultural).
    OSMnx will return a single Polygon for elements tagged in this way.
    For more information see:
    https://wiki.openstreetmap.org/wiki/One_feature,_one_OSM_element)

    Parameters
    ----------
    element : dict
        closed element type "way" from overpass response JSON
    polygon_features : dict
        dict of tag keys with associated values and blocklist/passlist

    Returns
    -------
    is_polygon : boolean
        True if the tags are for a polygon type geometry
    """
    is_polygon = False

    # get the element's tags
    element_tags = element.get('tags')

    # if the element doesn't have any tags it is a component of a multipolygon -> Polygon
    if element_tags is None:
        is_polygon = True
    
    # if the element has tags and is tagged 'area':'no' -> LineString
    elif element_tags.get('area') == 'no':
        is_polygon = False

    # if the element has tags and is not tagged 'area':'no'
    # compare its tags with the polygon_features dictionary
    else:
        # identify common keys in the element's tags and the polygon_features dictionary
        intersecting_keys = element_tags.keys() & polygon_features.keys()

        # if common keys are found
        if len(intersecting_keys) > 0:

            # for each key in the intersecting keys
            for key in intersecting_keys:
                # Get the key's value from the element's tags
                key_value = element_tags.get(key)
                # Determine if the key is for a blocklist or passlist in the polygon_features dictionary
                blocklist_or_passlist = polygon_features.get(key).get('polygon')
                # Get the values for the key from the polygon_features dictionary
                polygon_features_values = polygon_features.get(key).get('values')
                
                # if all features with that key should be polygons -> Polygon
                if blocklist_or_passlist == 'all':
                    is_polygon = True

                # if the key is for a blocklist
                elif blocklist_or_passlist == 'blocklist':
                    # if the value for that key in the element is in the blocklist -> LineString
                    # if the value for that key in the element is not in the blocklist -> Polygon
                    if key_value in polygon_features_values:
                        is_polygon = False
                    else:
                        is_polygon = True

                # if the key is for a passlist
                elif blocklist_or_passlist == 'passlist':
                    # if the value for that key in the element is in the passlist -> Polygon
                    # if the value for that key in the element is not in the passlist -> LineString
                    if key_value in polygon_features_values:
                        is_polygon = True
                    else:
                        is_polygon = False

        # the polygon_features dict is for determining which ways should become Polygons
        # therefore if no common keys are found the geometry should be a LineString by default
        else:
            is_polygon = False
        
    return is_polygon


def _parse_relation_to_multipolygon(element, geometries):
    """
    Parse multipolygon from OSM relation (type:MultiPolygon).

    See more information about relations from OSM documentation:
    http://wiki.openstreetmap.org/wiki/Relation

    Parameters
    ----------
    element : dict
        element type "relation" from overpass response JSON
    geometries : dict
        dictionary containing all linestrings and polygons generated from OSM ways

    Returns
    -------
    multipolygon : dictionary
        dictionary of tags and geometry for a single multipolygon
    """
    # Parse member 'way' ids
    member_way_refs = [
        member["ref"] for member in element["members"] if member["type"] == "way"
    ]
    # Extract the ways from the geometries dictionary using their unique id
    member_ways = [
        geometries[f"way/{member_way_ref}"] for member_way_ref in member_way_refs
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

    # Assemble MultiPolygon from component LineStrings and Polygons
    geometry = _assemble_multipolygon_geometry(element, geometries)

    multipolygon["geometry"] = geometry

    return multipolygon


def _assemble_multipolygon_geometry(element, geometries):
    """
    Assembles a MultiPolygon from its component LineStrings and Polygons.

    The OSM wiki suggests an algorithm for assembling multipolygon geometries
    https://wiki.openstreetmap.org/wiki/Relation:multipolygon/Algorithm.
    This method takes a simpler approach relying on the accurate tagging
    of component ways with 'inner' and 'outer' roles as required on this page
    https://wiki.openstreetmap.org/wiki/Relation:multipolygon.

    Parameters
    ----------
    element : dict
        element type "relation" from overpass response JSON
    geometries : dict
        dict containing all linestrings and polygons generated from OSM ways

    Returns
    -------
    geometry : MultiPolygon
        a single Shapely MultiPolygon
    """
    outer_polygons = []
    inner_polygons = []
    outer_linestrings = []
    inner_linestrings = []

    # get the linestrings and polygons that make up the multipolygon
    for member in element["members"]:
        if (member.get("type") == "way"):
            # get the member's geometry from linestrings_and_polygons
            linestring_or_polygon = geometries.get(f"way/{member['ref']}")
            # sort it into one of the lists according to its role and geometry
            if (member.get("role") == "outer") and (linestring_or_polygon["geometry"].geom_type == 'Polygon'):
                outer_polygons.append(linestring_or_polygon["geometry"])
            elif (member.get("role") == "inner") and (linestring_or_polygon["geometry"].geom_type == 'Polygon'):
                inner_polygons.append(linestring_or_polygon["geometry"])
            elif (member.get("role") == "outer") and (linestring_or_polygon["geometry"].geom_type == 'LineString'):
                outer_linestrings.append(linestring_or_polygon["geometry"])
            elif (member.get("role") == "inner") and (linestring_or_polygon["geometry"].geom_type == 'LineString'):
                inner_linestrings.append(linestring_or_polygon["geometry"])

    # Merge outer linestring fragments. Returns a single LineString or MultiLineString collection
    merged_outer_linestrings = linemerge(outer_linestrings)

    # polygonize each linestring separately and append to the list of outer polygons
    if merged_outer_linestrings.geom_type == 'LineString':
        outer_polygons += polygonize(merged_outer_linestrings)
    elif merged_outer_linestrings.geom_type == 'MultiLineString':
        for merged_outer_linestring in list(merged_outer_linestrings):
            outer_polygons += polygonize(merged_outer_linestring)

    # Merge inner linestring fragments. Returns a single LineString or MultiLineString collection
    merged_inner_linestrings = linemerge(inner_linestrings)

    # polygonize each linestring separately and append to the list of inner polygons
    if merged_inner_linestrings.geom_type == 'LineString':
        inner_polygons += polygonize(merged_inner_linestrings)
    elif merged_inner_linestrings.geom_type == 'MultiLineString':
        for merged_inner_linestring in list(merged_inner_linestrings):
            inner_polygons += polygonize(merged_inner_linestring)

    # create a new list to hold the outer polygons with the inner polygons subtracted
    outer_polygons_with_holes = []

    # loop through the outer polygons subtracting the inner polygons and appending to the list
    for outer_polygon in outer_polygons:
        for inner_polygon in inner_polygons:
            if inner_polygon.within(outer_polygon):
                try:
                    outer_polygon = outer_polygon.difference(inner_polygon)
                except TopologicalError as e:
                    print(e, f"\n MultiPolygon Relation OSM id {element['id']} difference failed, trying with geometries buffered by 0")
                    utils.log(f"relation https://www.openstreetmap.org/relation/{element['id']} caused a TopologicalError, trying with zero buffer")
                    outer_polygon = outer_polygon.buffer(0).difference(inner_polygon.buffer(0))

        # note: .buffer(0) can return either a Polygon or MultiPolygon
        # if it returns a MultiPolygon we need to extract the component
        # sub Polygons to add to outer_polygons_with_holes
        if outer_polygon.geom_type == 'Polygon':
            outer_polygons_with_holes.append(outer_polygon)
        elif outer_polygon.geom_type == 'MultiPolygon':
            outer_polygons_with_holes.extend(list(outer_polygon))

    # create a multipolygon from the list of outer polygons with holes
    geometry = MultiPolygon(outer_polygons_with_holes)

    return geometry


def _filter_final_gdf(gdf, polygon, tags):
    """
    Filters the final gdf to the requested tags and bounding polygon

    Filters the final gdf to the requested tags and bounding polygon. Removes
    columns of all NaNs (that held values only in rows removed by the filters).
    Resets the gdf index.

    Parameters
    ----------
    gdf : GeoDataFrame
        the GeoDataFrame to filter
    polygon : Polygon
        Shapely polygon defining the boundary of the requested area
    tags : dict
        the tags requested

    Returns
    -------
    gdf : GeoDataFrame
        final, filtered GeoDataFrame
    """
    if not gdf.empty:

        # filter retaining geometries within the bounding polygon using spatial index
        if polygon is not None:
                gdf_indices_in_polygon = utils_geo._intersect_index_quadrats(gdf.centroid, polygon)
                gdf = gdf[gdf.index.isin(gdf_indices_in_polygon)]

        # filter retaining geometries with the requested tags
        if tags is not None:
            # Intersect the tags and column names to be sure the tags are present in the columns
            tags_present_in_columns = set(tags.keys()).intersection(set(gdf.columns))
            # Select rows which have non-null values in any of these columns
            gdf = gdf[gdf[tags_present_in_columns].notna().any(axis=1)]

        # remove columns of all nulls (created by discarded component geometries)
        gdf.dropna(axis='columns', how='all', inplace=True)

        # reset the index.
        # Theoretically points, linestrings and polygons, multipolygons could share index numbers
        gdf.reset_index(drop=True, inplace=True)

    return gdf