"""
Download OpenStreetMap geospatial features' geometries and attributes.

Retrieve points of interest, building footprints, transit lines/stops, or any
other map features from OSM, including their geometries and attribute data,
then construct a GeoDataFrame of them. You can use this module to query for
nodes, ways, and relations (the latter of type "multipolygon" or "boundary"
only) by passing a dictionary of desired OSM tags.

For more details, see https://wiki.openstreetmap.org/wiki/Map_features and
https://wiki.openstreetmap.org/wiki/Elements

Refer to the Getting Started guide for usage limitations.
"""

import logging as lg
import warnings
from warnings import warn

import geopandas as gpd
import pandas as pd
from shapely.errors import GEOSException
from shapely.errors import TopologicalError
from shapely.geometry import LineString
from shapely.geometry import MultiPolygon
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.ops import linemerge
from shapely.ops import polygonize

from . import _overpass
from . import geocoder
from . import osm_xml
from . import settings
from . import utils
from . import utils_geo
from ._errors import CacheOnlyInterruptError
from ._errors import InsufficientResponseError

# dict of tags to determine if closed ways should be polygons, based on JSON
# from https://wiki.openstreetmap.org/wiki/Overpass_turbo/Polygon_Features
_POLYGON_FEATURES = {
    "building": {"polygon": "all"},
    "highway": {"polygon": "passlist", "values": ["services", "rest_area", "escape", "elevator"]},
    "natural": {
        "polygon": "blocklist",
        "values": ["coastline", "cliff", "ridge", "arete", "tree_row"],
    },
    "landuse": {"polygon": "all"},
    "waterway": {"polygon": "passlist", "values": ["riverbank", "dock", "boatyard", "dam"]},
    "amenity": {"polygon": "all"},
    "leisure": {"polygon": "all"},
    "barrier": {
        "polygon": "passlist",
        "values": ["city_wall", "ditch", "hedge", "retaining_wall", "spikes"],
    },
    "railway": {
        "polygon": "passlist",
        "values": ["station", "turntable", "roundhouse", "platform"],
    },
    "area": {"polygon": "all"},
    "boundary": {"polygon": "all"},
    "man_made": {"polygon": "blocklist", "values": ["cutline", "embankment", "pipeline"]},
    "power": {"polygon": "passlist", "values": ["plant", "substation", "generator", "transformer"]},
    "place": {"polygon": "all"},
    "shop": {"polygon": "all"},
    "aeroway": {"polygon": "blocklist", "values": ["taxiway"]},
    "tourism": {"polygon": "all"},
    "historic": {"polygon": "all"},
    "public_transport": {"polygon": "all"},
    "office": {"polygon": "all"},
    "building:part": {"polygon": "all"},
    "military": {"polygon": "all"},
    "ruins": {"polygon": "all"},
    "area:highway": {"polygon": "all"},
    "craft": {"polygon": "all"},
    "golf": {"polygon": "all"},
    "indoor": {"polygon": "all"},
}


def features_from_bbox(north=None, south=None, east=None, west=None, bbox=None, tags=None):
    """
    Create a GeoDataFrame of OSM features within a N, S, E, W bounding box.

    You can use the `settings` module to retrieve a snapshot of historical OSM
    data as of a certain date, or to configure the Overpass server timeout,
    memory allocation, and other custom settings.

    For more details, see: https://wiki.openstreetmap.org/wiki/Map_features

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
    tags : dict
        Dict of tags used for finding elements in the selected area. Results
        returned are the union, not intersection of each individual tag.
        Each result matches at least one given tag. The dict keys should be
        OSM tags, (e.g., `building`, `landuse`, `highway`, etc) and the dict
        values should be either `True` to retrieve all items with the given
        tag, or a string to get a single tag-value combination, or a list of
        strings to get multiple values for the given tag. For example,
        `tags = {'building': True}` would return all building footprints in
        the area. `tags = {'amenity':True, 'landuse':['retail','commercial'],
        'highway':'bus_stop'}` would return all amenities, landuse=retail,
        landuse=commercial, and highway=bus_stop.

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    if not (north is None and south is None and east is None and west is None):
        msg = (
            "The `north`, `south`, `east`, and `west` parameters are deprecated and "
            "will be removed in the v2.0.0 release. Use the `bbox` parameter instead. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
        )
        warn(msg, FutureWarning, stacklevel=2)
        bbox = (north, south, east, west)

    msg = (
        "The expected order of coordinates in `bbox` will change in the "
        "v2.0.0 release to `(left, bottom, right, top)`."
    )
    warn(msg, FutureWarning, stacklevel=2)

    # convert bounding box to a polygon
    polygon = utils_geo.bbox_to_poly(bbox=bbox)

    # create GeoDataFrame of features within this polygon
    return features_from_polygon(polygon, tags)


def features_from_point(center_point, tags, dist=1000):
    """
    Create GeoDataFrame of OSM features within some distance N, S, E, W of a point.

    You can use the `settings` module to retrieve a snapshot of historical OSM
    data as of a certain date, or to configure the Overpass server timeout,
    memory allocation, and other custom settings.

    For more details, see: https://wiki.openstreetmap.org/wiki/Map_features

    Parameters
    ----------
    center_point : tuple
        the (lat, lon) center point around which to get the features
    tags : dict
        Dict of tags used for finding elements in the selected area. Results
        returned are the union, not intersection of each individual tag.
        Each result matches at least one given tag. The dict keys should be
        OSM tags, (e.g., `building`, `landuse`, `highway`, etc) and the dict
        values should be either `True` to retrieve all items with the given
        tag, or a string to get a single tag-value combination, or a list of
        strings to get multiple values for the given tag. For example,
        `tags = {'building': True}` would return all building footprints in
        the area. `tags = {'amenity':True, 'landuse':['retail','commercial'],
        'highway':'bus_stop'}` would return all amenities, landuse=retail,
        landuse=commercial, and highway=bus_stop.
    dist : numeric
        distance in meters

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    # create bounding box from center point and distance in each direction
    bbox = utils_geo.bbox_from_point(center_point, dist)

    # convert the bounding box to a polygon
    polygon = utils_geo.bbox_to_poly(bbox=bbox)

    # create GeoDataFrame of features within this polygon
    return features_from_polygon(polygon, tags)


def features_from_address(address, tags, dist=1000):
    """
    Create GeoDataFrame of OSM features within some distance N, S, E, W of address.

    You can use the `settings` module to retrieve a snapshot of historical OSM
    data as of a certain date, or to configure the Overpass server timeout,
    memory allocation, and other custom settings.

    For more details, see: https://wiki.openstreetmap.org/wiki/Map_features

    Parameters
    ----------
    address : string
        the address to geocode and use as the central point around which to
        get the features
    tags : dict
        Dict of tags used for finding elements in the selected area. Results
        returned are the union, not intersection of each individual tag.
        Each result matches at least one given tag. The dict keys should be
        OSM tags, (e.g., `building`, `landuse`, `highway`, etc) and the dict
        values should be either `True` to retrieve all items with the given
        tag, or a string to get a single tag-value combination, or a list of
        strings to get multiple values for the given tag. For example,
        `tags = {'building': True}` would return all building footprints in
        the area. `tags = {'amenity':True, 'landuse':['retail','commercial'],
        'highway':'bus_stop'}` would return all amenities, landuse=retail,
        landuse=commercial, and highway=bus_stop.
    dist : numeric
        distance in meters

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    # geocode the address string to a (lat, lon) point
    center_point = geocoder.geocode(query=address)

    # create GeoDataFrame of features around this point
    return features_from_point(center_point, tags, dist=dist)


def features_from_place(query, tags, which_result=None, buffer_dist=None):
    """
    Create GeoDataFrame of OSM features within boundaries of some place(s).

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
    memory allocation, and other custom settings.

    For more details, see: https://wiki.openstreetmap.org/wiki/Map_features

    Parameters
    ----------
    query : string or dict or list
        the query or queries to geocode to get place boundary polygon(s)
    tags : dict
        Dict of tags used for finding elements in the selected area. Results
        returned are the union, not intersection of each individual tag.
        Each result matches at least one given tag. The dict keys should be
        OSM tags, (e.g., `building`, `landuse`, `highway`, etc) and the dict
        values should be either `True` to retrieve all items with the given
        tag, or a string to get a single tag-value combination, or a list of
        strings to get multiple values for the given tag. For example,
        `tags = {'building': True}` would return all building footprints in
        the area. `tags = {'amenity':True, 'landuse':['retail','commercial'],
        'highway':'bus_stop'}` would return all amenities, landuse=retail,
        landuse=commercial, and highway=bus_stop.
    which_result : int
        which geocoding result to use. if None, auto-select the first
        (Multi)Polygon or raise an error if OSM doesn't return one.
    buffer_dist : float
        deprecated, do not use

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    if buffer_dist is not None:
        warn(
            "The buffer_dist argument has been deprecated and will be removed "
            "in the v2.0.0 release. Buffer your query area directly, if desired. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123",
            FutureWarning,
            stacklevel=2,
        )

    # create a GeoDataFrame with the spatial boundaries of the place(s)
    if isinstance(query, (str, dict)):
        # if it is a string (place name) or dict (structured place query),
        # then it is a single place
        gdf_place = geocoder.geocode_to_gdf(
            query, which_result=which_result, buffer_dist=buffer_dist
        )
    elif isinstance(query, list):
        # if it is a list, it contains multiple places to get
        gdf_place = geocoder.geocode_to_gdf(query, buffer_dist=buffer_dist)
    else:  # pragma: no cover
        msg = "query must be dict, string, or list of strings"
        raise TypeError(msg)

    # extract the geometry from the GeoDataFrame to use in API query
    polygon = gdf_place["geometry"].unary_union
    utils.log("Constructed place geometry polygon(s) to query API")

    # create GeoDataFrame using this polygon(s) geometry
    return features_from_polygon(polygon, tags)


def features_from_polygon(polygon, tags):
    """
    Create GeoDataFrame of OSM features within boundaries of a (multi)polygon.

    You can use the `settings` module to retrieve a snapshot of historical OSM
    data as of a certain date, or to configure the Overpass server timeout,
    memory allocation, and other custom settings.

    For more details, see: https://wiki.openstreetmap.org/wiki/Map_features

    Parameters
    ----------
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        geographic boundaries to fetch features within
    tags : dict
        Dict of tags used for finding elements in the selected area. Results
        returned are the union, not intersection of each individual tag.
        Each result matches at least one given tag. The dict keys should be
        OSM tags, (e.g., `building`, `landuse`, `highway`, etc) and the dict
        values should be either `True` to retrieve all items with the given
        tag, or a string to get a single tag-value combination, or a list of
        strings to get multiple values for the given tag. For example,
        `tags = {'building': True}` would return all building footprints in
        the area. `tags = {'amenity':True, 'landuse':['retail','commercial'],
        'highway':'bus_stop'}` would return all amenities, landuse=retail,
        landuse=commercial, and highway=bus_stop.

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    # verify that the geometry is valid and a Polygon/MultiPolygon
    if not polygon.is_valid:
        msg = "The geometry of `polygon` is invalid"
        raise ValueError(msg)
    if not isinstance(polygon, (Polygon, MultiPolygon)):
        msg = (
            "Boundaries must be a shapely Polygon or MultiPolygon. If you "
            "requested features from place name, make sure your query resolves "
            "to a Polygon or MultiPolygon, and not some other geometry, like a "
            "Point. See OSMnx documentation for details."
        )
        raise TypeError(msg)

    # download the data from OSM
    response_jsons = _overpass._download_overpass_features(polygon, tags)

    # create GeoDataFrame from the downloaded data
    return _create_gdf(response_jsons, polygon, tags)


def features_from_xml(filepath, polygon=None, tags=None, encoding="utf-8"):
    """
    Create a GeoDataFrame of OSM features in an OSM-formatted XML file.

    Because this function creates a GeoDataFrame of features from an
    OSM-formatted XML file that has already been downloaded (i.e. no query is
    made to the Overpass API) the polygon and tags arguments are not required.
    If they are not supplied to the function, features_from_xml() will
    return features for all of the tagged elements in the file. If they are
    supplied they will be used to filter the final GeoDataFrame.

    For more details, see: https://wiki.openstreetmap.org/wiki/Map_features

    Parameters
    ----------
    filepath : string or pathlib.Path
        path to file containing OSM XML data
    polygon : shapely.geometry.Polygon
        optional geographic boundary to filter elements
    tags : dict
        optional dict of tags for filtering elements from the XML. Results
        returned are the union, not intersection of each individual tag.
        Each result matches at least one given tag. The dict keys should be
        OSM tags, (e.g., `building`, `landuse`, `highway`, etc) and the dict
        values should be either `True` to retrieve all items with the given
        tag, or a string to get a single tag-value combination, or a list of
        strings to get multiple values for the given tag. For example,
        `tags = {'building': True}` would return all building footprints in
        the area. `tags = {'amenity':True, 'landuse':['retail','commercial'],
        'highway':'bus_stop'}` would return all amenities, landuse=retail,
        landuse=commercial, and highway=bus_stop.
    encoding : string
        the XML file's character encoding

    Returns
    -------
    gdf : geopandas.GeoDataFrame
    """
    # transmogrify file of OSM XML data into JSON
    response_jsons = [osm_xml._overpass_json_from_file(filepath, encoding)]

    # create GeoDataFrame using this response JSON
    return _create_gdf(response_jsons, polygon=polygon, tags=tags)


def _create_gdf(response_jsons, polygon, tags):
    """
    Parse JSON responses from the Overpass API to a GeoDataFrame.

    Note: the `polygon` and `tags` arguments can both be `None` and the
    GeoDataFrame will still be created but it won't be filtered at the end
    i.e. the final GeoDataFrame will contain all tagged features in the
    `response_jsons`.

    Parameters
    ----------
    response_jsons : list
        list of JSON responses from from the Overpass API
    polygon : shapely.geometry.Polygon
        geographic boundary used for filtering the final GeoDataFrame
    tags : dict
        dict of tags used for filtering the final GeoDataFrame

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame of features and their associated tags
    """
    response_count = 0
    if settings.cache_only_mode:
        # if cache_only_mode, consume response_jsons then interrupt
        for _ in response_jsons:
            response_count += 1
        utils.log(f"Retrieved all data from API in {response_count} request(s)")
        msg = "Interrupted because `settings.cache_only_mode=True`"
        raise CacheOnlyInterruptError(msg)

    # Dictionaries to hold nodes and complete geometries
    coords = {}
    geometries = {}

    # Set to hold the unique IDs of elements that do not have tags
    untagged_element_ids = set()

    # identify which relation types to parse to (multi)polygons
    relation_types = {"boundary", "multipolygon"}

    # extract geometries from the downloaded osm data
    for response_json in response_jsons:
        response_count += 1

        # Parses the JSON of OSM nodes, ways and (multipolygon) relations
        # to dictionaries of coordinates, Shapely Points, LineStrings,
        # Polygons and MultiPolygons
        for element in response_json["elements"]:
            # id numbers are only unique within element types
            # create unique id from combination of type and id
            unique_id = f"{element['type']}/{element['id']}"

            # add elements that are not nodes and that are without tags or
            # with empty tags to the untagged_element_ids set (untagged
            # nodes are not added to the geometries dict at all)
            if (element["type"] != "node") and (("tags" not in element) or (not element["tags"])):
                untagged_element_ids.add(unique_id)

            if element["type"] == "node":
                # Parse all nodes to coords
                coords[element["id"]] = _parse_node_to_coords(element=element)

                # If the node has tags and the tags are not empty parse it
                # to a Point. Empty check is necessary for JSONs created
                # from XML where nodes without tags are assigned tags={}
                if "tags" in element and len(element["tags"]) > 0:
                    point = _parse_node_to_point(element=element)
                    geometries[unique_id] = point

            elif element["type"] == "way":
                # Parse all ways to linestrings or polygons
                linestring_or_polygon = _parse_way_to_linestring_or_polygon(
                    element=element, coords=coords
                )
                geometries[unique_id] = linestring_or_polygon

            elif (
                element["type"] == "relation" and element.get("tags").get("type") in relation_types
            ):
                # parse relations to (multi)polygons
                multipolygon = _parse_relation_to_multipolygon(
                    element=element, geometries=geometries
                )
                geometries[unique_id] = multipolygon

    utils.log(f"Retrieved all data from API in {response_count} request(s)")

    # ensure we got some node/way data back from the server request(s)
    if len(geometries) == 0:  # pragma: no cover
        msg = "No data elements in server response. Check log and query location/tags."
        raise InsufficientResponseError(msg)

    # remove untagged elements from the final dict of geometries
    utils.log(f"{len(geometries)} geometries created in the dict")
    for untagged_element_id in untagged_element_ids:
        geometries.pop(untagged_element_id, None)
    utils.log(f"{len(untagged_element_ids)} untagged features removed")

    # create GeoDataFrame, ensure it has geometry, then set crs
    gdf = gpd.GeoDataFrame.from_dict(geometries, orient="index")
    if "geometry" not in gdf.columns:
        # if there is no geometry column, create a null column
        gdf = gdf.set_geometry([None] * len(gdf))
    gdf = gdf.set_crs(settings.default_crs)

    # Apply .buffer(0) to any invalid geometries
    gdf = _buffer_invalid_geometries(gdf)

    # Filter final gdf to requested tags and query polygon
    gdf = _filter_gdf_by_polygon_and_tags(gdf, polygon=polygon, tags=tags)

    # bug in geopandas <0.9 raises a TypeError if trying to plot empty
    # geometries but missing geometries (gdf['geometry'] = None) cannot be
    # projected e.g. gdf.to_crs(). Remove rows with empty (e.g. Point())
    # or missing (e.g. None) geometry, and suppress gpd warning caused by
    # calling gdf["geometry"].isna() on GeoDataFrame with empty geometries
    if not gdf.empty:
        warnings.filterwarnings("ignore", "GeoSeries.isna", UserWarning)
        gdf = gdf[~(gdf["geometry"].is_empty | gdf["geometry"].isna())].copy()
        warnings.resetwarnings()

    utils.log(f"{len(gdf)} features in the final GeoDataFrame")
    return gdf


def _parse_node_to_coords(element):
    """
    Parse coordinates from a node in the overpass response.

    The coords are only used to create LineStrings and Polygons.

    Parameters
    ----------
    element : dict
        element type "node" from overpass response JSON

    Returns
    -------
    coords : dict
        dict of latitude/longitude coordinates
    """
    # return the coordinate of a single node element
    return {"lat": element["lat"], "lon": element["lon"]}


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

    point["geometry"] = Point(element["lon"], element["lat"])
    return point


def _parse_way_to_linestring_or_polygon(element, coords):
    """
    Parse open LineString, closed LineString or Polygon from OSM 'way'.

    Please see https://wiki.openstreetmap.org/wiki/Overpass_turbo/Polygon_Features
    for more information on which tags should be parsed to polygons

    Parameters
    ----------
    element : dict
        element type "way" from overpass response JSON
    coords : dict
        dict of node IDs and their latitude/longitude coordinates

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

    # un-nest individual tags
    if "tags" in element:
        for tag in element["tags"]:
            linestring_or_polygon[tag] = element["tags"][tag]

    # if the OSM element is an open way (i.e. first and last nodes are not the
    # same) the geometry should be a Shapely LineString
    if element["nodes"][0] != element["nodes"][-1]:
        try:
            geometry = LineString([(coords[node]["lon"], coords[node]["lat"]) for node in nodes])
        except KeyError as e:  # pragma: no cover
            # XMLs may include geometries that are incomplete, in which case
            # return an empty geometry
            utils.log(
                f"node/{e} was not found in `coords`.\n"
                f"https://www.openstreetmap.org/{element['type']}/{element['id']} was not created."
            )
            geometry = LineString()

    # if the OSM element is a closed way (i.e. first and last nodes are the
    # same) depending upon the tags the geometry could be a Shapely LineString
    # or Polygon
    elif element["nodes"][0] == element["nodes"][-1]:
        # determine if closed way represents LineString or Polygon
        if _is_closed_way_a_polygon(element):
            # if it is a Polygon
            try:
                geometry = Polygon([(coords[node]["lon"], coords[node]["lat"]) for node in nodes])
            except (GEOSException, ValueError) as e:
                # XMLs may include geometries that are incomplete, in which
                # case return an empty geometry
                utils.log(
                    f"{e} . The geometry for "
                    f"https://www.openstreetmap.org/{element['type']}/{element['id']} was not created."
                )
                geometry = Polygon()
        else:
            # if it is a LineString
            try:
                geometry = LineString(
                    [(coords[node]["lon"], coords[node]["lat"]) for node in nodes]
                )
            except (GEOSException, ValueError) as e:
                # XMLs may include geometries that are incomplete, in which
                # case return an empty geometry
                utils.log(
                    f"{e} . The geometry for "
                    f"https://www.openstreetmap.org/{element['type']}/{element['id']} was not created."
                )
                geometry = LineString()

    linestring_or_polygon["geometry"] = geometry
    return linestring_or_polygon


def _is_closed_way_a_polygon(element, polygon_features=_POLYGON_FEATURES):
    """
    Determine whether a closed OSM way represents a Polygon, not a LineString.

    Closed OSM ways may represent LineStrings (e.g. a roundabout or hedge
    round a field) or Polygons (e.g. a building footprint or land use area)
    depending on the tags applied to them.

    The starting assumption is that it is not a polygon, however any polygon
    type tagging will return a polygon unless explicitly tagged with area:no.

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
    is_polygon : bool
        True if the tags are for a polygon type geometry
    """
    # polygon_features dict is for determining which ways should become Polygons
    # therefore the starting assumption is that the geometry is a LineString
    is_polygon = False

    # get the element's tags
    element_tags = element.get("tags")

    # if the element doesn't have any tags leave it as a Linestring
    if element_tags is not None:
        # if the element is specifically tagged 'area':'no' -> LineString
        if element_tags.get("area") == "no":
            pass

        # if the element has tags and is not tagged 'area':'no'
        # compare its tags with the polygon_features dict
        else:
            # identify common keys in element's tags and polygon_features dict
            intersecting_keys = element_tags.keys() & polygon_features.keys()

            # for each key in the intersecting keys (if any found)
            for key in intersecting_keys:
                # Get the key's value from the element's tags
                key_value = element_tags.get(key)

                # Determine if the key is for a blocklist or passlist in
                # polygon_features dict
                blocklist_or_passlist = polygon_features.get(key).get("polygon")

                # Get values for the key from the polygon_features dict
                polygon_features_values = polygon_features.get(key).get("values")

                # if all features with that key should be polygons -> Polygon
                if blocklist_or_passlist == "all":
                    is_polygon = True

                # if the key is for a blocklist i.e. tags that should not
                # become Polygons
                elif blocklist_or_passlist == "blocklist":
                    # if the value for that key in the element is not in
                    # the blocklist -> Polygon
                    if key_value not in polygon_features_values:
                        is_polygon = True

                # if the key is for a passlist i.e. specific tags should
                # become Polygons, and if the value for that key in the
                # element is in the passlist -> Polygon
                elif (blocklist_or_passlist == "passlist") and (
                    key_value in polygon_features_values
                ):
                    is_polygon = True

    return is_polygon


def _parse_relation_to_multipolygon(element, geometries):
    """
    Parse multipolygon from OSM relation (type:MultiPolygon).

    See more information about relations from OSM documentation:
    https://wiki.openstreetmap.org/wiki/Relation

    Parameters
    ----------
    element : dict
        element type "relation" from overpass response JSON
    geometries : dict
        dict containing all linestrings and polygons generated from OSM ways

    Returns
    -------
    multipolygon : dict
        dict of tags and geometry for a single multipolygon
    """
    multipolygon = {}
    multipolygon["osmid"] = element["id"]
    multipolygon["element_type"] = "relation"

    # Parse member 'way' ids
    member_way_refs = [member["ref"] for member in element["members"] if member["type"] == "way"]
    multipolygon["ways"] = member_way_refs

    # Add the tags
    if "tags" in element:
        for tag in element["tags"]:
            multipolygon[tag] = element["tags"][tag]

    # Extract the ways from the geometries dict using their unique id.
    # XMLs exported from the openstreetmap.org homepage with a bounding box
    # may include the relation but not the ways outside the bounding box.
    try:
        member_ways = [geometries[f"way/{member_way_ref}"] for member_way_ref in member_way_refs]
    except KeyError as e:  # pragma: no cover
        utils.log(
            f"{e} was not found in `geometries`.\nThe geometry for "
            f"https://www.openstreetmap.org/{element['type']}/{element['id']} was not created."
        )
        multipolygon["geometry"] = MultiPolygon()
        return multipolygon

    # Extract the nodes of those ways
    member_nodes = [[member_way["nodes"] for member_way in member_ways]]
    multipolygon["nodes"] = member_nodes

    # Assemble MultiPolygon component polygons from component LineStrings and
    # Polygons
    outer_polygons, inner_polygons = _assemble_multipolygon_component_polygons(element, geometries)

    # Subtract inner polygons from outer polygons
    geometry = _subtract_inner_polygons_from_outer_polygons(element, outer_polygons, inner_polygons)

    multipolygon["geometry"] = geometry
    return multipolygon


def _assemble_multipolygon_component_polygons(element, geometries):
    """
    Assemble a MultiPolygon from its component LineStrings and Polygons.

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
    geometry : shapely.geometry.MultiPolygon
        a single MultiPolygon object
    """
    outer_polygons = []
    inner_polygons = []
    outer_linestrings = []
    inner_linestrings = []

    # get the linestrings and polygons that make up the multipolygon
    for member in element["members"]:
        if member.get("type") == "way":
            # get the member's geometry from linestrings_and_polygons
            linestring_or_polygon = geometries.get(f"way/{member['ref']}")
            # sort it into one of the lists according to its role and geometry
            if (member.get("role") == "outer") and (
                linestring_or_polygon["geometry"].geom_type == "Polygon"
            ):
                outer_polygons.append(linestring_or_polygon["geometry"])
            elif (member.get("role") == "inner") and (
                linestring_or_polygon["geometry"].geom_type == "Polygon"
            ):
                inner_polygons.append(linestring_or_polygon["geometry"])
            elif (member.get("role") == "outer") and (
                linestring_or_polygon["geometry"].geom_type == "LineString"
            ):
                outer_linestrings.append(linestring_or_polygon["geometry"])
            elif (member.get("role") == "inner") and (
                linestring_or_polygon["geometry"].geom_type == "LineString"
            ):
                inner_linestrings.append(linestring_or_polygon["geometry"])

    # Merge outer linestring fragments.
    # Returns a single LineString or MultiLineString collection
    merged_outer_linestrings = linemerge(outer_linestrings)

    # polygonize each linestring separately and append to list of outer polygons
    if merged_outer_linestrings.geom_type == "LineString":
        outer_polygons += polygonize(merged_outer_linestrings)
    elif merged_outer_linestrings.geom_type == "MultiLineString":
        for merged_outer_linestring in list(merged_outer_linestrings.geoms):
            outer_polygons += polygonize(merged_outer_linestring)

    # Merge inner linestring fragments.
    # Returns a single LineString or MultiLineString collection
    merged_inner_linestrings = linemerge(inner_linestrings)

    # polygonize each linestring separately and append to list of inner polygons
    if merged_inner_linestrings.geom_type == "LineString":
        inner_polygons += polygonize(merged_inner_linestrings)
    elif merged_inner_linestrings.geom_type == "MultiLineString":
        for merged_inner_linestring in merged_inner_linestrings.geoms:
            inner_polygons += polygonize(merged_inner_linestring)

    if not outer_polygons:
        utils.log(
            "No outer polygons were created for"
            f" https://www.openstreetmap.org/{element['type']}/{element['id']}"
        )

    return outer_polygons, inner_polygons


def _subtract_inner_polygons_from_outer_polygons(element, outer_polygons, inner_polygons):
    """
    Subtract inner polygons from outer polygons.

    Creates a Polygon or MultiPolygon with holes.

    Parameters
    ----------
    element : dict
        element type "relation" from overpass response JSON
    outer_polygons : list
        list of outer polygons that are part of a multipolygon
    inner_polygons : list
        list of inner polygons that are part of a multipolygon

    Returns
    -------
    geometry : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        a single Polygon or MultiPolygon
    """
    # create a new list to hold the outer polygons with the inner polygons
    # subtracted
    outer_polygons_with_holes = []

    # loop through the outer polygons subtracting the inner polygons and
    # appending to the list
    for outer_polygon in outer_polygons:
        outer_polygon_diff = outer_polygon
        for inner_polygon in inner_polygons:
            if inner_polygon.within(outer_polygon):
                try:
                    outer_polygon_diff = outer_polygon_diff.difference(inner_polygon)
                except TopologicalError:  # pragma: no cover
                    utils.log(
                        f"relation https://www.openstreetmap.org/relation/{element['id']} "
                        "caused a TopologicalError, trying with zero buffer."
                    )
                    outer_polygon_diff = outer_polygon.buffer(0).difference(inner_polygon.buffer(0))

        # note: .buffer(0) can return either a Polygon or MultiPolygon
        # if it returns a MultiPolygon we need to extract the component
        # sub Polygons to add to outer_polygons_with_holes
        if outer_polygon_diff.geom_type == "Polygon":
            outer_polygons_with_holes.append(outer_polygon_diff)
        elif outer_polygon_diff.geom_type == "MultiPolygon":
            outer_polygons_with_holes.extend(list(outer_polygon_diff.geoms))

    # if only one polygon with holes was created, return that single polygon
    if len(outer_polygons_with_holes) == 1:
        geometry = outer_polygons_with_holes[0]
    # otherwise create a multipolygon from list of outer polygons with holes
    else:
        geometry = MultiPolygon(outer_polygons_with_holes)

    return geometry


def _buffer_invalid_geometries(gdf):
    """
    Buffer any invalid geometries remaining in the GeoDataFrame.

    Invalid geometries in the GeoDataFrame (which may accurately reproduce
    invalid geometries in OpenStreetMap) can cause the filtering to the query
    polygon and other subsequent geometric operations to fail. This function
    logs the ids of the invalid geometries and applies a buffer of zero to try
    to make them valid.

    Note: the resulting geometries may differ from the originals - please
    check them against OpenStreetMap

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        a GeoDataFrame with possibly invalid geometries

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        the GeoDataFrame with .buffer(0) applied to invalid geometries
    """
    # only apply the filters if the GeoDataFrame is not empty
    if not gdf.empty:
        # create a filter for rows with invalid geometries
        invalid_geometry_filter = ~gdf["geometry"].is_valid

        # if there are invalid geometries
        if invalid_geometry_filter.any():
            # get their unique_ids from the index
            invalid_geometry_ids = gdf.loc[invalid_geometry_filter].index.to_list()

            # create a list of their urls and log them
            osm_url = "https://www.openstreetmap.org/"
            invalid_geom_urls = [osm_url + unique_id for unique_id in invalid_geometry_ids]
            utils.log(
                f"{len(invalid_geometry_ids)} invalid geometries"
                f".buffer(0) applied to {invalid_geom_urls}",
                level=lg.WARNING,
            )

            gdf.loc[invalid_geometry_filter, "geometry"] = gdf.loc[
                invalid_geometry_filter, "geometry"
            ].buffer(0)

    return gdf


def _filter_gdf_by_polygon_and_tags(gdf, polygon, tags):
    """
    Filter the GeoDataFrame to the requested bounding polygon and tags.

    Filters GeoDataFrame to query polygon and tags. Removes columns of all
    NaNs (that held values only in rows removed by the filters). Resets the
    index of GeoDataFrame, writing it into a new column called 'unique_id'.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        the GeoDataFrame to filter
    polygon : shapely.geometry.Polygon
        polygon defining the boundary of the requested area
    tags : dict
        the tags requested

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        final filtered GeoDataFrame
    """
    # only apply the filters if the GeoDataFrame is not empty
    if not gdf.empty:
        # create two filters, initially all True
        polygon_filter = pd.Series(True, index=gdf.index)
        combined_tag_filter = pd.Series(True, index=gdf.index)

        # if a polygon was supplied, create a filter that is True for
        # features that intersect with the polygon
        if polygon:
            # get set of index labels of features that intersect polygon
            gdf_indices_in_polygon = utils_geo._intersect_index_quadrats(gdf, polygon)
            # create boolean series, True for features whose index is in set
            polygon_filter = gdf.index.isin(gdf_indices_in_polygon)

            utils.log(f"{sum(~polygon_filter)} features removed by the polygon filter")

        # if tags were supplied, create filter that is True for features
        # that have at least one of the requested tags
        if tags:
            # Reset all values in the combined_tag_filter to False
            combined_tag_filter[:] = False

            # reduce the tags to those that are actually present in the
            # GeoDataFrame columns
            tags_in_columns = {key: tags[key] for key in tags if key in gdf.columns}

            for key, value in tags_in_columns.items():
                if value is True:
                    tag_filter = gdf[key].notna()
                elif isinstance(value, str):
                    tag_filter = gdf[key] == value
                elif isinstance(value, list):
                    tag_filter = gdf[key].isin(value)

                combined_tag_filter = combined_tag_filter | tag_filter

            utils.log(f"{sum(~combined_tag_filter)} features removed by the tag filter")

        # apply the filters
        gdf = gdf[polygon_filter & combined_tag_filter].copy()

        # remove columns of all nulls (created by discarded component features)
        gdf = gdf.dropna(axis="columns", how="all")

    # multi-index gdf by element_type and osmid then return
    idx_cols = ["element_type", "osmid"]
    if all(c in gdf.columns for c in idx_cols):
        gdf = gdf.set_index(idx_cols)
    return gdf
