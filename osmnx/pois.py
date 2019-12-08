################################################################################
# Module: pois.py
# Description: Download and plot points of interests (POIs) from OpenStreetMap
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

import geopandas as gpd
from shapely.geometry import MultiPolygon
from shapely.geometry import Point
from shapely.geometry import Polygon

from . import settings
from .core import bbox_from_point
from .core import gdf_from_place
from .core import overpass_request
from .geo_utils import geocode
from .utils import log


def create_overpass_query(north, south, east, west, amenities=None, tags=None,
                    timeout=180, maxsize=None):
    """
    Generate the Overpass QL query string given provided amenities and tags.

    Parameters
    ----------

    north : float
        Northernmost coordinate from bounding box of the search area.
    south : float
        Southernmost coordinate from bounding box of the search area.
    east : float
        Easternmost coordinate from bounding box of the search area.
    west : float
        Westernmost coordinate of the bounding box of the search area.
    amenities : list
        List of amenities that will be used for finding the POIs from the
        selected area. Results returned are the _union_, not _intersection_ of
        each individual amenity tag; each result matches at least one tag given.
        Equivalent to `tags={'amenity': amenities}`. Some POIs don't use the
        `amenity` key, in which case you should use the `tags` option
    tags : dict[str, Union[bool, str, List[str]]]
        Dict of tags used for finding POIs from the selected area. Results
        returned are the _union_, not _intersection_ of each individual tag;
        each result matches at least one tag given.

        The keys should be the OSM tag key, (i.e. `amenity`, `landuse`,
        `highway`, etc) and the values should be `True` to retrieve all objects
        with the given tag, a string to get a single tag-value combination, and
        a list of strings to get multiple values for the given tag. For example,

            tags = {
                'amenity':True,
                'landuse':['retail','commercial'],
                'highway':'bus_stop'
            }

        Would return all amenities, `landuse=retail` and `landuse=commercial`,
        and `highway=bus_stop`. If both `amenities` and `tags` are combined, any
        `amenity` tags will be combined. Providing `False` is not currently
        supported.
    timeout : int
        Timeout for the API request.
    """
    if amenities is None and tags is None:
        raise ValueError('Either amenities or tags must be provided')
    if amenities is not None and not isinstance(amenities, list):
        raise TypeError('amenities must be None or a list of str')
    if tags is not None and not isinstance(tags, dict):
        msg = 'tags must be None or a dict of bool, str, or list of str'
        raise TypeError(msg)

    # Merge `amenities` and `tags` parameters into a single dict (`_tags`).
    # `amenities` is kept for backwards compatibility.
    _tags = {}

    # If `tags` is provided, use that as _tags, making sure that every dict
    # value is either a list of strings or a boolean
    if tags is not None:
        for key, value in tags.items():
            # If value is boolean, set _tags value as boolean
            if isinstance(value, bool):
                _tags[key] = value
                continue

            # If value is string, set _tags value as a list containing that str
            if isinstance(value, str):
                _tags[key] = [value]
                continue

            # If value is list, set _tags value as that list
            # Make sure that every value of the list is a str, to prevent issues
            # later
            if isinstance(value, list):
                if not all(isinstance(s, str) for s in value):
                    msg = 'If a list is provided as a value of the `tags` dict, it must be a list of strings'
                    raise TypeError(msg)
                _tags[key] = value
                continue

            # If we've reached here, the value in the key-value combination is
            # neither a boolean, string, or list, so raise a TypeError
            raise TypeError('Values of tags dict must be bool, str, or list')

    # If `amenities` is provided, add those values to `tags['amenity']`
    # Note that if `tags['amenity']` is True, all values of tag `amenity` will
    # be retrieved, so don't try to add
    if amenities is not None:
        existing_val = _tags.get('amenity')
        # If the `amenity` key already exists and is a list, add `amenities` to
        # it
        if isinstance(existing_val, list):
            _tags['amenity'].extend(amenities)
        # If the `amenity` key already exists and is a boolean, don't add
        # `amenities`
        elif isinstance(existing_val, bool):
            pass
        # If the `amenity` key does not already exist, set it as `amenities`
        elif existing_val is None:
            _tags['amenity'] = amenities

    # Because Overpass QL doesn't have an OR operand, you can't do
    # `"amenity"="a|b"` to retrieve all amenities that have either `a` or `b` as
    # tags. Instead you have to do `"amenity"="a"` and then also `"amenity"="b"`
    #
    # To fit this model, here the _tags dictionary, which has either a boolean
    # or list of str as values, is transformed into a list of dictionaries, each
    # with a single value.
    # So
    # ```
    # {
    #     'amenity':True,
    #     'landuse':['retail','commercial'],
    #     'highway':'bus_stop'
    # }
    # ```
    # is transformed to
    # [
    #     {'amenity': True},
    #     {'landuse': 'retail'},
    #     {'landuse': 'commercial'},
    #     {'highway': 'bus_stop'}
    # ]
    _tags_list = []
    for key, value in _tags.items():
        if isinstance(value, bool):
            _tags_list.append({key: value})
            continue

        for value_item in value:
            _tags_list.append({key: value_item})

    # Create query template
    # Build query for each key-value pair for each of nodes, ways, and
    # relations.

    # Set query settings: Query settings must all be on the same line. The body
    # of the query can be over separate lines. See documentation here:
    # https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL#Settings
    # This sets the output format, server timeout, and bounding box
    query_settings = [
        '[out:json]',
        '[timeout:{t}]'.format(t=timeout),
        '[bbox:{s:.6f},{w:.6f},{n:.6f},{e:.6f}]'.format(s=south, w=west,
                                                        n=north, e=east)
    ]
    if maxsize is not None:
        query_settings.append('[maxsize:{m}]'.format(m=maxsize))

    # Start building query string
    # Keep it as list of str for now; then join at the end
    query = [
        ''.join(query_settings) + ';',
        '('
    ]

    # For each key, value pair in the list of dicts, build
    # (node["key"="value"];(._;>;);)
    # (way["key"="value"];(._;>;);)
    # (relation["key"="value"];(._;>;);)
    # lines in the query string
    for d in _tags_list:
        for key, value in d.items():
            tag_str = '["{key}"'.format(key=key)

            if not isinstance(value, bool):
                tag_str += '="{value}"'.format(value=value)
            # The `(._; >;);` joins the new result set into the previous result
            # set: https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL#Item
            tag_str += '];(._;>;);'

            for t in ['node', 'way', 'relation']:
                query.append('({t}{tag_str});'.format(t=t, tag_str=tag_str))

    query.append(');')
    query.append('out;')
    # So for the example of
    # {
    #     'amenity':True,
    #     'landuse':['retail','commercial'],
    #     'highway':'bus_stop'
    # }
    # the query string here is:
    # ['[out:json][timeout:180][bbox:south,west,north,east];',
    # '(',
    # '(node["amenity"];(._;>;););',
    # '(way["amenity"];(._;>;););',
    # '(relation["amenity"];(._;>;););',
    # '(node["landuse"="retail"];(._;>;););',
    # '(way["landuse"="retail"];(._;>;););',
    # '(relation["landuse"="retail"];(._;>;););',
    # '(node["landuse"="commercial"];(._;>;););',
    # '(way["landuse"="commercial"];(._;>;););',
    # '(relation["landuse"="commercial"];(._;>;););',
    # '(node["highway"="bus_stop"];(._;>;););',
    # '(way["highway"="bus_stop"];(._;>;););',
    # '(relation["highway"="bus_stop"];(._;>;););',
    # ');',
    # 'out;']
    return ''.join(query)


def osm_poi_download(polygon=None, amenities=None, tags=None, north=None,
                     south=None, east=None, west=None, timeout=180,
                     max_query_area_size=50*1000*50*1000):
    """
    Get points of interests (POIs) from OpenStreetMap based on selected amenity types.
    Note that if a polygon is passed-in, the query will be limited to its bounding box
    rather than to the shape of the polygon itself.

    Parameters
    ----------
    poly : shapely.geometry.Polygon
        Polygon that will be used to limit the POI search.
    north, south, east, west: float
        Bounding box used to limit the POI search.
    amenities : list
        List of amenities that will be used for finding the POIs from the
        selected area. Results returned are the _union_, not _intersection_ of
        each individual amenity tag; each result matches at least one tag given.
        Equivalent to `tags={'amenity': amenities}`. Some POIs don't use the
        `amenity` key, in which case you should use the `tags` option
    tags : dict[str, Union[bool, str, List[str]]]
        Dict of tags used for finding POIs from the selected area. Results
        returned are the _union_, not _intersection_ of each individual tag;
        each result matches at least one tag given.

        The keys should be the OSM tag key, (i.e. `amenity`, `landuse`,
        `highway`, etc) and the values should be `True` to retrieve all objects
        with the given tag, a string to get a single tag-value combination, and
        a list of strings to get multiple values for the given tag. For example,

            tags = {
                'amenity':True,
                'landuse':['retail','commercial'],
                'highway':'bus_stop'
            }

        Would return all amenities, `landuse=retail` and `landuse=commercial`,
        and `highway=bus_stop`. If both `amenities` and `tags` are combined, any
        `amenity` tags will be combined.
    timeout : int
        Timeout for the API request.

    Returns
    -------
    responses : dict
        JSON response with POIs from Overpass API server.
    """

    if polygon is not None:
        west, south, east, north = polygon.bounds

    # Parse the Overpass QL query
    query = create_overpass_query(amenities=amenities, tags=tags, west=west, south=south, east=east, north=north, maxsize=max_query_area_size)

    # Get the POIs
    responses = overpass_request(data={'data': query}, timeout=timeout)

    return responses


def parse_nodes_coords(osm_response):
    """
    Parse node coordinates from OSM response. Some nodes are standalone points
    of interest, others are vertices in polygonal (areal) POIs.

    Parameters
    ----------
    osm_response : string
        OSM response JSON string

    Returns
    -------
    coords : dict
        dict of node IDs and their lat, lon coordinates
    """

    coords = {}
    for result in osm_response['elements']:
        if 'type' in result and result['type'] == 'node':
            coords[result['id']] = {'lat': result['lat'],
                                    'lon': result['lon']}
    return coords


def parse_polygonal_poi(coords, response):
    """
    Parse areal POI way polygons from OSM node coords.

    Parameters
    ----------
    coords : dict
        dict of node IDs and their lat, lon coordinates

    Returns
    -------
    dict of POIs containing each's nodes, polygon geometry, and osmid
    """

    if 'type' in response and response['type'] == 'way':
        nodes = response['nodes']
        try:
            polygon = Polygon([(coords[node]['lon'], coords[node]['lat']) for node in nodes])

            poi = {'nodes': nodes,
                   'geometry': polygon,
                   'osmid': response['id']}

            if 'tags' in response:
                for tag in response['tags']:
                    poi[tag] = response['tags'][tag]
            return poi

        except Exception:
            log('Polygon has invalid geometry: {}'.format(nodes))

    return None


def parse_osm_node(response):
    """
    Parse points from OSM nodes.

    Parameters
    ----------
    response : JSON
        Nodes from OSM response.

    Returns
    -------
    Dict of vertex IDs and their lat, lon coordinates.
    """

    try:
        point = Point(response['lon'], response['lat'])

        poi = {
            'osmid': response['id'],
            'geometry': point
        }
        if 'tags' in response:
            for tag in response['tags']:
                poi[tag] = response['tags'][tag]

    except Exception:
        log('Point has invalid geometry: {}'.format(response['id']))

    return poi


def invalid_multipoly_handler(gdf, relation, way_ids):
    """
    Handles invalid multipolygon geometries when there exists e.g. a feature without
    geometry (geometry == NaN)

    Parameters
    ----------

    gdf : gpd.GeoDataFrame
        GeoDataFrame with Polygon geometries that should be converted into a MultiPolygon object.
    relation : dict
        OSM 'relation' dictionary
    way_ids : list
        A list of 'way' ids that should be converted into a MultiPolygon object.
    """

    try:
        gdf_clean = gdf.dropna(subset=['geometry'])
        multipoly = MultiPolygon(list(gdf_clean['geometry']))
        return multipoly

    except Exception:
        log("Invalid geometry at relation id %s.\nWay-ids of the invalid MultiPolygon:" % (
        relation['id'], str(way_ids)))
        return None


def parse_osm_relations(relations, osm_way_df):
    """
    Parses the osm relations (multipolygons) from osm
    ways and nodes. See more information about relations
    from OSM documentation: http://wiki.openstreetmap.org/wiki/Relation

    Parameters
    ----------
    relations : list
        OSM 'relation' items (dictionaries) in a list.
    osm_way_df : gpd.GeoDataFrame
        OSM 'way' features as a GeoDataFrame that contains all the
        'way' features that will constitute the multipolygon relations.

    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame with MultiPolygon representations of the
        relations and the attributes associated with them.
    """

    gdf_relations = gpd.GeoDataFrame()

    # Iterate over relations and extract the items
    for relation in relations:
        if relation['tags']['type'] == 'multipolygon':
            try:
                # Parse member 'way' ids
                member_way_ids = [member['ref'] for member in relation['members'] if member['type'] == 'way']
                # Extract the ways
                member_ways = osm_way_df.reindex(member_way_ids)
                # Extract the nodes of those ways
                member_nodes = list(member_ways['nodes'].values)
                try:
                    # Create MultiPolygon from geometries (exclude NaNs)
                    multipoly = MultiPolygon(list(member_ways['geometry']))
                except Exception:
                    multipoly = invalid_multipoly_handler(gdf=member_ways, relation=relation, way_ids=member_way_ids)

                if multipoly:
                    # Create GeoDataFrame with the tags and the MultiPolygon and its 'ways' (ids), and the 'nodes' of those ways
                    geo = gpd.GeoDataFrame(relation['tags'], index=[relation['id']])
                    # Initialize columns (needed for .loc inserts)
                    geo = geo.assign(geometry=None, ways=None, nodes=None, element_type=None, osmid=None)
                    # Add attributes
                    geo.loc[relation['id'], 'geometry'] = multipoly
                    geo.loc[relation['id'], 'ways'] = member_way_ids
                    geo.loc[relation['id'], 'nodes'] = member_nodes
                    geo.loc[relation['id'], 'element_type'] = 'relation'
                    geo.loc[relation['id'], 'osmid'] = relation['id']

                    # Append to relation GeoDataFrame
                    gdf_relations = gdf_relations.append(geo, sort=False)
                    # Remove such 'ways' from 'osm_way_df' that are part of the 'relation'
                    osm_way_df = osm_way_df.drop(member_way_ids)
            except Exception:
                log("Could not handle OSM 'relation': {}".format(relation['id']))

    # Merge 'osm_way_df' and the 'gdf_relations'
    osm_way_df = osm_way_df.append(gdf_relations, sort=False)
    return osm_way_df


def create_poi_gdf(polygon=None, amenities=None, tags=None, north=None, south=None, east=None, west=None):
    """
    Create GeoDataFrame from POI json returned by Overpass API.

    Parameters
    ----------
    polygon : shapely Polygon or MultiPolygon
        geographic shape to fetch the POIs within
    amenities: list
        List of amenities that will be used for finding the POIs from the selected area.
        See available amenities from: http://wiki.openstreetmap.org/wiki/Key:amenity

    north : float
        northern latitude of bounding box
    south : float
        southern latitude of bounding box
    east : float
        eastern longitude of bounding box
    west : float
        western longitude of bounding box

    Returns
    -------

    gdf : geopandas.GeoDataFrame
        Points of interest and the tags associated with them as geopandas GeoDataFrame.
    """

    responses = osm_poi_download(polygon=polygon, amenities=amenities, tags=tags, north=north, south=south, east=east, west=west)

    # Parse coordinates from all the nodes in the response
    coords = parse_nodes_coords(responses)

    # POI nodes
    poi_nodes = {}

    # POI ways
    poi_ways = {}

    # A list of POI relations
    relations = []

    for result in responses['elements']:
        if result['type'] == 'node' and 'tags' in result:
            poi = parse_osm_node(response=result)
            # Add element_type
            poi['element_type'] = 'node'
            # Add to 'pois'
            poi_nodes[result['id']] = poi
        elif result['type'] == 'way':
            # Parse POI area Polygon
            poi_area = parse_polygonal_poi(coords=coords, response=result)
            if poi_area:
                # Add element_type
                poi_area['element_type'] = 'way'
                # Add to 'poi_ways'
                poi_ways[result['id']] = poi_area

        elif result['type'] == 'relation':
            # Add relation to a relation list (needs to be parsed after all nodes and ways have been parsed)
            relations.append(result)

    # Create GeoDataFrames
    gdf_nodes = gpd.GeoDataFrame(poi_nodes).T
    gdf_nodes.crs = settings.default_crs

    gdf_ways = gpd.GeoDataFrame(poi_ways).T
    gdf_ways.crs = settings.default_crs

    # Parse relations (MultiPolygons) from 'ways'
    gdf_ways = parse_osm_relations(relations=relations, osm_way_df=gdf_ways)

    # Combine GeoDataFrames
    gdf = gdf_nodes.append(gdf_ways, sort=False)

    if polygon:
        gdf = gdf.loc[gdf['geometry'].centroid.within(polygon)==True]

    return gdf


def pois_from_point(point, distance=None, amenities=None):
    """
    Get point of interests (POIs) within some distance north, south, east, and west of
    a lat-long point.

    Parameters
    ----------
    point : tuple
        a lat-long point
    distance : numeric
        distance in meters
    amenities : list
        List of amenities that will be used for finding the POIs from the selected area.
        See available amenities from: http://wiki.openstreetmap.org/wiki/Key:amenity

    Returns
    -------
    GeoDataFrame
    """

    bbox = bbox_from_point(point=point, distance=distance)
    north, south, east, west = bbox
    return create_poi_gdf(amenities=amenities, north=north, south=south, east=east, west=west)


def pois_from_address(address, distance, amenities=None):
    """
    Get OSM points of Interests within some distance north, south, east, and west of
    an address.

    Parameters
    ----------
    address : string
        the address to geocode to a lat-long point
    distance : numeric
        distance in meters
    amenities : list
        List of amenities that will be used for finding the POIs from the selected area. See available
        amenities from: http://wiki.openstreetmap.org/wiki/Key:amenity

    Returns
    -------
    GeoDataFrame
    """

    # geocode the address string to a (lat, lon) point
    point = geocode(query=address)

    # get POIs within distance of this point
    return pois_from_point(point=point, amenities=amenities, distance=distance)


def pois_from_polygon(polygon, amenities=None):
    """
    Get OSM points of interest within some polygon.

    Parameters
    ----------
    polygon : Polygon
        Polygon where the POIs are search from.
    amenities : list
        List of amenities that will be used for finding the POIs from the selected area.
        See available amenities from: http://wiki.openstreetmap.org/wiki/Key:amenity

    Returns
    -------
    GeoDataFrame
    """

    return create_poi_gdf(polygon=polygon, amenities=amenities)


def pois_from_place(place, amenities=None, which_result=1):
    """
    Get points of interest (POIs) within the boundaries of some place.

    Parameters
    ----------
    place : string
        the query to geocode to get geojson boundary polygon.
    amenities : list
        List of amenities that will be used for finding the POIs from the selected area.
        See available amenities from: http://wiki.openstreetmap.org/wiki/Key:amenity
    which_result : int
        max number of results to return and which to process upon receipt

    Returns
    -------
    GeoDataFrame
    """

    city = gdf_from_place(place, which_result=which_result)
    polygon = city['geometry'].iloc[0]
    return create_poi_gdf(polygon=polygon, amenities=amenities)
