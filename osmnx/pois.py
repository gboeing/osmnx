################################################################################
# Module: pois.py
# Description: Download and plot points of interests (POIs) from OpenStreetMap
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

import geopandas as gpd
from shapely.geometry import box
from shapely.geometry import LineString
from shapely.geometry import MultiPolygon
from shapely.geometry import Point
from shapely.geometry import Polygon

from . import settings
from .core import bbox_from_point
from .core import gdf_from_place
from .core import overpass_request
from .geo_utils import geocode, bbox_to_poly
from .utils import log


def parse_poi_query(north, south, east, west, amenities=None, tags=None,
                    timeout=180, maxsize=None):
    """
    Parse the Overpass QL query based on the list of amenities.

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
        `amenity` tags will be combined.
    timeout : int
        Timeout for the API request.
    """
    # Convert input into dict of tags
    _tags = {}
    if amenities is not None:
        # Instatiate amenity key
        _tags['amenity'] = _tags.get('amenity', [])
        # Add amenities to _tags dict
        [_tags['amenity'].append(a) for a in amenities]

    if tags is not None:
        for k, v in tags.items():
            if type(v) == bool:
                if v:
                    # Use None to signal later that the query should be
                    # ["tag_key"]
                    # which gets all tags with key `tag_key`
                    _tags[k] = None

                continue

            if type(v) == list:
                # Instatiate key
                _tags[k] = _tags.get(k, [])
                [_tags[k].append(item) for item in v]

                continue

            if type(v) == str:
                # Instatiate key
                _tags[k] = _tags.get(k, [])
                _tags[k].append(v)

    # List of dictionaries with one key and one value.
    # Overpass QL doesn't have an OR, so you can't do "amenity"="a|b"
    # Instead you have to do "amenity"="a" and then also "amenity"="b"
    _tags_list = []
    for k, v in _tags.items():
        if v is None:
            _tags_list.append({k: None})
            continue

        for item in v:
            _tags_list.append({k: item})

    # Create query template
    # Build query for each tag for each of nodes, ways, and relations.

    # Set query settings: Query settings must all be on the same line. The body
    # of the query can be over separate lines. See documentation here:
    # https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL#Settings
    query_settings = [
        '[out:json]',
        '[timeout:{t}]'.format(t=timeout),
        '[bbox:{s:.6f},{w:.6f},{n:.6f},{e:.6f}]'.format(s=south, w=west,
                                                        n=north, e=east)
    ]
    if maxsize is not None:
        query_settings.append('[maxsize:{m}]'.format(m=maxsize))

    query = [
        ''.join(query_settings) + ';',
        '('
    ]

    for d in _tags_list:
        for k, v in d.items():
            tag_str = '["{k}"'.format(k=k)
            if v is not None:
                tag_str += '="{v}"'.format(v=v)
            # The `(._; >;);` joins the new result set into the previous result
            # set: https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL#Item
            tag_str += '];(._;>;);'

            for t in ['node', 'way', 'relation']:
                query.append('({t}{tag_str});'.format(t=t, tag_str=tag_str))

    query.append(');')
    query.append('out;')
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

    elif not all(isinstance(x, float) for x in [north, south, east, west]):
        msg = 'You must pass a polygon or north, south, east, and west'
        raise ValueError(msg)

    # Parse the Overpass QL query
    query = parse_poi_query(amenities=amenities, tags=tags, west=west, south=south, east=east, north=north, maxsize=max_query_area_size)

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
