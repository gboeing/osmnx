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
from .utils import bbox_to_poly
from .utils import geocode
from .utils import log

def parse_poi_query(north, south, east, west, tags=None, timeout=180, maxsize=''):
    """
    Construct an Overpass QL query to load POIs with certain tags.

    By default, queries all features with an amenity tag.

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
    tags : dict
        Dictionary of tag keys and values that will be used for finding POIs in the selected area.
        Keys may be strings or lists of strings.
        Values make be string, lists of strings, or None, if all values should be returned for a given key.
        By default, all POIs with an 'amenity' key of any value will be be returned.
    timeout : int
        Timeout for the API request.
    """

    # build default tags
    if not tags:
        tags = {'amenity':True}
    
    # define templates for objects and extents    
    object_template = '({object_type}[{{keys}}{{op}}"{{values}}"]{{extent}});'
    # object_template = '({object_type}[~"^({{keys}})$"{{op}}"{{values}}"]{{extent}});'
    
    re_keys_template = '~"^({keys})$"'
    single_key_template = '"{key}"'

    extent_template = '({south:.6f},{west:.6f},{north:.6f},{east:.6f});(._;>;);'
    extent = extent_template.format(south=south, west=west, north=north, east=east)
    
    # initate query string
    query_template = "[out:json][timeout:{timeout}]{maxsize};("
    query_str = query_template.format(timeout=timeout, maxsize=maxsize)
    
    # add statements for each object type
    # templates = [object_template.format(object_type=x) for x in ['node','way','relation']]
    templates = [object_template.format(object_type=x) for x in ['nwr']]

    for template in templates:
       
        # add statements for each key
        for keys, values in tags.items():

            # ensure keys is a list
            keys = [keys] if not isinstance(keys, list) else keys
           
            if values == True:
                # get features with any value for these keys
                # add positive statement with multiple keys and no specific values
                query_str += template.format(keys=re_keys_template.format(keys='|'.join(keys)), values='.*', extent=extent, op='~')

            elif values == False:
                # get features wihout these keys, not matter their values
                for key in keys:
                    # add negative statement with multiple keys and no specific values
                    # can only be added one at a time withough key regex
                    query_str += template.format(keys=single_key_template.format(key=key), values='.*', extent=extent, op='!~')

            else: 
                # get features with specified values for these keys
                # ensure values is a list
                values = [values] if not isinstance(values, list) else values
                # add positive statement with multiple keys in specific values
                query_str += template.format(keys='{}'.format('|'.join(keys)), values='|'.join(values), extent=extent, op='~')
                           
    # terminate query string
    query_str += ");out;"
   
    return query_str


def osm_poi_download(polygon=None, north=None, south=None, east=None, west=None,
                     timeout=180, max_query_area_size=50*1000*50*1000, tags=None,
                     return_query=False):
    """
    Get points of interests (POIs) from OpenStreetMap based on selected amenity types.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon
        Polygon that will be used to limit the POI search.
    north : float
    	Northern WGS84 extent (lat)
	south : float
		Southern WGS84 extent (lat)
	east : float
		Eastern WGS84 extent (lon)
	west : float
		Western WGS84 extent (lon)
    tags : dict
        Dictionary of tag keys and values that will be used for finding POIs in the selected area.
        Keys may be strings or lists of strings.
        Values make be string, lists of strings, or None, if all values should be returned for a given key.
        By default, all POIs with an 'amenity' key of any value will be be returned.
    return_query : bool (default=False)
    	If True, Overpass API query based on tags dictionary is returned as the second
    	object in a tuple

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        Points of interest and the tags associated with them as geopandas GeoDataFrame.
    query : str (optional; returned only if `return_query=True`)
    	Overpass API query constructed based on tags dictionary input
    """

    def set_query(west, south, east, north, tags):
        if isinstance(tags, dict):
            return parse_poi_query(west=west, south=south, east=east, north=north, tags=tags)
        elif isinstance(tags, str):
            return tags
        else:
            raise ValueError('tags must be defined with a dictionary or a pre-formatted query string')

    if polygon:
        # Bounds
        west, south, east, north = polygon.bounds

        # Parse the Overpass QL query
        query = set_query(west=west, south=south, east=east, north=north, tags=tags)
        
    elif not (north is None or south is None or east is None or west is None):
        # TODO: Add functionality for subdividing search area geometry based on max_query_area_size
        # Parse Polygon from bbox
        #polygon = bbox_to_poly(north=north, south=south, east=east, west=west)

        # Parse the Overpass QL query
        query = set_query(west=west, south=south, east=east, north=north, tags=tags)

    else:
        raise ValueError('You must pass a polygon or north, south, east, and west')

    # Get the POIs
    responses = overpass_request(data={'data': query}, timeout=timeout)

    if return_query:
        return responses, query
    else:
        return responses


def parse_nodes_coords(osm_response):
    """
    Parse node coordinates from OSM response. Some nodes are
    standalone points of interest, others are vertices in 
    polygonal (areal) POIs.
    
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


def create_poi_gdf(polygon=None, north=None, south=None, east=None, west=None, tags=None, return_query=False):
    """
    Parse GeoDataFrames from POI json that was returned by Overpass API.

    Parameters
    ----------
    polygon : shapely Polygon or MultiPolygon
        geographic shape to fetch the building footprints within
    north : float
        northern latitude of bounding box
    south : float
        southern latitude of bounding box
    east : float
        eastern longitude of bounding box
    west : float
        western longitude of bounding box
    tags : dict
        Dictionary of tag keys and values that will be used for finding POIs in the selected area.
        Keys may be strings or lists of strings.
        Values make be string, lists of strings, or None, if all values should be returned for a given key.
        By default, all POIs with an 'amenity' key of any value will be be returned.
        
    Returns
    -------
    Geopandas GeoDataFrame with POIs and the associated attributes. 
    """

    if return_query:
        responses, query = osm_poi_download(polygon=polygon, north=north, south=south, east=east, west=west, tags=tags, return_query=return_query)
    else:
        responses = osm_poi_download(polygon=polygon, north=north, south=south, east=east, west=west, tags=tags)

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

    if return_query:
        return gdf, query
    else:
        return gdf


def pois_from_point(point, distance=None, tags=None, return_query=False, **kwargs):
    """
    Get point of interests (POIs) within some distance north, south, east, and west of
    a lat-long point.

    Parameters
    ----------
    point : tuple
        a lat-long point
    distance : numeric
        distance in meters
    tags : dict
        Dictionary of tag keys and values that will be used for finding POIs in the selected area.
        Keys may be strings or lists of strings.
        Values make be string, lists of strings, or None, if all values should be returned for a given key.
        By default, all POIs with an 'amenity' key of any value will be be returned.

    Keyword Arguments
    -----------------
    amenities : list
        If provided, will override tags with the dictionary {'amenities':amenities}.
        Accommodates the depricatd amenities parameter.

    Returns
    -------
    GeoDataFrame 
    """
    if 'amenities' in locals():
        tags = {'amenity':amenities}
    if not tags:
        tags = {'amenity':True}

    bbox = bbox_from_point(point=point, distance=distance)
    north, south, east, west = bbox
    return create_poi_gdf(north=north, south=south, east=east, west=west, tags=tags)


def pois_from_address(address, distance, tags=None, return_query=False, **kwargs):
    """
    Get OSM points of Interests within some distance north, south, east, and west of
    an address.

    Parameters
    ----------
    address : string
        the address to geocode to a lat-long point
    distance : numeric
        distance in meters
    tags : dict
        Dictionary of tag keys and values that will be used for finding POIs in the selected area.
        Keys may be strings or lists of strings.
        Values make be string, lists of strings, or None, if all values should be returned for a given key.
        By default, all POIs with an 'amenity' key of any value will be be returned.

    Keyword Arguments
    -----------------
    amenities : list
        If provided, will override tags with the dictionary {'amenities':amenities}.
        Accommodates the depricatd amenities parameter.

    Returns
    -------
    GeoDataFrame
    """
    if 'amenities' in locals():
        tags = {'amenity':amenities}

    # geocode the address string to a (lat, lon) point
    point = geocode(query=address)

    # get buildings within distance of this point
    return pois_from_point(point=point, distance=distance, tags=tags, return_query=return_query)


def pois_from_polygon(polygon, tags=None, return_query=False, **kwargs):
    """
    Get OSM points of interest within some polygon.

    Parameters
    ----------
    polygon : Polygon
        Polygon where the POIs are search from. 
    tags : dict
        Dictionary of tag keys and values that will be used for finding POIs in the selected area.
        Keys may be strings or lists of strings.
        Values make be string, lists of strings, or None, if all values should be returned for a given key.
        By default, all POIs with an 'amenity' key of any value will be be returned.

    Keyword Arguments
    -----------------
    amenities : list
        If provided, will override tags with the dictionary {'amenities':amenities}.
        Accommodates the depricatd amenities parameter.

    Returns
    -------
    GeoDataFrame
    """
    if 'amenities' in locals():
        tags = {'amenity':amenities}
    if not tags:
        tags = {'amenity':True}
    return create_poi_gdf(polygon=polygon, tags=tags, return_query=return_query)


def pois_from_place(place, tags=None, return_query=False, **kwargs):
    """
    Get points of interest (POIs) within the boundaries of some place.

    Parameters
    ----------
    place : string
        the query to geocode to get geojson boundary polygon.
    tags : dict
        Dictionary of tag keys and values that will be used for finding POIs in the selected area.
        Keys may be strings or lists of strings.
        Values make be string, lists of strings, or None, if all values should be returned for a given key.
        By default, all POIs with an 'amenity' key of any value will be be returned.
    return_query : bool (default=False)
    	If True, Overpass API query based on tags dictionary is returned as the second
    	object in a tuple

    Keyword Arguments
    -----------------
    amenities : list
        If provided, will override tags with the dictionary {'amenities':amenities}.
        Accommodates the depricatd amenities parameter.
    
    Returns
    -------
    gdf : geopandas.GeoDataFrame
        Points of interest and the tags associated with them as geopandas GeoDataFrame.
    query : str (optional; returned only if `return_query=True`)
    	Overpass API query constructed based on tags dictionary input
    """
    if 'amenities' in locals():
        tags = {'amenity':amenities}
    if not tags:
        tags = {'amenity':True}

    city = gdf_from_place(place)
    polygon = city['geometry'].iloc[0]
    return create_poi_gdf(polygon=polygon, tags=tags, return_query=return_query)
