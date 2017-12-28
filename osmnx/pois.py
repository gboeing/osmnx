################################################################################
# Module: pois.py
# Description: Download and plot Points of Interests (POIs) from OpenStreetMap
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

from shapely.geometry import Point, LineString, MultiPolygon, Polygon, box
import geopandas as gpd
from .core import overpass_request, bbox_from_point, gdf_from_place
from .utils import log, geocode, bbox_to_poly

"""
# List of available amenities
available_amenities = ['fast_food', 'restaurant', 'food_court', 'pub', 'bar', 'nighclub', 'biergarten', 'cafe', 'bakery',
                       'hairdresser', 'beauty_shop', 'supermarket', 'kiosk', 'beverages', 'greengrocer', 'butcher', 'convenience', 'department_store', 'computer_shop', 'clothes', 'sports_shop', 'bicycle_shop', 'video_shop', 'furniture_shop', 'outdoor_shop', 'shoe_shop', 'bookshop', 'jeweller', 'gift_shop', 'mobile_phone_shop',
                       'car_dealership', 'chemist', 'florist', 'garden_centre', 'mall', 'sports_shop', 'toy_shop', 'vending_any', 'vending_machine',
                       'dentist', 'pharmacy', 'optician', 'hospital', 'doctors', 'nursing_home', 'veterinary',
                       'school', 'kindergarten', 'college', 'university',
                       'arts_centre', 'archaeological', 'artwork', 'attraction', 'battlefield', 'castle', 'fort', 'fountain', 'lighthouse', 'memorial', 'monument', 'museum', 'observation_tower', 'ruins', 'tower','town_hall', 'viewpoint', 'windmill',
                       'bank', 'atm', 'bicycle_rental', 'car_rental', 'car_sharing',
                       'camp_site', 'hotel', 'hostel', 'alpine_hut', 'caravan_site', 'chalet', 'guest_house',
                       'car_repair', 'car_wash', 'laundry', 'library', 'tourist_info', 'travel_agent',
                       'cinema', 'theatre', 'theme_park', 'zoo', 'dog_park', 'golf_course', 'hunting_stand', 'ice_rink', 'park', 'picnic_site', 'pitch', 'playground', 'sports_centre', 'stadium', 'swimming_pool', 'track', 'community_centre', 'courthouse', 'embassy', 'doityourself',
                       'fire_station', 'graveyard', 'police', 'post_box', 'post_office', 'prison', 'public_building', 'shelter', 'telephone',
                       'toilet', 'vending_parking', 'waste_basket', 'wastewater_plant', 'water_mill', 'water_tower', 'water_well', 'water_works', 'wayside_cross', 'wayside_shrine', 'bench', 'drinking_water',
                       'recycling', 'recycling_clothes', 'recycling_glass', 'recycling_metal', 'recycling_paper']
"""
def parse_poi_query(west, south, east, north, amenities=None, timeout=180, maxsize=''):
    """
    Parse the Overpass QL query based on the list of amenities.

    Parameters
    ----------
    west : float
        Westernmost coordinate of the bounding box of the search area.
    south : float
        Southernmost coordinate from bounding box of the search area.
    east : float
        Easternmost coordinate from bounding box of the search area.
    north : float
        Northernmost coordinate from bounding box of the search area.
    amenities : list
        List of amenities that will be used for finding the POIs from the selected area.
    timeout : int
        Timeout for the API request. 
    """
    if amenities:
        # Overpass QL template
        query_template = ('[out:json][timeout:{timeout}]{maxsize};((node["amenity"~"{amenities}"]({south:.8f},'
                          '{west:.8f},{north:.8f},{east:.8f});(._;>;););(way["amenity"~"{amenities}"]({south:.8f},'
                          '{west:.8f},{north:.8f},{east:.8f});(._;>;););(relation["amenity"~"{amenities}"]'
                          '({south:.8f},{west:.8f},{north:.8f},{east:.8f});(._;>;);));out;')

        # Parse amenties
        query_str = query_template.format(amenities="|".join(amenities), north=north, south=south, east=east, west=west,
                                          timeout=timeout, maxsize=maxsize)
    else:
        # Overpass QL template
        query_template = ('[out:json][timeout:{timeout}]{maxsize};((node["amenity"]({south:.8f},'
                          '{west:.8f},{north:.8f},{east:.8f});(._;>;););(way["amenity"]({south:.8f},'
                          '{west:.8f},{north:.8f},{east:.8f});(._;>;););(relation["amenity"]'
                          '({south:.8f},{west:.8f},{north:.8f},{east:.8f});(._;>;);));out;')

        # Parse amenties
        query_str = query_template.format(north=north, south=south, east=east, west=west,
                                          timeout=timeout, maxsize=maxsize)

    return query_str

def osm_poi_download(polygon=None, amenities=None, north=None, south=None, east=None, west=None,
                     timeout=180, max_query_area_size=50*1000*50*1000):
    """
    Get Point of interests (POIs) from OpenStreetMap based on selected amenity types.

    Parameters
    ----------
    poly : shapely.geometry.Polygon
        Polygon that will be used to limit the POI search. 
    amenities : list
        List of amenities that will be used for finding the POIs from the selected area.

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        Points of interests and the tags associated to them as Geopandas GeoDataFrame.
    """

    if polygon:
        # Bounds
        west, south, east, north = polygon.bounds

        # Parse the Overpass QL query
        query = parse_poi_query(amenities=amenities, west=west, south=south, east=east, north=north)

    elif not (north is None or south is None or east is None or west is None):
        # TODO: Add functionality for subdividing search area geometry based on max_query_area_size
        # Parse Polygon from bbox
        #polygon = bbox_to_poly(north=north, south=south, east=east, west=west)

        # Parse the Overpass QL query
        query = parse_poi_query(amenities=amenities, west=west, south=south, east=east, north=north)

    else:
        raise ValueError('You must pass a polygon or north, south, east, and west')

    # Get the POIs
    responses = overpass_request(data={'data': query}, timeout=timeout)

    return responses

def parse_vertice_nodes(osm_response):
    """
    Parse node vertices from OSM response.
     
    Parameters
    ==========
    osm_response : JSON 
        OSM response JSON 
    
    Returns
    =======
    Dict of vertice IDs and their lat, lon coordinates.
    """
    vertices = {}
    for result in osm_response['elements']:
        if 'type' in result and result['type'] == 'node':
            vertices[result['id']] = {'lat': result['lat'],
                                      'lon': result['lon']}
    return vertices

def parse_osm_way(vertices, response):
    """
    Parse ways (areas) from OSM node vertices.

    Parameters
    ==========
    vertices : Python dict 
        Node vertices parsed from OSM response.  

    Returns
    =======
    Dict of vertice IDs and their lat, lon coordinates.
    """
    if 'type' in response and response['type'] == 'way':
        nodes = response['nodes']
        try:
            polygon = Polygon([(vertices[node]['lon'], vertices[node]['lat']) for node in nodes])

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
    ==========
    response : JSON 
        Nodes from OSM response.  

    Returns
    =======
    Dict of vertice IDs and their lat, lon coordinates.
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
    Handles invalid multipolygon geometries when there exists e.g. a feature without geometry (geometry == NaN)

    Parameters
    ==========

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
    Parses the osm relations from osm ways and nodes.
         
    Parameters
    ==========
    relations : list
        OSM 'relation' items (dictionaries) in a list. 
    osm_way_df : gpd.GeoDataFrame
        OSM 'way' features as a GeoDataFrame that contains all the 'way' features that will constitute the multipolygon relations.
     
    Information
    ===========
    http://wiki.openstreetmap.org/wiki/Relation 
    
    Returns
    =======
    gpd.GeoDataFrame
        A GeoDataFrame with MultiPolygon representations of the relations and the attributes associated with them.   
    """
    gdf_relations = gpd.GeoDataFrame()

    # Iterate over relations and extract the items
    for relation in relations:
        if relation['tags']['type'] == 'multipolygon':
            try:
                # Parse member 'way' ids
                member_way_ids = [member['ref'] for member in relation['members'] if member['type'] == 'way']
                # Extract the ways
                member_ways = osm_way_df.loc[member_way_ids]
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
                    gdf_relations = gdf_relations.append(geo)
                    # Remove such 'ways' from 'osm_way_df' that are part of the 'relation'
                    osm_way_df = osm_way_df.drop(member_way_ids)
            except Exception:
                log("Could not handle OSM 'relation': {}".format(relation['id']))

    # Merge 'osm_way_df' and the 'gdf_relations'
    osm_way_df = osm_way_df.append(gdf_relations)
    return osm_way_df

def create_poi_gdf(polygon=None, amenities=None, north=None, south=None, east=None, west=None, retain_invalid=False):
    """
    Parse GeoDataFrames from POI json that was returned by Overpass API.

    Parameters
    ----------
    polygon : shapely Polygon or MultiPolygon
        geographic shape to fetch the building footprints within
    amenities: list
        List of amenities that will be used for finding the POIs from the selected area. See available amenities from: 
    north : float
        northern latitude of bounding box
    south : float
        southern latitude of bounding box
    east : float
        eastern longitude of bounding box
    west : float
        western longitude of bounding box
    retain_invalid : bool
        if False discard any building footprints with an invalid geometry
        
    Returns
    -------
    Geopandas GeoDataFrame with POIs and the associated attributes. 
    """
    responses = osm_poi_download(polygon=polygon, amenities=amenities, north=north, south=south, east=east, west=west)

    # Parse vertices
    vertices = parse_vertice_nodes(responses)

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
            poi_area = parse_osm_way(vertices=vertices, response=result)
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
    gdf_nodes.crs = {'init': 'epsg:4326'}

    gdf_ways = gpd.GeoDataFrame(poi_ways).T
    gdf_ways.crs = {'init': 'epsg:4326'}

    # Parse relations (MultiPolygons) from 'ways'
    gdf_ways = parse_osm_relations(relations=relations, osm_way_df=gdf_ways)

    # Combine GeoDataFrames
    gdf = gdf_nodes.append(gdf_ways)

    return gdf

def pois_from_point(point, distance=None, amenities=None, retain_invalid=False):
    """
    Get Point of interests (POIs) within some distance north, south, east, and west of
    a lat-long point.

    Parameters
    ----------
    point : tuple
        a lat-long point
    amenities : list
        List of amenities that will be used for finding the POIs from the selected area.
    distance : numeric
        distance in meters
    retain_invalid : bool
        if False discard any building footprints with an invalid geometry

    Returns
    -------
    GeoDataFrame 
    """
    bbox = bbox_from_point(point=point, distance=distance)
    north, south, east, west = bbox
    return create_poi_gdf(amenities=amenities, north=north, south=south, east=east, west=west)

def pois_from_address(address, distance, amenities=None, retain_invalid=False):
    """
    Get OSM Points of Interests within some distance north, south, east, and west of
    an address.

    Parameters
    ----------
    address : string
        the address to geocode to a lat-long point
    amenities : list
        List of amenities that will be used for finding the POIs from the selected area.
    distance : numeric
        distance in meters
    retain_invalid : bool
        if False discard any building footprints with an invalid geometry

    Returns
    -------
    GeoDataFrame
    """

    # geocode the address string to a (lat, lon) point
    point = geocode(query=address)

    # get buildings within distance of this point
    return pois_from_point(point=point, amenities=amenities, distance=distance)

def pois_from_polygon(polygon, amenities=None, retain_invalid=False):
    """
    Get OSM Points of Interests within some polygon.

    Parameters
    ----------
    polygon : Polygon
        Polygon where the POIs are search from. 
    amenities : list
        List of amenities that will be used for finding the POIs from the selected area.
    retain_invalid : bool
        if False discard any building footprints with an invalid geometry

    Returns
    -------
    GeoDataFrame
    """

    return create_poi_gdf(polygon=polygon, amenities=amenities)

def pois_from_place(place, amenities=None, retain_invalid=False):
    """
    Get building footprints within the boundaries of some place.

    Parameters
    ----------
    place : string
        the query to geocode to get geojson boundary polygon.
    amenities : list
        List of amenities that will be used for finding the POIs from the selected area.
    retain_invalid : bool
        if False discard any building footprints with an invalid geometry

    Returns
    -------
    GeoDataFrame
    """

    city = gdf_from_place(place)
    polygon = city['geometry'].iloc[0]
    return create_poi_gdf(polygon=polygon, amenities=amenities, retain_invalid=retain_invalid)

