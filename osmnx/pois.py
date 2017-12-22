################################################################################
# Module: buildings.py
# Description: Download and plot Points of Interests (POIs) from OpenStreetMap
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

from shapely.geometry import Point, LineString, box
import geopandas as gpd
from .core import overpass_request, bbox_from_point, gdf_from_place
from .utils import log, geocode, bbox_to_poly

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

def parse_poi_query(amenities, west, south, east, north, timeout=180, maxsize=''):
    """
    Parse the Overpass QL query based on the list of amenities.

    Parameters
    ----------
    amenities : list
        List of amenities that will be used for finding the POIs from the selected area.
    west : float
        Westernmost coordinate of the bounding box of the search area.
    south : float
        Southernmost coordinate from bounding box of the search area.
    east : float
        Easternmost coordinate from bounding box of the search area.
    north : float
        Northernmost coordinate from bounding box of the search area.
    timeout : int
        Timeout for the API request. 
    """

    # Overpass QL template
    query_template = ('[out:json][timeout:{timeout}]{maxsize};((node["amenity"~"{amenities}"]({south:.8f},'
                      '{west:.8f},{north:.8f},{east:.8f});(._;>;););(way["amenity"~"{amenities}"]({south:.8f},'
                      '{west:.8f},{north:.8f},{east:.8f});(._;>;););(relation["amenity"~"{amenities}"]'
                      '({south:.8f},{west:.8f},{north:.8f},{east:.8f});(._;>;);));out;')

    # Parse amenties
    query_str = query_template.format(amenities="|".join(amenities), north=north, south=south, east=east, west=west,
                                      timeout=timeout, maxsize=maxsize)
    return query_str

def validate_amenities(amenities=None):
    """
    Validate the list of amenities. Warn users about unrecognized amenities. 
    """
    # Count valid amenities
    valid_cnt = 0
    for amenity in amenities:
        # Check known amenities
        if amenity in available_amenities:
            valid_cnt += 1
        else:
            Warning("Amenity {} was not recognized.".format(amenity))
    if valid_cnt == 0:
        raise Warning("There were not any recognized amenities. You might not get any results. Check available amenities by running: ox.pois.available")

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

    # Validate amenities
    validate_amenities(amenities)

    if polygon:
        # Bounds
        west, south, east, north = polygon.bounds

        # Parse the Overpass QL query
        query = parse_poi_query(amenities, west, south, east, north)

    elif not (north is None or south is None or east is None or west is None):
        # TODO: Add functionality for subdividing search area geometry based on max_query_area_size
        # Parse Polygon from bbox
        #polygon = bbox_to_poly(north=north, south=south, east=east, west=west)

        # Parse the Overpass QL query
        query = parse_poi_query(amenities, west, south, east, north)

    else:
        raise ValueError('You must pass a polygon or north, south, east, and west')

    # Get the POIs
    responses = overpass_request(data={'data': query}, timeout=timeout)

    return responses

def create_poi_gdf(polygon=None, amenities=None, north=None, south=None, east=None, west=None, retain_invalid=False):
    """
    Parse GeoDataFrames from POI json that was returned by Overpass API.

    Parameters
    ----------
    poi_json : json
        Overpass API JSON object.

    Returns
    -------
    Geopandas GeoDataFrame with POIs and the associated attributes. 
    """
    responses = osm_poi_download(polygon=polygon, amenities=amenities, north=north, south=south, east=east, west=west)

    # Parse POIs
    pois = {}
    for result in responses['elements']:
        if result['type'] == 'node' and 'tags' in result:
            try:
                point = Point(result['lon'], result['lat'])

                poi = {
                    'osmid': result['id'],
                    'geometry': point
                }
                if 'tags' in result:
                    for tag in result['tags']:
                        poi[tag] = result['tags'][tag]
                pois[result['id']] = poi

            except Exception:
                log('Point has invalid geometry: {}'.format(result['id']))

        elif result['type'] == 'relation':
            # TODO: Add functionalities to parse 'relation' tags.
            pass
        elif result['type'] == 'way':
            # TODO: Add functionalities to parse 'way' tags as a separate GeoDataFrame.
            pass

    gdf = gpd.GeoDataFrame(pois).T
    gdf.crs = {'init': 'epsg:4326'}
    return gdf

def pois_from_point(point, amenities, distance=None, retain_invalid=False):
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

def pois_from_address(address, amenities, distance, retain_invalid=False):
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

def pois_from_polygon(polygon, amenities, retain_invalid=False):
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

def pois_from_place(place, amenities, retain_invalid=False):
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
    return create_poi_gdf(polygon, amenities, retain_invalid=retain_invalid)
