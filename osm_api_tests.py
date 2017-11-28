import sys, urllib.parse, http.client, json, logging

logging.basicConfig(level=logging.INFO, stream=sys.stdout, format='%(levelname)09s - %(message)s')

def get_nominatim_response(query, include_geojson):
    '''
    '''
    params = dict(q=query, format='json', dedupe=0, limit=1)
    
    if include_geojson:
        params.update(polygon_geojson=1)
    
    conn = http.client.HTTPConnection('nominatim.openstreetmap.org')
    conn.request('GET', '/search?'+urllib.parse.urlencode(params))
    
    return conn.getresponse()

def is_good_nominatim_response(body, expect_geojson):
    ''' Return True if JSON body contains expected Nominatim response.
    
        Sample response with GeoJSON included:
        
        [ {
            "boundingbox": [ "40.7702778", "40.8102778", "-73.9797222", "-73.9397222" ],
            "display_name": "Manhattan, NYC, New York, 83, United States of America",
            "geojson": { "coordinates": [ -73.9597222, 40.7902778 ], "type": "Point" },
            ...
        }, ... ]    
        
        Sample response with GeoJSON excluded:
        
        [ {
            "boundingbox": [ "40.6717833", "40.6725822", "-73.9852274", "-73.9840533" ],
            "display_name": "350, 5th Avenue, Park Slope, BK, Kings County, NYC, New York, 11215, United States of America",
            "lat": "40.67218135", "lon": "-73.9844516483024",
            ...
        }, ... ]    
    '''
    results = json.loads(body)
    
    if type(results) is not list:
        return False
    
    for result in results:
        if 'boundingbox' not in result:
            logging.warning('No boundingbox in Nominatim result')
            return False
        
        if type(result['boundingbox']) is not list:
            logging.warning('Bad boundingbox in Nominatim result')
            return False
        
        if len(result['boundingbox']) is not 4:
            logging.warning('Bad boundingbox in Nominatim result')
            return False
        
        if 'display_name' not in result:
            logging.warning('No display_name in Nominatim result')
            return False
        
        if expect_geojson:
            if 'geojson' not in result:
                logging.warning('No geojson in Nominatim result')
                return False
        
            geojson = result['geojson']
        
            if 'type' not in geojson:
                logging.warning('No type in Nominatim result GeoJSON')
                return False
        
            if geojson['type'] not in ('Point', 'Polygon', 'LineString',
                'MultiPoint', 'MultiPolygon', 'MultiLineString'):
                logging.warning('Unknown type in Nominatim result GeoJSON')
                return False
        
            if 'coordinates' not in geojson:
                logging.warning('No coordinates in Nominatim result GeoJSON')
                return False
        
            if type(geojson['coordinates']) is not list:
                logging.warning('Bad coordinates in Nominatim result GeoJSON')
                return False
        
        else:
            if 'lat' not in result or 'lon' not in result:
                logging.warning('No lon, lat in Nominatim result')
                return False
    
    return True

def main():
    '''
    '''
    logging.info(f'Querying Nominatim for Manhattan')
    resp1 = get_nominatim_response('Manhattan, New York City, New York, USA', True)
    assert resp1.status == 200, 'Should see an HTTP 200 response'

    body1 = resp1.read()
    assert is_good_nominatim_response(body1, True), 'Should be good Nominatim response'
    
    results1 = json.loads(body1)
    geometry1 = results1[0]['geojson']
    assert geometry1['type'] == 'Point', 'Should be a point geometry'
    
    lonlat1 = [round(n, 2) for n in geometry1['coordinates']]
    assert lonlat1 == [-73.96, 40.79], 'Should be located in Manhattan'

    #
    
    logging.info(f'Querying Nominatim for Piedmont')
    resp2 = get_nominatim_response('Piedmont, California, USA', True)
    assert resp2.status == 200, 'Should see an HTTP 200 response'

    body2 = resp2.read()
    assert is_good_nominatim_response(body2, True), 'Should be good Nominatim response'
    
    results2 = json.loads(body2)
    geometry2 = results2[0]['geojson']
    assert geometry2['type'] == 'Polygon', 'Should be a polygon geometry'
    
    lonlat2 = [round(n, 2) for n in geometry2['coordinates'][0][0]]
    assert lonlat2 == [-122.25, 37.82], 'Should be located in Piedmont'

    #
    
    logging.info(f'Querying Nominatim for 350 5th Ave')
    resp3 = get_nominatim_response('350 5th Ave, New York, NY', False)
    assert resp3.status == 200, 'Should see an HTTP 200 response'

    body3 = resp3.read()
    assert is_good_nominatim_response(body3, False), 'Should be good Nominatim response'
    
    results3 = json.loads(body3)
    lonlat3 = [round(float(results3[0]['lon']), 2), round(float(results3[0]['lat']), 2)]
    assert lonlat3 == [-73.98, 40.67], 'Should be located in New York'

def lambda_handler(event, context):
    return main()

if __name__ == '__main__':
    exit(main())
