import sys, urllib.parse, http.client, json, logging

logging.basicConfig(level=logging.INFO, stream=sys.stdout, format='%(levelname)09s - %(message)s')

def get_nominatim_response(query):
    '''
    '''
    params = {'q': query, 'format': 'json', 'limit': 1, 'dedupe': 0, 'polygon_geojson': 1}
    conn = http.client.HTTPConnection('nominatim.openstreetmap.org')
    conn.request('GET', '/search?'+urllib.parse.urlencode(params))
    
    return conn.getresponse()

def is_good_nominatim_response(body):
    '''
    
        Sample response:
        
        [ {
            "boundingbox": [ "40.7702778", "40.8102778", "-73.9797222", "-73.9397222" ],
            "display_name": "Manhattan, NYC, New York, 83, United States of America",
            "geojson": { "coordinates": [ -73.9597222, 40.7902778 ], "type": "Point" },
            ...
        }, ... ]    
    '''
    results = json.loads(body)
    
    if type(results) is not list:
        return False
    
    for result in results:
        if 'boundingbox' not in result:
            return False
        
        if type(result['boundingbox']) is not list:
            return False
        
        if len(result['boundingbox']) is not 4:
            return False
        
        if 'display_name' not in result:
            return False
        
        if 'geojson' not in result:
            return False
        
        geojson = result['geojson']
        
        if 'type' not in geojson:
            return False
        
        if geojson['type'] not in ('Point', 'Polygon', 'LineString',
            'MultiPoint', 'MultiPolygon', 'MultiLineString'):
            return False
        
        if 'coordinates' not in geojson:
            return False
        
        if type(geojson['coordinates']) is not list:
            return False
    
    return True

def main():
    '''
    '''
    logging.info(f'Querying Nominatim for Manhattan')
    resp1 = get_nominatim_response('Manhattan, New York City, New York, USA')
    assert resp1.status == 200, 'Should see an HTTP 200 response'

    body1 = resp1.read()
    assert is_good_nominatim_response(body1), 'Should be good Nominatim response'
    
    results1 = json.loads(body1)
    geometry1 = results1[0]['geojson']
    assert geometry1['type'] == 'Point', 'Should be a point geometry'
    
    lonlat1 = [round(n, 2) for n in geometry1['coordinates']]
    assert lonlat1 == [-73.96, 40.79], 'Should be located in Manhattan'

def lambda_handler(event, context):
    return main()

if __name__ == '__main__':
    exit(main())
