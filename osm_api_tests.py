import sys, urllib.parse, http.client, json, logging, time

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

def get_overpass_response(data, use_post):
    '''
    '''
    params = dict(data=data)
    conn = http.client.HTTPConnection('www.overpass-api.de')
    
    if use_post:
        conn.request('POST', '/api/interpreter', urllib.parse.urlencode(params))
    else:
        conn.request('GET', '/api/interpreter?'+urllib.parse.urlencode(params))
    
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
        logging.warning('Nominatim result is not a list')
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

def is_good_overpass_response(body):
    ''' Return True if JSON body contains expected Overpass API response.
    
        JSON response documentation: http://overpass-api.de/output_formats.html
    '''
    results = json.loads(body)
    
    if type(results) is not dict:
        logging.warning('Overpass response is not a dict')
        return False
    
    if 'elements' not in results:
        logging.warning('Overpass response has no elements')
        return False
    
    if type(results['elements']) is not list:
        logging.warning('Overpass elements is not a list')
        return False
    
    for element in results['elements']:
        if element['type'] == 'node':
            if 'lat' not in element or 'lon' not in element:
                logging.warning('No lon, lat in Overpass node element')
                return False
        
        elif element['type'] in ('way', 'relation'):
            pass
        
        else:
            logging.warning('Unknown Overpass element type')
            return False
    
    return True

def is_bounded_overpass_response(elements, minlat, maxlat, minlon, maxlon):
    ''' Return True if all Overpass API nodes are in bounds.
    '''
    node_lats = [float(el['lat']) for el in elements if el['type'] == 'node']
    node_lons = [float(el['lon']) for el in elements if el['type'] == 'node']
    
    if min(node_lats) < minlat:
        logging.warning('An Overpass node element is too far south')
        return False

    if max(node_lats) > maxlat:
        logging.warning('An Overpass node element is too far north')
        return False

    if min(node_lons) < minlon:
        logging.warning('An Overpass node element is too far west')
        return False

    if max(node_lons) > maxlon:
        logging.warning('An Overpass node element is too far east')
        return False
        
    return True

def is_correctly_tagged_overpass_response(elements, tag_name, expected_values, unwanted_values):
    ''' Return True if all way and relation tags match expectations.
    
        Look at one tag name, and check for presence of expected values
        and absence of unwanted values.
    '''
    all_values = set([el['tags'].get(tag_name) for el in elements
        if el['type'] in ('way', 'relation') and 'tags' in el])
    
    for expected_value in expected_values:
        if expected_value not in all_values:
            logging.warning('An expected tag is missing from Overpass way elements')
            return False

    for unwanted_value in unwanted_values:
        if unwanted_value in all_values:
            logging.warning('An unwanted tag is present from Overpass way elements')
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
    
    geometry1 = json.loads(body1)[0]['geojson']
    assert geometry1['type'] == 'Point', 'Should be a point geometry'
    
    lonlat1 = [round(n, 2) for n in geometry1['coordinates']]
    assert lonlat1 == [-73.96, 40.79], 'Should be located in Manhattan'

    #
    
    logging.info(f'Querying Nominatim for Piedmont')
    resp2 = get_nominatim_response('Piedmont, California, USA', True)
    assert resp2.status == 200, 'Should see an HTTP 200 response'

    body2 = resp2.read()
    assert is_good_nominatim_response(body2, True), 'Should be good Nominatim response'
    
    geometry2 = json.loads(body2)[0]['geojson']
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
    
    #
    
    logging.info(f'Querying Overpass for driving + service roads in a bounding box')
    data4 = '[out:json][timeout:180];(way["highway"]["area"!~"yes"]["highway"!~"cycleway|footway|path|pedestrian|steps|track|proposed|construction|bridleway|abandoned|platform|raceway"]["motor_vehicle"!~"no"]["motorcar"!~"no"]["access"!~"private"]["service"!~"parking|parking_aisle|private|emergency_access"](37.77549352,-122.43567860,37.79450647,-122.40432141);>;);out;'
    resp4 = get_overpass_response(data4, False)
    
    body4 = resp4.read()
    assert is_good_overpass_response(body4), 'Should be good Overpass response'
    
    elements4 = json.loads(body4)['elements']
    assert is_bounded_overpass_response(elements4, 37.7, 37.9, -122.5, -122.3)
    assert is_correctly_tagged_overpass_response(elements4, 'highway', ('secondary', 'service'), ('path', 'abandoned'))
    
    #
    
    time.sleep(1)
    logging.info(f'Querying Overpass for driving roads in a bounding box')
    data5 = '[out:json][timeout:180];(way["highway"]["area"!~"yes"]["highway"!~"cycleway|footway|path|pedestrian|steps|track|proposed|construction|bridleway|abandoned|platform|raceway|service"]["motor_vehicle"!~"no"]["motorcar"!~"no"]["access"!~"private"]["service"!~"parking|parking_aisle|driveway|private|emergency_access"](37.78016097,-122.42421511,37.80269301,-122.39582092);>;);out;'
    resp5 = get_overpass_response(data5, False)
    
    body5 = resp5.read()
    assert is_good_overpass_response(body5), 'Should be good Overpass response'
    
    elements5 = json.loads(body5)['elements']
    assert is_bounded_overpass_response(elements5, 37.7, 37.9, -122.5, -122.3)
    assert is_correctly_tagged_overpass_response(elements5, 'highway', ('secondary', ), ('service', 'path', 'abandoned'))
    
    #
    
    time.sleep(1)
    logging.info(f'Querying Overpass for all roads in a bounding polygon')
    data6 = '[out:json][timeout:180];(way["highway"]["area"!~"yes"]["highway"!~"proposed|construction|abandoned|platform|raceway"](poly:"37.8231128 -122.2550101 37.8231986 -122.2550266 37.8236514 -122.2550550 37.8241041 -122.2550260 37.8245522 -122.2549398 37.8249911 -122.2547972 37.8254165 -122.2545998 37.8258240 -122.2543495 37.8262095 -122.2540489 37.8265691 -122.2537009 37.8267985 -122.2534549 37.8269152 -122.2533304 37.8271762 -122.2532437 37.8273771 -122.2531777 37.8278002 -122.2530103 37.8282083 -122.2527913 37.8285973 -122.2525227 37.8289635 -122.2522072 37.8293033 -122.2518479 37.8296134 -122.2514482 37.8298909 -122.2510120 37.8301428 -122.2505734 37.8308034 -122.2494288 37.8308386 -122.2494069 37.8311972 -122.2491284 37.8312856 -122.2490523 37.8321234 -122.2484068 37.8323930 -122.2481996 37.8324077 -122.2481883 37.8324477 -122.2481573 37.8328059 -122.2478488 37.8331389 -122.2474983 37.8334437 -122.2471091 37.8337173 -122.2466847 37.8339574 -122.2462291 37.8341615 -122.2457465 37.8343279 -122.2452414 37.8344550 -122.2447184 37.8345417 -122.2441825 37.8345871 -122.2436385 37.8345889 -122.2433716 37.8346776 -122.2430385 37.8347305 -122.2428266 37.8347635 -122.2426846 37.8348903 -122.2419715 37.8349093 -122.2418195 37.8349546 -122.2413063 37.8349626 -122.2411523 37.8349718 -122.2408642 37.8349835 -122.2398543 37.8349880 -122.2397427 37.8350128 -122.2395611 37.8350401 -122.2392889 37.8350610 -122.2392269 37.8351885 -122.2387224 37.8352783 -122.2382051 37.8353140 -122.2378386 37.8353481 -122.2377173 37.8354480 -122.2372274 37.8354813 -122.2370275 37.8355202 -122.2367957 37.8355548 -122.2366028 37.8356908 -122.2358689 37.8357368 -122.2356199 37.8357384 -122.2356110 37.8359711 -122.2343422 37.8360408 -122.2339666 37.8364541 -122.2326301 37.8365616 -122.2322824 37.8368395 -122.2313839 37.8372227 -122.2301484 37.8373632 -122.2296209 37.8374620 -122.2290785 37.8375183 -122.2285266 37.8375314 -122.2279704 37.8375014 -122.2274152 37.8374285 -122.2268664 37.8373133 -122.2263292 37.8371570 -122.2258088 37.8369610 -122.2253101 37.8367273 -122.2248381 37.8364581 -122.2243971 37.8363208 -122.2241936 37.8359085 -122.2235811 37.8359015 -122.2235707 37.8358600 -122.2235093 37.8356553 -122.2230849 37.8355952 -122.2229848 37.8355864 -122.2229684 37.8355777 -122.2229521 37.8354337 -122.2226841 37.8354293 -122.2226760 37.8353813 -122.2225870 37.8353737 -122.2225730 37.8353222 -122.2224781 37.8352980 -122.2224335 37.8350627 -122.2219948 37.8350549 -122.2219803 37.8349738 -122.2218303 37.8347163 -122.2213989 37.8344274 -122.2210003 37.8341098 -122.2206380 37.8337663 -122.2203154 37.8334001 -122.2200354 37.8330145 -122.2198005 37.8328097 -122.2196906 37.8322717 -122.2194020 37.8319766 -122.2191834 37.8316135 -122.2189147 37.8306543 -122.2180950 37.8306465 -122.2180883 37.8305910 -122.2180411 37.8305154 -122.2179760 37.8304956 -122.2179590 37.8304115 -122.2178873 37.8301443 -122.2176586 37.8298068 -122.2173963 37.8297328 -122.2173443 37.8296309 -122.2172806 37.8295031 -122.2171815 37.8294300 -122.2171305 37.8289865 -122.2168595 37.8287005 -122.2167085 37.8286619 -122.2166884 37.8283788 -122.2165428 37.8275161 -122.2154528 37.8274323 -122.2153469 37.8271978 -122.2150505 37.8271903 -122.2150410 37.8270290 -122.2146739 37.8268701 -122.2143054 37.8267976 -122.2141549 37.8267402 -122.2140404 37.8267112 -122.2139569 37.8266442 -122.2137859 37.8264158 -122.2132719 37.8261464 -122.2127906 37.8258389 -122.2123469 37.8257489 -122.2122299 37.8253885 -122.2118098 37.8250070 -122.2114538 37.8250046 -122.2114493 37.8248915 -122.2112403 37.8247578 -122.2110060 37.8245017 -122.2105805 37.8237325 -122.2093009 37.8236484 -122.2091748 37.8235729 -122.2090403 37.8234918 -122.2088937 37.8232145 -122.2083889 37.8231457 -122.2082637 37.8230313 -122.2080474 37.8226193 -122.2072661 37.8224820 -122.2070055 37.8222254 -122.2065658 37.8219365 -122.2061591 37.8216180 -122.2057893 37.8212727 -122.2054597 37.8209040 -122.2051735 37.8205152 -122.2049331 37.8201099 -122.2047410 37.8196918 -122.2045987 37.8192648 -122.2045077 37.8188328 -122.2044688 37.8183999 -122.2044824 37.8179700 -122.2045483 37.8175471 -122.2046658 37.8171352 -122.2048340 37.8167379 -122.2050513 37.8163590 -122.2053157 37.8160020 -122.2056246 37.8156702 -122.2059753 37.8153665 -122.2063646 37.8149586 -122.2069406 37.8149297 -122.2069817 37.8146749 -122.2073486 37.8144811 -122.2073807 37.8143391 -122.2074098 37.8141774 -122.2074466 37.8140304 -122.2074836 37.8136337 -122.2076073 37.8134998 -122.2076573 37.8129183 -122.2079305 37.8129113 -122.2079345 37.8128273 -122.2079838 37.8126913 -122.2080658 37.8121319 -122.2084674 37.8119780 -122.2085971 37.8119663 -122.2086028 37.8115751 -122.2088208 37.8114831 -122.2088788 37.8112387 -122.2090450 37.8112127 -122.2090640 37.8111914 -122.2090796 37.8111587 -122.2091038 37.8111361 -122.2091202 37.8107766 -122.2094102 37.8107296 -122.2094522 37.8103449 -122.2098374 37.8102099 -122.2099884 37.8100346 -122.2101953 37.8099406 -122.2103123 37.8099220 -122.2103374 37.8096735 -122.2105369 37.8093272 -122.2108788 37.8090090 -122.2112617 37.8087219 -122.2116819 37.8084687 -122.2121353 37.8082517 -122.2126178 37.8080732 -122.2131245 37.8079347 -122.2136508 37.8078375 -122.2141916 37.8077827 -122.2147418 37.8077707 -122.2152961 37.8078016 -122.2158492 37.8078751 -122.2163959 37.8079906 -122.2169310 37.8081470 -122.2174493 37.8082402 -122.2177179 37.8084123 -122.2182155 37.8084985 -122.2184277 37.8085350 -122.2186465 37.8085364 -122.2186461 37.8085856 -122.2189431 37.8086188 -122.2191088 37.8086473 -122.2192442 37.8087398 -122.2196621 37.8087763 -122.2198502 37.8090868 -122.2212253 37.8091058 -122.2212984 37.8091814 -122.2216490 37.8092646 -122.2219914 37.8093176 -122.2221874 37.8093224 -122.2222052 37.8093894 -122.2224503 37.8094458 -122.2226780 37.8095655 -122.2234399 37.8095691 -122.2234625 37.8096171 -122.2237615 37.8097045 -122.2242108 37.8097675 -122.2244858 37.8099565 -122.2253108 37.8099606 -122.2253285 37.8100246 -122.2256045 37.8101102 -122.2259356 37.8101802 -122.2261806 37.8101842 -122.2261945 37.8103962 -122.2269304 37.8104115 -122.2269822 37.8105203 -122.2274651 37.8105218 -122.2274715 37.8109800 -122.2294944 37.8111036 -122.2300396 37.8111429 -122.2302504 37.8111852 -122.2304598 37.8112575 -122.2307906 37.8116883 -122.2329123 37.8116902 -122.2329213 37.8118381 -122.2336453 37.8120298 -122.2343745 37.8120978 -122.2345845 37.8121046 -122.2346055 37.8123101 -122.2352331 37.8123347 -122.2353083 37.8123649 -122.2354157 37.8125134 -122.2359627 37.8125694 -122.2361572 37.8126236 -122.2363351 37.8126715 -122.2365042 37.8126947 -122.2365901 37.8126990 -122.2366062 37.8127662 -122.2368526 37.8128712 -122.2374074 37.8129178 -122.2379612 37.8129194 -122.2379799 37.8129394 -122.2382099 37.8130063 -122.2387434 37.8131131 -122.2392665 37.8131975 -122.2396111 37.8134505 -122.2406431 37.8134535 -122.2406552 37.8134986 -122.2408377 37.8135174 -122.2409335 37.8136704 -122.2417134 37.8137213 -122.2419734 37.8137667 -122.2421868 37.8138577 -122.2425838 37.8141307 -122.2437747 37.8141335 -122.2437871 37.8142255 -122.2441850 37.8142776 -122.2443951 37.8143176 -122.2445461 37.8144376 -122.2449991 37.8144534 -122.2450501 37.8145799 -122.2454951 37.8146339 -122.2456621 37.8146414 -122.2456849 37.8148046 -122.2461834 37.8148578 -122.2463470 37.8150115 -122.2467692 37.8151921 -122.2471743 37.8152321 -122.2472560 37.8153542 -122.2475087 37.8153732 -122.2475474 37.8154056 -122.2476148 37.8155885 -122.2479663 37.8157916 -122.2482997 37.8159039 -122.2484707 37.8159255 -122.2485036 37.8160642 -122.2488061 37.8163701 -122.2493261 37.8164131 -122.2493911 37.8166433 -122.2497132 37.8166703 -122.2497482 37.8166746 -122.2497539 37.8167156 -122.2498069 37.8169757 -122.2501174 37.8171127 -122.2502684 37.8173429 -122.2505057 37.8174789 -122.2506367 37.8180255 -122.2510881 37.8181100 -122.2511473 37.8183931 -122.2514709 37.8184027 -122.2514819 37.8186676 -122.2518513 37.8195966 -122.2531463 37.8199015 -122.2535324 37.8202344 -122.2538799 37.8205921 -122.2541857 37.8209715 -122.2544468 37.8213689 -122.2546609 37.8217808 -122.2548261 37.8222033 -122.2549407 37.8226326 -122.2550038 37.8230647 -122.2550147 37.8231128 -122.2550101");>;);out;'
    resp6 = get_overpass_response(data6, True)
    
    body6 = resp6.read()
    assert is_good_overpass_response(body6), 'Should be good Overpass response'
    
    elements6 = json.loads(body6)['elements']
    assert is_bounded_overpass_response(elements6, 37.8, 38.0, -122.3, -122.1)
    assert is_correctly_tagged_overpass_response(elements6, 'highway', ('secondary', 'service', 'path'), ('abandoned', ))
    
    #
    
    time.sleep(1)
    logging.info(f'Querying Overpass for all buildings in a bounding box')
    data7 = '[out:json][timeout:180];((way["building"](37.79387218,-122.41182759,37.79927982,-122.40501281);(._;>;););(relation["building"](37.79387218,-122.41182759,37.79927982,-122.40501281);(._;>;);));out;'
    resp7 = get_overpass_response(data7, True)
    
    body7 = resp7.read()
    assert is_good_overpass_response(body7), 'Should be good Overpass response'
    
    elements7 = json.loads(body7)['elements']
    assert is_bounded_overpass_response(elements7, 37.7, 37.9, -122.5, -122.3)
    assert is_correctly_tagged_overpass_response(elements7, 'building', ('yes', ), ('no', 'false'))

def lambda_handler(event, context):
    return main()

if __name__ == '__main__':
    exit(main())
