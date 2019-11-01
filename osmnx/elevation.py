################################################################################
# Module: elevation.py
# Description: Get node elevations and edge grades from the Google Maps
#              Elevation API
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

import math
import networkx as nx
import pandas as pd
import rasterio
import requests
import time
from scipy.interpolate import interp2d

from .core import get_from_cache
from .core import save_to_cache
from .utils import log



def add_node_elevations_from_dem(G, dem_path, num_buffer=1, interp_kind='linear'):
    """
    Add elevation for each node in network using a Digital Elevation Model (DEM)
    and add to node as an attribute. Measurement unit is unchanged from the
    source DEM; most DEM files are in meters.

    DEM sources include SRTM, which has worldwide 1 arc-second (30m) coverage,
    and USGS 3DEP, which has 1/3 arc-second (10m) and 1 arc-second (30m)
    coverage for the continental United States.

    Parameters
    ----------
    G : networkx multidigraph
    dem_path : str
        path to DEM file with elevation data. Can be any file format readable by
        GDAL.
    num_buffer : int
        number of bordering cells around (lon, lat) to use when interpolating
    interp_kind : str
        kind of interpolation. Passed to scipy.interpolate.interp2d. Can be
        ['linear’, ‘cubic’, ‘quintic']

    Returns
    -------
    G : networkx multidigraph
    """

    # make a pandas series of all the nodes' coordinates as 'lat,lng'
    # round coorindates to 5 decimal places (approx 1 meter) to be able to fit
    # in more locations per API call
    node_points = pd.DataFrame([(node, data['x'], data['y'])
                                for node, data in G.nodes(data=True)],
                               columns=['node', 'lon', 'lat'])
    node_points = node_points.set_index('node')

    results = []
    for row in node_points.itertuples():
        value = query_dem(dem_path=dem_path, lon=row.lon, lat=row.lat)
        results.append(value)

    # sanity check that all our vectors have the same number of elements
    if not (len(results) == len(G.nodes()) == len(node_points)):
        raise Exception(
            'Graph has {} nodes but we received {} results from the elevation API.'
            .format(len(G.nodes()), len(results)))
    else:
        log('Graph has {} nodes and we received {} results from the elevation API.'
            .format(len(G.nodes()), len(results)))

    # add elevation as an attribute to the nodes
    node_points['elevation'] = results
    # round elevation to millimeter
    node_points['elevation'] = node_points['elevation'].round(3)
    nx.set_node_attributes(G,
                           name='elevation',
                           values=node_points['elevation'].to_dict())
    log('Added elevation data to all nodes.')

    return G


def query_dem(dem_path, lon, lat, num_buffer=1, interp_kind='linear'):
    """Query elevation data for given point

    Parameters
    ----------
    lon : float
        longitude
    lat : float
        latitude
    num_buffer : int
        number of bordering cells around (lon, lat) to use when interpolating
    interp_kind : str
        kind of interpolation. Passed to scipy.interpolate.interp2d. Can be
        ['linear’, ‘cubic’, ‘quintic']

    Returns
    -------
    value : float
        elevation in terms of the unit of the DEM (usually meters)
    """
    # Read metadata of file
    dataset = rasterio.open(dem_path)

    # Find x, y of elevation square inside raster
    x, y = dataset.index(lon, lat)

    # Make window include cells around it
    # The number of additional cells depends on the value of num_buffer
    # When num_buffer==1, an additional 8 cells will be loaded and
    # interpolated on;
    # When num_buffer==2, an additional 24 cells will be loaded and
    # interpolated on, etc.
    # When using kind='linear' interpolation, I'm not sure if having the
    # extra cells makes a difference; ie if it creates the plane based only
    # on the closest cells or from all. When using kind='cubic', it's
    # probably more accurate with more cells.

    minx = x - num_buffer if x >= num_buffer else x
    maxx = x + num_buffer if x + num_buffer <= dataset.width else x
    miny = y - num_buffer if y >= num_buffer else y
    maxy = y + num_buffer if y + num_buffer <= dataset.width else y

    # Add +1 to deal with range() not including end
    maxx += 1
    maxy += 1

    # Retrieve just this window of values from the DEM
    # Assumes that the elevation data is in band 1 of the file
    window = ([minx, maxx], [miny, maxy])
    val_arr = dataset.read(1, window=window)

    msg = 'array has too few or too many values'
    max_num = 2 * num_buffer + 1
    assert (1 <= val_arr.shape[0] <= max_num) and (1 <= val_arr.shape[1] <=
                                                   max_num), msg

    # Now linearly interpolate
    # Get actual lat/lons
    # Note that zipping together means that I get the diagonal, i.e. one of
    # each of x, y. Since these aren't projected coordinates, but rather the
    # original lat/lons, this is a regular grid and this is ok.
    lonlats = [
        dataset.xy(x, y) for x, y in zip(range(minx, maxx), range(miny, maxy))
    ]
    lons = [x[0] for x in lonlats]
    lats = [x[1] for x in lonlats]

    fun = interp2d(x=lons, y=lats, z=val_arr, kind=interp_kind)
    value = fun(lon, lat)
    return value[0]


def add_node_elevations(G,
                        api_key,
                        max_locations_per_batch=350,
                        pause_duration=0.02):  # pragma: no cover
    """
    Get the elevation (meters) of each node in the network and add it to the
    node as an attribute.

    Parameters
    ----------
    G : networkx multidigraph
    api_key : string
        your google maps elevation API key
    max_locations_per_batch : int
        max number of coordinate pairs to submit in each API call (if this is
        too high, the server will reject the request because its character
        limit exceeds the max)
    pause_duration : float
        time to pause between API calls

    Returns
    -------
    G : networkx multidigraph
    """

    # google maps elevation API endpoint
    url_template = 'https://maps.googleapis.com/maps/api/elevation/json?locations={}&key={}'

    # make a pandas series of all the nodes' coordinates as 'lat,lng'
    # round coorindates to 5 decimal places (approx 1 meter) to be able to fit
    # in more locations per API call
    node_points = pd.Series({node:'{:.5f},{:.5f}'.format(data['y'], data['x']) for node, data in G.nodes(data=True)})
    log('Requesting node elevations from the API in {} calls.'.format(math.ceil(len(node_points) / max_locations_per_batch)))

    # break the series of coordinates into chunks of size max_locations_per_batch
    # API format is locations=lat,lng|lat,lng|lat,lng|lat,lng...
    results = []
    for i in range(0, len(node_points), max_locations_per_batch):
        chunk = node_points.iloc[i : i + max_locations_per_batch]
        locations = '|'.join(chunk)
        url = url_template.format(locations, api_key)

        # check if this request is already in the cache (if global use_cache=True)
        cached_response_json = get_from_cache(url)
        if cached_response_json is not None:
            response_json = cached_response_json
        else:
            try:
                # request the elevations from the API
                log('Requesting node elevations: {}'.format(url))
                time.sleep(pause_duration)
                response = requests.get(url)
                response_json = response.json()
                save_to_cache(url, response_json)
            except Exception as e:
                log(e)
                log('Server responded with {}: {}'.format(response.status_code, response.reason))

        # append these elevation results to the list of all results
        results.extend(response_json['results'])

    # sanity check that all our vectors have the same number of elements
    if not (len(results) == len(G.nodes()) == len(node_points)):
        raise Exception('Graph has {} nodes but we received {} results from the elevation API.'.format(len(G.nodes()), len(results)))
    else:
        log('Graph has {} nodes and we received {} results from the elevation API.'.format(len(G.nodes()), len(results)))

    # add elevation as an attribute to the nodes
    df = pd.DataFrame(node_points, columns=['node_points'])
    df['elevation'] = [result['elevation'] for result in results]
    df['elevation'] = df['elevation'].round(3) # round to millimeter
    nx.set_node_attributes(G, name='elevation', values=df['elevation'].to_dict())
    log('Added elevation data to all nodes.')

    return G



def add_edge_grades(G, add_absolute=True): # pragma: no cover
    """
    Get the directed grade (ie, rise over run) for each edge in the network and
    add it to the edge as an attribute. Nodes must have elevation attributes to
    use this function.

    Parameters
    ----------
    G : networkx multidigraph
    add_absolute : bool
        if True, also add the absolute value of the grade as an edge attribute

    Returns
    -------
    G : networkx multidigraph
    """

    # for each edge, calculate the difference in elevation from origin to
    # destination, then divide by edge length
    for u, v, data in G.edges(keys=False, data=True):
        elevation_change = G.nodes[v]['elevation'] - G.nodes[u]['elevation']

        # round to ten-thousandths decimal place
        try:
            grade = round(elevation_change / data['length'], 4)
        except ZeroDivisionError:
            grade = None

        # add grade and (optionally) grade absolute value to the edge data
        data['grade'] = grade
        if add_absolute:
            data['grade_abs'] = abs(grade)

    log('Added grade data to all edges.')
    return G
