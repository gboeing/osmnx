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
import requests
import time

from .core import get_from_cache
from .core import save_to_cache
from .utils import log


def add_node_elevations(G, api_key, max_locations_per_batch=350,
                        pause_duration=0.02): # pragma: no cover
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
