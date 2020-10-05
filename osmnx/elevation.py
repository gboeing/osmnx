"""Get node elevations and calculate edge grades."""

import math
import time

import networkx as nx
import numpy as np
import pandas as pd
import requests

from . import downloader
from . import settings
from . import utils


def add_node_elevations(
    G, api_key, max_locations_per_batch=350, pause_duration=0.02, precision=3
):  # pragma: no cover
    """
    Add `elevation` (meters) attribute to each node.

    Uses the Google Maps Elevation API by default, but you can configure
    this to a different provider via ox.config()

    See also the `add_edge_grades` function.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    api_key : string
        your google maps elevation API key, or equivalent if using a different
        provider
    max_locations_per_batch : int
        max number of coordinate pairs to submit in each API call (if this is
        too high, the server will reject the request because its character
        limit exceeds the max)
    pause_duration : float
        time to pause between API calls
    precision : int
        decimal precision to round elevation

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with node elevation attributes
    """
    # different elevation API endpoints formatted ready for use
    endpoints = {
        "google": "https://maps.googleapis.com/maps/api/elevation/json?locations={}&key={}",
        "airmap": "https://api.airmap.com/elevation/v1/ele?points={}",
    }
    url_template = endpoints[settings.elevation_provider]

    # make a pandas series of all the nodes' coordinates as 'lat,lng'
    # round coordinates to 5 decimal places (approx 1 meter) to be able to fit
    # in more locations per API call
    node_points = pd.Series(
        {node: f'{data["y"]:.5f},{data["x"]:.5f}' for node, data in G.nodes(data=True)}
    )
    n_calls = math.ceil(len(node_points) / max_locations_per_batch)
    utils.log(f"Requesting node elevations from the API in {n_calls} calls")

    # break the series of coordinates into chunks of size max_locations_per_batch
    # API format is locations=lat,lng|lat,lng|lat,lng|lat,lng...
    results = []
    for i in range(0, len(node_points), max_locations_per_batch):
        chunk = node_points.iloc[i : i + max_locations_per_batch]

        if settings.elevation_provider == "google":
            locations = "|".join(chunk)
            url = url_template.format(locations, api_key)
        elif settings.elevation_provider == "airmap":
            locations = ",".join(chunk)
            url = url_template.format(locations)

        # check if this request is already in the cache (if global use_cache=True)
        cached_response_json = downloader._retrieve_from_cache(url)
        if cached_response_json is not None:
            response_json = cached_response_json
        else:
            try:
                # request the elevations from the API
                utils.log(f"Requesting node elevations: {url}")
                time.sleep(pause_duration)
                if settings.elevation_provider == "google":
                    response = requests.get(url)
                elif settings.elevation_provider == "airmap":
                    headers = {
                        "X-API-Key": api_key,
                        "Content-Type": "application/json; charset=utf-8",
                    }
                    response = requests.get(url, headers=headers)
                response_json = response.json()
                downloader._save_to_cache(url, response_json, response.status_code)
            except Exception as e:
                utils.log(e)
                utils.log(f"Server responded with {response.status_code}: {response.reason}")

        # append these elevation results to the list of all results
        if settings.elevation_provider == "google":
            results.extend(response_json["results"])
        elif settings.elevation_provider == "airmap":
            results.extend(response_json["data"])

    # sanity check that all our vectors have the same number of elements
    if not (len(results) == len(G) == len(node_points)):
        raise Exception(
            f"Graph has {len(G)} nodes but we received {len(results)} results from elevation API"
        )
    else:
        utils.log(
            f"Graph has {len(G)} nodes and we received {len(results)} results from elevation API"
        )

    # add elevation as an attribute to the nodes
    df = pd.DataFrame(node_points, columns=["node_points"])
    if settings.elevation_provider == "google":
        df["elevation"] = [result["elevation"] for result in results]
    elif settings.elevation_provider == "airmap":
        df["elevation"] = results
    df["elevation"] = df["elevation"].round(precision)
    nx.set_node_attributes(G, name="elevation", values=df["elevation"].to_dict())
    utils.log("Added elevation data to all nodes.")

    return G


def add_edge_grades(G, add_absolute=True, precision=3):
    """
    Add `grade` attribute to each graph edge.

    Get the directed grade (ie, rise over run) for each edge in the graph and
    add it to the edge as an attribute. Nodes must have `elevation` attributes
    to use this function.

    See also the `add_node_elevations` function.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    add_absolute : bool
        if True, also add absolute value of grade as `grade_abs` attribute
    precision : int
        decimal precision to round grade values

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with edge `grade` (and optionally `grade_abs`) attributes
    """
    for u, v, data in G.edges(keys=False, data=True):

        # for each edge, calculate elevation change from origin to destination
        elevation_change = G.nodes[v]["elevation"] - G.nodes[u]["elevation"]

        try:
            # divide by edge length then round
            edge_grade = float(elevation_change) / float(data["length"])
            data["grade"] = round(edge_grade, precision)
        except ZeroDivisionError:
            data["grade"] = np.nan

        # optionally add grade absolute value to the edge attributes
        if add_absolute:
            data["grade_abs"] = abs(data["grade"])

    utils.log("Added grade attributes to all edges.")
    return G
