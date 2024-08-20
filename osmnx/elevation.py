"""Add node elevations from raster files or web APIs, and calculate edge grades."""

from __future__ import annotations

import logging as lg
import multiprocessing as mp
import time
from hashlib import sha1
from pathlib import Path
from typing import TYPE_CHECKING
from typing import Any

import networkx as nx
import numpy as np
import pandas as pd
import requests

from . import _http
from . import convert
from . import settings
from . import utils
from ._errors import InsufficientResponseError

if TYPE_CHECKING:
    from collections.abc import Iterable

# rasterio and rio-vrt are optional dependencies for raster querying
try:
    import rasterio
    from rio_vrt import build_vrt
except ImportError:  # pragma: no cover
    rasterio = None
    build_vrt = None


def add_edge_grades(G: nx.MultiDiGraph, *, add_absolute: bool = True) -> nx.MultiDiGraph:
    """
    Calculate and add `grade` attributes to all graph edges.

    Vectorized function to calculate the directed grade (i.e., rise over run)
    for each edge in the graph and add it to the edge as an attribute. Nodes
    must already have `elevation` and `length` attributes before using this
    function.

    See also the `add_node_elevations_raster` and `add_node_elevations_google`
    functions.

    Parameters
    ----------
    G
        Graph with `elevation` node attributes.
    add_absolute
        If True, also add absolute value of grade as `grade_abs` attribute.

    Returns
    -------
    G
        Graph with `grade` (and optionally `grade_abs`) attributes on the
        edges.
    """
    elev_lookup = G.nodes(data="elevation")
    u, v, k, lengths = zip(*G.edges(keys=True, data="length"))
    uvk = tuple(zip(u, v, k))

    # calculate edges' elevation changes from u to v then divide by lengths
    elevs = np.array([(elev_lookup[u], elev_lookup[v]) for u, v, k in uvk])
    grades = (elevs[:, 1] - elevs[:, 0]) / np.array(lengths)
    nx.set_edge_attributes(G, dict(zip(uvk, grades)), name="grade")

    # optionally add grade absolute value to the edge attributes
    if add_absolute:
        nx.set_edge_attributes(G, dict(zip(uvk, np.abs(grades))), name="grade_abs")

    msg = "Added grade attributes to all edges"
    utils.log(msg, level=lg.INFO)
    return G


def _query_raster(
    nodes: pd.DataFrame,
    filepath: str | Path,
    band: int,
) -> Iterable[tuple[int, Any]]:
    """
    Query a raster file for values at coordinates in DataFrame x/y columns.

    Parameters
    ----------
    nodes
        DataFrame indexed by node ID and with two columns representing x and y
        coordinates.
    filepath
        Path to the raster file or VRT to query.
    band
        Which raster band to query.

    Returns
    -------
    nodes_values
        Zip of node IDs and corresponding raster values.
    """
    # must open raster file here: cannot pickle it to pass in multiprocessing
    with rasterio.open(filepath) as raster:
        values = np.array(tuple(raster.sample(nodes.to_numpy(), band)), dtype=float).squeeze()
        values[values == raster.nodata] = np.nan
        return zip(nodes.index, values)


def add_node_elevations_raster(
    G: nx.MultiDiGraph,
    filepath: str | Path | Iterable[str | Path],
    *,
    band: int = 1,
    cpus: int | None = None,
) -> nx.MultiDiGraph:
    """
    Add `elevation` attributes to all nodes from local raster file(s).

    If `filepath` is an iterable of paths, this will generate a virtual raster
    composed of the files at those paths as an intermediate step.

    See also the `add_edge_grades` function.

    Parameters
    ----------
    G
        Graph in same CRS as raster.
    filepath
        The path(s) to the raster file(s) to query.
    band
        Which raster band to query.
    cpus
        How many CPU cores to use. If None, use all available.

    Returns
    -------
    G
        Graph with `elevation` attributes on the nodes.
    """
    if rasterio is None or build_vrt is None:  # pragma: no cover
        msg = "rasterio and rio-vrt must be installed as optional dependencies to query rasters."
        raise ImportError(msg)

    if cpus is None:
        cpus = mp.cpu_count()
    cpus = min(cpus, mp.cpu_count())
    msg = f"Attaching elevations with {cpus} CPUs..."
    utils.log(msg, level=lg.INFO)

    # if multiple filepaths are passed in, compose them as a virtual raster
    # use the sha1 hash of the filepaths object as the vrt filename
    if not isinstance(filepath, (str, Path)):
        filepaths = [str(p) for p in filepath]
        checksum = sha1(str(filepaths).encode("utf-8")).hexdigest()  # noqa: S324
        filepath = f"./.osmnx_{checksum}.vrt"
        build_vrt(filepath, filepaths)

    nodes = convert.graph_to_gdfs(G, edges=False, node_geometry=False)[["x", "y"]]
    if cpus == 1:
        elevs = dict(_query_raster(nodes, filepath, band))
    else:
        # divide nodes into equal-sized chunks for multiprocessing
        size = int(np.ceil(len(nodes) / cpus))
        args = ((nodes.iloc[i : i + size], filepath, band) for i in range(0, len(nodes), size))
        with mp.get_context("spawn").Pool(cpus) as pool:
            results = pool.starmap_async(_query_raster, args).get()
        elevs = {k: v for kv in results for k, v in kv}

    nx.set_node_attributes(G, elevs, name="elevation")
    msg = "Added elevation data from raster to all nodes"
    utils.log(msg, level=lg.INFO)
    return G


def add_node_elevations_google(
    G: nx.MultiDiGraph,
    *,
    api_key: str | None = None,
    batch_size: int = 512,
    pause: float = 0,
) -> nx.MultiDiGraph:
    """
    Add `elevation` (meters) attributes to all nodes using a web API.

    By default this uses the Google Maps Elevation API, but you could instead
    use any equivalent API with the same interface and response format (such
    as the Open Topo Data API or the Open-Elevation API) via the `settings`
    module's `elevation_url_template`. Adjust the `batch_size` and `pause`
    arguments as needed for the provider. The Google Maps Elevation API
    requires an API key but other providers may not. You can find more
    information about the Google Maps Elevation API interface and format at:
    https://developers.google.com/maps/documentation/elevation

    For a free local alternative see the `add_node_elevations_raster`
    function. See also the `add_edge_grades` function.

    Parameters
    ----------
    G
        Graph to add elevation data to.
    api_key
        A valid API key. Can be None if the API does not require a key.
    batch_size
        Max number of coordinate pairs to submit in each request (depends on
        provider's limits). Google's limit is 512.
    pause
        How long to pause in seconds between API calls, which can be increased
        if you get rate limited.

    Returns
    -------
    G
        Graph with `elevation` attributes on the nodes.
    """
    # make a pandas series of all the nodes' coordinates as "lat,lon" and
    # round coordinates to 6 decimal places (approx 5 to 10 cm resolution)
    node_points = pd.Series({n: f"{d['y']:.6f},{d['x']:.6f}" for n, d in G.nodes(data=True)})
    n_calls = int(np.ceil(len(node_points) / batch_size))
    domain = _http._hostname_from_url(settings.elevation_url_template)

    msg = f"Requesting node elevations from {domain!r} in {n_calls} request(s)"
    utils.log(msg, level=lg.INFO)

    # break the series of coordinates into chunks of batch_size
    # API format is locations=lat,lon|lat,lon|lat,lon|lat,lon...
    results = []
    for i in range(0, len(node_points), batch_size):
        chunk = node_points.iloc[i : i + batch_size]
        locations = "|".join(chunk)
        url = settings.elevation_url_template.format(locations=locations, key=api_key)

        # download and append these elevation results to list of all results
        response_json = _elevation_request(url, pause)
        if "results" in response_json and len(response_json["results"]) > 0:
            results.extend(response_json["results"])
        else:
            raise InsufficientResponseError(str(response_json))

    # sanity check that all our vectors have the same number of elements
    msg = f"Graph has {len(G):,} nodes and we received {len(results):,} results from {domain!r}"
    utils.log(msg, level=lg.INFO)
    if not (len(results) == len(G) == len(node_points)):  # pragma: no cover
        err_msg = f"{msg}\n{response_json}"
        raise InsufficientResponseError(err_msg)

    # add elevation as an attribute to the nodes
    df_elev = pd.DataFrame(node_points, columns=["node_points"])
    df_elev["elevation"] = [result["elevation"] for result in results]
    nx.set_node_attributes(G, name="elevation", values=df_elev["elevation"].to_dict())
    msg = f"Added elevation data from {domain!r} to all nodes."
    utils.log(msg, level=lg.INFO)

    return G


def _elevation_request(url: str, pause: float) -> dict[str, Any]:
    """
    Send a HTTP GET request to a Google Maps-style Elevation API.

    Parameters
    ----------
    url
        URL of API endpoint, populated with request data.
    pause
        How long to pause in seconds before request.

    Returns
    -------
    response_json
    """
    # check if request already exists in cache
    cached_response_json = _http._retrieve_from_cache(url)
    if isinstance(cached_response_json, dict):
        return cached_response_json

    # pause then request this URL
    domain = _http._hostname_from_url(url)
    msg = f"Pausing {pause} second(s) before making HTTP GET request to {domain!r}"
    utils.log(msg, level=lg.INFO)
    time.sleep(pause)

    # transmit the HTTP GET request
    msg = f"Get {url} with timeout={settings.requests_timeout}"
    utils.log(msg, level=lg.INFO)
    response = requests.get(
        url,
        timeout=settings.requests_timeout,
        headers=_http._get_http_headers(),
        **settings.requests_kwargs,
    )

    response_json = _http._parse_response(response)
    if not isinstance(response_json, dict):  # pragma: no cover
        msg = "Elevation API did not return a dict of results."
        raise InsufficientResponseError(msg)
    _http._save_to_cache(url, response_json, response.ok)
    return response_json
