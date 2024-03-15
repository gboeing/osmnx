"""Add node elevations from raster files or web APIs, and calculate edge grades."""

import multiprocessing as mp
import time
from hashlib import sha1
from pathlib import Path
from warnings import warn

import networkx as nx
import numpy as np
import pandas as pd
import requests

from . import _downloader
from . import convert
from . import settings
from . import utils
from ._errors import InsufficientResponseError

# rasterio and gdal are optional dependencies for raster querying
try:
    import rasterio
    from osgeo import gdal
except ImportError:  # pragma: no cover
    rasterio = gdal = None


def add_edge_grades(G, add_absolute=True, precision=None):
    """
    Add `grade` attribute to each graph edge.

    Vectorized function to calculate the directed grade (ie, rise over run)
    for each edge in the graph and add it to the edge as an attribute. Nodes
    must already have `elevation` attributes to use this function.

    See also the `add_node_elevations_raster` and `add_node_elevations_google`
    functions.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph with `elevation` node attribute
    add_absolute : bool
        if True, also add absolute value of grade as `grade_abs` attribute
    precision : int
        deprecated, do not use

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with edge `grade` (and optionally `grade_abs`) attributes
    """
    if precision is None:
        precision = 3
    else:
        warn(
            "The `precision` parameter is deprecated and will be removed in the v2.0.0 release. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123",
            FutureWarning,
            stacklevel=2,
        )

    elev_lookup = G.nodes(data="elevation")
    u, v, k, lengths = zip(*G.edges(keys=True, data="length"))
    uvk = tuple(zip(u, v, k))

    # calculate edges' elevation changes from u to v then divide by lengths
    elevs = np.array([(elev_lookup[u], elev_lookup[v]) for u, v, k in uvk])
    grades = ((elevs[:, 1] - elevs[:, 0]) / np.array(lengths)).round(precision)
    nx.set_edge_attributes(G, dict(zip(uvk, grades)), name="grade")

    # optionally add grade absolute value to the edge attributes
    if add_absolute:
        nx.set_edge_attributes(G, dict(zip(uvk, np.abs(grades))), name="grade_abs")

    utils.log("Added grade attributes to all edges.")
    return G


def _query_raster(nodes, filepath, band):
    """
    Query a raster for values at coordinates in a DataFrame's x/y columns.

    Parameters
    ----------
    nodes : pandas.DataFrame
        DataFrame indexed by node ID and with two columns: x and y
    filepath : string or pathlib.Path
        path to the raster file or VRT to query
    band : int
        which raster band to query

    Returns
    -------
    nodes_values : zip
        zipped node IDs and corresponding raster values
    """
    # must open raster file here: cannot pickle it to pass in multiprocessing
    with rasterio.open(filepath) as raster:
        values = np.array(tuple(raster.sample(nodes.to_numpy(), band)), dtype=float).squeeze()
        values[values == raster.nodata] = np.nan
        return zip(nodes.index, values)


def add_node_elevations_raster(G, filepath, band=1, cpus=None):
    """
    Add `elevation` attribute to each node from local raster file(s).

    If `filepath` is a list of paths, this will generate a virtual raster
    composed of the files at those paths as an intermediate step.

    See also the `add_edge_grades` function.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph, in same CRS as raster
    filepath : string or pathlib.Path or list of strings/Paths
        path (or list of paths) to the raster file(s) to query
    band : int
        which raster band to query
    cpus : int
        how many CPU cores to use; if None, use all available

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with node elevation attributes
    """
    if rasterio is None or gdal is None:  # pragma: no cover
        msg = "gdal and rasterio must be installed to query raster files"
        raise ImportError(msg)

    if cpus is None:
        cpus = mp.cpu_count()
    cpus = min(cpus, mp.cpu_count())
    utils.log(f"Attaching elevations with {cpus} CPUs...")

    # if a list of filepaths is passed, compose them all as a virtual raster
    # use the sha1 hash of the filepaths list as the vrt filename
    if not isinstance(filepath, (str, Path)):
        filepaths = [str(p) for p in filepath]
        sha = sha1(str(filepaths).encode("utf-8")).hexdigest()
        filepath = f"./.osmnx_{sha}.vrt"
        gdal.UseExceptions()
        gdal.BuildVRT(filepath, filepaths).FlushCache()

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

    assert len(G) == len(elevs)
    nx.set_node_attributes(G, elevs, name="elevation")
    utils.log("Added elevation data from raster to all nodes.")
    return G


def add_node_elevations_google(
    G,
    api_key=None,
    batch_size=350,
    pause=0,
    max_locations_per_batch=None,
    precision=None,
    url_template=None,
):
    """
    Add an `elevation` (meters) attribute to each node using a web service.

    By default, this uses the Google Maps Elevation API but you can optionally
    use an equivalent API with the same interface and response format, such as
    Open Topo Data, via the `settings` module's `elevation_url_template`. The
    Google Maps Elevation API requires an API key but other providers may not.

    For a free local alternative see the `add_node_elevations_raster`
    function. See also the `add_edge_grades` function.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    api_key : string
        a valid API key, can be None if the API does not require a key
    batch_size : int
        max number of coordinate pairs to submit in each API call (if this is
        too high, the server will reject the request because its character
        limit exceeds the max allowed)
    pause : float
        time to pause between API calls, which can be increased if you get
        rate limited
    max_locations_per_batch : int
        deprecated, do not use
    precision : int
        deprecated, do not use
    url_template : string
        deprecated, do not use

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with node elevation attributes
    """
    if max_locations_per_batch is None:
        max_locations_per_batch = batch_size
    else:
        warn(
            "The `max_locations_per_batch` parameter is deprecated and will be "
            "removed the v2.0.0 release, use the `batch_size` parameter instead. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123",
            FutureWarning,
            stacklevel=2,
        )

    if precision is None:
        precision = 3
    else:
        warn(
            "The `precision` parameter is deprecated and will be removed in the v2.0.0 release. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123",
            FutureWarning,
            stacklevel=2,
        )

    if url_template is None:
        url_template = settings.elevation_url_template
    else:
        warn(
            "The `url_template` parameter is deprecated and will be removed "
            "in the v2.0.0 release. Configure the `settings` module's "
            "`elevation_url_template` instead. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123",
            FutureWarning,
            stacklevel=2,
        )

    # make a pandas series of all the nodes' coordinates as 'lat,lon'
    # round coordinates to 5 decimal places (approx 1 meter) to be able to fit
    # in more locations per API call
    node_points = pd.Series(
        {node: f'{data["y"]:.5f},{data["x"]:.5f}' for node, data in G.nodes(data=True)}
    )
    n_calls = int(np.ceil(len(node_points) / max_locations_per_batch))
    domain = _downloader._hostname_from_url(url_template)
    utils.log(f"Requesting node elevations from {domain!r} in {n_calls} request(s)")

    # break the series of coordinates into chunks of max_locations_per_batch
    # API format is locations=lat,lon|lat,lon|lat,lon|lat,lon...
    results = []
    for i in range(0, len(node_points), max_locations_per_batch):
        chunk = node_points.iloc[i : i + max_locations_per_batch]
        locations = "|".join(chunk)
        url = url_template.format(locations=locations, key=api_key)

        # download and append these elevation results to list of all results
        response_json = _elevation_request(url, pause)
        if "results" in response_json and len(response_json["results"]) > 0:
            results.extend(response_json["results"])
        else:
            raise InsufficientResponseError(str(response_json))

    # sanity check that all our vectors have the same number of elements
    msg = f"Graph has {len(G):,} nodes and we received {len(results):,} results from {domain!r}"
    utils.log(msg)
    if not (len(results) == len(G) == len(node_points)):  # pragma: no cover
        err_msg = f"{msg}\n{response_json}"
        raise InsufficientResponseError(err_msg)

    # add elevation as an attribute to the nodes
    df_elev = pd.DataFrame(node_points, columns=["node_points"])
    df_elev["elevation"] = [result["elevation"] for result in results]
    df_elev["elevation"] = df_elev["elevation"].round(precision)
    nx.set_node_attributes(G, name="elevation", values=df_elev["elevation"].to_dict())
    utils.log(f"Added elevation data from {domain!r} to all nodes.")

    return G


def _elevation_request(url, pause):
    """
    Send a HTTP GET request to Google Maps-style Elevation API.

    Parameters
    ----------
    url : string
        URL for API endpoint populated with request data
    pause : float
        how long to pause in seconds before request

    Returns
    -------
    response_json : dict
    """
    if settings.timeout is None:
        timeout = settings.requests_timeout
    else:
        timeout = settings.timeout
        msg = (
            "`settings.timeout` is deprecated and will be removed in the v2.0.0 "
            "release: use `settings.requests_timeout` instead. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
        )
        warn(msg, FutureWarning, stacklevel=2)

    # check if request already exists in cache
    cached_response_json = _downloader._retrieve_from_cache(url)
    if cached_response_json is not None:
        return cached_response_json

    # pause then request this URL
    domain = _downloader._hostname_from_url(url)
    utils.log(f"Pausing {pause} second(s) before making HTTP GET request to {domain!r}")
    time.sleep(pause)

    # transmit the HTTP GET request
    utils.log(f"Get {url} with timeout={timeout}")
    response = requests.get(
        url,
        timeout=timeout,
        headers=_downloader._get_http_headers(),
        **settings.requests_kwargs,
    )

    response_json = _downloader._parse_response(response)
    _downloader._save_to_cache(url, response_json, response.status_code)
    return response_json
