"""Tools to work with the Nominatim API."""

from __future__ import annotations

import logging as lg
import time
from collections import OrderedDict
from typing import Any

import requests

from . import _downloader
from . import settings
from . import utils


def _download_nominatim_element(
    query: str | dict[str, str],
    by_osmid: bool = False,
    limit: int = 1,
    polygon_geojson: bool = True,
) -> list[dict[str, Any]]:
    """
    Retrieve an OSM element from the Nominatim API.

    Parameters
    ----------
    query : string or dict
        query string or structured query dict
    by_osmid : bool
        if True, treat query as an OSM ID lookup rather than text search
    limit : int
        max number of results to return
    polygon_geojson : bool
        whether to retrieve the place's geometry from the API

    Returns
    -------
    response_json : list
        JSON response from the Nominatim server
    """
    # define the parameters
    params: OrderedDict[str, int | str] = OrderedDict()
    params["format"] = "json"
    params["polygon_geojson"] = int(polygon_geojson)  # bool -> int

    if by_osmid:
        # if querying by OSM ID, use the lookup endpoint
        if not isinstance(query, str):
            msg = "`query` must be a string if `by_osmid` is True"
            raise TypeError(msg)
        request_type = "lookup"
        params["osm_ids"] = query

    else:
        # if not querying by OSM ID, use the search endpoint
        request_type = "search"

        # prevent OSM from deduping so we get precise number of results
        params["dedupe"] = 0
        params["limit"] = limit

        if isinstance(query, str):
            params["q"] = query
        elif isinstance(query, dict):
            # add query keys in alphabetical order so URL is the same string
            # each time, for caching purposes
            for key in sorted(query):
                params[key] = query[key]
        else:  # pragma: no cover
            msg = "query must be a dict or a string"
            raise TypeError(msg)

    # request the URL, return the JSON
    return _nominatim_request(params=params, request_type=request_type)


def _nominatim_request(
    params: OrderedDict[str, int | str],
    request_type: str = "search",
    pause: float = 1,
    error_pause: float = 60,
) -> list[dict[str, Any]]:
    """
    Send a HTTP GET request to the Nominatim API and return response.

    Parameters
    ----------
    params : OrderedDict
        key-value pairs of parameters
    request_type : string {"search", "reverse", "lookup"}
        which Nominatim API endpoint to query
    pause : float
        how long to pause before request, in seconds. per the nominatim usage
        policy: "an absolute maximum of 1 request per second" is allowed
    error_pause : float
        how long to pause in seconds before re-trying request if error

    Returns
    -------
    response_json : list
    """
    if request_type not in {"search", "reverse", "lookup"}:  # pragma: no cover
        msg = 'Nominatim request_type must be "search", "reverse", or "lookup"'
        raise ValueError(msg)

    # add nominatim API key to params if one has been provided in settings
    if settings.nominatim_key is not None:
        params["key"] = settings.nominatim_key

    # prepare Nominatim API URL and see if request already exists in cache
    url = settings.nominatim_endpoint.rstrip("/") + "/" + request_type
    prepared_url = str(requests.Request("GET", url, params=params).prepare().url)

    cached_response_json: list[dict[str, Any]] | None = _downloader._retrieve_from_cache(
        prepared_url
    )  # type: ignore[assignment]
    if cached_response_json is not None:
        return cached_response_json

    # pause then request this URL
    domain = _downloader._hostname_from_url(url)
    utils.log(f"Pausing {pause} second(s) before making HTTP GET request to {domain!r}")
    time.sleep(pause)

    # transmit the HTTP GET request
    utils.log(f"Get {prepared_url} with timeout={settings.timeout}")
    response = requests.get(
        url,
        params=params,
        timeout=settings.timeout,
        headers=_downloader._get_http_headers(),
        **settings.requests_kwargs,
    )

    # handle 429 and 504 errors by pausing then recursively re-trying request
    if response.status_code in {429, 504}:  # pragma: no cover
        msg = (
            f"{domain!r} responded {response.status_code} {response.reason}: "
            f"we'll retry in {error_pause} secs"
        )
        utils.log(msg, level=lg.WARNING)
        time.sleep(error_pause)
        return _nominatim_request(params, request_type, pause, error_pause)

    response_json: list[dict[str, Any]] = _downloader._parse_response(response)  # type: ignore[assignment]
    _downloader._save_to_cache(prepared_url, response_json, response.ok)
    if not isinstance(response_json, list):
        msg = "Nominatim API did not return a list of results."
        raise TypeError(msg)
    return response_json
