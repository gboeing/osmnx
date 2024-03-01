"""Tools to work with the Nominatim API."""

from __future__ import annotations

import logging as lg
import time
from collections import OrderedDict
from typing import Any

import requests

from . import _http
from . import settings
from . import utils
from ._errors import InsufficientResponseError


def _download_nominatim_element(
    query: str | dict[str, str],
    *,
    by_osmid: bool = False,
    limit: int = 1,
    polygon_geojson: bool = True,
) -> list[dict[str, Any]]:
    """
    Retrieve an OSM element from the Nominatim API.

    Parameters
    ----------
    query
        Query string or structured query dict.
    by_osmid
        If True, treat `query` as an OSM ID lookup rather than text search.
    limit
        Max number of results to return.
    polygon_geojson
        Whether to retrieve the place's geometry from the API.

    Returns
    -------
    response_json
    """
    # define the parameters
    params: OrderedDict[str, int | str] = OrderedDict()
    params["format"] = "json"
    params["polygon_geojson"] = int(polygon_geojson)  # bool -> int

    if by_osmid:
        # if querying by OSM ID, use the lookup endpoint
        if not isinstance(query, str):
            msg = "`query` must be a string if `by_osmid` is True."
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
            msg = "Each query must be a dict or a string."  # type: ignore[unreachable]
            raise TypeError(msg)

    # request the URL, return the JSON
    return _nominatim_request(params=params, request_type=request_type)


def _nominatim_request(
    params: OrderedDict[str, int | str],
    *,
    request_type: str = "search",
    pause: float = 1,
    error_pause: float = 60,
) -> list[dict[str, Any]]:
    """
    Send a HTTP GET request to the Nominatim API and return response.

    Parameters
    ----------
    params
        Key-value pairs of parameters.
    request_type
        {"search", "reverse", "lookup"}
        Which Nominatim API endpoint to query.
    pause
        How long to pause before request, in seconds. Per the Nominatim usage
        policy: "an absolute maximum of 1 request per second" is allowed.
    error_pause
        How long to pause in seconds before re-trying request if error.

    Returns
    -------
    response_json
    """
    if request_type not in {"search", "reverse", "lookup"}:  # pragma: no cover
        msg = "Nominatim `request_type` must be 'search', 'reverse', or 'lookup'."
        raise ValueError(msg)

    # add nominatim API key to params if one has been provided in settings
    if settings.nominatim_key is not None:
        params["key"] = settings.nominatim_key

    # prepare Nominatim API URL and see if request already exists in cache
    url = settings.nominatim_url.rstrip("/") + "/" + request_type
    prepared_url = str(requests.Request("GET", url, params=params).prepare().url)
    cached_response_json = _http._retrieve_from_cache(prepared_url)
    if isinstance(cached_response_json, list):
        return cached_response_json

    # pause then request this URL
    domain = _http._hostname_from_url(url)
    msg = f"Pausing {pause} second(s) before making HTTP GET request to {domain!r}"
    utils.log(msg, level=lg.INFO)
    time.sleep(pause)

    # transmit the HTTP GET request
    msg = f"Get {prepared_url} with timeout={settings.requests_timeout}"
    utils.log(msg, level=lg.INFO)
    response = requests.get(
        url,
        params=params,
        timeout=settings.requests_timeout,
        headers=_http._get_http_headers(),
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
        return _nominatim_request(
            params,
            request_type=request_type,
            pause=pause,
            error_pause=error_pause,
        )

    response_json = _http._parse_response(response)
    if not isinstance(response_json, list):
        msg = "Nominatim API did not return a list of results."
        raise InsufficientResponseError(msg)
    _http._save_to_cache(prepared_url, response_json, response.ok)
    return response_json
