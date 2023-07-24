"""Interact with web APIs."""

import json
import logging as lg
import socket
import time
from collections import OrderedDict
from hashlib import sha1
from pathlib import Path
from urllib.parse import urlparse

import requests
from requests.exceptions import JSONDecodeError

from . import _overpass
from . import settings
from . import utils
from ._errors import InsufficientResponseError
from ._errors import ResponseStatusCodeError

# capture getaddrinfo function to use original later after mutating it
_original_getaddrinfo = socket.getaddrinfo


def _save_to_cache(url, response_json, ok):
    """
    Save a HTTP response JSON object to a file in the cache folder.

    Function calculates the checksum of url to generate the cache file's name.
    If the request was sent to server via POST instead of GET, then URL should
    be a GET-style representation of request. Response is only saved to a
    cache file if settings.use_cache is True, response_json is not None, and
    ok is True.

    Users should always pass OrderedDicts instead of dicts of parameters into
    request functions, so the parameters remain in the same order each time,
    producing the same URL string, and thus the same hash. Otherwise the cache
    will eventually contain multiple saved responses for the same request
    because the URL's parameters appeared in a different order each time.

    Parameters
    ----------
    url : string
        the URL of the request
    response_json : dict
        the JSON response
    ok : bool
        requests response.ok value

    Returns
    -------
    None
    """
    if settings.use_cache:
        if not ok:  # pragma: no cover
            utils.log("Did not save to cache because response status_code is not OK")

        elif response_json is None:  # pragma: no cover
            utils.log("Did not save to cache because response_json is None")

        else:
            # create the folder on the disk if it doesn't already exist
            cache_folder = Path(settings.cache_folder)
            cache_folder.mkdir(parents=True, exist_ok=True)

            # hash the url to make the filename succinct but unique
            # sha1 digest is 160 bits = 20 bytes = 40 hexadecimal characters
            filename = sha1(url.encode("utf-8")).hexdigest() + ".json"
            cache_filepath = cache_folder / filename

            # dump to json, and save to file
            cache_filepath.write_text(json.dumps(response_json), encoding="utf-8")
            utils.log(f"Saved response to cache file {str(cache_filepath)!r}")


def _url_in_cache(url):
    """
    Determine if a URL's response exists in the cache.

    Calculates the checksum of url to determine the cache file's name.

    Parameters
    ----------
    url : string
        the URL to look for in the cache

    Returns
    -------
    filepath : pathlib.Path
        path to cached response for url if it exists, otherwise None
    """
    # hash the url to generate the cache filename
    filename = sha1(url.encode("utf-8")).hexdigest() + ".json"
    filepath = Path(settings.cache_folder) / filename

    # if this file exists in the cache, return its full path
    return filepath if filepath.is_file() else None


def _retrieve_from_cache(url, check_remark=True):
    """
    Retrieve a HTTP response JSON object from the cache, if it exists.

    Parameters
    ----------
    url : string
        the URL of the request
    check_remark : string
        if True, only return filepath if cached response does not have a
        remark key indicating a server warning

    Returns
    -------
    response_json : dict
        cached response for url if it exists in the cache, otherwise None
    """
    # if the tool is configured to use the cache
    if settings.use_cache:
        # return cached response for this url if exists, otherwise return None
        cache_filepath = _url_in_cache(url)
        if cache_filepath is not None:
            response_json = json.loads(cache_filepath.read_text(encoding="utf-8"))

            # return None if check_remark is True and there is a server
            # remark in the cached response
            if check_remark and "remark" in response_json:  # pragma: no cover
                utils.log(
                    f"Ignoring cache file {str(cache_filepath)!r} because "
                    f"it contains a remark: {response_json['remark']!r}"
                )
                return None

            utils.log(f"Retrieved response from cache file {str(cache_filepath)!r}")
            return response_json
    return None


def _get_http_headers(user_agent=None, referer=None, accept_language=None):
    """
    Update the default requests HTTP headers with OSMnx info.

    Parameters
    ----------
    user_agent : string
        the user agent string, if None will set with OSMnx default
    referer : string
        the referer string, if None will set with OSMnx default
    accept_language : string
        make accept-language explicit e.g. for consistent nominatim result
        sorting

    Returns
    -------
    headers : dict
    """
    if user_agent is None:
        user_agent = settings.default_user_agent
    if referer is None:
        referer = settings.default_referer
    if accept_language is None:
        accept_language = settings.default_accept_language

    headers = requests.utils.default_headers()
    headers.update(
        {"User-Agent": user_agent, "referer": referer, "Accept-Language": accept_language}
    )
    return headers


def _resolve_host_via_doh(hostname):
    """
    Resolve hostname to IP address via Google's public DNS-over-HTTPS API.

    Necessary fallback as socket.gethostbyname will not always work when using
    a proxy. See https://developers.google.com/speed/public-dns/docs/doh/json
    If the user has set `settings.doh_url_template=None` or if resolution
    fails (e.g., due to local network blocking DNS-over-HTTPS) the hostname
    itself will be returned instead. Note that this means that server slot
    management may be violated: see `_config_dns` documentation for details.

    Parameters
    ----------
    hostname : string
        the hostname to consistently resolve the IP address of

    Returns
    -------
    ip_address : string
        resolved IP address of host, or hostname itself if resolution failed
    """
    if settings.doh_url_template is None:
        # if user has set the url template to None, return hostname itself
        utils.log("User set `doh_url_template=None`, requesting host by name", level=lg.WARNING)
        return hostname

    err_msg = f"Failed to resolve {hostname!r} IP via DoH, requesting host by name"
    try:
        url = settings.doh_url_template.format(hostname=hostname)
        response = requests.get(url, timeout=settings.timeout)
        data = response.json()
        if response.ok and data["Status"] == 0:
            # status 0 means NOERROR, so return the IP address
            return data["Answer"][0]["data"]
        else:  # pragma: no cover
            # if we cannot reach DoH server or cannot resolve host, return hostname itself
            utils.log(err_msg, level=lg.ERROR)
            return hostname

    # if we cannot reach DoH server or cannot resolve host, return hostname itself
    except requests.exceptions.RequestException:  # pragma: no cover
        utils.log(err_msg, level=lg.ERROR)
        return hostname


def _config_dns(url):
    """
    Force socket.getaddrinfo to use IP address instead of hostname.

    Resolves the URL's domain to an IP address so that we use the same server
    for both 1) checking the necessary pause duration and 2) sending the query
    itself even if there is round-robin redirecting among multiple server
    machines on the server-side. Mutates the getaddrinfo function so it uses
    the same IP address everytime it finds the hostname in the URL.

    For example, the server overpass-api.de just redirects to one of the other
    servers (currently gall.openstreetmap.de and lambert.openstreetmap.de). So
    if we check the status endpoint of overpass-api.de, we may see results for
    server gall, but when we submit the query itself it gets redirected to
    server lambert. This could result in violating server lambert's slot
    management timing.

    Parameters
    ----------
    url : string
        the URL to consistently resolve the IP address of

    Returns
    -------
    None
    """
    hostname = _hostname_from_url(url)
    try:
        ip = socket.gethostbyname(hostname)
    except socket.gaierror:  # pragma: no cover
        # may occur when using a proxy, so instead resolve IP address via DoH
        utils.log(
            f"Encountered gaierror while trying to resolve {hostname!r}, trying again via DoH...",
            level=lg.ERROR,
        )
        ip = _resolve_host_via_doh(hostname)

    # mutate socket.getaddrinfo to map hostname -> IP address
    def _getaddrinfo(*args, **kwargs):
        if args[0] == hostname:
            utils.log(f"Resolved {hostname!r} to {ip!r}")
            return _original_getaddrinfo(ip, *args[1:], **kwargs)
        else:
            return _original_getaddrinfo(*args, **kwargs)

    socket.getaddrinfo = _getaddrinfo


def _hostname_from_url(url):
    """
    Extract the hostname (domain) from a URL.

    Parameters
    ----------
    url : string
        the url from which to extract the hostname

    Returns
    -------
    hostname : string
        the extracted hostname (domain)
    """
    return urlparse(url).netloc.split(":")[0]


def _parse_response(response):
    """
    Parse JSON from a requests response and log the details.

    Parameters
    ----------
    response : requests.response
        the response object

    Returns
    -------
    response_json : dict
    """
    # log the response size and domain
    domain = _hostname_from_url(response.url)
    size_kb = len(response.content) / 1000
    utils.log(f"Downloaded {size_kb:,.1f}kB from {domain!r} with code {response.status_code}")

    # parse the response to JSON and log/raise exceptions
    try:
        response_json = response.json()
    except JSONDecodeError as e:  # pragma: no cover
        msg = f"{domain!r} responded: {response.status_code} {response.reason} {response.text}"
        utils.log(msg, level=lg.ERROR)
        if response.ok:
            raise InsufficientResponseError(msg) from e
        raise ResponseStatusCodeError(msg) from e

    # log any remarks if they exist
    if "remark" in response_json:  # pragma: no cover
        utils.log(f'{domain!r} remarked: {response_json["remark"]!r}', level=lg.WARNING)

    return response_json


def _retrieve_nominatim_element(query, by_osmid=False, limit=1, polygon_geojson=1):
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
    polygon_geojson : int
        retrieve the place's geometry from the API, 0=no, 1=yes

    Returns
    -------
    response_json : dict
        JSON response from the Nominatim server
    """
    # define the parameters
    params = OrderedDict()
    params["format"] = "json"
    params["polygon_geojson"] = polygon_geojson

    if by_osmid:
        # if querying by OSM ID, use the lookup endpoint
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


def _nominatim_request(params, request_type="search", pause=1, error_pause=60):
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
    response_json : dict
    """
    if request_type not in {"search", "reverse", "lookup"}:  # pragma: no cover
        msg = 'Nominatim request_type must be "search", "reverse", or "lookup"'
        raise ValueError(msg)

    # prepare Nominatim API URL and see if request already exists in cache
    url = settings.nominatim_endpoint.rstrip("/") + "/" + request_type
    params["key"] = settings.nominatim_key
    prepared_url = requests.Request("GET", url, params=params).prepare().url
    cached_response_json = _retrieve_from_cache(prepared_url)
    if cached_response_json is not None:
        return cached_response_json

    # pause then request this URL
    domain = _hostname_from_url(url)
    utils.log(f"Pausing {pause} second(s) before making HTTP GET request to {domain!r}")
    time.sleep(pause)

    # transmit the HTTP GET request
    utils.log(f"Get {prepared_url} with timeout={settings.timeout}")
    response = requests.get(
        url,
        params=params,
        timeout=settings.timeout,
        headers=_get_http_headers(),
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

    response_json = _parse_response(response)
    _save_to_cache(prepared_url, response_json, response.status_code)
    return response_json


def _overpass_request(data, pause=None, error_pause=60):
    """
    Send a HTTP POST request to the Overpass API and return response.

    Parameters
    ----------
    data : OrderedDict
        key-value pairs of parameters
    pause : float
        how long to pause in seconds before request, if None, will query API
        status endpoint to find when next slot is available
    error_pause : float
        how long to pause in seconds (in addition to `pause`) before re-trying
        request if error

    Returns
    -------
    response_json : dict
    """
    # resolve url to same IP even if there is server round-robin redirecting
    _config_dns(settings.overpass_endpoint)

    # prepare the Overpass API URL and see if request already exists in cache
    url = settings.overpass_endpoint.rstrip("/") + "/interpreter"
    prepared_url = requests.Request("GET", url, params=data).prepare().url
    cached_response_json = _retrieve_from_cache(prepared_url)
    if cached_response_json is not None:
        return cached_response_json

    # pause then request this URL
    if pause is None:
        this_pause = _overpass._get_overpass_pause(settings.overpass_endpoint)
    domain = _hostname_from_url(url)
    utils.log(f"Pausing {this_pause} second(s) before making HTTP POST request to {domain!r}")
    time.sleep(this_pause)

    # transmit the HTTP POST request
    utils.log(f"Post {prepared_url} with timeout={settings.timeout}")
    response = requests.post(
        url,
        data=data,
        timeout=settings.timeout,
        headers=_get_http_headers(),
        **settings.requests_kwargs,
    )

    # handle 429 and 504 errors by pausing then recursively re-trying request
    if response.status_code in {429, 504}:  # pragma: no cover
        this_pause = error_pause + _overpass._get_overpass_pause(settings.overpass_endpoint)
        msg = (
            f"{domain!r} responded {response.status_code} {response.reason}: "
            f"we'll retry in {this_pause} secs"
        )
        utils.log(msg, level=lg.WARNING)
        time.sleep(this_pause)
        return _overpass_request(data, pause, error_pause)

    response_json = _parse_response(response)
    _save_to_cache(prepared_url, response_json, response.status_code)
    return response_json


def _google_request(url, pause):
    """
    Send a HTTP GET request to Google Maps Elevation API and return response.

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
    # check if request already exists in cache
    cached_response_json = _retrieve_from_cache(url)
    if cached_response_json is not None:
        return cached_response_json

    # pause then request this URL
    domain = _hostname_from_url(url)
    utils.log(f"Pausing {pause} second(s) before making HTTP GET request to {domain!r}")
    time.sleep(pause)

    # transmit the HTTP GET request
    utils.log(f"Get {url} with timeout={settings.timeout}")
    response = requests.get(
        url, timeout=settings.timeout, headers=_get_http_headers(), **settings.requests_kwargs
    )

    response_json = _parse_response(response)
    _save_to_cache(url, response_json, response.status_code)
    return response_json
