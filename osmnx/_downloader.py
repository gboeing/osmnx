"""Handle HTTP requests to web APIs."""

import json
import logging as lg
import socket
from hashlib import sha1
from pathlib import Path
from urllib.parse import urlparse
from warnings import warn

import requests
from requests.exceptions import JSONDecodeError

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
    if settings.default_accept_language is None:
        default_accept_language = settings.http_accept_language
    else:
        default_accept_language = settings.default_accept_language
        msg = (
            "`settings.default_accept_language` is deprecated and will be removed "
            "in the v2.0.0 release: use `settings.http_accept_language` instead. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
        )
        warn(msg, FutureWarning, stacklevel=2)

    if settings.default_referer is None:
        default_referer = settings.http_referer
    else:
        default_referer = settings.default_referer
        msg = (
            "`settings.default_referer` is deprecated and will be removed in the "
            "v2.0.0 release: use `settings.http_referer` instead. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
        )
        warn(msg, FutureWarning, stacklevel=2)

    if settings.default_user_agent is None:
        default_user_agent = settings.http_user_agent
    else:
        default_user_agent = settings.default_user_agent
        msg = (
            "`settings.default_user_agent` is deprecated and will be removed in "
            "the v2.0.0 release: use `settings.http_user_agent` instead. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
        )
        warn(msg, FutureWarning, stacklevel=2)

    if user_agent is None:
        user_agent = default_user_agent
    if referer is None:
        referer = default_referer
    if accept_language is None:
        accept_language = default_accept_language

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

    if settings.doh_url_template is None:
        # if user has set the url template to None, return hostname itself
        utils.log("User set `doh_url_template=None`, requesting host by name", level=lg.WARNING)
        return hostname

    err_msg = f"Failed to resolve {hostname!r} IP via DoH, requesting host by name"
    try:
        url = settings.doh_url_template.format(hostname=hostname)
        response = requests.get(url, timeout=timeout)
        data = response.json()

    # if we cannot reach DoH server or resolve host, return hostname itself
    except requests.exceptions.RequestException:  # pragma: no cover
        utils.log(err_msg, level=lg.ERROR)
        return hostname

    # if there were no exceptions, return
    else:
        if response.ok and data["Status"] == 0:
            # status 0 means NOERROR, so return the IP address
            return data["Answer"][0]["data"]

        # otherwise, if we cannot reach DoH server or cannot resolve host
        # just return the hostname itself
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

        # otherwise
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
    utils.log(f"Downloaded {size_kb:,.1f}kB from {domain!r} with status {response.status_code}")

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
