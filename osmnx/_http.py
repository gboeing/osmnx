"""Handle HTTP requests to web APIs."""

from __future__ import annotations

import json
import logging as lg
import socket
from hashlib import sha1
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import requests
from requests.exceptions import JSONDecodeError

from . import settings
from . import utils
from ._errors import InsufficientResponseError
from ._errors import ResponseStatusCodeError

# capture getaddrinfo function to use original later after mutating it
_original_getaddrinfo = socket.getaddrinfo


def _save_to_cache(
    url: str,
    response_json: dict[str, Any] | list[dict[str, Any]],
    ok: bool,  # noqa: FBT001
) -> None:
    """
    Save a HTTP response JSON object to a file in the cache folder.

    This calculates the checksum of `url` to generate the cache file name. If
    the request was sent to server via POST instead of GET, then `url` should
    be a GET-style representation of the request. Response is only saved to a
    cache file if `settings.use_cache` is True, `response_json` is not None,
    and `ok` is True.

    Users should always pass OrderedDicts instead of dicts of parameters into
    request functions, so the parameters remain in the same order each time,
    producing the same URL string, and thus the same hash. Otherwise the cache
    will eventually contain multiple saved responses for the same request
    because the URL's parameters appeared in a different order each time.

    Parameters
    ----------
    url
        The URL of the request.
    response_json
        The JSON response from the server.
    ok
        A `requests.response.ok` value.

    Returns
    -------
    None
    """
    if settings.use_cache:
        if not ok:  # pragma: no cover
            msg = "Did not save to cache because HTTP status code is not OK"
            utils.log(msg, level=lg.WARNING)
        else:
            # create the folder on the disk if it doesn't already exist
            cache_folder = Path(settings.cache_folder)
            cache_folder.mkdir(parents=True, exist_ok=True)

            # hash the url to make the filename succinct but unique
            # sha1 digest is 160 bits = 20 bytes = 40 hexadecimal characters
            checksum = sha1(url.encode("utf-8")).hexdigest()  # noqa: S324
            cache_filepath = cache_folder / f"{checksum}.json"

            # dump to json, and save to file
            cache_filepath.write_text(json.dumps(response_json), encoding="utf-8")
            msg = f"Saved response to cache file {str(cache_filepath)!r}"
            utils.log(msg, level=lg.INFO)


def _url_in_cache(url: str) -> Path | None:
    """
    Determine if a URL's response exists in the cache.

    Calculates the checksum of `url` to determine the cache file's name.
    Returns None if it cannot be found in the cache.

    Parameters
    ----------
    url
        The URL to look for in the cache.

    Returns
    -------
    cache_filepath
        Path to cached response for `url` if it exists, otherwise None.
    """
    # hash the url to generate the cache filename
    checksum = sha1(url.encode("utf-8")).hexdigest()  # noqa: S324
    cache_filepath = Path(settings.cache_folder) / f"{checksum}.json"

    # if this file exists in the cache, return its full path
    return cache_filepath if cache_filepath.is_file() else None


def _retrieve_from_cache(url: str) -> dict[str, Any] | list[dict[str, Any]] | None:
    """
    Retrieve a HTTP response JSON object from the cache if it exists.

    Returns None if there is a server remark in the cached response.

    Parameters
    ----------
    url
        The URL of the request.

    Returns
    -------
    response_json
        Cached response for `url` if it exists in the cache and does not
        contain a server remark, otherwise None.
    """
    # if the tool is configured to use the cache
    if settings.use_cache:
        # return cached response for this url if exists, otherwise return None
        cache_filepath = _url_in_cache(url)
        if cache_filepath is not None:
            response_json: dict[str, Any] | list[dict[str, Any]] = json.loads(
                cache_filepath.read_text(encoding="utf-8"),
            )

            # return None if there is a server remark in the cached response
            if isinstance(response_json, dict) and ("remark" in response_json):  # pragma: no cover
                msg = (
                    f"Ignoring cache file {str(cache_filepath)!r} because "
                    f"it contains a remark: {response_json['remark']!r}"
                )
                utils.log(msg, lg.WARNING)
                return None

            msg = f"Retrieved response from cache file {str(cache_filepath)!r}"
            utils.log(msg, lg.INFO)
            return response_json

    return None


def _get_http_headers(
    *,
    user_agent: str | None = None,
    referer: str | None = None,
    accept_language: str | None = None,
) -> dict[str, str]:
    """
    Update the default requests HTTP headers with OSMnx information.

    Parameters
    ----------
    user_agent
        The user agent. If None, use `settings.http_user_agent` value.
    referer
        The referer. If None, use `settings.http_referer` value.
    accept_language
        The accept language. If None, use `settings.http_accept_language`
        value.

    Returns
    -------
    headers
    """
    if user_agent is None:
        user_agent = settings.http_user_agent
    if referer is None:
        referer = settings.http_referer
    if accept_language is None:
        accept_language = settings.http_accept_language

    info = {"User-Agent": user_agent, "referer": referer, "Accept-Language": accept_language}
    headers = dict(requests.utils.default_headers())
    headers.update(info)
    return headers


def _resolve_host_via_doh(hostname: str) -> str:
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
    hostname
        The hostname to consistently resolve the IP address of.

    Returns
    -------
    ip_address
        Resolved IP address of host, or hostname itself if resolution failed.
    """
    if settings.doh_url_template is None:
        # if user has set the url template to None, return hostname itself
        msg = "User set `doh_url_template=None`, requesting host by name"
        utils.log(msg, level=lg.WARNING)
        return hostname

    err_msg = f"Failed to resolve {hostname!r} IP via DoH, requesting host by name"
    try:
        url = settings.doh_url_template.format(hostname=hostname)
        response = requests.get(url, timeout=settings.requests_timeout)
        data = response.json()

    # if we cannot reach DoH server or resolve host, return hostname itself
    except requests.exceptions.RequestException:  # pragma: no cover
        utils.log(err_msg, level=lg.ERROR)
        return hostname

    # if there were no request exceptions, return
    else:
        if response.ok and data["Status"] == 0:
            # status 0 means NOERROR, so return the IP address
            ip_address: str = data["Answer"][0]["data"]
            return ip_address

        # otherwise, if we cannot reach DoH server or cannot resolve host
        # just return the hostname itself
        utils.log(err_msg, level=lg.ERROR)
        return hostname


def _config_dns(url: str) -> None:
    """
    Force socket.getaddrinfo to use IP address instead of hostname.

    Resolves URL's hostname to an IP address so that we use the same server
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
    url
        The URL to consistently resolve the IP address of.

    Returns
    -------
    None
    """
    hostname = _hostname_from_url(url)
    try:
        ip = socket.gethostbyname(hostname)
    except socket.gaierror:  # pragma: no cover
        # may occur when using a proxy, so instead resolve IP address via DoH
        msg = f"Encountered gaierror while trying to resolve {hostname!r}, trying again via DoH..."
        utils.log(msg, level=lg.ERROR)
        ip = _resolve_host_via_doh(hostname)

    # mutate socket.getaddrinfo to map hostname -> IP address
    def _getaddrinfo(*args: Any, **kwargs: Any) -> Any:  # noqa: ANN401
        if args[0] == hostname:
            msg = f"Resolved {hostname!r} to {ip!r}"
            utils.log(msg, level=lg.INFO)
            return _original_getaddrinfo(ip, *args[1:], **kwargs)

        # otherwise
        return _original_getaddrinfo(*args, **kwargs)

    socket.getaddrinfo = _getaddrinfo


def _hostname_from_url(url: str) -> str:
    """
    Extract the hostname (domain) from a URL.

    Parameters
    ----------
    url
        The url from which to extract the hostname.

    Returns
    -------
    hostname
        The extracted hostname (domain).
    """
    return urlparse(url).netloc.split(":")[0]


def _parse_response(response: requests.Response) -> dict[str, Any] | list[dict[str, Any]]:
    """
    Parse JSON from a requests response and log the details.

    Parameters
    ----------
    response
        The response object.

    Returns
    -------
    response_json
        Value will be a dict if the response is from the Google or Overpass
        APIs, and a list if the response is from the Nominatim API.
    """
    # log the response size and hostname
    hostname = _hostname_from_url(response.url)
    size_kb = len(response.content) / 1000
    msg = f"Downloaded {size_kb:,.1f}kB from {hostname!r} with status {response.status_code}"
    utils.log(msg, level=lg.INFO)

    # parse the response to JSON and log/raise exceptions
    try:
        response_json: dict[str, Any] | list[dict[str, Any]] = response.json()
    except JSONDecodeError as e:  # pragma: no cover
        msg = f"{hostname!r} responded: {response.status_code} {response.reason} {response.text}"
        utils.log(msg, level=lg.ERROR)
        if response.ok:
            raise InsufficientResponseError(msg) from e
        raise ResponseStatusCodeError(msg) from e

    # log any remarks if they exist
    if isinstance(response_json, dict) and "remark" in response_json:  # pragma: no cover
        msg = f"{hostname!r} remarked: {response_json['remark']!r}"
        utils.log(msg, level=lg.WARNING)

    # log if the response status_code is not OK
    if not response.ok:
        msg = f"{hostname!r} returned HTTP status code {response.status_code}"
        utils.log(msg, level=lg.WARNING)

    return response_json
