"""
Global settings that can be configured by the user.

all_oneway : bool
    Only use if subsequently saving graph to an OSM XML file via the
    `save_graph_xml` function. If True, forces all ways to be added as one-way
    ways, preserving the original order of the nodes in the OSM way. This also
    retains the original OSM way's oneway tag's string value as edge attribute
    values, rather than converting them to True/False bool values. Default is
    `False`.
bidirectional_network_types : list[str]
    Network types for which a fully bidirectional graph will be created.
    Default is `["walk"]`.
cache_folder : str | Path
    Path to folder to save/load HTTP response cache files, if the `use_cache`
    setting is True. Default is `"./cache"`.
cache_only_mode : bool
    If True, download network data from Overpass then raise a
    `CacheOnlyModeInterrupt` error for user to catch. This prevents graph
    building from taking place and instead just saves Overpass response to
    cache. Useful for sequentially caching lots of raw data (as you can
    only query Overpass one request at a time) then using the local cache to
    quickly build many graphs simultaneously with multiprocessing. Default is
    `False`.
data_folder : str | Path
    Path to folder to save/load graph files by default. Default is `"./data"`.
default_access : str
    Filter for the OSM "access" tag. Default is `'["access"!~"private"]'`.
    Note that also filtering out "access=no" ways prevents including
    transit-only bridges (e.g., Tilikum Crossing) from appearing in drivable
    road network (e.g., `'["access"!~"private|no"]'`). However, some drivable
    tollroads have "access=no" plus a "access:conditional" tag to clarify when
    it is accessible, so we can't filter out all "access=no" ways by default.
    Best to be permissive here then remove complicated combinations of tags
    programatically after the full graph is downloaded and constructed.
default_crs : str
    Default coordinate reference system to set when creating graphs. Default
    is `"epsg:4326"`.
doh_url_template : str | None
    Endpoint to resolve DNS-over-HTTPS if local DNS resolution fails. Set to
    None to disable DoH, but see `downloader._config_dns` documentation for
    caveats. Default is: `"https://8.8.8.8/resolve?name={hostname}"`
elevation_url_template : str
    Endpoint of the Google Maps Elevation API (or equivalent), containing
    exactly two parameters: `locations` and `key`. Default is:
    `"https://maps.googleapis.com/maps/api/elevation/json?locations={locations}&key={key}"`
    One example of an alternative equivalent would be Open Topo Data:
    `"https://api.opentopodata.org/v1/aster30m?locations={locations}&key={key}"`
http_accept_language : str
    HTTP header accept-language. Default is `"en"`. Note that Nominatim's
    default language is "en" and it may sort its results' importance scores
    differently if a different language is specified.
http_referer : str
    HTTP header referer. Default is
    `"OSMnx Python package (https://github.com/gboeing/osmnx)"`.
http_user_agent : str
    HTTP header user-agent. Default is
    `"OSMnx Python package (https://github.com/gboeing/osmnx)"`.
imgs_folder : str | Path
    Path to folder in which to save plotted images by default. Default is
    `"./images"`.
log_file : bool
    If True, save log output to a file in `logs_folder`. Default is `False`.
log_filename : str
    Name of the log file, without file extension. Default is `"osmnx"`.
log_console : bool
    If True, print log output to the console (terminal window). Default is
    `False`.
log_level : int
    One of Python's `logger.level` constants. Default is `logging.INFO`.
log_name : str
    Name of the logger. Default is `"OSMnx"`.
logs_folder : str | Path
    Path to folder in which to save log files. Default is `"./logs"`.
max_query_area_size : float
    Maximum area for any part of the geometry in meters: any polygon bigger
    than this will get divided up for multiple queries to the API. Default is
    `2500000000`.
nominatim_key : str | None
    Your Nominatim API key, if you are using an API instance that requires
    one. Default is `None`.
nominatim_url : str
    The base API url to use for Nominatim queries. Default is
    `"https://nominatim.openstreetmap.org/"`.
overpass_memory : int | None
    Overpass server memory allocation size for the query, in bytes. If
    None, server will choose its default allocation size. Use with caution.
    Default is `None`.
overpass_rate_limit : bool
    If True, check the Overpass server status endpoint for how long to
    pause before making request. Necessary if server uses slot management,
    but can be set to False if you are running your own Overpass instance
    without rate limiting. Default is `True`.
overpass_settings : str
    Settings string for Overpass queries. Default is
    `"[out:json][timeout:{timeout}]{maxsize}"`. By default, the {timeout} and
    {maxsize} values are set dynamically by OSMnx when used.
    To query, for example, historical OSM data as of a certain date:
    `'[out:json][timeout:90][date:"2019-10-28T19:20:00Z"]'`. Use with caution.
overpass_url : str
    The base API url to use for Overpass queries. Default is
    `"https://overpass-api.de/api"`.
requests_kwargs : dict[str, Any]
    Optional keyword args to pass to the requests package when connecting
    to APIs, for example to configure authentication or provide a path to
    a local certificate file. More info on options such as auth, cert,
    verify, and proxies can be found in the requests package advanced docs.
    Default is `{}`.
requests_timeout : int
    The timeout interval in seconds for HTTP requests, and (when applicable)
    for Overpass server to use for executing the query. Default is `180`.
use_cache : bool
    If True, cache HTTP responses locally in `cache_folder` instead of calling
    API repeatedly for the same request. Default is `True`.
useful_tags_node : list[str]
    OSM "node" tags to add as graph node attributes, when present in the data
    retrieved from OSM. Default is `["highway", "junction", "railway", "ref"]`.
useful_tags_way : list[str]
    OSM "way" tags to add as graph edge attributes, when present in the data
    retrieved from OSM. Default is `["access", "area", "bridge", "est_width",
    "highway", "junction", "landuse", "lanes", "maxspeed", "name", "oneway",
    "ref", "service", "tunnel", "width"]`.
implicit_maxspeed_values: dict
    Implicit maxspeed values mapping that is used by `add_edge_speeds`
    function. Default values are obtained from wiki page
    https://wiki.openstreetmap.org/wiki/Key:maxspeed.
"""

from __future__ import annotations

import logging as lg
from typing import TYPE_CHECKING
from typing import Any

if TYPE_CHECKING:
    from pathlib import Path

all_oneway: bool = False
bidirectional_network_types: list[str] = ["walk"]
cache_folder: str | Path = "./cache"
cache_only_mode: bool = False
data_folder: str | Path = "./data"
default_access: str = '["access"!~"private"]'
default_crs: str = "epsg:4326"
doh_url_template: str | None = "https://8.8.8.8/resolve?name={hostname}"
elevation_url_template: str = (
    "https://maps.googleapis.com/maps/api/elevation/json?locations={locations}&key={key}"
)
http_accept_language: str = "en"
http_referer: str = "OSMnx Python package (https://github.com/gboeing/osmnx)"
http_user_agent: str = "OSMnx Python package (https://github.com/gboeing/osmnx)"
imgs_folder: str | Path = "./images"
log_console: bool = False
log_file: bool = False
log_filename: str = "osmnx"
log_level: int = lg.INFO
log_name: str = "OSMnx"
logs_folder: str | Path = "./logs"
max_query_area_size: float = 50 * 1000 * 50 * 1000
nominatim_key: str | None = None
nominatim_url: str = "https://nominatim.openstreetmap.org/"
overpass_memory: int | None = None
overpass_rate_limit: bool = True
overpass_settings: str = "[out:json][timeout:{timeout}]{maxsize}"
overpass_url: str = "https://overpass-api.de/api"
requests_kwargs: dict[str, Any] = {}
requests_timeout: float = 180
use_cache: bool = True
useful_tags_node: list[str] = ["highway", "junction", "railway", "ref"]
useful_tags_way: list[str] = [
    "access",
    "area",
    "bridge",
    "est_width",
    "highway",
    "junction",
    "landuse",
    "lanes",
    "maxspeed",
    "name",
    "oneway",
    "ref",
    "service",
    "tunnel",
    "width",
]
implicit_maxspeed_values = {
    "AR:rural": 110.0,
    "AR:urban": 40.0,
    "AR:urban:primary": 60.0,
    "AR:urban:secondary": 60.0,
    "AT:bicycle_road": 30.0,
    "AT:motorway": 130.0,
    "AT:rural": 100.0,
    "AT:trunk": 100.0,
    "AT:urban": 50.0,
    "BE-BRU:rural": 70.0,
    "BE-BRU:urban": 30.0,
    "BE-VLG:rural": 70.0,
    "BE-VLG:urban": 50.0,
    "BE-WAL:rural": 90.0,
    "BE-WAL:urban": 50.0,
    "BE:cyclestreet": 30.0,
    "BE:living_street": 20.0,
    "BE:motorway": 120.0,
    "BE:trunk": 120.0,
    "BE:zone30": 30.0,
    "BG:living_street": 20.0,
    "BG:motorway": 140.0,
    "BG:rural": 90.0,
    "BG:trunk": 120.0,
    "BG:urban": 50.0,
    "BY:living_street": 20.0,
    "BY:motorway": 110.0,
    "BY:rural": 90.0,
    "BY:urban": 60.0,
    "CA-AB:rural": 90.0,
    "CA-AB:urban": 65.0,
    "CA-BC:rural": 80.0,
    "CA-BC:urban": 50.0,
    "CA-MB:rural": 90.0,
    "CA-MB:urban": 50.0,
    "CA-ON:rural": 80.0,
    "CA-ON:urban": 50.0,
    "CA-QC:motorway": 100.0,
    "CA-QC:rural": 75.0,
    "CA-QC:urban": 50.0,
    "CA-SK:nsl": 80.0,
    "CH:motorway": 120.0,
    "CH:rural": 80.0,
    "CH:trunk": 100.0,
    "CH:urban": 50.0,
    "CZ:living_street": 20.0,
    "CZ:motorway": 130.0,
    "CZ:pedestrian_zone": 20.0,
    "CZ:rural": 90.0,
    "CZ:trunk": 110.0,
    "CZ:urban": 50.0,
    "CZ:urban_motorway": 80.0,
    "CZ:urban_trunk": 80.0,
    "DE:bicycle_road": 30.0,
    "DE:living_street": 15.0,
    "DE:motorway": 120.0,
    "DE:rural": 80.0,
    "DE:urban": 50.0,
    "DK:motorway": 130.0,
    "DK:rural": 80.0,
    "DK:urban": 50.0,
    "EE:rural": 90.0,
    "EE:urban": 50.0,
    "ES:living_street": 20.0,
    "ES:motorway": 120.0,
    "ES:rural": 90.0,
    "ES:trunk": 90.0,
    "ES:urban": 50.0,
    "ES:zone30": 30.0,
    "FI:motorway": 120.0,
    "FI:rural": 80.0,
    "FI:trunk": 100.0,
    "FI:urban": 50.0,
    "FR:motorway": 120.0,
    "FR:rural": 80.0,
    "FR:urban": 50.0,
    "FR:zone30": 30.0,
    "GB:nsl_restricted": 48.28,
    "GR:motorway": 130.0,
    "GR:rural": 90.0,
    "GR:trunk": 110.0,
    "GR:urban": 50.0,
    "HU:living_street": 20.0,
    "HU:motorway": 130.0,
    "HU:rural": 90.0,
    "HU:trunk": 110.0,
    "HU:urban": 50.0,
    "IT:motorway": 130.0,
    "IT:rural": 90.0,
    "IT:trunk": 110.0,
    "IT:urban": 50.0,
    "JP:express": 100.0,
    "JP:nsl": 60.0,
    "LT:rural": 90.0,
    "LT:urban": 50.0,
    "NO:rural": 80.0,
    "NO:urban": 50.0,
    "PH:express": 100.0,
    "PH:rural": 80.0,
    "PH:urban": 30.0,
    "PT:motorway": 120.0,
    "PT:rural": 90.0,
    "PT:trunk": 100.0,
    "PT:urban": 50.0,
    "RO:motorway": 130.0,
    "RO:rural": 90.0,
    "RO:trunk": 100.0,
    "RO:urban": 50.0,
    "RS:living_street": 10.0,
    "RS:motorway": 130.0,
    "RS:rural": 80.0,
    "RS:trunk": 100.0,
    "RS:urban": 50.0,
    "RU:living_street": 20.0,
    "RU:motorway": 110.0,
    "RU:rural": 90.0,
    "RU:urban": 60.0,
    "SE:rural": 70.0,
    "SE:urban": 50.0,
    "SI:motorway": 130.0,
    "SI:rural": 90.0,
    "SI:trunk": 110.0,
    "SI:urban": 50.0,
    "SK:living_street": 20.0,
    "SK:motorway": 130.0,
    "SK:motorway_urban": 90.0,
    "SK:rural": 90.0,
    "SK:trunk": 90.0,
    "SK:urban": 50.0,
    "TR:living_street": 20.0,
    "TR:motorway": 130.0,
    "TR:rural": 90.0,
    "TR:trunk": 110.0,
    "TR:urban": 50.0,
    "TR:zone30": 30.0,
    "UA:living_street": 20.0,
    "UA:motorway": 130.0,
    "UA:rural": 90.0,
    "UA:trunk": 110.0,
    "UA:urban": 50.0,
    "UK:motorway": 112.65,
    "UK:nsl_dual": 112.65,
    "UK:nsl_single": 96.56,
    "UZ:living_street": 30.0,
    "UZ:motorway": 110.0,
    "UZ:rural": 100.0,
    "UZ:urban": 70.0,
}
