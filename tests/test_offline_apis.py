# ruff: noqa: D103, PLR2004, S101
# numpydoc ignore=GL08,PR01,RT01
"""Offline tests for cached and mocked API workflows."""

from __future__ import annotations

import json
import socket
import time as stdlib_time
from collections import OrderedDict
from pathlib import Path
from typing import cast

import pandas as pd
import pytest
import requests
from shapely import Point
from shapely import Polygon
from typeguard import suppress_type_checks

import osmnx as ox
from osmnx import _nominatim
from tests.helpers import ADDRESS
from tests.helpers import LOCATION_POINT
from tests.helpers import PLACE
from tests.helpers import TAGS
from tests.helpers import Response
from tests.helpers import drive_graph


@pytest.mark.unit
def test_cache_fixture_manifest(http_cache: Path) -> None:
    manifest = json.loads((http_cache / "manifest.json").read_text(encoding="utf-8"))

    assert manifest["expected"]["network_nodes"] == 5
    assert manifest["expected"]["network_ways"] == 4
    assert manifest["expected"]["feature_elements"] == 8
    assert len(manifest["records"]) >= 10
    assert all((http_cache / record["cache_file"]).is_file() for record in manifest["records"])
    assert all(record["endpoint"].startswith("https://") for record in manifest["records"])
    assert all(record["query"] for record in manifest["records"])


@pytest.mark.integration
def test_geocoder_uses_committed_cache(
    http_cache: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    ox.settings.cache_folder = http_cache

    city_by_id = ox.geocode_to_gdf("R2999176", by_osmid=True)
    city_by_query = ox.geocode_to_gdf(PLACE, which_result=1)
    city_list = ox.geocode_to_gdf([PLACE, PLACE], which_result=[1, 1])
    first_polygon = {"geojson": {"type": "Polygon"}}

    nonpolygon = {
        "boundingbox": ["0", "1", "2", "3"],
        "geojson": {"type": "Point", "coordinates": [2, 0]},
        "importance": 1,
        "lat": "0",
        "lon": "2",
        "name": "Fixture Point",
    }

    def _fake_download_nominatim_element(
        *_args: object,
        **_kwargs: object,
    ) -> list[dict[str, object]]:
        return [nonpolygon]

    with monkeypatch.context() as m:
        m.setattr(
            _nominatim,
            "_download_nominatim_element",
            _fake_download_nominatim_element,
        )
        point_result = ox.geocoder._geocode_query_to_gdf(
            "Fixture Point",
            which_result=1,
            by_osmid=False,
        )
    city_projected = ox.projection.project_gdf(city_by_query, to_crs="epsg:3395")

    assert city_by_id.geometry.iloc[0].equals_exact(city_by_query.geometry.iloc[0], tolerance=0)
    assert len(city_list) == 2
    assert point_result.geometry.iloc[0].geom_type == "Point"
    assert (
        ox.geocoder._get_first_polygon([{"geojson": {"type": "Point"}}, first_polygon])
        == first_polygon
    )
    assert city_by_query.loc[0, "name"] == "Piedmont Fixture"
    assert city_projected.crs == "epsg:3395"
    assert ox.geocode(ADDRESS) == LOCATION_POINT

    with pytest.raises(ox._errors.InsufficientResponseError):
        ox.geocode("!@#$%^&*")
    with pytest.raises(ox._errors.InsufficientResponseError):
        ox.geocode_to_gdf(query="AAAZZZ")
    with pytest.raises(TypeError, match="did not geocode"):
        ox.geocode_to_gdf("Bunker Hill, Los Angeles, California, USA")


@pytest.mark.integration
def test_graph_downloaders_use_committed_cache(http_cache: Path) -> None:
    ox.settings.cache_folder = http_cache

    bbox = ox.utils_geo.bbox_from_point(LOCATION_POINT, dist=500)
    polygon = ox.geocode_to_gdf(PLACE, which_result=1).geometry.iloc[0]

    graphs = [
        ox.graph_from_bbox(bbox, network_type="drive", simplify=False, retain_all=True),
        ox.graph_from_point(LOCATION_POINT, dist=500, network_type="drive", retain_all=True),
        ox.graph_from_address(ADDRESS, dist=500, network_type="drive", retain_all=True),
        ox.graph_from_place(PLACE, network_type="all", which_result=1, retain_all=True),
        ox.graph_from_polygon(
            polygon,
            network_type="walk",
            simplify=False,
            retain_all=True,
            truncate_by_edge=True,
        ),
    ]

    for G in graphs:
        ox.convert.validate_graph(G)
        assert set(G.nodes) == {101, 102, 103, 104, 105}
        assert G.graph["crs"] == ox.settings.default_crs
        assert "Fixture Street" in {
            data.get("name") for _, _, _, data in G.edges(keys=True, data=True)
        }

    G_drive = graphs[0]
    assert len(G_drive.edges) == 14
    assert dict(G_drive.nodes(data="street_count")) == {101: 2, 102: 4, 103: 2, 104: 4, 105: 4}


@pytest.mark.integration
def test_features_downloaders_use_committed_cache(http_cache: Path) -> None:
    ox.settings.cache_folder = http_cache

    bbox = ox.utils_geo.bbox_from_point(LOCATION_POINT, dist=500)
    polygon = ox.geocode_to_gdf(PLACE, which_result=1).geometry.iloc[0]
    gdfs = [
        ox.features_from_bbox(bbox, tags=TAGS),
        ox.features_from_point(LOCATION_POINT, tags=TAGS, dist=500),
        ox.features_from_address(ADDRESS, tags=TAGS, dist=500),
        ox.features_from_place(PLACE, tags=TAGS, which_result=1),
        ox.features_from_polygon(polygon, tags=TAGS),
    ]

    expected_index = [("node", 301), ("relation", 304), ("way", 302), ("way", 303)]
    for gdf in gdfs:
        ox.convert.validate_features_gdf(gdf)
        assert gdf.index.to_list() == expected_index
        assert gdf.crs == ox.settings.default_crs
        assert gdf.geometry.geom_type.to_list() == ["Point", "Polygon", "Polygon", "LineString"]
        assert set(gdf["name"]) == {
            "Fixture Cafe",
            "Fixture Retail",
            "Fixture Building",
            "Fixture Stop",
        }

    with pytest.raises(ValueError, match="The geometry of `polygon` is invalid"):
        ox.features.features_from_polygon(Polygon(((0, 0), (0, 0), (0, 0), (0, 0))), tags={})
    with suppress_type_checks(), pytest.raises(TypeError):
        ox.features.features_from_polygon(Point(0, 0), tags={})


@pytest.mark.integration
def test_elevation_from_cache_and_raster(
    http_cache: Path,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    ox.settings.cache_folder = http_cache
    G = drive_graph()

    with monkeypatch.context() as m:

        def _empty_elevation_response(_url: str, _pause: float) -> dict[str, list[object]]:
            return {"results": []}

        m.setattr(ox.elevation, "_elevation_request", _empty_elevation_response)
        with pytest.raises(ox._errors.InsufficientResponseError):
            ox.elevation.add_node_elevations_google(G.copy(), api_key="", batch_size=350)

    ox.settings.elevation_url_template = (
        "https://api.opentopodata.org/v1/aster30m?locations={locations}&key={key}"
    )
    G = ox.elevation.add_node_elevations_google(G, batch_size=100, pause=0)
    assert dict(G.nodes(data="elevation")) == {
        101: 10.0,
        102: 11.0,
        103: 12.0,
        104: 13.0,
        105: 14.0,
    }
    G_grades = ox.add_edge_grades(G.copy(), add_absolute=True)
    assert all("grade_abs" in data for _, _, data in G_grades.edges(data=True))

    ox.settings.cache_folder = tmp_path / "raster-cache"
    rasters = sorted(Path("tests/input_data").glob("elevation*.tif"))
    G = ox.elevation.add_node_elevations_raster(G, rasters[0], cpus=1)
    assert pd.notna(pd.Series(dict(G.nodes(data="elevation")))).any()
    G = ox.elevation.add_node_elevations_raster(G, rasters[0], cpus=2)
    assert pd.notna(pd.Series(dict(G.nodes(data="elevation")))).any()
    G = ox.elevation.add_node_elevations_raster(G, rasters, cpus=1)
    assert pd.notna(pd.Series(dict(G.nodes(data="elevation")))).sum() >= 2


@pytest.mark.unit
def test_http_parse_and_request_validation() -> None:
    response = cast("requests.Response", Response({"status": "ok"}))
    assert ox._http._parse_response(response) == {"status": "ok"}
    response = cast("requests.Response", Response([{"status": "ok"}]))
    assert ox._http._parse_response(response) == [{"status": "ok"}]

    params: OrderedDict[str, int | str] = OrderedDict()
    params["format"] = "json"
    params["address_details"] = 0
    with pytest.raises(ValueError, match="Nominatim `request_type` must be"):
        ox._nominatim._nominatim_request(params=params, request_type="xyz")
    with pytest.raises(TypeError, match="`query` must be a string"):
        ox.geocode_to_gdf(query={"City": "Boston"}, by_osmid=True)

    assert ox._overpass._get_network_filter("drive").startswith('["highway"]')
    with pytest.raises(ValueError, match="Unrecognized network_type"):
        ox._overpass._get_network_filter("not-real")
    ox.settings.overpass_memory = 123
    assert "[maxsize:123]" in ox._overpass._make_overpass_settings()


@pytest.mark.unit
def test_uncached_nominatim_request_paths(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    ox.settings.use_cache = False
    ox.settings.cache_folder = tmp_path

    def _sleep(_seconds: float) -> None:
        return None

    monkeypatch.setattr(stdlib_time, "sleep", _sleep)

    nominatim_responses = iter(
        [
            Response([], ok=False, status_code=429),
            Response([{"place_id": 1}]),
        ],
    )

    def _fake_nominatim_get(*_args: object, **_kwargs: object) -> Response:
        return next(nominatim_responses)

    ox.settings.nominatim_key = "fixture-key"
    params: OrderedDict[str, int | str] = OrderedDict()
    params["format"] = "json"
    params["q"] = "fixture"
    monkeypatch.setattr(requests, "get", _fake_nominatim_get)
    assert ox._nominatim._nominatim_request(params=params) == [{"place_id": 1}]
    assert params["key"] == "fixture-key"

    def _fake_nominatim_dict(*_args: object, **_kwargs: object) -> Response:
        return Response({"not": "a-list"})

    monkeypatch.setattr(requests, "get", _fake_nominatim_dict)
    with pytest.raises(ox._errors.InsufficientResponseError, match="did not return a list"):
        ox._nominatim._nominatim_request(params=OrderedDict(format="json", q="fixture"))


@pytest.mark.unit
def test_uncached_overpass_and_elevation_request_paths(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    ox.settings.use_cache = False
    ox.settings.cache_folder = tmp_path

    def _sleep(_seconds: float) -> None:
        return None

    monkeypatch.setattr(stdlib_time, "sleep", _sleep)

    config_dns_calls: list[str] = []

    def _fake_config_dns(url: str) -> None:
        config_dns_calls.append(url)

    overpass_responses = iter(
        [
            Response({"remark": "retry"}, ok=False, status_code=429),
            Response({"elements": []}),
        ],
    )

    def _fake_overpass_post(*_args: object, **_kwargs: object) -> Response:
        return next(overpass_responses)

    ox.settings.overpass_rate_limit = False
    monkeypatch.setattr(ox._http, "_config_dns", _fake_config_dns)
    monkeypatch.setattr(requests, "post", _fake_overpass_post)
    assert ox._overpass._overpass_request(OrderedDict(data="node(1);out;")) == {"elements": []}
    assert config_dns_calls == [ox.settings.overpass_url, ox.settings.overpass_url]

    def _fake_overpass_list(*_args: object, **_kwargs: object) -> Response:
        return Response([{"not": "a-dict"}])

    monkeypatch.setattr(requests, "post", _fake_overpass_list)
    with pytest.raises(ox._errors.InsufficientResponseError, match="did not return a dict"):
        ox._overpass._overpass_request(OrderedDict(data="node(2);out;"))

    elevation_responses = iter(
        [
            Response({"results": [{"elevation": 10.0}]}),
            Response([{"not": "a-dict"}]),
        ],
    )

    def _fake_elevation_get(*_args: object, **_kwargs: object) -> Response:
        return next(elevation_responses)

    monkeypatch.setattr(requests, "get", _fake_elevation_get)
    assert ox.elevation._elevation_request("https://example.com/elevation", pause=0) == {
        "results": [{"elevation": 10.0}],
    }
    with pytest.raises(ox._errors.InsufficientResponseError, match="did not return a dict"):
        ox.elevation._elevation_request("https://example.com/elevation?second", pause=0)


@pytest.mark.unit
def test_http_cache_and_dns_helper_branches(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    ox.settings.cache_folder = tmp_path

    assert ox._http._retrieve_from_cache("https://not-in-cache.example") is None
    assert ox._http._parse_response(Response({"status": "error"}, ok=False, status_code=500)) == {
        "status": "error",
    }

    real_config_dns = ox._http._config_dns
    original_getaddrinfo = socket.getaddrinfo

    def _fake_gethostbyname(_hostname: str) -> str:
        return "127.0.0.1"

    def _fake_original_getaddrinfo(*_args: object, **_kwargs: object) -> list[object]:
        return []

    monkeypatch.setattr(socket, "gethostbyname", _fake_gethostbyname)
    monkeypatch.setattr(ox._http, "_original_getaddrinfo", _fake_original_getaddrinfo)
    try:
        real_config_dns("https://example.com/api")
        assert socket.getaddrinfo("example.com", 443) == []
        assert socket.getaddrinfo("other.example", 443) == []
    finally:
        monkeypatch.setattr(socket, "getaddrinfo", original_getaddrinfo)


@pytest.mark.unit
def test_overpass_query_and_pause_branches(monkeypatch: pytest.MonkeyPatch) -> None:

    def _sleep(_seconds: float) -> None:
        return None

    class _StatusResponse:
        def __init__(self, text: str) -> None:
            self.text = text

    status_responses = iter(
        [
            _StatusResponse("bad"),
            _StatusResponse("\n\n\n\nMysterious server status\n"),
            _StatusResponse("\n\n\n\nSlot available after: 2000-01-01T00:00:00Z,\n"),
            _StatusResponse("\n\n\n\nCurrently running a query\n"),
            _StatusResponse("\n\n\n\n2 slots available now.\n"),
        ],
    )

    def _fake_status_get(*_args: object, **_kwargs: object) -> _StatusResponse:
        return next(status_responses)

    ox.settings.overpass_rate_limit = True
    monkeypatch.setattr(stdlib_time, "sleep", _sleep)
    monkeypatch.setattr(requests, "get", _fake_status_get)
    assert ox._overpass._get_overpass_pause("https://overpass.example/api", default_pause=7) == 7
    assert ox._overpass._get_overpass_pause("https://overpass.example/api", default_pause=8) == 8
    assert ox._overpass._get_overpass_pause("https://overpass.example/api") == 1
    assert ox._overpass._get_overpass_pause("https://overpass.example/api") == 0

    query = ox._overpass._create_overpass_features_query(
        "0 0 0 1 1 1 0 0",
        {"amenity": "cafe", "shop": ["books", "toys"]},
    )
    assert "amenity" in query
    assert "books" in query
    bad_tags = cast("dict[str, bool | str | list[str]]", {"amenity": [1]})
    with suppress_type_checks(), pytest.raises(TypeError, match="`tags` must be a dict"):
        ox._overpass._create_overpass_features_query("0 0", bad_tags)

    payloads: list[str] = []

    def _fake_overpass_request(data: OrderedDict[str, object]) -> dict[str, list[object]]:
        payloads.append(str(data["data"]))
        return {"elements": []}

    polygon = Polygon([(-1, -1), (1, -1), (1, 1), (-1, 1)])
    monkeypatch.setattr(ox._overpass, "_overpass_request", _fake_overpass_request)
    responses_str = list(ox._overpass._download_overpass_network(polygon, "all", '["highway"]'))
    responses_list = list(
        ox._overpass._download_overpass_network(
            polygon,
            "all",
            ['["railway"]', '["power"]'],
        ),
    )
    assert all(response == {"elements": []} for response in responses_str)
    assert all(response == {"elements": []} for response in responses_list)
    assert len(responses_list) == 2 * len(responses_str)
    assert any("railway" in payload for payload in payloads)
