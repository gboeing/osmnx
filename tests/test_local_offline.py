#!/usr/bin/env python
"""Offline CI tests for OSMnx."""

# ruff: noqa: D103, PLR2004, S101

from __future__ import annotations

import bz2
import gzip
import json
import logging as lg
import socket
import time as stdlib_time
from collections import OrderedDict
from pathlib import Path
from typing import cast

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
import pytest
import requests
from conftest import ADDRESS
from conftest import LOCATION_POINT
from conftest import PLACE
from conftest import TAGS
from conftest import _drive_graph
from conftest import _Response
from conftest import _toy_graph
from lxml import etree
from shapely import LineString
from shapely import Point
from shapely import Polygon
from shapely import wkt
from typeguard import suppress_type_checks

import osmnx as ox
from osmnx import _nominatim

# Cached API workflows


@pytest.mark.offline
def test_cache_fixture_manifest(http_cache: Path) -> None:
    manifest = json.loads((http_cache / "manifest.json").read_text(encoding="utf-8"))

    assert manifest["expected"]["network_nodes"] == 5
    assert manifest["expected"]["network_ways"] == 4
    assert manifest["expected"]["feature_elements"] == 8
    assert len(manifest["records"]) >= 10
    assert all((http_cache / record["cache_file"]).is_file() for record in manifest["records"])
    assert all(record["endpoint"].startswith("https://") for record in manifest["records"])
    assert all(record["query"] for record in manifest["records"])


@pytest.mark.offline
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


@pytest.mark.offline
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


@pytest.mark.offline
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


@pytest.mark.offline
def test_elevation_from_cache_and_raster(
    http_cache: Path,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    ox.settings.cache_folder = http_cache
    G = _drive_graph()

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


@pytest.mark.offline
def test_http_parse_and_request_validation() -> None:
    response = cast("requests.Response", _Response({"status": "ok"}))
    assert ox._http._parse_response(response) == {"status": "ok"}
    response = cast("requests.Response", _Response([{"status": "ok"}]))
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


@pytest.mark.offline
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
            _Response([], ok=False, status_code=429),
            _Response([{"place_id": 1}]),
        ],
    )

    def _fake_nominatim_get(*_args: object, **_kwargs: object) -> _Response:
        return next(nominatim_responses)

    ox.settings.nominatim_key = "fixture-key"
    params: OrderedDict[str, int | str] = OrderedDict()
    params["format"] = "json"
    params["q"] = "fixture"
    monkeypatch.setattr(requests, "get", _fake_nominatim_get)
    assert ox._nominatim._nominatim_request(params=params) == [{"place_id": 1}]
    assert params["key"] == "fixture-key"

    def _fake_nominatim_dict(*_args: object, **_kwargs: object) -> _Response:
        return _Response({"not": "a-list"})

    monkeypatch.setattr(requests, "get", _fake_nominatim_dict)
    with pytest.raises(ox._errors.InsufficientResponseError, match="did not return a list"):
        ox._nominatim._nominatim_request(params=OrderedDict(format="json", q="fixture"))


@pytest.mark.offline
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
            _Response({"remark": "retry"}, ok=False, status_code=429),
            _Response({"elements": []}),
        ],
    )

    def _fake_overpass_post(*_args: object, **_kwargs: object) -> _Response:
        return next(overpass_responses)

    ox.settings.overpass_rate_limit = False
    monkeypatch.setattr(ox._http, "_config_dns", _fake_config_dns)
    monkeypatch.setattr(requests, "post", _fake_overpass_post)
    assert ox._overpass._overpass_request(OrderedDict(data="node(1);out;")) == {"elements": []}
    assert config_dns_calls == [ox.settings.overpass_url, ox.settings.overpass_url]

    def _fake_overpass_list(*_args: object, **_kwargs: object) -> _Response:
        return _Response([{"not": "a-dict"}])

    monkeypatch.setattr(requests, "post", _fake_overpass_list)
    with pytest.raises(ox._errors.InsufficientResponseError, match="did not return a dict"):
        ox._overpass._overpass_request(OrderedDict(data="node(2);out;"))

    elevation_responses = iter(
        [
            _Response({"results": [{"elevation": 10.0}]}),
            _Response([{"not": "a-dict"}]),
        ],
    )

    def _fake_elevation_get(*_args: object, **_kwargs: object) -> _Response:
        return next(elevation_responses)

    monkeypatch.setattr(requests, "get", _fake_elevation_get)
    assert ox.elevation._elevation_request("https://example.com/elevation", pause=0) == {
        "results": [{"elevation": 10.0}],
    }
    with pytest.raises(ox._errors.InsufficientResponseError, match="did not return a dict"):
        ox.elevation._elevation_request("https://example.com/elevation?second", pause=0)


@pytest.mark.offline
def test_http_cache_and_dns_helpers_are_deterministic(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    ox.settings.cache_folder = tmp_path

    assert ox._http._retrieve_from_cache("https://not-in-cache.example") is None
    assert ox._http._parse_response(_Response({"status": "error"}, ok=False, status_code=500)) == {
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


@pytest.mark.offline
def test_overpass_status_and_query_helpers_parse_responses(monkeypatch: pytest.MonkeyPatch) -> None:

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


# Graph, GeoDataFrame, and file IO workflows


@pytest.mark.offline
def test_stats_simplification_and_conversion(http_cache: Path) -> None:
    ox.settings.cache_folder = http_cache
    G = _drive_graph()
    G_proj = ox.project_graph(G)
    G_proj = ox.project_graph(G_proj)

    stats = ox.basic_stats(G)
    assert stats["n"] == 5
    assert stats["m"] == 14
    assert stats["streets_per_node_counts"] == {0: 0, 1: 0, 2: 2, 3: 0, 4: 3}
    assert stats["edge_length_total"] == pytest.approx(1996.7794675)

    G_clean = ox.consolidate_intersections(G_proj, tolerance=10, rebuild_graph=True)
    G_reconnected = ox.consolidate_intersections(
        G_proj,
        tolerance=10,
        rebuild_graph=True,
        reconnect_edges=False,
    )
    centroids = ox.consolidate_intersections(G_proj, tolerance=10, rebuild_graph=False)

    assert len(G_clean) <= len(G_proj)
    assert len(G_reconnected) <= len(G_proj)
    assert len(centroids) <= len(G_proj)

    gdf_nodes, gdf_edges = ox.graph_to_gdfs(G)
    G2 = ox.graph_from_gdfs(gdf_nodes, gdf_edges, graph_attrs=G.graph)
    ox.convert.validate_graph(G2)
    assert set(G2.nodes) == set(G.nodes)
    assert set(G2.edges) == set(G.edges)

    D = ox.convert.to_digraph(G)
    Gu = ox.convert.to_undirected(G)
    assert isinstance(D, nx.DiGraph)
    assert isinstance(Gu, nx.MultiGraph)


@pytest.mark.offline
def test_save_load_graph_files(http_cache: Path) -> None:
    ox.settings.cache_folder = http_cache
    G = _drive_graph()
    ox.convert.validate_graph(G)

    ox.save_graph_geopackage(G, directed=False)
    fp = Path(ox.settings.data_folder) / "graph-dir.gpkg"
    ox.save_graph_geopackage(G, filepath=fp, directed=True)
    gdf_nodes1 = gpd.read_file(fp, layer="nodes").set_index("osmid")
    gdf_edges1 = gpd.read_file(fp, layer="edges").set_index(["u", "v", "key"])
    G2 = ox.graph_from_gdfs(gdf_nodes1, gdf_edges1, graph_attrs=G.graph)
    gdf_nodes2, gdf_edges2 = ox.graph_to_gdfs(G2)

    assert set(gdf_nodes1.index) == set(gdf_nodes2.index) == set(G.nodes)
    assert set(gdf_edges1.index) == set(gdf_edges2.index) == set(G.edges)

    with pytest.raises(ValueError, match="You must request nodes or edges or both"):
        ox.graph_to_gdfs(G2, nodes=False, edges=False)
    with pytest.raises(ValueError, match="Invalid literal for boolean"):
        ox.io._convert_bool_string("T")

    attr_name = "test_bool"
    G.graph[attr_name] = False
    nx.set_node_attributes(G, {n: bool(i % 2) for i, n in enumerate(G.nodes)}, attr_name)
    nx.set_edge_attributes(G, {e: bool(i % 2) for i, e in enumerate(G.edges)}, attr_name)

    graphml_fp = Path(ox.settings.data_folder) / "graph.graphml"
    ox.save_graphml(G, gephi=True)
    ox.save_graphml(G, gephi=False, filepath=graphml_fp)
    G2 = ox.load_graphml(
        graphml_fp,
        graph_dtypes={attr_name: ox.io._convert_bool_string},
        node_dtypes={attr_name: ox.io._convert_bool_string},
        edge_dtypes={attr_name: ox.io._convert_bool_string},
    )
    ox.convert.validate_graph(G2)
    assert set(G2.nodes) == set(G.nodes)
    assert set(G2.edges) == set(G.edges)

    graphml = Path("tests/input_data/short.graphml").read_text(encoding="utf-8")
    G3 = ox.load_graphml(graphml_str=graphml, node_dtypes={"osmid": str})
    assert len(G3) > 0


@pytest.mark.offline
def test_osm_xml_read_write(http_cache: Path, tmp_path: Path) -> None:
    ox.settings.cache_folder = http_cache
    node_id = 53098262
    neighbor_ids = 53092170, 53060438, 53027353, 667744075

    path_bz2 = Path("tests/input_data/West-Oakland.osm.bz2")
    file_contents = bz2.decompress(path_bz2.read_bytes())
    path_osm_temp = tmp_path / "West-Oakland.osm"
    path_gz_temp = tmp_path / "West-Oakland.osm.gz"
    path_osm_temp.write_bytes(file_contents)
    path_gz_temp.write_bytes(gzip.compress(file_contents))

    for filepath in (path_bz2, path_gz_temp, path_osm_temp):
        G = ox.graph_from_xml(filepath)
        ox.convert.validate_graph(G, strict=False)
        assert node_id in G.nodes
        for neighbor_id in neighbor_ids:
            edge_key = (node_id, neighbor_id, 0)
            assert neighbor_id in G.nodes
            assert edge_key in G.edges
            assert G.edges[edge_key]["name"] in {"8th Street", "Willow Street"}

    default_all_oneway = ox.settings.all_oneway
    ox.settings.all_oneway = True
    G = _drive_graph()
    fp = Path(ox.settings.data_folder) / "graph.osm"
    ox.io.save_graph_xml(G, filepath=fp, way_tag_aggs={"lanes": "sum"})

    parser = etree.XMLParser(schema=etree.XMLSchema(file="./tests/input_data/osm_schema.xsd"))
    _ = etree.parse(fp, parser=parser)

    with pytest.raises(ox._errors.GraphSimplificationError):
        ox.io.save_graph_xml(ox.simplification.simplify_graph(G))
    ox.settings.all_oneway = default_all_oneway


@pytest.mark.offline
def test_graph_creation_validates_inputs(http_cache: Path) -> None:
    ox.settings.cache_folder = http_cache
    G_network = ox.graph_from_point(
        LOCATION_POINT,
        dist=500,
        dist_type="network",
        network_type="drive",
        retain_all=True,
    )
    assert len(G_network) > 0

    bbox = ox.utils_geo.bbox_from_point(LOCATION_POINT, dist=500)
    G_largest = ox.graph_from_bbox(bbox, network_type="drive", simplify=False, retain_all=False)
    ox.convert.validate_graph(G_largest)

    with pytest.raises(ValueError, match="`dist_type` must be"):
        ox.graph_from_point(LOCATION_POINT, dist=500, dist_type="bad")
    with suppress_type_checks(), pytest.raises(TypeError, match="Geometry must be"):
        ox.graph_from_polygon(Point(0, 0))
    with pytest.raises(ValueError, match="geometry of `polygon` is invalid"):
        ox.graph_from_polygon(Polygon(((0, 0), (0, 0), (0, 0), (0, 0))))

    ox.settings.cache_only_mode = True
    with pytest.raises(ox._errors.CacheOnlyInterruptError, match="Interrupted because"):
        ox.graph._create_graph([{"elements": []}], bidirectional=False)
    ox.settings.cache_only_mode = False
    with pytest.raises(ox._errors.InsufficientResponseError, match="No data elements"):
        ox.graph._create_graph([{"elements": []}], bidirectional=False)

    assert ox.graph._is_path_one_way(
        {"junction": "roundabout"},
        bidirectional=False,
        oneway_values=set(),
    )

    G = nx.MultiDiGraph()
    G.add_nodes_from([1, 2])
    path = {"osmid": 1, "nodes": [1, 2], "oneway": "-1", "highway": "residential"}
    ox.graph._add_paths(G, [path], bidirectional=False)
    assert (2, 1, 0) in G.edges


@pytest.mark.offline
def test_conversion_and_graphml_validate_inputs(tmp_path: Path) -> None:
    empty = nx.MultiDiGraph(crs="epsg:4326")
    with pytest.raises(ValueError, match="contains no nodes"):
        ox.graph_to_gdfs(empty, nodes=True, edges=False)
    with pytest.raises(ValueError, match="contains no edges"):
        ox.graph_to_gdfs(empty, nodes=False, edges=True)

    gdf_nodes = gpd.GeoDataFrame({"x": [0.0, 1.0], "y": [0.0, 1.0]}, index=[1, 2])
    gdf_edges = gpd.GeoDataFrame(
        {"osmid": [1], "length": [1.0], "geometry": [LineString([(0, 0), (1, 1)])]},
        index=pd.MultiIndex.from_tuples([(1, 2, 0)], names=["u", "v", "key"]),
        crs="epsg:4326",
    )
    G = ox.graph_from_gdfs(gdf_nodes, gdf_edges)
    assert G.graph["crs"] == "epsg:4326"

    with pytest.raises(ValueError, match="one and only one"):
        ox.load_graphml()
    with pytest.raises(ValueError, match="one and only one"):
        ox.load_graphml(tmp_path / "missing.graphml", graphml_str="<graphml />")
    assert ox.io._convert_bool_string(value=True)

    G_types = nx.MultiDiGraph(consolidated="False", simplified="True")
    G_types.add_node(1, x="0", y="0", osmid="1", street_count="1", tags="[1, 2]")
    G_types.add_node(2, x="1", y="1", osmid="2", street_count="1")
    G_types.add_edge(
        1,
        2,
        osmid="[10, 11]",
        length="1.5",
        oneway="False",
        reversed="True",
        geometry="LINESTRING (0 0, 1 1)",
    )
    ox.io._convert_graph_attr_types(G_types, {"consolidated": ox.io._convert_bool_string})
    ox.io._convert_node_attr_types(G_types, {"osmid": int, "street_count": int, "x": float})
    ox.io._convert_edge_attr_types(
        G_types,
        {"osmid": int, "length": float, "oneway": ox.io._convert_bool_string},
    )
    assert G_types.edges[1, 2, 0]["osmid"] == [10, 11]

    G_parallel = _toy_graph()
    G_parallel.add_edge(
        2,
        1,
        osmid=10,
        length=1.0,
        highway="residential",
        geometry=LineString([(1, 0), (0, 0)]),
    )
    G_parallel.add_edge(
        1,
        2,
        key=1,
        osmid=12,
        length=1.0,
        highway="service",
        geometry=LineString([(0, 0), (0.5, 0.2), (1, 0)]),
    )
    G_parallel.add_edge(
        2,
        1,
        key=1,
        osmid=12,
        length=1.0,
        highway="service",
        geometry=LineString([(1, 0), (0.5, -0.2), (0, 0)]),
    )
    Gu = ox.convert.to_undirected(G_parallel)
    assert isinstance(Gu, nx.MultiGraph)

    G_duplicate = nx.MultiDiGraph(crs="epsg:4326")
    G_duplicate.add_node(1, x=0, y=0)
    G_duplicate.add_node(2, x=1, y=1)
    G_duplicate.add_edge(1, 2, key=0, osmid=1)
    G_duplicate.add_edge(2, 1, key=1, osmid=1)
    Gu_duplicate = ox.convert.to_undirected(G_duplicate)
    assert len(Gu_duplicate.edges) == 1

    same_geom = LineString([(0, 0), (1, 1)])
    assert ox.convert._is_duplicate_edge(
        {"osmid": [1, 2], "geometry": same_geom},
        {"osmid": [2, 1], "geometry": LineString([(1, 1), (0, 0)])},
    )
    assert ox.convert._is_duplicate_edge({"osmid": 1}, {"osmid": 1})
    assert not ox.convert._is_duplicate_edge(
        {"osmid": 1, "geometry": LineString([(0, 0), (1, 1)])},
        {"osmid": 1},
    )


@pytest.mark.offline
def test_osm_xml_roundtrip_warns_for_projected_graphs(tmp_path: Path) -> None:
    G = nx.MultiDiGraph(crs="epsg:3857")
    G.add_node(1, x=0.0, y=0.0, street_count=1)
    G.add_node(2, x=1.0, y=0.0, street_count=2)
    G.add_node(3, x=2.0, y=0.0, street_count=1)
    G.add_node(1, x=0.0, y=0.0, street_count=1, uid=None)
    G.add_edge(1, 2, osmid=10, length=1.0, highway="residential", uid=None)
    G.add_edge(2, 3, osmid=11, length=1.0, highway="primary")
    fp = tmp_path / "projected.osm"
    with pytest.warns(UserWarning, match="all_oneway=True|unprojected"):
        ox.io.save_graph_xml(G, filepath=fp)
    with pytest.warns(UserWarning, match="generated by OSMnx"):
        G_loaded = ox.graph_from_xml(fp, simplify=False, retain_all=True)
    assert len(G_loaded) > 0

    G_source_cycle = nx.MultiDiGraph([(1, 2), (2, 3), (3, 1), (1, 4)])
    assert set(ox._osm_xml._sort_nodes(G_source_cycle, osmid=1)) >= {1, 2, 3, 4}

    G_target_cycle = nx.MultiDiGraph([(1, 2), (2, 3), (3, 1), (4, 1)])
    assert set(ox._osm_xml._sort_nodes(G_target_cycle, osmid=2)) >= {1, 2, 3, 4}

    G_multi_cycle = nx.MultiDiGraph([(1, 2), (2, 3), (3, 1), (3, 4), (4, 2)])
    assert len(ox._osm_xml._sort_nodes(G_multi_cycle, osmid=3)) > 0
    G_simple_cycle = nx.MultiDiGraph([(1, 2), (2, 3), (3, 1)])
    assert len(ox._osm_xml._sort_nodes(G_simple_cycle, osmid=4)) > 0


# Analysis, routing, distance, and plotting workflows


@pytest.mark.offline
def test_logging_and_utils(tmp_path: Path) -> None:
    ox.settings.log_console = True
    ox.settings.log_file = True
    ox.settings.logs_folder = tmp_path

    ox.utils.log("test a fake default message")
    ox.utils.log("test a fake debug", level=lg.DEBUG)
    ox.utils.log("test a fake info", level=lg.INFO)
    ox.utils.log("test a fake warning", level=lg.WARNING)
    ox.utils.log("test a fake error", level=lg.ERROR)

    ox.utils.citation(style="apa")
    ox.utils.citation(style="bibtex")
    ox.utils.citation(style="ieee")

    assert "T" in ox.utils.ts(style="iso8601")
    assert len(ox.utils.ts(style="date")) == 10
    assert len(ox.utils.ts(style="time")) == 8
    assert ox.settings.logs_folder == tmp_path


@pytest.mark.offline
def test_bearings_routing_and_nearest(http_cache: Path) -> None:
    ox.settings.cache_folder = http_cache
    G = _drive_graph()

    assert ox.bearing.calculate_bearing(0, 0, 1, 1) == pytest.approx(44.99563646)

    G = ox.add_edge_bearings(G)
    with pytest.warns(UserWarning, match="directional"):
        bearings, weights = ox.bearing._extract_edge_bearings(G, min_length=0, weight=None)
    assert len(bearings) == len(weights) == len(G.edges)
    Gu = ox.convert.to_undirected(G)
    entropy = ox.bearing.orientation_entropy(Gu, weight="length")
    assert entropy == pytest.approx(1.777310166)

    _, ax = ox.plot.plot_orientation(Gu, area=True, title="Title")
    _, _ = ox.plot.plot_orientation(Gu, ax=ax, area=False, title="Title")

    G = ox.add_edge_speeds(G)
    G = ox.add_edge_speeds(G, hwy_speeds={"motorway": 100})
    G = ox.add_edge_travel_times(G)
    route = ox.shortest_path(G, 101, 103, weight="travel_time")
    assert route == [101, 102, 103]
    assert ox.routing.route_to_gdf(G, route, weight="travel_time").index.to_list() == [
        (101, 102, 0),
        (102, 103, 0),
    ]

    assert ox.routing._clean_maxspeed("100,2") == 100.2
    assert ox.routing._clean_maxspeed("100 mph") == pytest.approx(160.934)
    assert ox.routing._clean_maxspeed("signal") is None
    assert ox.routing._collapse_multiple_maxspeed_values(["25 mph", "30 mph"], np.mean) == 44.25685

    paths1 = ox.shortest_path(G, [101, 102], [103, 104], weight="length", cpus=1)
    paths2 = ox.shortest_path(G, [101, 102], [103, 104], weight="length", cpus=2)
    assert paths1 == paths2 == [[101, 102, 103], [102, 105, 104]]

    routes = list(ox.routing.k_shortest_paths(G, 101, 103, k=2, weight="travel_time"))
    assert routes[0] == [101, 102, 103]

    assert ox.distance.great_circle(0, 0, 1, 1) == pytest.approx(157249.6034105)
    assert ox.distance.euclidean(0, 0, 1, 1) == pytest.approx(1.4142135)
    assert ox.distance.nearest_nodes(G, -122.4100, 37.79125) == 105

    Gp = ox.project_graph(G)
    points = ox.utils_geo.sample_points(ox.convert.to_undirected(Gp), 5)
    X = points.x.to_numpy()
    Y = points.y.to_numpy()
    assert len(ox.distance.nearest_nodes(Gp, X, Y, return_dist=True)[0]) == 5
    nearest_edges = ox.distance.nearest_edges(Gp, X, Y, return_dist=False)
    assert len(nearest_edges) == 5
    assert len(nearest_edges[0]) == 3


@pytest.mark.offline
def test_plots_and_colors(http_cache: Path) -> None:
    ox.settings.cache_folder = http_cache
    G = _drive_graph()
    Gp = ox.project_graph(G)

    colors = ox.plot.get_colors(n=5, cmap="plasma", start=0.1, stop=0.9, alpha=0.5)
    node_colors = ox.plot.get_node_colors_by_attr(G, "x")
    edge_colors = ox.plot.get_edge_colors_by_attr(G, "length", num_bins=5)

    assert len(colors) == 5
    assert len(node_colors) == len(G)
    assert len(edge_colors) == len(G.edges)

    filepath = Path(ox.settings.data_folder) / "test.svg"
    _, ax = ox.plot_graph(G, show=False, save=True, close=True, filepath=filepath)
    _, ax = ox.plot_graph(Gp, edge_linewidth=0, figsize=(5, 5), bgcolor="y")
    _, _ = ox.plot_graph(
        Gp,
        ax=ax,
        dpi=180,
        node_color="k",
        node_size=5,
        node_alpha=0.1,
        node_edgecolor="b",
        node_zorder=5,
        edge_color="r",
        edge_linewidth=2,
        edge_alpha=0.1,
        show=False,
        save=True,
        close=True,
    )
    _, _ = ox.plot_figure_ground(G=G)
    _, _ = ox.plot_graph_route(G, [101, 102, 103], save=True)
    _, _ = ox.plot_graph_routes(G, [[101, 102, 103], [101, 102, 105, 104]])

    assert filepath.is_file()


@pytest.mark.offline
def test_utils_geo_projection_and_http_helpers(tmp_path: Path) -> None:
    geom = ox.utils_geo.buffer_geometry(Point(-122.41, 37.79), dist=100)
    bbox = ox.utils_geo.bbox_from_point(LOCATION_POINT, dist=500, project_utm=True)
    bbox_with_crs = ox.utils_geo.bbox_from_point(
        LOCATION_POINT,
        dist=500,
        project_utm=True,
        return_crs=True,
    )
    poly = ox.utils_geo.bbox_to_poly(ox.utils_geo.bbox_from_point(LOCATION_POINT, dist=500))
    subdivided = ox.utils_geo._consolidate_subdivide_geometry(poly)

    assert geom.area > 0
    assert len(bbox) == 4
    assert len(bbox_with_crs) == 2
    assert len(subdivided.geoms) == 1
    assert list(ox.utils_geo.interpolate_points(LineString([(0, 0), (1, 0)]), 0.25)) == [
        (0.0, 0.0),
        (0.25, 0.0),
        (0.5, 0.0),
        (0.75, 0.0),
        (1.0, 0.0),
    ]

    assert ox.projection.is_projected("epsg:3857")
    assert not ox.projection.is_projected("epsg:4326")
    geom_proj, crs_proj = ox.projection.project_geometry(poly)
    assert ox.projection.is_projected(crs_proj)
    geom_latlon, crs_latlon = ox.projection.project_geometry(
        geom_proj,
        crs=crs_proj,
        to_latlong=True,
    )
    assert crs_latlon == ox.settings.default_crs
    assert geom_latlon.area == pytest.approx(poly.area, rel=0.01)

    ox.settings.cache_folder = tmp_path
    ox._http._save_to_cache("https://example.com/test", {"ok": True}, ok=True)
    assert ox._http._retrieve_from_cache("https://example.com/test") == {"ok": True}
    assert ox._http._hostname_from_url("https://example.com:443/api") == "example.com"
    assert "User-Agent" in ox._http._get_http_headers(user_agent="test-agent")


@pytest.mark.offline
def test_routing_and_http_helpers_validate_inputs_and_statuses(
    http_cache: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    ox.settings.cache_folder = http_cache
    G = _drive_graph()

    G_dist = ox.truncate.truncate_graph_dist(G, 101, 150)
    assert set(G_dist.nodes).issubset(G.nodes)
    assert len(G_dist) < len(G)

    bbox = (-122.412, 37.790, -122.410, 37.793)
    G_bbox = ox.truncate.truncate_graph_bbox(G, bbox, truncate_by_edge=True)
    assert 101 in G_bbox.nodes
    with pytest.raises(ValueError, match="Found no graph nodes"):
        ox.truncate.truncate_graph_bbox(G, (0, 0, 1, 1))

    G_disconnected = G.copy()
    G_disconnected.add_node(999, x=0, y=0, street_count=0)
    G_largest = ox.truncate.largest_component(G_disconnected)
    assert 999 not in G_largest.nodes

    with pytest.raises(TypeError, match="must either both be iterable"):
        ox.shortest_path(G, 101, [103])  # type: ignore[call-overload]
    with pytest.raises(ValueError, match="must be of equal length"):
        ox.shortest_path(G, [101], [102, 103])

    G_no_speed = G.copy()
    for _, _, _, data in G_no_speed.edges(keys=True, data=True):
        data.pop("maxspeed", None)
    with pytest.raises(ValueError, match="must pass `hwy_speeds` or `fallback`"):
        ox.routing.add_edge_speeds(G_no_speed)

    ox.settings.doh_url_template = None
    assert ox._http._resolve_host_via_doh("example.com") == "example.com"
    ox.settings.doh_url_template = "https://dns.example.test/{hostname}"

    doh_response = _Response({"Status": 0, "Answer": [{"data": "192.0.2.1"}]})

    def _fake_doh_get(*_args: object, **_kwargs: object) -> _Response:
        return doh_response

    monkeypatch.setattr(requests, "get", _fake_doh_get)
    assert ox._http._resolve_host_via_doh("example.com") == "192.0.2.1"

    doh_response = _Response({"Status": 1})
    assert ox._http._resolve_host_via_doh("example.com") == "example.com"

    class _StatusResponse:
        text = "\n\n\n\n2 slots available now.\n"

    def _fake_status_get(*_args: object, **_kwargs: object) -> _StatusResponse:
        return _StatusResponse()

    ox.settings.overpass_rate_limit = True
    monkeypatch.setattr(requests, "get", _fake_status_get)
    assert ox._overpass._get_overpass_pause("https://overpass.example.test/api") == 0


@pytest.mark.offline
def test_projection_stats_distance_and_geometry_behaviors() -> None:
    gdf_north = gpd.GeoDataFrame(geometry=[Point(0, 85.1)], crs="epsg:4326")
    gdf_south = gpd.GeoDataFrame(geometry=[Point(0, -84.1)], crs="epsg:4326")
    assert ox.projection.project_gdf(gdf_north).crs == "epsg:32661"
    assert ox.projection.project_gdf(gdf_south).crs == "epsg:32761"
    with pytest.raises(ValueError, match="valid CRS"):
        ox.projection.project_gdf(gpd.GeoDataFrame(geometry=[]))

    G_simple = nx.MultiDiGraph(crs="epsg:4326")
    G_simple.add_node(1, x=0.0, y=0.0, street_count=1)
    G_simple.add_node(2, x=0.01, y=0.0, street_count=2)
    G_simple.add_node(3, x=0.02, y=0.0, street_count=1)
    G_simple.add_edge(1, 2, osmid=1, length=1.0)
    G_simple.add_edge(2, 3, osmid=2, length=1.0)
    G_simple = ox.simplification.simplify_graph(G_simple)
    G_projected = ox.project_graph(G_simple)
    assert ox.project_graph(G_projected, to_latlong=True).graph["crs"] == ox.settings.default_crs

    G_stats = _toy_graph(crs="epsg:3857")
    G_missing_count = G_stats.copy()
    del G_missing_count.nodes[1]["street_count"]
    assert set(ox.stats.streets_per_node(G_missing_count)) != set(G_missing_count.nodes)
    Gu = ox.convert.to_undirected(G_stats)
    assert ox.stats.circuity_avg(Gu) == pytest.approx(1.0)
    assert ox.stats.count_streets_per_node(G_stats) == {1: 1, 2: 2, 3: 1}
    stats = ox.basic_stats(G_stats, area=1000, clean_int_tol=1)
    assert "clean_intersection_density_km" in stats

    G_bad = nx.MultiDiGraph(crs="epsg:4326")
    G_bad.add_nodes_from([(1, {"x": np.nan, "y": 0}), (2, {"x": 1, "y": 1})])
    G_bad.add_edge(1, 2)
    with pytest.raises(ValueError, match="possibly due to input data clipping"):
        ox.distance.add_edge_lengths(G_bad)
    with pytest.raises(ValueError, match="cannot contain nulls"):
        ox.distance.nearest_nodes(G_stats, [np.nan], [0])
    with pytest.raises(ValueError, match="cannot contain nulls"):
        ox.distance.nearest_edges(G_stats, [0], [np.nan])
    assert ox.distance.nearest_nodes(G_stats, 0, 0, return_dist=True)[0] == 1
    assert ox.distance.nearest_nodes(G_stats, [0], [0], return_dist=False)[0] == 1
    assert ox.distance.nearest_edges(G_stats, 0, 0, return_dist=True)[0] == (1, 2, 0)
    assert ox.distance.nearest_edges(G_stats, 0, 0, return_dist=False) == (1, 2, 0)
    assert len(ox.distance.nearest_edges(G_stats, [0], [0], return_dist=True)[0]) == 1

    with pytest.warns(UserWarning, match="undirected"):
        assert len(ox.utils_geo.sample_points(G_stats, 1)) == 1
    with suppress_type_checks(), pytest.raises(TypeError, match="LineString"):
        list(ox.utils_geo.interpolate_points(Point(0, 0), 1))
    with suppress_type_checks(), pytest.raises(TypeError, match="Geometry must be"):
        ox.utils_geo._consolidate_subdivide_geometry(Point(0, 0))

    assert ox.simplification._build_path(nx.MultiDiGraph([(1, 2), (2, 1)]), 1, 2, {1}) == [
        1,
        2,
    ]
    empty_graph = nx.MultiDiGraph(crs="epsg:3857")
    assert len(ox.consolidate_intersections(empty_graph, rebuild_graph=True)) == 0
    assert len(ox.consolidate_intersections(empty_graph, rebuild_graph=False)) == 0
    assert (
        len(
            ox.consolidate_intersections(
                G_stats,
                tolerance={1: 0.5, 2: 0.5},
                rebuild_graph=False,
                dead_ends=True,
            ),
        )
        > 0
    )


@pytest.mark.offline
def test_truncation_plotting_and_routing_errors(tmp_path: Path) -> None:
    G_trunc = nx.MultiDiGraph(crs="epsg:4326")
    G_trunc.add_node(1, x=0.0, y=0.0)
    G_trunc.add_node(2, x=2.0, y=0.0)
    G_trunc.add_node(3, x=4.0, y=0.0)
    G_trunc.add_edges_from([(1, 2), (2, 3)])
    polygon = Polygon([(-1, -1), (1, -1), (1, 1), (-1, 1)])
    G_truncated = ox.truncate.truncate_graph_polygon(G_trunc, polygon, truncate_by_edge=True)
    assert set(G_truncated.nodes) == {1, 2}

    G_strong = nx.MultiDiGraph(crs="epsg:4326")
    G_strong.add_edges_from([(1, 2), (2, 1), (3, 4)])
    assert set(ox.truncate.largest_component(G_strong, strongly=True).nodes) == {1, 2}

    G_plot = _toy_graph()
    G_plot.add_node(99, x=0.0, y=1.0, street_count=0)
    _, _ = ox.plot_graph_route(G_plot, [1, 2, 3], show=False, close=True)
    _, _ = ox.plot_figure_ground(G_plot, show=False, close=True)
    footprints = gpd.GeoDataFrame(
        geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
        crs=ox.settings.default_crs,
    )
    _, _ = ox.plot_footprints(
        footprints,
        show=False,
        close=True,
        save=True,
        filepath=tmp_path / "footprints.png",
    )
    color_values = pd.Series([1.0, 2.0, 3.0, 4.0])
    assert (
        len(
            ox.plot._get_colors_by_value(
                color_values,
                num_bins=2,
                cmap="viridis",
                start=0,
                stop=1,
                na_color="none",
                equal_size=True,
            ),
        )
        == 4
    )
    with pytest.raises(ValueError, match="no attribute values"):
        ox.plot._get_colors_by_value(
            pd.Series(dtype=float),
            num_bins=None,
            cmap="viridis",
            start=0,
            stop=1,
            na_color="none",
            equal_size=False,
        )

    G_route = _toy_graph()
    G_route.add_node(4, x=10.0, y=10.0, street_count=0)
    assert ox.shortest_path(G_route, 1, 4) is None
    assert ox.shortest_path(G_route, [1], [3], cpus=None) == [[1, 2, 3]]
    assert ox.routing._collapse_multiple_maxspeed_values(["10"], lambda _values: "bad") is None
    assert ox.routing._collapse_multiple_maxspeed_values(["mph", "kph"], np.mean) is None
    assert ox.routing._collapse_multiple_maxspeed_values(["signal"], np.mean) is None

    with pytest.raises(ValueError, match="Invalid citation style"):
        ox.utils.citation(style="bad")
    with pytest.raises(ValueError, match="Invalid timestamp style"):
        ox.utils.ts(style="bad")


# Geometry, validation, and simplification workflows


@pytest.mark.offline
def test_validation_errors() -> None:
    G = nx.MultiDiGraph()
    G.add_edge(0, 1)
    with pytest.raises(ox._errors.ValidationError):
        ox._validate._verify_numeric_edge_attribute(G, "length", strict=True)

    with pytest.raises(ox._errors.ValidationError):
        ox.convert.validate_features_gdf(gpd.GeoDataFrame(index=[0, 0]))

    gdf_nodes = pd.DataFrame(index=[0, 0])
    gdf_edges = pd.DataFrame()
    with suppress_type_checks(), pytest.raises(ox._errors.ValidationError):
        ox.convert.validate_node_edge_gdfs(gdf_nodes, gdf_edges)

    gdf_nodes = gpd.GeoDataFrame(geometry=[Polygon(), Polygon()])
    gdf_edges = gpd.GeoDataFrame()
    with pytest.raises(ox._errors.ValidationError):
        ox.convert.validate_node_edge_gdfs(gdf_nodes, gdf_edges)

    data = {"x": [0, 1], "y": [2, 3]}
    gdf_nodes = gpd.GeoDataFrame(data=data, geometry=[Point((6, 7)), Point((8, 9))])
    with pytest.raises(ox._errors.ValidationError):
        ox.convert.validate_node_edge_gdfs(gdf_nodes, gdf_edges)

    G = nx.Graph()
    with suppress_type_checks(), pytest.raises(ox._errors.ValidationError):
        ox.convert.validate_graph(G)

    G = nx.MultiDiGraph()
    del G.graph
    G.add_edge("0", "1")
    with pytest.raises(ox._errors.ValidationError):
        ox.convert.validate_graph(G)

    G = nx.MultiDiGraph()
    G.graph["crs"] = "epsg:4326"
    G.add_edge(0, 1)
    with pytest.raises(ox._errors.ValidationError):
        ox.convert.validate_graph(G)

    nx.set_node_attributes(G, values=0, name="x")
    nx.set_node_attributes(G, values=0, name="y")
    nx.set_node_attributes(G, values=0, name="street_count")
    nx.set_edge_attributes(G, values=[0], name="osmid")
    nx.set_edge_attributes(G, values=1.5, name="length")
    ox.convert.validate_graph(G)


@pytest.mark.offline
def test_validation_reports_bad_graph_attributes() -> None:
    G = nx.MultiDiGraph(crs="epsg:4326")
    G.add_node("a", x="bad", y=0)
    G.add_node("b", x=1, y="bad")
    G.add_edge("a", "b", osmid="bad", length="bad")

    with pytest.raises(ox._errors.ValidationError, match="contains non-numeric"):
        ox._validate._verify_numeric_edge_attribute(G, "length", strict=True)

    with pytest.raises(ox._errors.ValidationError, match="Node 'x' and 'y'"):
        ox.convert.validate_graph(G)

    with pytest.warns(UserWarning, match="Node 'x' and 'y'"):
        ox.convert.validate_graph(G, strict=False)

    G_nan = nx.MultiDiGraph(crs="epsg:4326")
    G_nan.add_node(1, x=0, y=0, street_count=1)
    G_nan.add_node(2, x=1, y=1, street_count=1)
    G_nan.add_edge(1, 2, osmid=1, length=np.nan)
    with pytest.raises(ox._errors.ValidationError, match="missing or null"):
        ox._validate._verify_numeric_edge_attribute(G_nan, "length", strict=True)
    with pytest.warns(UserWarning, match="missing or null"):
        ox._validate._verify_numeric_edge_attribute(G_nan, "length", strict=False)

    G_bad_crs = nx.MultiDiGraph(crs="epsg:999999")
    G_bad_crs.add_node(1, x=0, y=0, street_count=1)
    G_bad_crs.add_node(2, x=1, y=1, street_count=1)
    G_bad_crs.add_edge(1, 2, osmid=1, length=1.0)
    with pytest.raises(ox._errors.ValidationError, match="valid CRS"):
        ox.convert.validate_graph(G_bad_crs)

    gdf_nodes = gpd.GeoDataFrame(
        {"x": [0.0, 1.0], "y": [0.0, 1.0]},
        geometry=[Point(0, 0), Point(1, 1)],
        index=[1, 2],
        crs="epsg:4326",
    )
    gdf_edges = gpd.GeoDataFrame(
        {"osmid": [1], "length": [1.0], "geometry": [LineString([(0, 0), (1, 1)])]},
        index=pd.MultiIndex.from_tuples([(1, 2, 0)], names=["u", "v", "key"]),
        crs="epsg:4326",
    )
    ox.convert.validate_node_edge_gdfs(gdf_nodes, gdf_edges)


@pytest.mark.offline
def test_features_from_xml_and_polygon_holes() -> None:
    gdf = ox.features_from_xml("tests/input_data/planet_10.068,48.135_10.071,48.137.osm")
    assert len(gdf) > 0

    gdf = ox.features_from_xml("tests/input_data/West-Oakland.osm.bz2")
    assert "Willow Street" in gdf["name"].to_numpy()

    outer1 = Polygon(((0, 0), (4, 0), (4, 4), (0, 4)))
    inner1 = Polygon(((1, 1), (2, 1), (2, 3), (1, 3)))
    inner2 = Polygon(((2, 1), (3, 1), (3, 3), (2, 3)))
    outer2 = Polygon(((1.5, 1.5), (2.5, 1.5), (2.5, 2.5), (1.5, 2.5)))
    result = ox.features._remove_polygon_holes([outer1, outer2], [inner1, inner2])
    geom_wkt = (
        "MULTIPOLYGON (((4 4, 4 0, 0 0, 0 4, 4 4), "
        "(3 1, 3 3, 2 3, 1 3, 1 1, 2 1, 3 1)), "
        "((2.5 2.5, 2.5 1.5, 1.5 1.5, 1.5 2.5, 2.5 2.5)))"
    )
    assert result.equals(wkt.loads(geom_wkt))


@pytest.mark.offline
def test_simplification_preserves_merged_edge_attrs(monkeypatch: pytest.MonkeyPatch) -> None:
    G = nx.MultiDiGraph(crs="epsg:4326")
    G.add_node(1, x=0.0, y=0.0, street_count=1)
    G.add_node(2, x=1.0, y=0.0, street_count=2)
    G.add_node(3, x=2.0, y=0.0, street_count=1)
    G.add_edge(1, 2, osmid=1, length=1.0, travel_time=1.0, name="a")
    G.add_edge(2, 3, osmid=2, length=1.0, travel_time=1.0, name="b")
    G.add_edge(2, 3, osmid=3, length=2.0, travel_time=2.0, name="c")

    def _fake_paths(
        _G: nx.MultiDiGraph,
        _node_attrs_include: object,
        _edge_attrs_differ: object,
    ) -> list[list[int]]:
        return [[1, 2, 3]]

    monkeypatch.setattr(ox.simplification, "_get_paths_to_simplify", _fake_paths)
    Gs = ox.simplification.simplify_graph(G, track_merged=True)
    edge_data = next(iter(Gs.edges(data=True)))[2]
    assert edge_data["length"] == 2.0
    assert edge_data["travel_time"] == 2.0
    assert edge_data["merged_edges"] == [(1, 2), (2, 3)]

    Gs.graph["simplified"] = True
    with pytest.raises(ox._errors.GraphSimplificationError, match="already been simplified"):
        ox.simplification.simplify_graph(Gs)

    H = nx.MultiDiGraph()
    H.add_node(1)
    H.add_node(2, highway="traffic_signals")
    H.add_node(3)
    H.add_edges_from([(1, 2, {"osmid": 1}), (2, 1, {"osmid": 1})])
    H.add_edges_from([(2, 3, {"osmid": 2}), (3, 2, {"osmid": 3})])
    assert ox.simplification._is_endpoint(H, 2, ["highway"], None)
    assert ox.simplification._is_endpoint(H, 2, None, ["osmid"])

    B = nx.MultiDiGraph([(1, 2), (2, 3), (3, 1)])
    assert ox.simplification._build_path(B, 1, 2, {1}) == [1, 2, 3, 1]
    C = nx.MultiDiGraph([(1, 2), (2, 3)])
    assert ox.simplification._build_path(C, 1, 2, {1}) == [1, 2, 3]
    D = nx.MultiDiGraph([(1, 2), (2, 3), (3, 4), (3, 5)])
    with pytest.raises(ox._errors.GraphSimplificationError, match="Impossible"):
        ox.simplification._build_path(D, 1, 2, {1, 4, 5})

    R = nx.MultiDiGraph(crs="epsg:4326")
    R.add_edges_from([(1, 2), (2, 3), (3, 1)])
    assert len(ox.simplification._remove_rings(R, None, None)) == 0


@pytest.mark.offline
def test_consolidation_merges_nearby_intersections() -> None:
    Gm = nx.MultiDiGraph(crs="epsg:3857")
    Gm.add_node(1, x=0.0, y=0.0, street_count=2, elevation=1.0, color="red")
    Gm.add_node(2, x=1.0, y=0.0, street_count=2, elevation=3.0, color="blue")
    Gm.add_node(3, x=5.0, y=0.0, street_count=1, elevation=5.0, color="green")
    Gm.add_edge(1, 2, osmid=1, length=1.0, geometry=LineString([(0, 0), (1, 0)]))
    Gm.add_edge(2, 3, osmid=2, length=4.0, geometry=LineString([(1, 0), (5, 0)]))
    Gc = ox.consolidate_intersections(
        Gm,
        tolerance=2,
        dead_ends=True,
        node_attr_aggs={"elevation": "mean"},
    )
    merged_nodes = [
        data for _, data in Gc.nodes(data=True) if isinstance(data["osmid_original"], list)
    ]
    assert merged_nodes[0]["elevation"] == 2.0
    assert set(merged_nodes[0]["color"]) == {"red", "blue"}
    assert any(data["length"] > 4.0 for _, _, data in Gc.edges(data=True))

    Gm.graph["consolidated"] = True
    with pytest.raises(ox._errors.GraphSimplificationError, match="already been consolidated"):
        ox.consolidate_intersections(Gm, tolerance=2, dead_ends=True)

    Gsplit = nx.MultiDiGraph(crs="epsg:3857")
    Gsplit.add_node(1, x=0.0, y=0.0, street_count=2)
    Gsplit.add_node(2, x=1.0, y=0.0, street_count=2)
    Gsplit.add_node(3, x=10.0, y=0.0, street_count=1)
    Gsplit.add_node(4, x=11.0, y=0.0, street_count=1)
    Gsplit.add_edge(1, 3, osmid=13, length=10.0, geometry=LineString([(0, 0), (10, 0)]))
    Gsplit.add_edge(2, 4, osmid=24, length=10.0, geometry=LineString([(1, 0), (11, 0)]))
    Gc_split = ox.consolidate_intersections(Gsplit, tolerance=2, dead_ends=True)
    assert len(Gc_split) == 4


@pytest.mark.offline
def test_feature_processing_filters_tags_and_geometry() -> None:
    outer_line = LineString([(0, 0), (4, 0), (4, 4), (0, 4), (0, 0)])
    inner_line = LineString([(1, 1), (2, 1), (2, 2), (1, 2), (1, 1)])
    inner_poly = Polygon([(2.5, 2.5), (3, 2.5), (3, 3), (2.5, 3)])
    relation_geom = ox.features._build_relation_geometry(
        [
            {"type": "way", "ref": 1, "role": "outer"},
            {"type": "way", "ref": 2, "role": "inner"},
            {"type": "way", "ref": 3, "role": "inner"},
        ],
        {1: outer_line, 2: inner_line, 3: inner_poly},
    )
    assert relation_geom.area > 0
    assert ox.features._build_relation_geometry(
        [{"type": "way", "ref": 999, "role": "outer"}],
        {},
    ).is_empty
    assert (
        ox.features._remove_polygon_holes([Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])], []).area == 1
    )
    assert ox.features._build_way_geometry(1, [1, 2], {}, {}).is_empty

    idx = pd.MultiIndex.from_tuples(
        [("node", 1), ("node", 2), ("node", 3)],
        names=["element", "id"],
    )
    gdf = gpd.GeoDataFrame(
        {
            "amenity": ["cafe", "school", None],
            "landuse": [None, "retail", "industrial"],
            "geometry": [Point(0.5, 0.5), Point(2, 2), Point(5, 5)],
        },
        index=idx,
        crs=ox.settings.default_crs,
    )
    polygon = Polygon([(0, 0), (3, 0), (3, 3), (0, 3)])
    filtered = ox.features._filter_features(
        gdf,
        polygon,
        {"amenity": "cafe", "landuse": ["retail"]},
    )
    assert filtered.index.to_list() == [("node", 1), ("node", 2)]
    with (
        suppress_type_checks(),
        pytest.raises(
            ox._errors.InsufficientResponseError,
            match="No matching features",
        ),
    ):
        ox.features._filter_features(gdf, polygon, {"amenity": "library"})

    ox.settings.cache_only_mode = True
    with pytest.raises(ox._errors.CacheOnlyInterruptError, match="Interrupted because"):
        ox.features._create_gdf([{"elements": []}], Polygon(), {})
    ox.settings.cache_only_mode = False
    with pytest.raises(ox._errors.InsufficientResponseError, match="No matching features"):
        ox.features._create_gdf(
            [{"elements": [{"type": "node", "id": 1, "lat": 0, "lon": 0, "tags": {}}]}],
            Polygon(),
            {"amenity": True},
        )
