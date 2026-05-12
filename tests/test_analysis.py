"""Offline tests for analysis, routing, distance, and plotting helpers."""

# ruff: noqa: D103, PLR2004, S101

from __future__ import annotations

import logging as lg
from pathlib import Path

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
import pytest
import requests
from shapely import LineString
from shapely import Point
from shapely import Polygon
from typeguard import suppress_type_checks

import osmnx as ox
from tests.conftest import LOCATION_POINT
from tests.conftest import _drive_graph
from tests.conftest import _Response
from tests.conftest import _toy_graph


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
