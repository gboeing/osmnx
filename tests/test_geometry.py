"""Offline tests for validation, features, simplification, and geometry workflows."""

# ruff: noqa: D103, PLR2004, S101

from __future__ import annotations

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
import pytest
from shapely import LineString
from shapely import Point
from shapely import Polygon
from shapely import wkt
from typeguard import suppress_type_checks

import osmnx as ox


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
