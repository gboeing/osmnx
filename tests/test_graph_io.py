# ruff: noqa: D103, PLR2004, S101
# numpydoc ignore=GL08,PR01,RT01
"""Offline tests for graph creation, conversion, and file IO."""

from __future__ import annotations

import bz2
import gzip
from pathlib import Path

import geopandas as gpd
import networkx as nx
import pandas as pd
import pytest
from lxml import etree
from shapely import LineString
from shapely import Point
from shapely import Polygon
from typeguard import suppress_type_checks

import osmnx as ox
from tests.helpers import LOCATION_POINT
from tests.helpers import drive_graph
from tests.helpers import toy_graph


@pytest.mark.integration
def test_stats_simplification_and_conversion(http_cache: Path) -> None:
    ox.settings.cache_folder = http_cache
    G = drive_graph()
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


@pytest.mark.integration
def test_save_load_graph_files(http_cache: Path) -> None:
    ox.settings.cache_folder = http_cache
    G = drive_graph()
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


@pytest.mark.integration
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
    G = drive_graph()
    fp = Path(ox.settings.data_folder) / "graph.osm"
    ox.io.save_graph_xml(G, filepath=fp, way_tag_aggs={"lanes": "sum"})

    parser = etree.XMLParser(schema=etree.XMLSchema(file="./tests/input_data/osm_schema.xsd"))
    _ = etree.parse(fp, parser=parser)

    with pytest.raises(ox._errors.GraphSimplificationError):
        ox.io.save_graph_xml(ox.simplification.simplify_graph(G))
    ox.settings.all_oneway = default_all_oneway


@pytest.mark.integration
def test_graph_creation_edge_cases(http_cache: Path) -> None:
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


@pytest.mark.unit
def test_convert_and_io_edge_cases(tmp_path: Path) -> None:
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

    G_parallel = toy_graph()
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


@pytest.mark.unit
def test_osm_xml_warning_and_sort_branches(tmp_path: Path) -> None:
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
