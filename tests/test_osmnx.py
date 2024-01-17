# ruff: noqa: E402 F841 PLR2004
"""Test suite for the package."""

from __future__ import annotations

# use agg backend so you don't need a display on CI
# do this first before pyplot is imported by anything
import matplotlib as mpl

mpl.use("Agg")

import bz2
import logging as lg
import os
import tempfile
from collections import OrderedDict
from pathlib import Path
from xml.etree import ElementTree as etree

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
import pytest
from requests.exceptions import ConnectionError
from shapely import wkt
from shapely.geometry import Point
from shapely.geometry import Polygon
from typeguard import TypeCheckError

import osmnx as ox

ox.settings.log_console = True
ox.settings.log_file = True
ox.settings.use_cache = True
ox.settings.data_folder = ".temp/data"
ox.settings.logs_folder = ".temp/logs"
ox.settings.imgs_folder = ".temp/imgs"
ox.settings.cache_folder = ".temp/cache"

# define queries to use throughout tests
location_point = (37.791427, -122.410018)
address = "600 Montgomery St, San Francisco, California, USA"
place1 = {"city": "Piedmont", "state": "California", "country": "USA"}
place2 = "SoHo, New York, NY"
p = (
    "POLYGON ((-122.262 37.869, -122.255 37.869, -122.255 37.874,"
    "-122.262 37.874, -122.262 37.869))"
)
polygon = wkt.loads(p)


def test_logging() -> None:
    """Test the logger."""
    ox.log("test a fake default message")
    ox.log("test a fake debug", level=lg.DEBUG)
    ox.log("test a fake info", level=lg.INFO)
    ox.log("test a fake warning", level=lg.WARNING)
    ox.log("test a fake error", level=lg.ERROR)

    ox.citation(style="apa")
    ox.citation(style="bibtex")
    ox.citation(style="ieee")
    ox.ts(style="date")
    ox.ts(style="time")


def test_exceptions() -> None:
    """Test the custom errors."""
    message = "testing exception"

    with pytest.raises(ox._errors.ResponseStatusCodeError):
        raise ox._errors.ResponseStatusCodeError(message)

    with pytest.raises(ox._errors.CacheOnlyInterruptError):
        raise ox._errors.CacheOnlyInterruptError(message)

    with pytest.raises(ox._errors.InsufficientResponseError):
        raise ox._errors.InsufficientResponseError(message)

    with pytest.raises(ox._errors.GraphSimplificationError):
        raise ox._errors.GraphSimplificationError(message)


def test_geocoder() -> None:
    """Test retrieving elements by place name and OSM ID."""
    city = ox.geocode_to_gdf("R2999176", by_osmid=True)
    city = ox.geocode_to_gdf(place1, which_result=1)
    city = ox.geocode_to_gdf(place2)
    city_projected = ox.project_gdf(city, to_crs="epsg:3395")

    # test geocoding a bad query: should raise exception
    with pytest.raises(ox._errors.InsufficientResponseError):
        _ = ox.geocode("!@#$%^&*")

    with pytest.raises(ox._errors.InsufficientResponseError):
        ox.geocode_to_gdf(query="AAAZZZ")


def test_stats() -> None:
    """Test generating graph stats."""
    # create graph, add a new node, add bearings, project it
    G = ox.graph_from_place(place1, network_type="all")
    G.add_node(0, x=location_point[1], y=location_point[0], street_count=0)
    _ = ox.bearing.calculate_bearing(0, 0, 1, 1)
    G = ox.add_edge_bearings(G)
    G_proj = ox.project_graph(G)
    G_proj = ox.distance.add_edge_lengths(G_proj, edges=tuple(G_proj.edges)[0:3])

    # calculate stats
    cspn = ox.stats.count_streets_per_node(G)
    stats = ox.basic_stats(G)
    stats = ox.basic_stats(G, area=1000)
    stats = ox.basic_stats(G_proj, area=1000, clean_int_tol=15)

    # calculate entropy
    Gu = ox.get_undirected(G)
    entropy = ox.bearing.orientation_entropy(Gu, weight="length")

    fig, ax = ox.plot.plot_orientation(Gu, area=True, title="Title")
    fig, ax = ox.plot.plot_orientation(Gu, ax=ax, area=False, title="Title")

    # test cleaning and rebuilding graph
    G_clean = ox.consolidate_intersections(G_proj, tolerance=10, rebuild_graph=True, dead_ends=True)
    G_clean = ox.consolidate_intersections(
        G_proj, tolerance=10, rebuild_graph=True, reconnect_edges=False
    )
    G_clean = ox.consolidate_intersections(G_proj, tolerance=10, rebuild_graph=False)

    # try consolidating an empty graph
    G = nx.MultiDiGraph(crs="epsg:4326")
    G_clean = ox.consolidate_intersections(G, rebuild_graph=True)
    G_clean = ox.consolidate_intersections(G, rebuild_graph=False)


def test_osm_xml() -> None:
    """Test working with .osm XML data."""
    # test loading a graph from a local .osm xml file
    node_id = 53098262
    neighbor_ids = 53092170, 53060438, 53027353, 667744075

    with bz2.BZ2File("tests/input_data/West-Oakland.osm.bz2") as f:
        handle, temp_filename = tempfile.mkstemp(suffix=".osm")
        os.write(handle, f.read())
        os.close(handle)

    for filename in ("tests/input_data/West-Oakland.osm.bz2", temp_filename):
        G = ox.graph_from_xml(filename)
        assert node_id in G.nodes

        for neighbor_id in neighbor_ids:
            edge_key = (node_id, neighbor_id, 0)
            assert neighbor_id in G.nodes
            assert edge_key in G.edges
            assert G.edges[edge_key]["name"] in {"8th Street", "Willow Street"}

    Path.unlink(Path(temp_filename))

    # test .osm xml saving
    G = ox.graph_from_point(location_point, dist=500, network_type="drive")
    fp = Path(ox.settings.data_folder) / "graph.osm"
    ox.save_graph_xml(G, merge_edges=False, filepath=fp)  # issues UserWarning
    G = ox.graph_from_xml(fp)  # issues UserWarning

    # test osm xml output merge edges
    default_all_oneway = ox.settings.all_oneway
    ox.settings.all_oneway = True
    ox.io.save_graph_xml(G, merge_edges=True, edge_tag_aggs=[("length", "sum")], precision=5)

    # test osm xml output from gdfs
    nodes, edges = ox.graph_to_gdfs(G)
    ox.io.save_graph_xml((nodes, edges))

    # test ordered nodes from way
    df_uv = pd.DataFrame({"u": [54, 2, 5, 3, 10, 19, 20], "v": [76, 3, 8, 10, 5, 20, 15]})
    ordered_nodes = ox._osm_xml._get_unique_nodes_ordered_from_way(df_uv)
    assert ordered_nodes == [2, 3, 10, 5, 8]

    # test roundabout handling
    default_overpass_settings = ox.settings.overpass_settings
    ox.settings.overpass_settings += '[date:"2023-04-01T00:00:00Z"]'
    point = (39.0290346, -84.4696884)
    G = ox.graph_from_point(point, dist=500, dist_type="bbox", network_type="drive", simplify=False)
    gdf_edges = ox.graph_to_gdfs(G, nodes=False)
    gdf_way = gdf_edges[gdf_edges["osmid"] == 570883705]  # roundabout
    first = gdf_way.iloc[0].dropna().astype(str)
    root = etree.Element("osm", attrib={"version": "0.6", "generator": "OSMnx"})
    edge = etree.SubElement(root, "way")
    ox._osm_xml._append_nodes_as_edge_attrs(edge, first.to_dict(), gdf_way)

    # restore settings
    ox.settings.overpass_settings = default_overpass_settings
    ox.settings.all_oneway = default_all_oneway


def test_elevation() -> None:
    """Test working with elevation data."""
    G = ox.graph_from_address(address=address, dist=500, dist_type="bbox", network_type="bike")

    # add node elevations from Google (fails without API key)
    with pytest.raises(ox._errors.InsufficientResponseError):
        _ = ox.elevation.add_node_elevations_google(G, api_key="", batch_size=350)

    # add node elevations from Open Topo Data (works without API key)
    ox.settings.elevation_url_template = (
        "https://api.opentopodata.org/v1/aster30m?locations={locations}&key={key}"
    )
    _ = ox.elevation.add_node_elevations_google(G, batch_size=100, pause=0.01)

    # same thing again, to hit the cache
    _ = ox.elevation.add_node_elevations_google(G, batch_size=100, pause=0.01)

    # add node elevations from a single raster file (some nodes will be null)
    rasters = list(Path("tests/input_data").glob("elevation*.tif"))
    G = ox.elevation.add_node_elevations_raster(G, rasters[0], cpus=1)
    assert pd.notna(pd.Series(dict(G.nodes(data="elevation")))).any()

    # add node elevations from multiple raster files (no nodes should be null)
    G = ox.elevation.add_node_elevations_raster(G, rasters)
    assert pd.notna(pd.Series(dict(G.nodes(data="elevation")))).all()

    # add edge grades and their absolute values
    G = ox.add_edge_grades(G, add_absolute=True)


def test_routing() -> None:
    """Test working with speed, travel time, and routing."""
    G = ox.graph_from_address(address=address, dist=500, dist_type="bbox", network_type="bike")

    # give each edge speed and travel time attributes
    G = ox.add_edge_speeds(G)
    G = ox.add_edge_speeds(G, hwy_speeds={"motorway": 100})
    G = ox.add_edge_travel_times(G)

    # test value cleaning
    assert ox.speed._clean_maxspeed("100,2") == 100.2
    assert ox.speed._clean_maxspeed("100.2") == 100.2
    assert ox.speed._clean_maxspeed("100 km/h") == 100.0
    assert ox.speed._clean_maxspeed("100 mph") == pytest.approx(160.934)
    assert ox.speed._clean_maxspeed("60|100") == 80
    assert ox.speed._clean_maxspeed("60|100 mph") == pytest.approx(128.7472)
    assert ox.speed._clean_maxspeed("signal") is None
    assert ox.speed._clean_maxspeed("100;70") is None

    # test collapsing multiple mph values to single kph value
    assert ox.speed._collapse_multiple_maxspeed_values(["25 mph", "30 mph"], np.mean) == 44.25685

    # test collapsing invalid values: should return None
    assert ox.speed._collapse_multiple_maxspeed_values(["mph", "kph"], np.mean) is None

    orig_x = np.array([-122.404771])
    dest_x = np.array([-122.401429])
    orig_y = np.array([37.794302])
    dest_y = np.array([37.794987])
    orig_node = int(ox.distance.nearest_nodes(G, orig_x, orig_y)[0])
    dest_node = int(ox.distance.nearest_nodes(G, dest_x, dest_y)[0])

    # test non-numeric weight, should raise ValueError
    with pytest.raises(ValueError, match="contains non-numeric values"):
        route1 = ox.shortest_path(G, orig_node, dest_node, weight="highway")

    # mismatch iterable and non-iterable orig/dest, should raise TypeError
    msg = "orig and dest must either both be iterable or neither must be iterable"
    with pytest.raises(TypeError, match=msg):
        route2 = ox.shortest_path(G, orig_node, [dest_node])  # type: ignore[call-overload]

    # mismatch lengths of orig/dest, should raise ValueError
    msg = "orig and dest must be of equal length"
    with pytest.raises(ValueError, match=msg):
        route2 = ox.shortest_path(G, [orig_node] * 2, [dest_node] * 3)

    # test missing weight (should raise warning)
    route3 = ox.shortest_path(G, orig_node, dest_node, weight="time")
    # test good weight
    route4 = ox.routing.shortest_path(G, orig_node, dest_node, weight="travel_time")
    route5 = ox.shortest_path(G, orig_node, dest_node, weight="travel_time")
    assert route5 is not None

    route_edges = ox.utils_graph.route_to_gdf(G, route5, "travel_time")

    fig, ax = ox.plot_graph_route(G, route5, save=True)

    # test multiple origins-destinations
    n = 5
    nodes = np.array(G.nodes)
    origs = [int(x) for x in np.random.default_rng().choice(nodes, size=n, replace=True)]
    dests = [int(x) for x in np.random.default_rng().choice(nodes, size=n, replace=True)]
    paths1 = ox.shortest_path(G, origs, dests, weight="length", cpus=1)
    paths2 = ox.shortest_path(G, origs, dests, weight="length", cpus=2)
    paths3 = ox.shortest_path(G, origs, dests, weight="length", cpus=None)
    assert paths1 == paths2 == paths3

    # test k shortest paths
    routes = ox.routing.k_shortest_paths(G, orig_node, dest_node, k=2, weight="travel_time")
    fig, ax = ox.plot_graph_routes(G, list(routes))

    # test great circle and euclidean distance calculators
    assert ox.distance.great_circle(0, 0, 1, 1) == pytest.approx(157249.6034105)
    assert ox.distance.euclidean(0, 0, 1, 1) == pytest.approx(1.4142135)


def test_plots() -> None:
    """Test visualization methods."""
    G = ox.graph_from_point(location_point, dist=500, network_type="drive")
    Gp = ox.project_graph(G)
    G = ox.project_graph(G, to_latlong=True)

    # test getting colors
    co1 = ox.plot.get_colors(n=5, cmap="plasma", start=0.1, stop=0.9, alpha=0.5)
    co2 = ox.plot.get_colors(n=5, cmap="plasma", start=0.1, stop=0.9, alpha=None)
    nc = ox.plot.get_node_colors_by_attr(G, "x")
    ec = ox.plot.get_edge_colors_by_attr(G, "length", num_bins=5)

    # plot and save to disk
    filepath = Path(ox.settings.data_folder) / "test.svg"
    fig, ax = ox.plot_graph(G, show=False, save=True, close=True, filepath=filepath)
    fig, ax = ox.plot_graph(Gp, edge_linewidth=0, figsize=(5, 5), bgcolor="y")
    fig, ax = ox.plot_graph(
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

    # figure-ground plots
    fig, ax = ox.plot_figure_ground(G=G)
    fig, ax = ox.plot_figure_ground(point=location_point, dist=500, network_type="drive")
    fig, ax = ox.plot_figure_ground(address=address, dist=500, network_type="bike")


def test_find_nearest() -> None:
    """Test nearest node/edge searching."""
    # get graph and x/y coords to search
    G = ox.graph_from_point(location_point, dist=500, network_type="drive", simplify=False)
    Gp = ox.project_graph(G)
    points = ox.utils_geo.sample_points(ox.get_undirected(Gp), 5)
    X = points.x.to_numpy()
    Y = points.y.to_numpy()

    # get nearest nodes
    _ = ox.distance.nearest_nodes(G, X, Y, return_dist=True)
    _ = ox.distance.nearest_nodes(G, X, Y, return_dist=False)
    nn0, dist0 = ox.distance.nearest_nodes(G, X[0], Y[0], return_dist=True)
    nn1 = ox.distance.nearest_nodes(Gp, X[0], Y[0], return_dist=False)

    # get nearest edge
    _ = ox.distance.nearest_edges(Gp, X, Y, return_dist=False)
    _ = ox.distance.nearest_edges(Gp, X, Y, return_dist=True)
    _ = ox.distance.nearest_edges(Gp, X[0], Y[0], return_dist=False)
    _ = ox.distance.nearest_edges(Gp, X[0], Y[0], return_dist=True)


def test_api_endpoints() -> None:
    """Test different API endpoints."""
    default_timeout = ox.settings.timeout
    default_key = ox.settings.nominatim_key
    default_nominatim_endpoint = ox.settings.nominatim_endpoint
    default_overpass_endpoint = ox.settings.overpass_endpoint
    default_overpass_rate_limit = ox.settings.overpass_rate_limit

    # test good and bad DNS resolution
    ox.settings.timeout = 1
    ip = ox._downloader._resolve_host_via_doh("overpass-api.de")
    ip = ox._downloader._resolve_host_via_doh("AAAAAAAAAAA")
    _doh_url_template_default = ox.settings.doh_url_template
    ox.settings.doh_url_template = "http://aaaaaa.hostdoesntexist.org/nothinguseful"
    ip = ox._downloader._resolve_host_via_doh("overpass-api.de")
    ox.settings.doh_url_template = None
    ip = ox._downloader._resolve_host_via_doh("overpass-api.de")
    ox.settings.doh_url_template = _doh_url_template_default

    # Test changing the Overpass endpoint.
    # This should fail because we didn't provide a valid endpoint
    ox.settings.overpass_rate_limit = False
    ox.settings.overpass_endpoint = "http://NOT_A_VALID_ENDPOINT/api/"
    with pytest.raises(ConnectionError, match="Max retries exceeded with url"):
        G = ox.graph_from_place(place1, network_type="all")

    ox.settings.overpass_rate_limit = default_overpass_rate_limit
    ox.settings.timeout = default_timeout

    params: OrderedDict[str, int | str] = OrderedDict()
    params["format"] = "json"
    params["address_details"] = 0

    # Bad Address - should return an empty response
    params["q"] = "AAAAAAAAAAA"
    response_json = ox._nominatim._nominatim_request(params=params, request_type="search")

    # Good Address - should return a valid response with a valid osm_id
    params["q"] = "Newcastle A186 Westgate Rd"
    response_json = ox._nominatim._nominatim_request(params=params, request_type="search")

    # Lookup
    params = OrderedDict()
    params["format"] = "json"
    params["address_details"] = 0
    params["osm_ids"] = "W68876073"

    # good call
    response_json = ox._nominatim._nominatim_request(params=params, request_type="lookup")

    # bad call
    with pytest.raises(
        ox._errors.InsufficientResponseError, match="Nominatim API did not return a list of results"
    ):
        response_json = ox._nominatim._nominatim_request(params=params, request_type="search")

    # query must be a str if by_osmid=True
    with pytest.raises(TypeError, match="`query` must be a string if `by_osmid` is True"):
        ox.geocode_to_gdf(query={"City": "Boston"}, by_osmid=True)

    # Invalid nominatim query type
    with pytest.raises(ValueError, match="Nominatim request_type must be"):
        response_json = ox._nominatim._nominatim_request(params=params, request_type="xyz")

    # Searching on public nominatim should work even if a (bad) key was provided
    ox.settings.nominatim_key = "NOT_A_KEY"
    response_json = ox._nominatim._nominatim_request(params=params, request_type="lookup")

    ox.settings.nominatim_key = default_key
    ox.settings.nominatim_endpoint = default_nominatim_endpoint
    ox.settings.overpass_endpoint = default_overpass_endpoint


def test_graph_save_load() -> None:
    """Test saving/loading graphs to/from disk."""
    G = ox.graph_from_point(location_point, dist=500, network_type="drive")

    # save/load geopackage and convert graph to/from node/edge GeoDataFrames
    ox.save_graph_geopackage(G, directed=False)
    fp = ".temp/data/graph-dir.gpkg"
    ox.save_graph_geopackage(G, filepath=fp, directed=True)
    gdf_nodes1 = gpd.read_file(fp, layer="nodes").set_index("osmid")
    gdf_edges1 = gpd.read_file(fp, layer="edges").set_index(["u", "v", "key"])
    G2 = ox.graph_from_gdfs(gdf_nodes1, gdf_edges1)
    G2 = ox.graph_from_gdfs(gdf_nodes1, gdf_edges1, graph_attrs=G.graph)
    gdf_nodes2, gdf_edges2 = ox.graph_to_gdfs(G2)
    _ = list(ox.utils_geo.interpolate_points(gdf_edges2["geometry"].iloc[0], 0.001))
    assert set(gdf_nodes1.index) == set(gdf_nodes2.index) == set(G.nodes) == set(G2.nodes)
    assert set(gdf_edges1.index) == set(gdf_edges2.index) == set(G.edges) == set(G2.edges)

    # test code branches that should raise exceptions
    with pytest.raises(ValueError, match="you must request nodes or edges or both"):
        ox.graph_to_gdfs(G2, nodes=False, edges=False)
    with pytest.raises(ValueError, match="invalid literal for boolean"):
        ox.io._convert_bool_string("T")

    # create random boolean graph/node/edge attributes
    attr_name = "test_bool"
    G.graph[attr_name] = False
    bools = np.random.default_rng().integers(low=0, high=2, size=len(G.nodes))
    node_attrs = {n: bool(b) for n, b in zip(G.nodes, bools)}
    nx.set_node_attributes(G, node_attrs, attr_name)
    bools = np.random.default_rng().integers(low=0, high=2, size=len(G.edges))
    edge_attrs = {n: bool(b) for n, b in zip(G.edges, bools)}
    nx.set_edge_attributes(G, edge_attrs, attr_name)

    # create list, set, and dict attributes for nodes and edges
    rand_ints_nodes = np.random.default_rng().integers(low=0, high=10, size=len(G.nodes))
    rand_ints_edges = np.random.default_rng().integers(low=0, high=10, size=len(G.edges))
    list_node_attrs = {n: [n, r] for n, r in zip(G.nodes, rand_ints_nodes)}
    nx.set_node_attributes(G, list_node_attrs, "test_list")
    list_edge_attrs = {e: [e, r] for e, r in zip(G.edges, rand_ints_edges)}
    nx.set_edge_attributes(G, list_edge_attrs, "test_list")
    set_node_attrs = {n: {n, r} for n, r in zip(G.nodes, rand_ints_nodes)}
    nx.set_node_attributes(G, set_node_attrs, "test_set")
    set_edge_attrs = {e: {e, r} for e, r in zip(G.edges, rand_ints_edges)}
    nx.set_edge_attributes(G, set_edge_attrs, "test_set")
    dict_node_attrs = {n: {n: r} for n, r in zip(G.nodes, rand_ints_nodes)}
    nx.set_node_attributes(G, dict_node_attrs, "test_dict")
    dict_edge_attrs = {e: {e: r} for e, r in zip(G.edges, rand_ints_edges)}
    nx.set_edge_attributes(G, dict_edge_attrs, "test_dict")

    # save/load graph as graphml file
    ox.save_graphml(G, gephi=True)
    ox.save_graphml(G, gephi=False)
    ox.save_graphml(G, gephi=False, filepath=fp)
    G2 = ox.load_graphml(
        fp,
        graph_dtypes={attr_name: ox.io._convert_bool_string},
        node_dtypes={attr_name: ox.io._convert_bool_string},
        edge_dtypes={attr_name: ox.io._convert_bool_string},
    )

    # verify everything in G is equivalent in G2
    assert tuple(G.graph.keys()) == tuple(G2.graph.keys())
    assert tuple(G.graph.values()) == tuple(G2.graph.values())
    z = zip(G.nodes(data=True), G2.nodes(data=True))
    for (n1, d1), (n2, d2) in z:
        assert n1 == n2
        assert tuple(d1.keys()) == tuple(d2.keys())
        assert tuple(d1.values()) == tuple(d2.values())
    z = zip(G.edges(keys=True, data=True), G2.edges(keys=True, data=True))
    for (u1, v1, k1, d1), (u2, v2, k2, d2) in z:
        assert u1 == u2
        assert v1 == v2
        assert k1 == k2
        assert tuple(d1.keys()) == tuple(d2.keys())
        assert tuple(d1.values()) == tuple(d2.values())

    # test custom data types
    nd = {"osmid": str}
    ed = {"length": str, "osmid": float}
    G2 = ox.load_graphml(fp, node_dtypes=nd, edge_dtypes=ed)

    # test loading graphml from a file stream
    file_bytes = Path.open(Path("tests/input_data/short.graphml"), "rb").read()
    data = str(file_bytes.decode())
    G = ox.load_graphml(graphml_str=data, node_dtypes=nd, edge_dtypes=ed)


def test_graph_from_functions() -> None:
    """Test downloading graphs from Overpass."""
    # test subdividing a large geometry (raises a UserWarning)
    bbox = ox.utils_geo.bbox_from_point((0, 0), dist=1e5, project_utm=True)
    poly = ox.utils_geo.bbox_to_poly(bbox)
    _ = ox.utils_geo._consolidate_subdivide_geometry(poly)

    # graph from bounding box
    _ = ox.utils_geo.bbox_from_point(location_point, project_utm=True, return_crs=True)
    bbox = ox.utils_geo.bbox_from_point(location_point, dist=500)
    G = ox.graph_from_bbox(bbox, network_type="drive")
    G = ox.graph_from_bbox(bbox, network_type="drive_service", truncate_by_edge=True)

    # truncate graph by bounding box
    bbox = ox.utils_geo.bbox_from_point(location_point, dist=400)
    G = ox.truncate.truncate_graph_bbox(G, bbox)
    G = ox.utils_graph.get_largest_component(G, strongly=True)

    # graph from address
    G = ox.graph_from_address(address=address, dist=500, dist_type="bbox", network_type="bike")

    # graph from list of places
    G = ox.graph_from_place([place1], which_result=[None], network_type="all")

    # graph from polygon
    G = ox.graph_from_polygon(polygon, network_type="walk", truncate_by_edge=True, simplify=False)
    G = ox.simplify_graph(G, strict=False, remove_rings=False, track_merged=True)

    # test custom query filter
    cf = (
        '["highway"]'
        '["area"!~"yes"]'
        '["highway"!~"motor|proposed|construction|abandoned|platform|raceway"]'
        '["foot"!~"no"]'
        '["service"!~"private"]'
        '["access"!~"private"]'
    )
    G = ox.graph_from_point(
        location_point, dist=500, custom_filter=cf, dist_type="bbox", network_type="all"
    )

    ox.settings.memory = 1073741824
    G = ox.graph_from_point(
        location_point,
        dist=500,
        dist_type="network",
        network_type="all_private",
    )


def test_features() -> None:
    """Test downloading features from Overpass."""
    bbox = ox.utils_geo.bbox_from_point(location_point, dist=500)
    tags1: dict[str, bool | str | list[str]] = {"landuse": True, "building": True, "highway": True}

    with pytest.raises(ValueError, match="The geometry of `polygon` is invalid"):
        ox.features.features_from_polygon(Polygon(((0, 0), (0, 0), (0, 0), (0, 0))), tags={})
    with pytest.raises(TypeCheckError):
        ox.features.features_from_polygon(Point(0, 0), tags={})

    # test cache_only_mode
    ox.settings.cache_only_mode = True
    with pytest.raises(ox._errors.CacheOnlyInterruptError, match="Interrupted because"):
        _ = ox.features_from_bbox(bbox, tags=tags1)
    ox.settings.cache_only_mode = False

    # features_from_bbox - bounding box query to return no data
    with pytest.raises(ox._errors.InsufficientResponseError):
        gdf = ox.features_from_bbox(bbox=(0.009, -0.009, 0.009, -0.009), tags={"building": True})

    # features_from_bbox - successful
    gdf = ox.features_from_bbox(bbox, tags=tags1)
    fig, ax = ox.plot_footprints(gdf)
    fig, ax = ox.plot_footprints(gdf, ax=ax, bbox=(10, 0, 10, 0))

    # features_from_point - tests multipolygon creation
    gdf = ox.utils_geo.bbox_from_point(location_point, dist=500)

    # features_from_place - includes test of list of places
    tags2: dict[str, bool | str | list[str]] = {
        "amenity": True,
        "landuse": ["retail", "commercial"],
        "highway": "bus_stop",
    }
    gdf = ox.features_from_place(place1, tags=tags2)
    gdf = ox.features_from_place([place1], which_result=[None], tags=tags2)

    # features_from_polygon
    polygon = ox.geocode_to_gdf(place1).geometry.iloc[0]
    ox.features_from_polygon(polygon, tags2)

    # features_from_address - includes testing overpass settings and snapshot from 2019
    ox.settings.overpass_settings = '[out:json][timeout:200][date:"2019-10-28T19:20:00Z"]'
    gdf = ox.features_from_address(address, tags=tags2)

    # features_from_xml - tests error handling of clipped XMLs with incomplete geometry
    gdf = ox.features_from_xml("tests/input_data/planet_10.068,48.135_10.071,48.137.osm")

    # test loading a geodataframe from a local .osm xml file
    with bz2.BZ2File("tests/input_data/West-Oakland.osm.bz2") as f:
        handle, temp_filename = tempfile.mkstemp(suffix=".osm")
        os.write(handle, f.read())
        os.close(handle)
    for filename in ("tests/input_data/West-Oakland.osm.bz2", temp_filename):
        gdf = ox.features_from_xml(filename)
        assert "Willow Street" in gdf["name"].to_numpy()
    Path.unlink(Path(temp_filename))
