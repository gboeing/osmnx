"""Unit tests for the package."""

# use agg backend so you don't need a display on ci
# do this first before pyplot is imported by anything
import matplotlib as mpl

from osmnx.speed import _clean_maxspeed

mpl.use("Agg")

import bz2
import logging as lg
import os
import tempfile
from collections import OrderedDict
from pathlib import Path

import folium
import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
import pytest
from shapely import wkt
from shapely.geometry import LineString
from shapely.geometry import MultiLineString
from shapely.geometry import MultiPoint
from shapely.geometry import MultiPolygon
from shapely.geometry import Point
from shapely.geometry import Polygon

import osmnx as ox


ox.config(log_console=True)
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


def test_logging():
    # test OSMnx's logger
    ox.log("test a fake default message")
    ox.log("test a fake debug", level=lg.DEBUG)
    ox.log("test a fake info", level=lg.INFO)
    ox.log("test a fake warning", level=lg.WARNING)
    ox.log("test a fake error", level=lg.ERROR)

    ox.citation()
    ox.ts(style="date")
    ox.ts(style="time")


def test_coords_rounding():
    # test the rounding of geometry coordinates
    precision = 3

    shape1 = Point(1.123456, 2.123456)
    shape2 = ox.utils_geo.round_geometry_coords(shape1, precision)

    shape1 = MultiPoint([(1.123456, 2.123456), (3.123456, 4.123456)])
    shape2 = ox.utils_geo.round_geometry_coords(shape1, precision)

    shape1 = LineString([(1.123456, 2.123456), (3.123456, 4.123456)])
    shape2 = ox.utils_geo.round_geometry_coords(shape1, precision)

    shape1 = MultiLineString(
        [
            [(1.123456, 2.123456), (3.123456, 4.123456)],
            [(11.123456, 12.123456), (13.123456, 14.123456)],
        ]
    )

    shape2 = ox.utils_geo.round_geometry_coords(shape1, precision)

    shape1 = Polygon([(1.123456, 2.123456), (3.123456, 4.123456), (6.123456, 5.123456)])
    shape2 = ox.utils_geo.round_geometry_coords(shape1, precision)

    shape1 = MultiPolygon(
        [
            Polygon([(1.123456, 2.123456), (3.123456, 4.123456), (6.123456, 5.123456)]),
            Polygon([(16.123456, 15.123456), (13.123456, 14.123456), (12.123456, 11.123456)]),
        ]
    )
    shape2 = ox.utils_geo.round_geometry_coords(shape1, precision)


def test_geocode_to_gdf():
    # test loading spatial boundaries and plotting
    city = ox.geocode_to_gdf("R2999176", by_osmid=True)
    city = ox.geocode_to_gdf(place1, which_result=1, buffer_dist=100)
    city = ox.geocode_to_gdf(place2)
    city_projected = ox.project_gdf(city, to_crs="epsg:3395")


def test_stats():
    # create graph, add a new node, add bearings, project it
    G = ox.graph_from_place(place1, network_type="all")
    G.add_node(0, x=location_point[1], y=location_point[0])
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
    fig, ax = ox.bearing.plot_orientation(Gu, area=True, title="Title")
    fig, ax = ox.plot_orientation(Gu, area=True, title="Title")
    fig, ax = ox.plot_orientation(Gu, ax=ax, area=False, title="Title")

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


def test_osm_xml():
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
            assert G.edges[edge_key]["name"] in ("8th Street", "Willow Street")

    os.remove(temp_filename)

    # test .osm xml saving
    default_all_oneway = ox.settings.all_oneway
    ox.settings.all_oneway = True
    G = ox.graph_from_point(location_point, dist=500, network_type="drive")
    ox.save_graph_xml(G, merge_edges=False, filepath=Path(ox.settings.data_folder) / "graph.osm")

    # test osm xml output merge edges
    ox.save_graph_xml(G, merge_edges=True, edge_tag_aggs=[("length", "sum")])

    # test osm xml output from gdfs
    nodes, edges = ox.graph_to_gdfs(G)
    ox.save_graph_xml([nodes, edges])

    # test ordered nodes from way
    df = pd.DataFrame({"u": [54, 2, 5, 3, 10, 19, 20], "v": [76, 3, 8, 10, 5, 20, 15]})
    ordered_nodes = ox.osm_xml._get_unique_nodes_ordered_from_way(df)
    assert ordered_nodes == [2, 3, 10, 5, 8]

    ox.settings.all_oneway = default_all_oneway


def test_elevation():
    G = ox.graph_from_address(address=address, dist=500, dist_type="bbox", network_type="bike")
    rasters = list(Path("tests/input_data").glob("elevation*.tif"))

    # add node elevations from a single raster file (some nodes will be null)
    G = ox.elevation.add_node_elevations_raster(G, rasters[0], cpus=1)

    # add node elevations from multiple raster files
    G = ox.elevation.add_node_elevations_raster(G, rasters)
    assert pd.notnull(pd.Series(dict(G.nodes(data="elevation")))).all()

    # add edge grades and their absolute values
    G = ox.add_edge_grades(G, add_absolute=True)


def test_routing():
    G = ox.graph_from_address(address=address, dist=500, dist_type="bbox", network_type="bike")

    # give each edge speed and travel time attributes
    G = ox.add_edge_speeds(G)
    G = ox.add_edge_speeds(G, hwy_speeds={"motorway": 100})
    G = ox.add_edge_travel_times(G)

    # test value cleaning
    assert _clean_maxspeed("100,2") == 100.2
    assert _clean_maxspeed("100.2") == 100.2
    assert _clean_maxspeed("100 km/h") == 100.0
    assert _clean_maxspeed("100 mph") == pytest.approx(160.934)
    assert _clean_maxspeed("60|100") == 80
    assert _clean_maxspeed("60|100 mph") == pytest.approx(128.7472)
    assert _clean_maxspeed("signal") is None
    assert _clean_maxspeed("100;70") is None

    orig_x = np.array([-122.404771])
    dest_x = np.array([-122.401429])
    orig_y = np.array([37.794302])
    dest_y = np.array([37.794987])
    orig_node = ox.distance.nearest_nodes(G, orig_x, orig_y)[0]
    dest_node = ox.distance.nearest_nodes(G, dest_x, dest_y)[0]

    route = ox.shortest_path(G, orig_node, dest_node, weight="travel_time")
    route_edges = ox.utils_graph.route_to_gdf(G, route, "travel_time")
    attributes = ox.utils_graph.get_route_edge_attributes(G, route)
    attributes = ox.utils_graph.get_route_edge_attributes(G, route, "travel_time")

    fig, ax = ox.plot_graph_route(G, route, save=True)

    # test multiple origins-destinations
    n = 5
    nodes = np.array(G.nodes)
    origs = np.random.choice(nodes, size=n, replace=True)
    dests = np.random.choice(nodes, size=n, replace=True)
    paths1 = ox.shortest_path(G, origs, dests, weight="length", cpus=1)
    paths2 = ox.shortest_path(G, origs, dests, weight="length", cpus=2)
    paths3 = ox.shortest_path(G, origs, dests, weight="length", cpus=None)
    assert paths1 == paths2 == paths3

    # test k shortest paths
    routes = ox.k_shortest_paths(G, orig_node, dest_node, k=2, weight="travel_time")
    fig, ax = ox.plot_graph_routes(G, list(routes))

    # test folium with keyword arguments to pass to folium.PolyLine
    gm = ox.plot_graph_folium(G, popup_attribute="name", color="#333333", weight=5, opacity=0.7)
    rm = ox.plot_route_folium(G, route, color="#cc0000", weight=5, opacity=0.7)

    # test calling folium plotters with FeatureGroup instead of Map, and extra kwargs
    fg = folium.FeatureGroup(name="legend name", show=True)
    gm = ox.plot_graph_folium(G, graph_map=fg)
    assert isinstance(gm, folium.FeatureGroup)

    rm = ox.plot_route_folium(G, route, route_map=fg, tooltip="x")
    assert isinstance(rm, folium.FeatureGroup)


def test_plots():
    G = ox.graph_from_point(location_point, dist=500, network_type="drive")
    Gp = ox.project_graph(G)

    # test getting colors
    co = ox.plot.get_colors(n=5, return_hex=True)
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


def test_find_nearest():
    # get graph and x/y coords to search
    G = ox.graph_from_point(location_point, dist=500, network_type="drive", simplify=False)
    Gp = ox.project_graph(G)
    points = ox.utils_geo.sample_points(ox.get_undirected(Gp), 5)
    X = points.x.values
    Y = points.y.values

    # get nearest nodes
    nn0, dist0 = ox.distance.nearest_nodes(G, X[0], Y[0], return_dist=True)
    nn1, dist1 = ox.distance.nearest_nodes(Gp, X[0], Y[0], return_dist=True)

    # get nearest edge
    ne0 = ox.distance.nearest_edges(Gp, X[0], Y[0], interpolate=None)
    ne1 = ox.distance.nearest_edges(Gp, X[0], Y[0], interpolate=50)
    ne2 = ox.distance.nearest_edges(G, X[0], Y[0], interpolate=50, return_dist=True)


def test_api_endpoints():
    ip = ox.downloader._resolve_host_via_doh("overpass-api.de")
    ip = ox.downloader._resolve_host_via_doh("AAAAAAAAAAA")

    _doh_url_template_default = ox.settings.doh_url_template
    ox.settings.doh_url_template = "http://aaaaaa.hostdoesntexist.org/nothinguseful"
    ip = ox.downloader._resolve_host_via_doh("overpass-api.de")
    ox.settings.doh_url_template = None
    ip = ox.downloader._resolve_host_via_doh("overpass-api.de")
    ox.settings.doh_url_template = _doh_url_template_default

    params = OrderedDict()
    params["format"] = "json"
    params["address_details"] = 0

    # Bad Address - should return an empty response
    params["q"] = "AAAAAAAAAAA"
    response_json = ox.downloader.nominatim_request(params=params, request_type="search")

    # Good Address - should return a valid response with a valid osm_id
    params["q"] = "Newcastle A186 Westgate Rd"
    response_json = ox.downloader.nominatim_request(params=params, request_type="search")

    # Lookup
    params = OrderedDict()
    params["format"] = "json"
    params["address_details"] = 0
    params["osm_ids"] = "W68876073"

    response_json = ox.downloader.nominatim_request(params=params, request_type="lookup")

    # Invalid nominatim query type
    with pytest.raises(ValueError):
        response_json = ox.downloader.nominatim_request(params=params, request_type="xyz")

    default_key = ox.settings.nominatim_key
    default_nominatim_endpoint = ox.settings.nominatim_endpoint
    default_overpass_endpoint = ox.settings.overpass_endpoint

    # Searching on public nominatim should work even if a key was provided
    ox.settings.nominatim_key = "NOT_A_KEY"
    response_json = ox.downloader.nominatim_request(params=params, request_type="search")

    # Test changing the endpoint.
    # This should fail because we didn't provide a valid endpoint
    ox.settings.overpass_endpoint = "http://NOT_A_VALID_ENDPOINT/api/"
    with pytest.raises(Exception) as ex:
        G = ox.graph_from_place(place1, network_type="all")

    ox.settings.nominatim_key = default_key
    ox.settings.nominatim_endpoint = default_nominatim_endpoint
    ox.settings.overpass_endpoint = default_overpass_endpoint


def test_graph_save_load():
    # save graph as shapefile and geopackage
    G = ox.graph_from_point(location_point, dist=500, network_type="drive")
    ox.save_graph_shapefile(G, directed=True)
    ox.save_graph_shapefile(G, filepath=Path(ox.settings.data_folder) / "graph_shapefile")
    ox.save_graph_geopackage(G, directed=False)

    # save/load geopackage and convert graph to/from node/edge GeoDataFrames
    fp = ".temp/data/graph-dir.gpkg"
    ox.save_graph_geopackage(G, filepath=fp, directed=True)
    gdf_nodes1 = gpd.read_file(fp, layer="nodes").set_index("osmid")
    gdf_edges1 = gpd.read_file(fp, layer="edges").set_index(["u", "v", "key"])
    G2 = ox.graph_from_gdfs(gdf_nodes1, gdf_edges1)
    G2 = ox.graph_from_gdfs(gdf_nodes1, gdf_edges1, graph_attrs=G.graph)
    gdf_nodes2, gdf_edges2 = ox.graph_to_gdfs(G2)
    assert set(gdf_nodes1.index) == set(gdf_nodes2.index) == set(G.nodes) == set(G2.nodes)
    assert set(gdf_edges1.index) == set(gdf_edges2.index) == set(G.edges) == set(G2.edges)

    # create random boolean graph/node/edge attributes
    attr_name = "test_bool"
    G.graph[attr_name] = False
    bools = np.random.randint(0, 2, len(G.nodes))
    node_attrs = {n: bool(b) for n, b in zip(G.nodes, bools)}
    nx.set_node_attributes(G, node_attrs, attr_name)
    bools = np.random.randint(0, 2, len(G.edges))
    edge_attrs = {n: bool(b) for n, b in zip(G.edges, bools)}
    nx.set_edge_attributes(G, edge_attrs, attr_name)

    # save/load graph as graphml file
    ox.save_graphml(G, gephi=True)
    ox.save_graphml(G, gephi=False)
    ox.save_graphml(G, gephi=False, filepath=Path(ox.settings.data_folder) / "graph.graphml")
    filepath = Path(ox.settings.data_folder) / "graph.graphml"
    G2 = ox.load_graphml(
        filepath,
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
    G2 = ox.load_graphml(filepath, node_dtypes=nd, edge_dtypes=ed)

    # test loading graphml from a file stream
    file_bytes = open("tests/input_data/short.graphml", "rb").read()
    data = str(file_bytes.decode())
    G = ox.load_graphml(graphml_str=data, node_dtypes=nd, edge_dtypes=ed)


def test_graph_from_functions():
    # graph from bounding box
    _ = ox.utils_geo.bbox_from_point(location_point, project_utm=True, return_crs=True)
    north, south, east, west = ox.utils_geo.bbox_from_point(location_point, dist=500)
    G = ox.graph_from_bbox(north, south, east, west, network_type="drive")
    G = ox.graph_from_bbox(
        north, south, east, west, network_type="drive_service", truncate_by_edge=True
    )

    # truncate graph by bounding box
    north, south, east, west = ox.utils_geo.bbox_from_point(location_point, dist=400)
    G = ox.truncate.truncate_graph_bbox(G, north, south, east, west)
    G = ox.utils_graph.get_largest_component(G, strongly=True)

    # graph from address
    G = ox.graph_from_address(address=address, dist=500, dist_type="bbox", network_type="bike")

    # graph from list of places
    G = ox.graph_from_place([place1], network_type="all", clean_periphery=False)

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

    ox.settings.memory = "1073741824"
    G = ox.graph_from_point(
        location_point,
        dist=500,
        dist_type="network",
        network_type="all_private",
    )


def test_geometries():
    # geometries_from_bbox - bounding box query to return no data
    try:
        gdf = ox.geometries_from_bbox(0.009, -0.009, 0.009, -0.009, tags={"building": True})
    except ox._errors.EmptyOverpassResponse:
        pass

    # geometries_from_bbox - successful
    north, south, east, west = ox.utils_geo.bbox_from_point(location_point, dist=500)
    tags = {"landuse": True, "building": True, "highway": True}
    gdf = ox.geometries_from_bbox(north, south, east, west, tags=tags)
    fig, ax = ox.plot_footprints(gdf)
    fig, ax = ox.plot_footprints(gdf, ax=ax, bbox=(10, 0, 10, 0))

    # geometries_from_point - tests multipolygon creation
    gdf = ox.geometries_from_point((48.15, 10.02), tags={"landuse": True}, dist=2000)

    # geometries_from_place - includes test of list of places
    tags = {"amenity": True, "landuse": ["retail", "commercial"], "highway": "bus_stop"}
    gdf = ox.geometries_from_place(place1, tags=tags)
    gdf = ox.geometries_from_place([place1], tags=tags)

    # geometries_from_address - includes testing overpass settings and snapshot from 2019
    ox.settings.overpass_settings = '[out:json][timeout:200][date:"2019-10-28T19:20:00Z"]'
    gdf = ox.geometries_from_address(address, tags=tags)

    # geometries_from_xml - tests error handling of clipped XMLs with incomplete geometry
    gdf = ox.geometries_from_xml("tests/input_data/planet_10.068,48.135_10.071,48.137.osm")

    # test loading a geodataframe from a local .osm xml file
    with bz2.BZ2File("tests/input_data/West-Oakland.osm.bz2") as f:
        handle, temp_filename = tempfile.mkstemp(suffix=".osm")
        os.write(handle, f.read())
        os.close(handle)
    for filename in ("tests/input_data/West-Oakland.osm.bz2", temp_filename):
        gdf = ox.geometries_from_xml(filename)
        assert "Willow Street" in gdf["name"].values
    os.remove(temp_filename)
