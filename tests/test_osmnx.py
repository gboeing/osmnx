"""
OSMnx tests
-----------
"""

import matplotlib as mpl
mpl.use('Agg') #use agg backend so you don't need a display on travis-ci

import os, shutil, bz2, tempfile
if os.path.exists('.temp'):
    shutil.rmtree('.temp')

import osmnx as ox, logging as lg
ox.config(log_console=True, log_file=True, use_cache=True,
          data_folder='.temp/data', logs_folder='.temp/logs', imgs_folder='.temp/imgs', cache_folder='.temp/cache')

ox.log('test debug', level=lg.DEBUG)
ox.log('test info', level=lg.INFO)
ox.log('test warning', level=lg.WARNING)
ox.log('test error', level=lg.ERROR)


def test_imports():

    import json, math, sys, os, io, ast, unicodedata, hashlib, re, random, time, warnings, datetime as dt, logging as lg
    from collections import OrderedDict, Counter
    from itertools import groupby, chain
    from dateutil import parser as date_parser
    import requests, numpy as np, pandas as pd, geopandas as gpd, networkx as nx, matplotlib.pyplot as plt, matplotlib.cm as cm
    from matplotlib.collections import LineCollection
    from shapely.geometry import Point, LineString, Polygon, MultiPolygon
    from shapely import wkt
    from shapely.ops import unary_union
    from descartes import PolygonPatch
    from rtree.index import Index as RTreeIndex


def test_gdf_shapefiles():

    city = ox.gdf_from_place('Manhattan, New York City, New York, USA')
    city_projected = ox.project_gdf(city, to_crs={'init':'epsg:3395'})
    ox.save_gdf_shapefile(city_projected)

    city = ox.gdf_from_place('Manhattan, New York City, New York, USA', buffer_dist=100)
    ox.plot_shape(city)


def test_graph_from_file():

    node_id = 53098262
    neighbor_ids = 53092170, 53060438, 53027353, 667744075

    with bz2.BZ2File('tests/input_data/West-Oakland.osm.bz2') as input:
        handle, temp_filename = tempfile.mkstemp(suffix='.osm')
        os.write(handle, input.read())
        os.close(handle)

    for filename in ('tests/input_data/West-Oakland.osm.bz2', temp_filename):
        G = ox.graph_from_file(filename)
        assert node_id in G.nodes

        for neighbor_id in neighbor_ids:
            edge_key = (node_id, neighbor_id, 0)
            assert neighbor_id in G.nodes
            assert edge_key in G.edges
            assert G.edges[edge_key]['name'] in ('8th Street', 'Willow Street')

    os.remove(temp_filename)


def test_network_saving_loading():

    G = ox.graph_from_place('Piedmont, California, USA')
    G_projected = ox.project_graph(G)
    ox.save_graph_shapefile(G_projected)
    ox.save_graphml(G_projected)
    G2 = ox.load_graphml('graph.graphml')

    gdf_edges = ox.graph_to_gdfs(G, nodes=False, edges=True, fill_edge_geometry=False)
    gdf_nodes, gdf_edges = ox.graph_to_gdfs(G, nodes=True, edges=True, node_geometry=True, fill_edge_geometry=True)
    G3 = ox.gdfs_to_graph(gdf_nodes, gdf_edges)


def test_get_network_methods():

    import geopandas as gpd
    north, south, east, west = 37.79, 37.78, -122.41, -122.43
    G1 = ox.graph_from_bbox(north, south, east, west, network_type='drive_service')
    G1 = ox.graph_from_bbox(north, south, east, west, network_type='drive_service', truncate_by_edge=True)

    location_point = (37.791427, -122.410018)
    bbox = ox.bbox_from_point(location_point, project_utm=True)
    G2 = ox.graph_from_point(location_point, distance=750, distance_type='bbox', network_type='drive')
    G3 = ox.graph_from_point(location_point, distance=500, distance_type='network')

    G4 = ox.graph_from_address(address='350 5th Ave, New York, NY', distance=1000, distance_type='network', network_type='bike')

    places = ['Los Altos, California, USA', {'city':'Los Altos Hills', 'state':'California'}, 'Loyola, California']
    G5 = ox.graph_from_place(places, network_type='all', clean_periphery=False)

    calif = gpd.read_file('tests/input_data/ZillowNeighborhoods-CA')
    mission_district = calif[(calif['CITY']=='San Francisco') & (calif['NAME']=='Mission')]
    polygon = mission_district['geometry'].iloc[0]
    G6 = ox.graph_from_polygon(polygon, network_type='walk')


def test_stats():

    location_point = (37.791427, -122.410018)
    G = ox.graph_from_point(location_point, distance=500, distance_type='network')
    G = ox.add_edge_bearings(G)
    G_proj = ox.project_graph(G)
    stats1 = ox.basic_stats(G)
    stats2 = ox.basic_stats(G, area=1000)
    stats3 = ox.basic_stats(G_proj, area=1000, clean_intersects=True, tolerance=15, circuity_dist='euclidean')
    stats4 = ox.extended_stats(G, connectivity=True, anc=True, ecc=True, bc=True, cc=True)


def test_plots():

    G = ox.graph_from_place('Piedmont, California, USA', network_type='drive', simplify=False)
    G2 = ox.simplify_graph(G, strict=False)
    nc = ox.get_node_colors_by_attr(G2, 'osmid')
    ec = ox.get_edge_colors_by_attr(G2, 'length')

    fig, ax = ox.plot_graph(G, save=True, file_format='png')

    G_simplified = ox.simplify_graph(G)
    fig, ax = ox.plot_graph(G_simplified, show=False, save=True, close=True, file_format='svg')

    G_projected = ox.project_graph(G_simplified)
    fig, ax = ox.plot_graph(G_projected)

    fig, ax = ox.plot_graph(G_projected, fig_height=5, fig_width=5, margin=0.05, axis_off=False, bgcolor='y',
                            file_format='png', filename='x', dpi=180, annotate=True, node_color='k', node_size=5,
                            node_alpha=0.1, node_edgecolor='b', node_zorder=5, edge_color='r', edge_linewidth=2,
                            edge_alpha=0.1, use_geom=False, show=False, save=True, close=True)

    fig, ax = ox.plot_figure_ground(G=G_simplified, file_format='png')
    fig, ax = ox.plot_figure_ground(point=(33.694981, -117.841375), file_format='png')
    fig, ax = ox.plot_figure_ground(address='Denver, Colorado, USA', file_format='png')


def test_routing_folium():

    import networkx as nx
    G = ox.graph_from_address('398 N. Sicily Pl., Chandler, Arizona', distance=800, network_type='drive')
    origin = (33.307792, -111.894940)
    destination = (33.312994, -111.894998)
    origin_node = ox.get_nearest_node(G, origin)
    destination_node = ox.get_nearest_node(G, destination)
    route = nx.shortest_path(G, origin_node, destination_node)

    attributes = ox.get_route_edge_attributes(G, route, 'length')

    fig, ax = ox.plot_graph_route(G, route, save=True, filename='route', file_format='png')
    fig, ax = ox.plot_graph_route(G, route, origin_point=origin, destination_point=destination,
                                  save=True, filename='route', file_format='png')

    graph_map = ox.plot_graph_folium(G, popup_attribute='name')
    route_map = ox.plot_route_folium(G, route)


def test_buildings():

    gdf = ox.buildings_from_place(place='Piedmont, California, USA')
    gdf = ox.buildings_from_address(address='260 Stockton Street, San Francisco, California, USA', distance=300)
    fig, ax = ox.plot_buildings(gdf)
