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

import httmock, gzip, io

try:
    from urllib.parse import parse_qsl
except ImportError:
    from urlparse import parse_qsl


def get_mock_response_content(overpass_filename=None):
    ''' Return a mock HTTP handler suitable for use by httmock.HTTMock().
    '''
    tests_dir = os.path.join(os.path.dirname(__file__), 'input_data')
    Nominatim_Searches = {}

    # Load all Nominatim searches into a dictionary
    with io.open(os.path.join(tests_dir, 'nominatim-responses.txt'), 'r', encoding='utf8') as file:
        while True:
            try:
                key, value, _ = next(file).strip(), next(file), next(file)
            except StopIteration:
                break
            else:
                Nominatim_Searches[key] = value

    def response_content(url, request):
        # Return with a stock Overpass status that doesn't trigger a wait
        if (url.netloc, url.path) == ('overpass-api.de', '/api/status'):
            with open(os.path.join(tests_dir, 'overpass-response-0.txt'), 'rb') as file:
                return httmock.response(200, file.read())

        # For larger Overpass responses, look to the given filename for data
        if (url.netloc, url.path) == ('overpass-api.de', '/api/interpreter'):
            with gzip.open(os.path.join(tests_dir, overpass_filename), 'r') as file:
                return httmock.response(200, file.read())

        # Look for Nominatim responses in the dictionary loaded earlier
        query_dict = dict(parse_qsl(url.query)) if url.query else {}
        body_dict = dict(parse_qsl(request.body)) if request.body else {}

        if (url.netloc, url.path) == ('nominatim.openstreetmap.org', '/search'):
            key1 = query_dict.get('q') #try this method for any single string queries
            key2 = '{c}, {s}'.format(c=query_dict.get('city'), s=query_dict.get('state')) #try this method for any dict queries
            for key in (key1, key2):
                if key in Nominatim_Searches:
                    return httmock.response(200, Nominatim_Searches[key])

            raise ValueError('Unknown Nominatim query: {}'.format(url.query))

        raise ValueError('Unknown URL: {}'.format(url))

    return response_content


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

    with httmock.HTTMock(get_mock_response_content()):
        city = ox.gdf_from_place('Manhattan, New York City, New York, USA')
    city_projected = ox.project_gdf(city, to_crs={'init':'epsg:3395'})
    ox.save_gdf_shapefile(city_projected)

    with httmock.HTTMock(get_mock_response_content()):
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

    with httmock.HTTMock(get_mock_response_content('overpass-response-7.json.gz')):
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
    with httmock.HTTMock(get_mock_response_content('overpass-response-6.json.gz')):
        G1 = ox.graph_from_bbox(north, south, east, west, network_type='drive_service')

    with httmock.HTTMock(get_mock_response_content()):
        G1 = ox.graph_from_bbox(north, south, east, west, network_type='drive_service', truncate_by_edge=True)

    location_point = (37.791427, -122.410018)
    bbox = ox.bbox_from_point(location_point, project_utm=True)
    with httmock.HTTMock(get_mock_response_content('overpass-response-8.json.gz')):
        G2 = ox.graph_from_point(location_point, distance=750, distance_type='bbox', network_type='drive')

    with httmock.HTTMock(get_mock_response_content('overpass-response-9.json.gz')):
        G3 = ox.graph_from_point(location_point, distance=500, distance_type='network')

    with httmock.HTTMock(get_mock_response_content('overpass-response-10.json.gz')):
        G4 = ox.graph_from_address(address='350 5th Ave, New York, NY', distance=1000, distance_type='network', network_type='bike')

    places = ['Los Altos, California, USA', {'city':'Los Altos Hills', 'state':'California'}, 'Loyola, California']
    with httmock.HTTMock(get_mock_response_content('overpass-response-11.json.gz')):
        G5 = ox.graph_from_place(places, network_type='all', clean_periphery=False)

    calif = gpd.read_file('tests/input_data/ZillowNeighborhoods-CA')
    mission_district = calif[(calif['CITY']=='San Francisco') & (calif['NAME']=='Mission')]
    polygon = mission_district['geometry'].iloc[0]
    with httmock.HTTMock(get_mock_response_content('overpass-response-12.json.gz')):
        G6 = ox.graph_from_polygon(polygon, network_type='walk')


def test_stats():

    location_point = (37.791427, -122.410018)
    with httmock.HTTMock(get_mock_response_content('overpass-response-5.json.gz')):
        G = ox.graph_from_point(location_point, distance=500, distance_type='network')
    G = ox.add_edge_bearings(G)
    G_proj = ox.project_graph(G)
    stats1 = ox.basic_stats(G)
    stats2 = ox.basic_stats(G, area=1000)
    stats3 = ox.basic_stats(G_proj, area=1000, clean_intersects=True, tolerance=15, circuity_dist='euclidean')
    stats4 = ox.extended_stats(G, connectivity=True, anc=True, ecc=True, bc=True, cc=True)


def test_plots():

    with httmock.HTTMock(get_mock_response_content('overpass-response-4.json.gz')):
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
    with httmock.HTTMock(get_mock_response_content('overpass-response-3.json.gz')):
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

    with httmock.HTTMock(get_mock_response_content('overpass-response-1.json.gz')):
        gdf = ox.buildings_from_place(place='Piedmont, California, USA')

    with httmock.HTTMock(get_mock_response_content('overpass-response-2.json.gz')):
        gdf = ox.buildings_from_address(address='260 Stockton Street, San Francisco, California, USA', distance=300)

    fig, ax = ox.plot_buildings(gdf)
