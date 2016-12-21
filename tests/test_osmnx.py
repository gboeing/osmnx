"""
OSMnx tests
-------
"""

import osmnx as ox
ox.config(log_console=True, log_file=True, use_cache=True, 
          data_folder='.temp/data', logs_folder='.temp/logs', imgs_folder='.temp/imgs', cache_folder='.temp/cache')
          
def test_imports():
    import json, math, sys, os, io, ast, unicodedata, hashlib, re, random, time, datetime as dt, logging as lg
    from collections import OrderedDict, Counter
    from itertools import groupby, chain
    from dateutil import parser as date_parser
    import requests, numpy as np, pandas as pd, geopandas as gpd, networkx as nx, matplotlib.pyplot as plt, matplotlib.cm as cm
    from matplotlib.collections import LineCollection
    from shapely.geometry import Point, LineString, Polygon, MultiPolygon
    from shapely import wkt
    from shapely.ops import unary_union
    from descartes import PolygonPatch
    from geopy.distance import great_circle, vincenty
    from geopy.geocoders import Nominatim
    from rtree.index import Index as RTreeIndex
    
def test_gdf_shapefiles():
    city = ox.gdf_from_place('Manhattan, New York City, New York, USA')
    city_projected = ox.project_gdf(city)
    ox.save_gdf_shapefile(city_projected)
    
def test_network_saving_loading():
    G = ox.graph_from_place('Piedmont, California, USA')
    G_projected = ox.project_graph(G)
    ox.save_graph_shapefile(G_projected)
    ox.save_graphml(G_projected)
    G2 = ox.load_graphml('graph.graphml')
    
def test_get_network_methods():
    
    import geopandas as gpd
    north, south, east, west = 37.79, 37.78, -122.41, -122.43
    G1 = ox.graph_from_bbox(north, south, east, west, network_type='drive_service')

    location_point = (37.791427, -122.410018)
    G2 = ox.graph_from_point(location_point, distance=750, distance_type='bbox', network_type='drive')
    G3 = ox.graph_from_point(location_point, distance=500, distance_type='network')

    G4 = ox.graph_from_address(address='350 5th Ave, New York, NY', distance=1000, distance_type='network', network_type='bike')

    places = ['Los Altos, California, USA', {'city':'Los Altos Hills', 'state':'California'}, 'Loyola, California']
    G5 = ox.graph_from_place(places, network_type='all', clean_periphery=False)

    calif = gpd.read_file('examples/input_data/ZillowNeighborhoods-CA')
    mission_district = calif[(calif['CITY']=='San Francisco') & (calif['NAME']=='Mission')]
    polygon = mission_district['geometry'].iloc[0]
    G6 = ox.graph_from_polygon(polygon, network_type='walk')
    
def test_stats():
    location_point = (37.791427, -122.410018)
    G = ox.graph_from_point(location_point, distance=500, distance_type='network')
    stats1 = ox.basic_stats(G)
    stats2 = ox.extended_stats(G, connectivity=True, anc=True, ecc=True, bc=True, cc=True)
    
def test_plots():
    
    G = ox.graph_from_place('Piedmont, California, USA', network_type='drive', simplify=False)
    fig, ax = ox.plot_graph(G, show=False, save=True, close=True)

    G_simplified = ox.simplify_graph(G)
    fig, ax = ox.plot_graph(G_simplified, show=False, save=True, close=True)

    G_projected = ox.project_graph(G_simplified)
    fig, ax = ox.plot_graph(G_projected, show=False, save=True, close=True)

    fig, ax = ox.plot_graph(G_projected, fig_height=5, fig_width=5, margin=0.05, axis_off=False, bgcolor='y',
                            file_format='png', filename='x', dpi=180, annotate=True, node_color='k', node_size=5,
                            node_alpha=0.1, node_edgecolor='b', node_zorder=5, edge_color='r', edge_linewidth=2,
                            edge_alpha=0.1, use_geom=False, show=False, save=True, close=True)
                            
def test_routing():
    import networkx as nx
    G = ox.graph_from_address('N. Sicily Pl., Chandler, Arizona', distance=800, network_type='drive')
    origin = (33.307792, -111.894940)
    destination = (33.312994, -111.894998)
    origin_node = ox.get_nearest_node(G, origin)
    destination_node = ox.get_nearest_node(G, destination)
    route = nx.shortest_path(G, origin_node, destination_node)
    fig, ax = ox.plot_graph_route(G, route, save=True, filename='route', show=False, close=True)
    
    
