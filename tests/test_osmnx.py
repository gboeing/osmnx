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

def test_get_network():
    G = ox.graph_from_place('Piedmont, California, USA')
    G_projected = ox.project_graph(G)
    ox.save_graph_shapefile(G_projected)
    
