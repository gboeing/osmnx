################################################################################
# test_osmnx.py
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

import matplotlib as mpl
mpl.use('Agg') #use agg backend so you don't need a display on travis-ci

# remove the .temp folder if it already exists so we start fresh with tests
import os, shutil
if os.path.exists('.temp'):
    shutil.rmtree('.temp')

import osmnx as ox

# configure OSMnx
ox.config(log_console=True, log_file=True, use_cache=True,
          data_folder='.temp/data', logs_folder='.temp/logs',
          imgs_folder='.temp/imgs', cache_folder='.temp/cache')


def test_imports():

    # test all of OSMnx's module imports
    import ast
    import datetime
    import geopandas
    import hashlib
    import io
    import json
    import logging
    import math
    import matplotlib.cm
    import matplotlib.pyplot
    import networkx
    import numpy
    import os
    import pandas
    import random
    import requests
    import re
    import sys
    import time
    import unicodedata
    import warnings
    from collections import Counter
    from collections import OrderedDict
    from dateutil import parser
    from descartes import PolygonPatch
    from itertools import chain
    from itertools import groupby
    from matplotlib.collections import LineCollection
    from rtree.index import Index
    from shapely import wkt
    from shapely.geometry import Point
    from shapely.geometry import MultiPoint
    from shapely.geometry import LineString
    from shapely.geometry import MultiLineString
    from shapely.geometry import Polygon
    from shapely.geometry import MultiPolygon
    from shapely.ops import unary_union


def test_logging():

    # test OSMnx's logger
    import logging as lg
    ox.log('test a fake debug', level=lg.DEBUG)
    ox.log('test a fake info', level=lg.INFO)
    ox.log('test a fake warning', level=lg.WARNING)
    ox.log('test a fake error', level=lg.ERROR)

    ox.citation()


def test_geometry_coords_rounding():

    # test the rounding of geometry coordinates
    from shapely.geometry import Point, MultiPoint, LineString, MultiLineString, Polygon, MultiPolygon

    precision = 3

    shape1 = Point(1.123456, 2.123456)
    shape2 = ox.round_shape_coords(shape1, precision)

    shape1 = MultiPoint([(1.123456, 2.123456), (3.123456, 4.123456)])
    shape2 = ox.round_shape_coords(shape1, precision)

    shape1 = LineString([(1.123456, 2.123456), (3.123456, 4.123456)])
    shape2 = ox.round_shape_coords(shape1, precision)

    shape1 = MultiLineString([[(1.123456, 2.123456), (3.123456, 4.123456)],
                              [(11.123456, 12.123456), (13.123456, 14.123456)]])
    shape2 = ox.round_shape_coords(shape1, precision)

    shape1 = Polygon([(1.123456, 2.123456), (3.123456, 4.123456), (6.123456, 5.123456)])
    shape2 = ox.round_shape_coords(shape1, precision)

    shape1 = MultiPolygon([Polygon([(1.123456, 2.123456), (3.123456, 4.123456), (6.123456, 5.123456)]),
                           Polygon([(16.123456, 15.123456), (13.123456, 14.123456), (12.123456, 11.123456)])])
    shape2 = ox.round_shape_coords(shape1, precision)


def test_gdf_shapefiles():

    # test loading spatial boundaries, saving as shapefile, and plotting
    city = ox.gdf_from_place('Manhattan, New York City, New York, USA')
    city_projected = ox.project_gdf(city, to_crs={'init':'epsg:3395'})
    ox.save_gdf_shapefile(city_projected)

    city = ox.gdf_from_place('Manhattan, New York City, New York, USA', buffer_dist=100)
    ox.plot_shape(city)


def test_graph_from_file():

    # test loading a graph from a local .osm file
    import bz2, tempfile

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

    # save/load graph as shapefile and graphml file
    G = ox.graph_from_place('Piedmont, California, USA')
    G_projected = ox.project_graph(G)
    ox.save_graph_shapefile(G_projected)
    ox.save_graphml(G_projected)
    ox.save_graphml(G_projected, filename='gephi.graphml', gephi=True)
    G2 = ox.load_graphml('graph.graphml')
    G3 = ox.load_graphml('graph.graphml', node_type=str)

    # convert graph to node/edge GeoDataFrames and back again
    gdf_edges = ox.graph_to_gdfs(G, nodes=False, edges=True, fill_edge_geometry=False)
    gdf_nodes, gdf_edges = ox.graph_to_gdfs(G, nodes=True, edges=True, node_geometry=True, fill_edge_geometry=True)
    G4 = ox.gdfs_to_graph(gdf_nodes, gdf_edges)

    # find graph nodes nearest to some set of points
    X = gdf_nodes['x'].head()
    Y = gdf_nodes['y'].head()
    nn1 = ox.get_nearest_nodes(G, X, Y)
    nn2 = ox.get_nearest_nodes(G, X, Y, method='kdtree')
    nn3 = ox.get_nearest_nodes(G, X, Y, method='balltree')

    # find graph edges nearest to some set of points
    ne1 = ox.get_nearest_edges(G, X, Y)
    ne2 = ox.get_nearest_edges(G, X, Y, method='kdtree')
    ne3 = ox.get_nearest_edges(G, X, Y, method='kdtree', dist=50)


def test_get_network_methods():

    from shapely import wkt

    # graph from bounding box
    north, south, east, west = 37.79, 37.78, -122.41, -122.43
    G1 = ox.graph_from_bbox(north, south, east, west, network_type='drive_service')
    G1 = ox.graph_from_bbox(north, south, east, west, network_type='drive_service', truncate_by_edge=True)

    # graph from point
    location_point = (37.791427, -122.410018)
    bbox = ox.bbox_from_point(location_point, project_utm=True)
    G2 = ox.graph_from_point(location_point, distance=750, distance_type='bbox', network_type='drive')
    G3 = ox.graph_from_point(location_point, distance=500, distance_type='network')

    # graph from address
    G4 = ox.graph_from_address(address='350 5th Ave, New York, NY', distance=1000, distance_type='network', network_type='bike')

    # graph from list of places
    places = ['Los Altos, California, USA', {'city':'Los Altos Hills', 'state':'California'}, 'Loyola, California']
    G5 = ox.graph_from_place(places, network_type='all', clean_periphery=False)

    # graph from polygon
    polygon = wkt.loads('POLYGON ((-122.418083 37.754154, -122.418082 37.766028, -122.410909 37.766028, -122.410908 37.754154, -122.418083 37.754154))')
    G6 = ox.graph_from_polygon(polygon, network_type='walk')

    # test custom query filter
    filtr = ('["area"!~"yes"]'
             '["highway"!~"motor|proposed|construction|abandoned|platform|raceway"]'
             '["foot"!~"no"]'
             '["service"!~"private"]'
             '["access"!~"private"]')
    G = ox.graph_from_point(location_point, network_type='walk', custom_filter=filtr)


def test_stats():

    # create graph, add bearings, project it
    location_point = (37.791427, -122.410018)
    G = ox.graph_from_point(location_point, distance=500, distance_type='network')
    G = ox.add_edge_bearings(G)
    G_proj = ox.project_graph(G)

    # calculate stats
    stats1 = ox.basic_stats(G)
    stats2 = ox.basic_stats(G, area=1000)
    stats3 = ox.basic_stats(G_proj, area=1000, clean_intersects=True, tolerance=15, circuity_dist='euclidean')
    stats4 = ox.extended_stats(G, connectivity=True, anc=True, ecc=True, bc=True, cc=True)


def test_plots():

    G = ox.graph_from_place('Piedmont, California, USA', network_type='drive', simplify=False)
    G2 = ox.simplify_graph(G, strict=False)

    # test getting colors
    co = ox.get_colors(n=5, return_hex=True)
    nc = ox.get_node_colors_by_attr(G2, 'osmid')
    ec = ox.get_edge_colors_by_attr(G2, 'length')

    # save a plot to disk as png
    fig, ax = ox.plot_graph(G, save=True, file_format='png')

    # save a plot to disk as svg
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

    # calculate shortest path and plot as static image and leaflet web map
    import networkx as nx
    G = ox.graph_from_address('398 N. Sicily Pl., Chandler, Arizona', distance=800, network_type='drive')
    origin = (33.307792, -111.894940)
    destination = (33.312994, -111.894998)
    origin_node = ox.get_nearest_node(G, origin)
    destination_node = ox.get_nearest_node(G, destination, method='euclidean')
    route = nx.shortest_path(G, origin_node, destination_node)

    attributes = ox.get_route_edge_attributes(G, route, 'length')

    fig, ax = ox.plot_graph_route(G, route, save=True, filename='route', file_format='png')
    fig, ax = ox.plot_graph_route(G, route, origin_point=origin, destination_point=destination,
                                  save=True, filename='route', file_format='png')

    # test multiple routes
    fig, ax = ox.plot_graph_routes(G, [route, route])

    graph_map = ox.plot_graph_folium(G, popup_attribute='name')
    route_map = ox.plot_route_folium(G, route)


def test_nearest_edge():

    # test in closest edge section
    sheik_sayed_dubai = [25.09, 25.06, 55.16, 55.11]
    location_coordinates = (25.071764, 55.138978)
    G = ox.graph_from_bbox(*sheik_sayed_dubai, simplify=False, retain_all=True, network_type='drive')
    geometry, u, v = ox.get_nearest_edge(G, location_coordinates)


def test_nearest_edges():

    from pyproj import Proj

    # test in closest edge section
    sheik_sayed_dubai = [25.09, 25.06, 55.16, 55.11]
    location_coordinates = (25.071764, 55.138978)
    G = ox.graph_from_bbox(*sheik_sayed_dubai, simplify=False, retain_all=True, network_type='drive')

    # Unprojected
    ne1 = ox.get_nearest_edges(G, X=[location_coordinates[1], location_coordinates[1]],
                                  Y=[location_coordinates[0], location_coordinates[0]], method='balltree', dist=0.0001)

    # Projected
    G2 = ox.project_graph(G)
    crs = Proj(G2.graph['crs'])

    projected_point = crs(location_coordinates[1], location_coordinates[0])
    ne2 = ox.get_nearest_edges(G2, X=[projected_point[0], projected_point[0]],
                                   Y=[projected_point[1], projected_point[1]], method='kdtree', dist=10)
    assert (ne1 == ne2).all()


def test_footprints():
    import json
    import pytest
    from shapely.geometry import Polygon

    # download footprints and plot them
    gdf = ox.footprints_from_place(place='Emeryville, California, USA')
    gdf = ox.footprints_from_polygon(Polygon([(17.574, -4.145), (17.575, -4.144), (17.576, -4.145)]))
    gdf = ox.footprints_from_address(address='600 Montgomery St, San Francisco, California, USA', distance=300)
    fig, ax = ox.plot_footprints(gdf)

    # new_river_head.json contains a relation with 1 outer closed way and 2 inner closed ways
    # inner way 665593284 is directly tagged as a building and should create its own polygon
    with open("tests/input_data/new_river_head.json", "r") as read_file:
        new_river_head_responses = [json.load(read_file)]
    new_river_head_gdf = ox.create_footprints_gdf(responses=new_river_head_responses)
    assert 665593284 in new_river_head_gdf.index
    assert new_river_head_gdf.loc[9246394]['geometry'].type=='Polygon'
    assert len(new_river_head_gdf.loc[9246394,'geometry'].interiors)==2

    # clapham_common.json contains a relation with 5 outer rings and 1 inner ring. One of the outer rings is a chain of open ways
    with open("tests/input_data/clapham_common.json", "r") as read_file:
        clapham_common_responses = [json.load(read_file)]
    clapham_common_gdf = ox.create_footprints_gdf(footprint_type='leisure', responses=clapham_common_responses)
    assert clapham_common_gdf.loc[1290065]['geometry'].type=='MultiPolygon'

    # relation_no_outer.json contains a relation with 0 outer rings and 1 inner ring
    with open("tests/input_data/relation_no_outer.json", "r") as read_file:
        relation_no_outer_responses = [json.load(read_file)]
    ox.create_footprints_gdf(responses=relation_no_outer_responses)

    # inner_chain.json contains a relation with 1 outer rings and several inner rings one of which is a chain of open ways
    with open("tests/input_data/inner_chain.json", "r") as read_file:
        inner_chain_responses = [json.load(read_file)]
    ox.create_footprints_gdf(responses=inner_chain_responses)

    # mis_tagged_bus_route.json contains a relation with out 'inner' or 'inner' rings
    with open("tests/input_data/mis_tagged_bus_route.json", "r") as read_file:
        mis_tagged_bus_route_responses = [json.load(read_file)]
    ox.create_footprints_gdf(responses=mis_tagged_bus_route_responses)

    # test plotting multipolygon
    fig, ax = ox.plot_footprints(clapham_common_gdf)

    # should raise an exception
    # polygon or -north, south, east, west- should be provided
    with pytest.raises(ValueError):
        ox.create_footprints_gdf(polygon=None, north=None, south=None, east=None, west=None)

    gdf = ox.footprints_from_place(place='kusatsu, shiga, japan', which_result=2)

def test_pois():

    import pytest
    # download all points of interests from place
    gdf = ox.pois_from_place(place='Kamppi, Helsinki, Finland')

    # get all restaurants and schools from place
    restaurants = ox.pois_from_place(place='Emeryville, California, USA', amenities=['restaurant'])
    schools = ox.pois_from_place(place='Emeryville, California, USA', amenities=['school'])

    # get all universities from Boston area (with 2 km buffer to cover also Cambridge)
    boston_q = "Boston, Massachusetts, United States of America"
    boston_poly = ox.gdf_from_place(boston_q, buffer_dist=2000)
    universities = ox.pois_from_polygon(boston_poly.geometry.values[0], amenities=['university'])

    # by point and by address
    restaurants = ox.pois_from_point(point=(42.344490, -71.070570), distance=1000, amenities=['restaurant'])
    restaurants = ox.pois_from_address(address='Emeryville, California, USA', distance=1000, amenities=['restaurant'])

    # should raise an exception
    # polygon or -north, south, east, west- should be provided
    with pytest.raises(ValueError):
        ox.create_poi_gdf(polygon=None, north=None, south=None, east=None, west=None)

    gdf = ox.pois_from_place(place='kusatsu, shiga, japan', which_result=2)


def test_nominatim():

    import pytest
    from collections import OrderedDict

    params = OrderedDict()
    params['format'] = "json"
    params['address_details'] = 0

    # Bad Address - should return an empty response
    params['q'] = "AAAAAAAAAAA"
    response_json = ox.nominatim_request(params = params,
                                         type = "search")

    # Good Address - should return a valid response with a valid osm_id
    params['q'] = "Newcastle A186 Westgate Rd"
    response_json = ox.nominatim_request(params = params,
                                         type = "search")

    # Lookup
    params = OrderedDict()
    params['format'] = "json"
    params['address_details'] = 0
    params['osm_ids'] = "W68876073"

    response_json = ox.nominatim_request(params = params,
                                         type = "lookup")

    # Invalid nominatim query type
    with pytest.raises(ValueError):
        response_json = ox.nominatim_request(
                            params = params,
                            type = "transfer")

    # Searching on public nominatim should work even if a key was provided
    ox.config(
        nominatim_key="NOT_A_KEY"
    )
    response_json = ox.nominatim_request(params = params,
                                         type = "search")

    # Test changing the endpoint. It should fail because we didn't provide a valid key
    ox.config(
        nominatim_endpoint="http://open.mapquestapi.com/nominatim/v1/"
    )
    with pytest.raises(Exception):
        response_json = ox.nominatim_request(params=params,
                                             type="search")

    ox.config(log_console=True, log_file=True, use_cache=True,
              data_folder='.temp/data', logs_folder='.temp/logs',
              imgs_folder='.temp/imgs', cache_folder='.temp/cache')


def test_osm_xml_output():
    G = ox.graph_from_place('Piedmont, California, USA')
    ox.save_graph_osm(G)
