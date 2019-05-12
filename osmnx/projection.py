################################################################################
# Module: projection.py
# Description: Project spatial geometries and street networks to/from UTM
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

import time
import math
import numpy as np
import geopandas as gpd
import networkx as nx
from shapely.geometry import Point

from . import settings
from .utils import log
from .utils import make_str


def project_geometry(geometry, crs=None, to_crs=None, to_latlong=False):
    """
    Project a shapely Polygon or MultiPolygon from lat-long to UTM, or
    vice-versa

    Parameters
    ----------
    geometry : shapely Polygon or MultiPolygon
        the geometry to project
    crs : dict
        the starting coordinate reference system of the passed-in geometry,
        default value (None) will set settings.default_crs as the CRS
    to_crs : dict
        if not None, just project to this CRS instead of to UTM
    to_latlong : bool
        if True, project from crs to lat-long, if False, project from crs to
        local UTM zone

    Returns
    -------
    tuple
        (geometry_proj, crs), the projected shapely geometry and the crs of the
        projected geometry
    """

    if crs is None:
        crs = settings.default_crs

    gdf = gpd.GeoDataFrame()
    gdf.crs = crs
    gdf.gdf_name = 'geometry to project'
    gdf['geometry'] = None
    gdf.loc[0, 'geometry'] = geometry
    gdf_proj = project_gdf(gdf, to_crs=to_crs, to_latlong=to_latlong)
    geometry_proj = gdf_proj['geometry'].iloc[0]
    return geometry_proj, gdf_proj.crs


def project_gdf(gdf, to_crs=None, to_latlong=False):
    """
    Project a GeoDataFrame to the UTM zone appropriate for its geometries'
    centroid.

    The simple calculation in this function works well for most latitudes, but
    won't work for some far northern locations like Svalbard and parts of far
    northern Norway.

    Parameters
    ----------
    gdf : GeoDataFrame
        the gdf to be projected
    to_crs : dict
        if not None, just project to this CRS instead of to UTM
    to_latlong : bool
        if True, projects to latlong instead of to UTM

    Returns
    -------
    GeoDataFrame
    """
    assert len(gdf) > 0, 'You cannot project an empty GeoDataFrame.'
    start_time = time.time()

    # if gdf has no gdf_name attribute, create one now
    if not hasattr(gdf, 'gdf_name'):
        gdf.gdf_name = 'unnamed'

    # if to_crs was passed-in, use this value to project the gdf
    if to_crs is not None:
        projected_gdf = gdf.to_crs(to_crs)

    # if to_crs was not passed-in, calculate the centroid of the geometry to
    # determine UTM zone
    else:
        if to_latlong:
            # if to_latlong is True, project the gdf to latlong
            latlong_crs = settings.default_crs
            projected_gdf = gdf.to_crs(latlong_crs)
            log('Projected the GeoDataFrame "{}" to default_crs in {:,.2f} seconds'.format(gdf.gdf_name, time.time()-start_time))
        else:
            # else, project the gdf to UTM
            # if GeoDataFrame is already in UTM, just return it
            if (gdf.crs is not None) and ('+proj=utm ' in gdf.crs):
                return gdf

            # calculate the centroid of the union of all the geometries in the
            # GeoDataFrame
            avg_longitude = gdf['geometry'].unary_union.centroid.x

            # calculate the UTM zone from this avg longitude and define the UTM
            # CRS to project
            utm_zone = int(math.floor((avg_longitude + 180) / 6.) + 1)
            utm_crs = '+proj=utm +zone={} +ellps=WGS84 +datum=WGS84 +units=m +no_defs'.format(utm_zone)

            # project the GeoDataFrame to the UTM CRS
            projected_gdf = gdf.to_crs(utm_crs)
            log('Projected the GeoDataFrame "{}" to UTM-{} in {:,.2f} seconds'.format(gdf.gdf_name, utm_zone, time.time()-start_time))

    projected_gdf.gdf_name = gdf.gdf_name
    return projected_gdf


def project_graph(G, to_crs=None):
    """
    Project a graph from lat-long to the UTM zone appropriate for its geographic
    location.

    Parameters
    ----------
    G : networkx multidigraph
        the networkx graph to be projected
    to_crs : dict
        if not None, just project to this CRS instead of to UTM

    Returns
    -------
    networkx multidigraph
    """

    G_proj = G.copy()
    start_time = time.time()

    # create a GeoDataFrame of the nodes, name it, convert osmid to str
    nodes, data = zip(*G_proj.nodes(data=True))
    gdf_nodes = gpd.GeoDataFrame(list(data), index=nodes)
    gdf_nodes.crs = G_proj.graph['crs']
    gdf_nodes.gdf_name = '{}_nodes'.format(G_proj.name)

    # create new lat/lon columns just to save that data for later, and create a
    # geometry column from x/y
    gdf_nodes['lon'] = gdf_nodes['x']
    gdf_nodes['lat'] = gdf_nodes['y']
    gdf_nodes['geometry'] = gdf_nodes.apply(lambda row: Point(row['x'], row['y']), axis=1)
    log('Created a GeoDataFrame from graph in {:,.2f} seconds'.format(time.time()-start_time))

    # project the nodes GeoDataFrame to UTM
    gdf_nodes_utm = project_gdf(gdf_nodes, to_crs=to_crs)

    # extract data for all edges that have geometry attribute
    edges_with_geom = []
    for u, v, key, data in G_proj.edges(keys=True, data=True):
        if 'geometry' in data:
            edges_with_geom.append({'u':u, 'v':v, 'key':key, 'geometry':data['geometry']})

    # create an edges GeoDataFrame and project to UTM, if there were any edges
    # with a geometry attribute. geom attr only exists if graph has been
    # simplified, otherwise you don't have to project anything for the edges
    # because the nodes still contain all spatial data
    if len(edges_with_geom) > 0:
        gdf_edges = gpd.GeoDataFrame(edges_with_geom)
        gdf_edges.crs = G_proj.graph['crs']
        gdf_edges.gdf_name = '{}_edges'.format(G_proj.name)
        gdf_edges_utm = project_gdf(gdf_edges, to_crs=to_crs)

    # extract projected x and y values from the nodes' geometry column
    start_time = time.time()
    gdf_nodes_utm['x'] = gdf_nodes_utm['geometry'].map(lambda point: point.x)
    gdf_nodes_utm['y'] = gdf_nodes_utm['geometry'].map(lambda point: point.y)
    gdf_nodes_utm = gdf_nodes_utm.drop('geometry', axis=1)
    log('Extracted projected node geometries from GeoDataFrame in {:,.2f} seconds'.format(time.time()-start_time))

    # clear the graph to make it a blank slate for the projected data
    start_time = time.time()
    edges = list(G_proj.edges(keys=True, data=True))
    graph_name = G_proj.graph['name']
    G_proj.clear()

    # add the projected nodes and all their attributes to the graph
    G_proj.add_nodes_from(gdf_nodes_utm.index)
    attributes = gdf_nodes_utm.to_dict()
    for label in gdf_nodes_utm.columns:
        nx.set_node_attributes(G_proj, name=label, values=attributes[label])

    # add the edges and all their attributes (including reconstructed geometry,
    # when it exists) to the graph
    for u, v, key, attributes in edges:
        if 'geometry' in attributes:
            row = gdf_edges_utm[(gdf_edges_utm['u']==u) & (gdf_edges_utm['v']==v) & (gdf_edges_utm['key']==key)]
            attributes['geometry'] = row['geometry'].iloc[0]

        # attributes dict contains key, so we don't need to explicitly pass it here
        G_proj.add_edge(u, v, **attributes)

    # set the graph's CRS attribute to the new, projected CRS and return the
    # projected graph
    G_proj.graph['crs'] = gdf_nodes_utm.crs
    G_proj.graph['name'] = '{}_UTM'.format(graph_name)
    if 'streets_per_node' in G.graph:
        G_proj.graph['streets_per_node'] = G.graph['streets_per_node']
    log('Rebuilt projected graph in {:,.2f} seconds'.format(time.time()-start_time))
    return G_proj
