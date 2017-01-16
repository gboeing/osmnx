import re
import time
import os
import ast
import numpy as np
import pandas as pd
import geopandas as gpd
import networkx as nx
from shapely.geometry import Point, LineString
from shapely import wkt
from . import globals
from .utils import log, make_str, get_undirected


def save_graph_shapefile(G, filename='graph', folder=None, encoding='utf-8'):
    """
    Save graph nodes and edges as ESRI shapefiles to disk.
    
    Parameters
    ----------
    G : graph
    filename : string, the name of the shapefiles (not including file extensions)
    folder : string, the folder to contain the shapefiles, if None, use default data folder
    encoding : string, the character encoding for the saved shapefiles
    
    Returns
    -------
    None
    """
    
    start_time = time.time()
    if folder is None:
        folder = globals.global_data_folder
    
    # convert directed graph G to an undirected graph for saving as a shapefile
    G_save = G.copy()
    G_save = get_undirected(G_save)
    
    # create a GeoDataFrame of the nodes and set CRS
    nodes = {node:data for node, data in G_save.nodes(data=True)}
    gdf_nodes = gpd.GeoDataFrame(nodes).T
    gdf_nodes.crs = G_save.graph['crs']
    
    # create a geometry column then drop the x and y columns
    gdf_nodes['geometry'] = gdf_nodes.apply(lambda row: Point(row['x'], row['y']), axis=1)
    gdf_nodes = gdf_nodes.drop(['x', 'y'], axis=1)

    # make everything but geometry column a string
    gdf_nodes['osmid'] = gdf_nodes['osmid'].astype(np.int64)
    for col in [c for c in gdf_nodes.columns if not c == 'geometry']:
        gdf_nodes[col] = gdf_nodes[col].fillna('').map(make_str)
        
    # create a list to hold our edges, then loop through each edge in the graph
    edges = []
    for u, v, key, data in G_save.edges(keys=True, data=True):

        # for each edge, add key and all attributes in data dict to the edge_details
        edge_details = {'key':key}
        for attr_key in data:
            edge_details[attr_key] = data[attr_key]

        # if edge doesn't already have a geometry attribute, create one now
        if not 'geometry' in data:
            point_u = Point((G_save.node[u]['x'], G_save.node[u]['y']))
            point_v = Point((G_save.node[v]['x'], G_save.node[v]['y']))
            edge_details['geometry'] = LineString([point_u, point_v])
        
        edges.append(edge_details)
    
    # create a geodataframe from the list of edges and set the CRS
    gdf_edges = gpd.GeoDataFrame(edges)
    gdf_edges.crs = G_save.graph['crs']

    # make everything but geometry column a string
    for col in [c for c in gdf_edges.columns if not c == 'geometry']:
        gdf_edges[col] = gdf_edges[col].fillna('').map(make_str)
        
    # if the save folder does not already exist, create it with a filename subfolder
    folder = '{}/{}'.format(folder, filename)
    if not os.path.exists(folder):
        os.makedirs(folder)
    
    # save the nodes and edges as separate ESRI shapefiles
    gdf_nodes.to_file('{}/nodes'.format(folder), encoding=encoding)
    gdf_edges.to_file('{}/edges'.format(folder), encoding=encoding)
    log('Saved graph "{}" to disk as shapefiles at "{}" in {:,.2f} seconds'.format(G_save.name, folder, time.time()-start_time))
    
    
def save_graphml(G, filename='graph.graphml', folder=None):
    """
    Save graph as GraphML file to disk.
    
    Parameters
    ----------
    G : graph
    filename : string, the name of the graphml file (including file extension)
    folder : string, the folder to contain the file, if None, use default data folder
    
    Returns
    -------
    None
    """
    
    start_time = time.time()
    if folder is None:
        folder = globals.global_data_folder
    
    # create a copy and convert all the node/edge attribute values to string or it won't save
    G_save = G.copy()
    for dict_key in G_save.graph:
        # convert all the graph attribute values to strings
        G_save.graph[dict_key] = make_str(G_save.graph[dict_key])
    for node, data in G_save.nodes(data=True):
        for dict_key in data:
            # convert all the node attribute values to strings
            data[dict_key] = make_str(data[dict_key])
    for u, v, key, data in G_save.edges(keys=True, data=True):
        for dict_key in data:
            # convert all the edge attribute values to strings
            data[dict_key] = make_str(data[dict_key])
                
    if not os.path.exists(folder):
        os.makedirs(folder)
    
    nx.write_graphml(G_save, '{}/{}'.format(folder, filename))
    log('Saved graph "{}" to disk as GraphML at "{}/{}" in {:,.2f} seconds'.format(G_save.name, folder, filename, time.time()-start_time))
    
    
def load_graphml(filename, folder=None):
    """
    Load a GraphML file from disk and convert the node/edge attributes to correct data types.
    
    Parameters
    ----------
    filename : string, the name of the graphml file (including file extension)
    folder : string, the folder containing the file, if None, use default data folder
    
    Returns
    -------
    G : graph
    """
    start_time = time.time()
    
    # read the graph from disk
    if folder is None:
        folder = globals.global_data_folder
    path = '{}/{}'.format(folder, filename)
    G = nx.MultiDiGraph(nx.read_graphml(path, node_type=int))
    
    # convert graph crs attribute from saved string to correct dict data type
    G.graph['crs'] = ast.literal_eval(G.graph['crs'])
     
    if 'streets_per_node' in G.graph:
        G.graph['streets_per_node'] = ast.literal_eval(G.graph['streets_per_node'])
        
    # convert numeric node tags from string to numeric data types
    log('Converting node and edge attribute data types')
    for node, data in G.nodes(data=True):
        data['osmid'] = int(data['osmid'])
        data['x'] = float(data['x'])
        data['y'] = float(data['y'])
    
    # convert numeric, bool, and list node tags from string to correct data types
    for u, v, key, data in G.edges(keys=True, data=True):

        # first parse oneway to bool and length to float - they should always have only 1 value each
        data['oneway'] = bool(data['oneway'])
        data['length'] = float(data['length'])

        # these attributes might have a single value, or a list if edge's topology was simplified
        for attr in ['highway', 'name', 'bridge', 'tunnel', 'lanes', 'ref', 'maxspeed', 'service', 'access', 'area', 'landuse', 'width', 'est_width']:
            # if this edge has this attribute, and it starts with '[' and ends with ']', then it's a list to be parsed
            if attr in data and data[attr][0] == '[' and data[attr][-1] == ']':
                # convert the string list to a list type, else leave as single-value string
                data[attr] = ast.literal_eval(data[attr])

        # osmid might have a single value or a list, but if single value, then parse int
        if 'osmid' in data:
            if data['osmid'][0] == '[' and data['osmid'][-1] == ']':
                data['osmid'] = ast.literal_eval(data['osmid'])
            else:
                data['osmid'] = int(data['osmid'])

        # if geometry attribute exists, load the string as well-known text to shapely LineString
        if 'geometry' in data:
            data['geometry'] = wkt.loads(data['geometry'])
    
    log('Loaded graph with {:,} nodes and {:,} edges in {:,.2f} seconds from "{}"'.format(len(G.nodes()),
                                                                                          len(G.edges()),
                                                                                          time.time()-start_time,
                                                                                          path))
    return G
    
    

    
    

    


    
    
def graph_to_gdfs(G, nodes=True, edges=True, node_geometry=True, fill_edge_geometry=True):
    """
    Convert a graph into node and/or edge GeoDataFrames
    
    Parameters
    ----------
    G : graph
    nodes : bool, if True, convert graph nodes to a GeoDataFrame and return it
    edges : bool, if True, convert graph edges to a GeoDataFrame and return it
    node_geometry : bool, if True, create a geometry column from node x and y data
    fill_edge_geometry : bool, if True, fill in missing edge geometry fields using origin and destination nodes
    
    Returns
    -------
    gdf_nodes : GeoDataFrame (optional)
    gdf_edges : GeoDataFrame (optional)
    """

    if not (nodes or edges):
        raise ValueError('You must request nodes or edges, or both.')
    
    to_return = []

    if nodes:
        
        start_time = time.time()
        
        nodes = {node:data for node, data in G.nodes(data=True)}
        gdf_nodes = gpd.GeoDataFrame(nodes).T
        if node_geometry:
            gdf_nodes['geometry'] = gdf_nodes.apply(lambda row: Point(row['x'], row['y']), axis=1)
        gdf_nodes.crs = G.graph['crs']
        gdf_nodes.gdf_name = '{}_nodes'.format(G.graph['name'])
        gdf_nodes['osmid'] = gdf_nodes['osmid'].astype(np.int64).map(make_str)
        
        to_return.append(gdf_nodes)
        log('Created GeoDataFrame "{}" from graph in {:,.2f} seconds'.format(gdf_nodes.gdf_name, time.time()-start_time))

    if edges:
        
        start_time = time.time()
        
        # create a list to hold our edges, then loop through each edge in the graph
        edges = []
        for u, v, key, data in G.edges(keys=True, data=True):

            # for each edge, add key and all attributes in data dict to the edge_details
            edge_details = {'u':u, 'v':v, 'key':key}
            for attr_key in data:
                edge_details[attr_key] = data[attr_key]

            # if edge doesn't already have a geometry attribute, create one now if fill_edge_geometry==True
            if not 'geometry' in data:
                if fill_edge_geometry:
                    point_u = Point((G.node[u]['x'], G.node[u]['y']))
                    point_v = Point((G.node[v]['x'], G.node[v]['y']))
                    edge_details['geometry'] = LineString([point_u, point_v])
                else:
                    edge_details['geometry'] = np.nan
            
            edges.append(edge_details)
        
        # create a GeoDataFrame from the list of edges and set the CRS
        gdf_edges = gpd.GeoDataFrame(edges)
        gdf_edges.crs = G.graph['crs']
        gdf_edges.gdf_name = '{}_edges'.format(G.graph['name'])
        
        to_return.append(gdf_edges)
        log('Created GeoDataFrame "{}" from graph in {:,.2f} seconds'.format(gdf_edges.gdf_name, time.time()-start_time))
        
    if len(to_return) > 1:
        return tuple(to_return)
    else:
        return to_return[0]
    

def gdfs_to_graph(gdf_nodes, gdf_edges):
    """
    Convert node and edge GeoDataFrames into a graph
    
    Parameters
    ----------
    gdf_nodes : GeoDataFrame
    gdf_edges : GeoDataFrame
    
    Returns
    -------
    G : graph
    """
    
    G = nx.MultiDiGraph()
    G.graph['crs'] = gdf_nodes.crs
    G.graph['name'] = gdf_nodes.gdf_name.rstrip('_nodes')
    
    # add the nodes and their attributes to the graph
    G.add_nodes_from(gdf_nodes.index)
    attributes = gdf_nodes.to_dict()
    for attribute_name in gdf_nodes.columns:
        # only add this attribute to nodes which have a non-null value for it
        attribute_values = {k:v for k, v in attributes[attribute_name].items() if pd.notnull(v)}
        nx.set_node_attributes(G, attribute_name, attribute_values)
    
    # add the edges and attributes that are not u, v, key (as they're added separately) or null
    for _, row in gdf_edges.iterrows():
        attrs = {}
        for label, value in row.iteritems():
            if (label not in ['u', 'v', 'key']) and (isinstance(value, list) or pd.notnull(value)):
                attrs[label] = value
        G.add_edge(u=row['u'], v=row['v'], key=row['key'], attr_dict=attrs)
    
    return G    

        
        
        
def make_shp_filename(place_name):
    """
    Create a filename string in a consistent format from a place name string.
    
    Parameters
    ----------
    place_name : string, place name to convert into a filename
    
    Returns
    -------
    filename : string
    """
    name_pieces = list(reversed(place_name.split(', ')))
    filename = '-'.join(name_pieces).lower().replace(' ','_')
    filename = re.sub('[^0-9a-zA-Z_-]+', '', filename)
    return filename


def save_gdf_shapefile(gdf, filename=None, folder=None):
    """
    Save GeoDataFrame as an ESRI shapefile.
    
    Parameters
    ----------
    gdf : GeoDataFrame, the gdf to be saved
    filename : string, what to call the shapefile (file extensions are added automatically)
    folder : string, where to save the shapefile, if none, then default folder
    
    Returns
    -------
    None
    """
    if folder is None:
        folder = globals.global_data_folder
    
    if filename is None:
        filename = make_shp_filename(gdf.name)
    
    # give the save folder a filename subfolder to make the full path to the files
    folder_path = '{}/{}'.format(folder, filename)
    
    # if the save folder does not already exist, create it with a filename subfolder
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    gdf.to_file(folder_path)
    if not hasattr(gdf, 'name'):
        gdf.name = 'unnamed'
    log('Saved the GeoDataFrame "{}" as shapefile "{}"'.format(gdf.name, folder_path))        
    

def convert_shp2graph(p, make_G_bidi = True, name='unamed'):
    """
    Converts shapefile to routable networkx graph.
    
    Parameters
    ----------
    p : str, File path - allowed formats geojson and ESRI Shapefile and other formats Fiona can read and write
    make_G_bidi : bool, if True, assumes linestrings are bidirectional
    name : str, Optional name of graph
    
    Returns
    -------
    G : graph
    """
    
    # Load shapefile into GeoDataFrame
    gdf = gp.read_file(p)
    
    # shapefile needs to include minimal: geometry linestring and the length computed (e.g. in QGIS)
    if 'length' not in gdf.columns:
        raise Exception('Shapefile is invalid: length not in attributes:\n{}'.format(gdf.columns))

    if  not gdf.geometry.map(lambda x: type(x) ==  LineString).all():
        s_invalid_geo = gdf.geometry[gdf.geometry.map(lambda x: type(x) ==  LineString)]
        raise Exception('Shapefile is invalid: geometry not all linestring \n{}'.format(s_invalid_geo))
        
    # Compute the start- and end-position based on linestring 
    gdf['Start_pos'] = gdf.geometry.apply(lambda x: x.coords[0])
    gdf['End_pos'] = gdf.geometry.apply(lambda x: x.coords[-1])
    
    # Create Series of unique nodes and their associated position
    s_points = gdf.Start_pos.append(gdf.End_pos).reset_index(drop=True)
    s_points = s_points.drop_duplicates()   
    log('GeoDataFrame has {} elements (linestrings) and {} unique nodes'.format(len(gdf),len(s_points)))
    
    # Add index of start and end node of linestring to geopandas DataFrame
    df_points = pd.DataFrame(s_points, columns=['Start_pos'])
    df_points['FNODE_'] = df_points.index
    gdf = pd.merge(gdf, df_points, on='Start_pos', how='inner')

    df_points = pd.DataFrame(s_points, columns=['End_pos'])
    df_points['TNODE_'] = df_points.index
    gdf = pd.merge(gdf, df_points, on='End_pos', how='inner')
    
    # Bring nodes and their position in form needed for osmnx (give arbitrary osmid (index) despite not osm file)
    df_points.columns = ['pos', 'osmid'] 
    df_points[['x', 'y']] = df_points['pos'].apply(pd.Series)
    df_node_xy = df_points.drop('pos', 1)
    
    # Create Graph Object
    G = nx.MultiDiGraph(name=name, crs=gdf.crs)
    
    # Add nodes to graph
    for node, data in df_node_xy.T.to_dict().items():
        G.add_node(node, **data)
        
    # Add edges to graph
    for i, row  in gdf.iterrows():
        dict_row  = row.to_dict()
        if 'geometry' in dict_row: del dict_row['geometry']
        G.add_edge(u=dict_row['FNODE_'], v=dict_row['TNODE_'], **dict_row)
        
    if make_G_bidi:
        gdf.rename(columns={'Start_pos': 'End_pos',
                   'End_pos': 'Start_pos',
                   'FNODE_': 'TNODE_', 
                   'TNODE_': 'FNODE_', }, inplace=True)
        
        # Add edges to graph
        for i, row  in gdf.iterrows():
            dict_row  = row.to_dict()
            if 'geometry' in dict_row: del dict_row['geometry']
            G.add_edge(u=dict_row['FNODE_'], v=dict_row['TNODE_'], **dict_row)
#         Alternative implementation        
#         G = G.to_undirected() # Some function in osmnx do not work anymore
        
    # Log information
    log('Graph has been successfully generated /n {}'.format(nx.info(G)))
    log('Show graph data structure EDGE'.format(G.get_edge_data(*list(G.edges())[0])))
    log('Show graph data structure NODE'.format(list(G.nodes())[0]))
    return G
