################################################################################
# Module: save_load.py
# Description: Save and load networks to/from disk
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

import ast
import geopandas as gpd
import os
import networkx as nx
import numpy as np
import pandas as pd
import re
from shapely import wkt
from shapely.geometry import LineString
from shapely.geometry import Point
from xml.etree import ElementTree as etree
from . import settings
from . import utils
from . import utils_graph



def save_graph_geopackage(G, filepath=None, encoding='utf-8'):
    """
    Save graph nodes and edges as to disk as layers in a GeoPackage file.

    Parameters
    ----------
    G : networkx.MultiDiGraph
    filepath : string
        path to the GeoPackage file including extension. if None, use
        default data folder + graph.gpkg
    encoding : string
        the character encoding for the saved file

    Returns
    -------
    None
    """

    # default filepath if none was provided
    if filepath is None:
        filepath = os.path.join(settings.data_folder, 'graph.gpkg')

    # if save folder does not already exist, create it
    folder, filename = os.path.split(filepath)
    if not folder == '' and not os.path.exists(folder):
        os.makedirs(folder)

    # convert undirected graph to geodataframes
    gdf_nodes, gdf_edges = utils_graph.graph_to_gdfs(get_undirected(G))

    # make every non-numeric edge attribute (besides geometry) a string
    for col in [c for c in gdf_edges.columns if not c == 'geometry']:
        if not pd.api.types.is_numeric_dtype(gdf_edges[col]):
            gdf_edges[col] = gdf_edges[col].fillna('').astype(str)

    # save the nodes and edges as GeoPackage layers
    gdf_nodes.to_file(filepath, layer='nodes', driver='GPKG', encoding=encoding)
    gdf_edges.to_file(filepath, layer='edges', driver='GPKG', encoding=encoding)
    utils.log(f'Saved graph as GeoPackage at "{filepath}"')



def save_graph_shapefile(G, filepath=None, encoding='utf-8'):
    """
    Save graph nodes and edges to disk as ESRI shapefiles.

    Parameters
    ----------
    G : networkx.MultiDiGraph
    filepath : string
        path to the shapefiles folder (no file extension). if None,
        use default data folder
    encoding : string
        the character encoding for the saved files

    Returns
    -------
    None
    """

    # default filepath if none was provided
    if filepath is None:
        filepath = settings.data_folder

    # if save folder does not already exist, create it (shapefiles
    # get saved as set of files)
    if not filepath == '' and not os.path.exists(filepath):
        os.makedirs(filepath)
    filepath_nodes = os.path.join(filepath, 'nodes.shp')
    filepath_edges = os.path.join(filepath, 'edges.shp')

    # convert undirected graph to geodataframes
    gdf_nodes, gdf_edges = utils_graph.graph_to_gdfs(get_undirected(G))

    # make every non-numeric edge attribute (besides geometry) a string
    for col in [c for c in gdf_edges.columns if not c == 'geometry']:
        if not pd.api.types.is_numeric_dtype(gdf_edges[col]):
            gdf_edges[col] = gdf_edges[col].fillna('').astype(str)

    # make every non-numeric node attribute (besides geometry) a string
    for col in [c for c in gdf_nodes.columns if not c == 'geometry']:
        if not pd.api.types.is_numeric_dtype(gdf_nodes[col]):
            gdf_nodes[col] = gdf_nodes[col].fillna('').astype(str)

    # save the nodes and edges as separate ESRI shapefiles
    gdf_nodes.to_file(filepath_nodes, encoding=encoding)
    gdf_edges.to_file(filepath_edges, encoding=encoding)
    utils.log(f'Saved graph as shapefiles at "{filepath}"')



def save_graph_osm(data, filepath=None,
                   node_tags=settings.osm_xml_node_tags,
                   node_attrs=settings.osm_xml_node_attrs,
                   edge_tags=settings.osm_xml_way_tags,
                   edge_attrs=settings.osm_xml_way_attrs,
                   oneway=False, merge_edges=True, edge_tag_aggs=None):
    """
    Save a graph as a .osm XML formatted file. Note: for very large
    networks this function can take a long time to finish.

    Parameters
    ----------
    data : networkx multi(di)graph OR a length 2 iterable of nodes/edges
        geopandas.GeoDataFrames
    filepath : string
        path to the .osm file including extension
    node_tags : list
        osm node tags to include in output OSM XML
    node_attrs: list
        osm node attributes to include in output OSM XML
    edge_tags : list
        osm way tags to include in output OSM XML
    edge_attrs : list
        osm way attributes to include in output OSM XML
    oneway : bool
        the default oneway value used to fill this tag where missing
    merge_edges : bool
        if True merges graph edges such that each OSM way has one entry
        and one entry only in the OSM XML. Otherwise, every OSM way
        will have a separate entry for each node pair it contains.
    edge_tag_aggs : list of length-2 string tuples
        useful only if merge_edges is True, this argument allows the user
        to specify edge attributes to aggregate such that the merged
        OSM way entry tags accurately represent the sum total of
        their component edge attributes. For example, if the user
        wants the OSM way to have a "length" attribute, the user must
        specify `edge_tag_aggs=[('length', 'sum')]` in order to tell
        this method to aggregate the lengths of the individual
        component edges. Otherwise, the length attribute will simply
        reflect the length of the first edge associated with the way.

    Returns
    -------
    None
    """

    # default filepath if none was provided
    if filepath is None:
        filepath = os.path.join(settings.data_folder, 'graph.osm')

    # if save folder does not already exist, create it
    folder, filename = os.path.split(filepath)
    if not folder == '' and not os.path.exists(folder):
        os.makedirs(folder)

    try:
        assert settings.all_oneway
    except AssertionError:
        raise UserWarning('In order for save_graph_osm to behave properly '
                          'the graph must have been created with the '
                          '`all_oneway` setting set to True.')

    try:
        gdf_nodes, gdf_edges = data
    except ValueError:
        gdf_nodes, gdf_edges = utils_graph.graph_to_gdfs(
            data, node_geometry=False, fill_edge_geometry=False)

    # rename columns per osm specification
    gdf_nodes.rename(
        columns={'osmid': 'id', 'x': 'lon', 'y': 'lat'}, inplace=True)
    if 'id' in gdf_edges.columns:
        gdf_edges = gdf_edges[[col for col in gdf_edges if col != 'id']]
    if 'uniqueid' in gdf_edges.columns:
        gdf_edges = gdf_edges.rename(columns={'uniqueid': 'id'})
    else:
        gdf_edges = gdf_edges.reset_index().rename(columns={'index': 'id'})

    # add default values for required attributes
    for table in [gdf_nodes, gdf_edges]:
        table['uid'] = '1'
        table['user'] = 'osmnx'
        table['version'] = '1'
        table['changeset'] = '1'
        table['timestamp'] = '2017-01-01T00:00:00Z'

    # convert all datatypes to str
    nodes = gdf_nodes.applymap(str)
    edges = gdf_edges.applymap(str)

    # misc. string replacements to meet OSM XML spec
    if 'oneway' in edges.columns:

        # fill blank oneway tags with default (False)
        edges.loc[pd.isnull(edges['oneway']), 'oneway'] = oneway
        edges.loc[:, 'oneway'] = edges['oneway'].astype(str)
        edges.loc[:, 'oneway'] = edges['oneway'].str.replace(
            'False', 'no').replace('True', 'yes')

    # initialize XML tree with an OSM root element
    root = etree.Element('osm', attrib={'version': '1', 'generator': 'OSMnx'})

    # append nodes to the XML tree
    for i, row in nodes.iterrows():
        node = etree.SubElement(
            root, 'node', attrib=row[node_attrs].dropna().to_dict())
        for tag in node_tags:
            if tag in nodes.columns:
                etree.SubElement(
                    node, 'tag', attrib={'k': tag, 'v': row[tag]})

    # append edges to the XML tree
    if merge_edges:
        for e in edges['id'].unique():
            all_way_edges = edges[edges['id'] == e]
            first = all_way_edges.iloc[0]
            edge = etree.SubElement(
                root, 'way', attrib=first[edge_attrs].dropna().to_dict())

            if len(all_way_edges) == 1:

                etree.SubElement(edge, 'nd', attrib={'ref': first['u']})
                etree.SubElement(edge, 'nd', attrib={'ref': first['v']})

            else:

                # topological sort
                ordered_nodes = get_unique_nodes_ordered_from_way(all_way_edges)

                for node in ordered_nodes:
                    etree.SubElement(edge, 'nd', attrib={'ref': node})

            if edge_tag_aggs is None:
                for tag in edge_tags:
                    if tag in all_way_edges.columns:
                        etree.SubElement(
                            edge, 'tag', attrib={'k': tag, 'v': first[tag]})
            else:
                for tag in edge_tags:
                    if tag in all_way_edges.columns:
                        if tag not in [t for t, agg in edge_tag_aggs]:
                            etree.SubElement(
                                edge, 'tag',
                                attrib={'k': tag, 'v': first[tag]})

                for tag, agg in edge_tag_aggs:
                    if tag in all_way_edges.columns:
                        etree.SubElement(edge, 'tag', attrib={
                            'k': tag, 'v': all_way_edges[tag].aggregate(agg)})

    else:

        # NOTE: this will generate separate OSM ways for each network edge,
        # even if the edges are all part of the same original OSM way. As
        # such, each way will be comprised of two nodes, and there will be
        # many ways with the same OSM id. This does not conform to the
        # OSM XML schema standard, however, the data will still comprise a
        # valid network and will be readable by *most* OSM tools.
        for i, row in edges.iterrows():
            edge = etree.SubElement(
                root, 'way', attrib=row[edge_attrs].dropna().to_dict())
            etree.SubElement(edge, 'nd', attrib={'ref': row['u']})
            etree.SubElement(edge, 'nd', attrib={'ref': row['v']})
            for tag in edge_tags:
                if tag in edges.columns:
                    etree.SubElement(
                        edge, 'tag', attrib={'k': tag, 'v': row[tag]})

    et = etree.ElementTree(root)
    et.write(filepath)
    utils.log(f'Saved graph as .osm file at "{filepath}"')



def get_unique_nodes_ordered_from_way(way_edges_df):
    """
    Function to recover the original order of nodes from a dataframe
    of edges associated with a single OSM way.

    Parameters
    ----------
    way_edges_df : pandas.DataFrame
        Dataframe containing columns 'u' and 'v' corresponding to
        origin/desitination nodes.

    Returns
    -------
    unique_ordered_nodes : list
        An ordered list of unique node IDs

    NOTE: If the edges do not all connect (e.g. [(1, 2), (2,3),
    (10, 11), (11, 12), (12, 13)]), then this method will return
    only those nodes associated with the largest component of
    connected edges, even if subsequent connected chunks are contain
    more total nodes. This is done to ensure a proper topological
    representation of nodes in the XML way records because if there
    are unconnected components, the sorting algorithm cannot recover
    their original order. I don't believe that we would ever encounter
    this kind of disconnected structure of nodes within a given way,
    but as best I could tell it is not explicitly forbidden in the
    OSM XML design schema.
    """

    G = nx.MultiDiGraph()
    all_nodes = list(way_edges_df['u'].values) + \
        list(way_edges_df['v'].values)

    G.add_nodes_from(all_nodes)
    G.add_edges_from(way_edges_df[['u', 'v']].values)

    # copy nodes into new graph
    H = utils_graph.get_largest_component(G, strongly=False)
    unique_ordered_nodes = list(nx.topological_sort(H))
    num_unique_nodes = len(np.unique(all_nodes))

    if len(unique_ordered_nodes) < num_unique_nodes:
        utils.log(f'Recovered order for {len(unique_ordered_nodes)} of {num_unique_nodes} nodes')

    return unique_ordered_nodes



def save_graphml(G, filepath=None, gephi=False, encoding='utf-8'):
    """
    Save graph to disk as GraphML file.

    Parameters
    ----------
    G : networkx.MultiDiGraph
    filepath : string
        path to the GraphML file including extension
    gephi : bool
        if True, give each edge a unique key to work around Gephi's
        restrictive interpretation of the GraphML specification
    encoding : string
        the character encoding for the saved file

    Returns
    -------
    None
    """

    # default filepath if none was provided
    if filepath is None:
        filepath = os.path.join(settings.data_folder, 'graph.graphml')

    # if save folder does not already exist, create it
    folder, filename = os.path.split(filepath)
    if not folder == '' and not os.path.exists(folder):
        os.makedirs(folder)

    # create a copy to convert all the node/edge attribute values to string
    G_save = G.copy()

    if gephi:

        gdf_nodes, gdf_edges = utils_graph.graph_to_gdfs(G_save,
                                                         nodes=True,
                                                         edges=True,
                                                         node_geometry=True,
                                                         fill_edge_geometry=True)

        # turn each edge's key into a unique ID for Gephi compatibility
        gdf_edges['key'] = range(len(gdf_edges))

        # gephi doesn't handle node attrs named x and y well, so rename
        gdf_nodes['xcoord'] = gdf_nodes['x']
        gdf_nodes['ycoord'] = gdf_nodes['y']
        G_save = utils_graph.gdfs_to_graph(gdf_nodes, gdf_edges)

        # remove graph attributes as Gephi only accepts node and edge attrs
        G_save.graph = {}

    else:
        # if not gephi, keep graph attrs and stringify them
        for dict_key in G_save.graph:
            # convert all the graph attribute values to strings
            G_save.graph[dict_key] = str(G_save.graph[dict_key])

    # stringify node and edge attributes
    for _, data in G_save.nodes(data=True):
        for dict_key in data:
            if gephi and dict_key in ['xcoord', 'ycoord']:
                # don't convert x y values to string if saving for gephi
                continue
            else:
                # convert all the node attribute values to strings
                data[dict_key] = str(data[dict_key])

    for _, _, data in G_save.edges(keys=False, data=True):
        for dict_key in data:
            # convert all the edge attribute values to strings
            data[dict_key] = str(data[dict_key])

    nx.write_graphml(G_save, path=filepath, encoding=encoding)
    utils.log(f'Saved graph as GraphML file at "{filepath}"')



def load_graphml(filepath, node_type=int):
    """
    Load an OSMnx-saved GraphML file from disk and convert the node/edge
    attributes to appropriate data types.

    Parameters
    ----------
    filepath : string
        the name of the graphml file (including file extension)
    folder : string
        the folder containing the file, if None, use default data folder
    node_type : type
        convert node ids to this type

    Returns
    -------
    networkx.MultiDiGraph
    """

    # read the graph from disk
    G = nx.MultiDiGraph(nx.read_graphml(filepath, node_type=node_type))

    # convert graph crs attribute from saved string to correct dict data type
    # if it is a stringified dict rather than a proj4 string
    if 'crs' in G.graph and G.graph['crs'].startswith('{') and G.graph['crs'].endswith('}'):
        G.graph['crs'] = ast.literal_eval(G.graph['crs'])

    if 'streets_per_node' in G.graph:
        G.graph['streets_per_node'] = ast.literal_eval(G.graph['streets_per_node'])

    # convert numeric node tags from string to numeric data types
    utils.log('Converting node and edge attribute data types')
    for _, data in G.nodes(data=True):
        data['osmid'] = node_type(data['osmid'])
        data['x'] = float(data['x'])
        data['y'] = float(data['y'])
        if 'elevation' in data:
            data['elevation'] = float(data['elevation'])
        if 'elevation_res' in data:
            data['elevation_res'] = float(data['elevation_res'])

    # convert numeric, bool, and list edge attributes from string to correct data types
    for _, _, data in G.edges(data=True, keys=False):

        # first parse oneway to bool and length to float - they should always
        # have only 1 value each
        data['oneway'] = ast.literal_eval(data['oneway'])
        data['length'] = float(data['length'])
        if 'grade' in data:
            data['grade'] = float(data['grade'])
        if 'grade_abs' in data:
            data['grade_abs'] = float(data['grade_abs'])

        # these attributes might have a single value, or a list if edge's
        # topology was simplified
        for attr in ['highway', 'name', 'bridge', 'tunnel', 'lanes', 'ref', 'maxspeed',
                     'service', 'access', 'area', 'landuse', 'width', 'est_width']:
            # if this edge has this attribute, and it starts with '[' and ends
            # with ']', then it's a list to be parsed
            if attr in data and data[attr].startswith('[') and data[attr].endswith(']'):
                # try to convert the string list to a list type, else leave as
                # single-value string (and leave as string if error)
                try:
                    data[attr] = ast.literal_eval(data[attr])
                except:
                    pass

        # osmid might have a single value or a list
        if 'osmid' in data:
            if data['osmid'][0] == '[' and data['osmid'][-1] == ']':
                # if it's a list, eval the list then convert each element to node_type
                data['osmid'] = [node_type(i) for i in ast.literal_eval(data['osmid'])]
            else:
                # if it's not a list, convert it to the node_type
                data['osmid'] = node_type(data['osmid'])

        # if geometry attribute exists, load the string as well-known text to
        # shapely LineString
        if 'geometry' in data:
            data['geometry'] = wkt.loads(data['geometry'])

    # remove node_default and edge_default metadata keys if they exist
    if 'node_default' in G.graph:
        del G.graph['node_default']
    if 'edge_default' in G.graph:
        del G.graph['edge_default']

    utils.log(f'Loaded graph with {len(G)} nodes and {len(G.edges())} edges from "{filepath}"')
    return G



def is_duplicate_edge(data, data_other):
    """
    Check if two edge data dictionaries are the same based on OSM ID and
    geometry.

    Parameters
    ----------
    data : dict
        the first edge's data
    data_other : dict
        the second edge's data

    Returns
    -------
    is_dupe : bool
    """

    is_dupe = False

    # if either edge's OSM ID contains multiple values (due to simplification), we want
    # to compare as sets so they are order-invariant, otherwise uv does not match vu
    osmid = set(data['osmid']) if isinstance(data['osmid'], list) else data['osmid']
    osmid_other = set(data_other['osmid']) if isinstance(data_other['osmid'], list) else data_other['osmid']

    if osmid == osmid_other:
        # if they contain the same OSM ID or set of OSM IDs (due to simplification)
        if ('geometry' in data) and ('geometry' in data_other):
            # if both edges have a geometry attribute
            if is_same_geometry(data['geometry'], data_other['geometry']):
                # if their edge geometries have the same coordinates
                is_dupe = True
        elif ('geometry' in data) and ('geometry' in data_other):
            # if neither edge has a geometry attribute
            is_dupe = True
        else:
            # if one edge has geometry attribute but the other doesn't, keep it
            pass

    return is_dupe



def is_same_geometry(ls1, ls2):
    """
    Check if LineString geometries in two edges are the same, in
    normal or reversed order of points.

    Parameters
    ----------
    ls1 : LineString
        the first edge's geometry
    ls2 : LineString
        the second edge's geometry

    Returns
    -------
    bool
    """

    # extract geometries from each edge data dict
    geom1 = [list(coords) for coords in ls1.xy]
    geom2 = [list(coords) for coords in ls2.xy]

    # reverse the first edge's list of x's and y's to look for a match in
    # either order
    geom1_r = [list(reversed(list(coords))) for coords in ls1.xy]

    # if the edge's geometry matches its reverse's geometry in either order,
    # return True
    return (geom1 == geom2 or geom1_r == geom2)



def update_edge_keys(G):
    """
    Update the keys of edges that share a u, v with another edge but differ in
    geometry. For example, two one-way streets from u to v that bow away from
    each other as separate streets, rather than opposite direction edges of a
    single street.

    Parameters
    ----------
    G : networkx multidigraph

    Returns
    -------
    networkx multigraph
    """

    # identify all the edges that are duplicates based on a sorted combination
    # of their origin, destination, and key. that is, edge uv will match edge vu
    # as a duplicate, but only if they have the same key
    edges = utils_graph.graph_to_gdfs(G, nodes=False, fill_edge_geometry=False)
    edges['uvk'] = edges.apply(lambda row: '_'.join(sorted([str(row['u']), str(row['v'])]) + [str(row['key'])]), axis=1)
    edges['dupe'] = edges['uvk'].duplicated(keep=False)
    dupes = edges[edges['dupe']==True].dropna(subset=['geometry'])

    different_streets = []
    groups = dupes[['geometry', 'uvk', 'u', 'v', 'key', 'dupe']].groupby('uvk')

    # for each set of duplicate edges
    for label, group in groups:

        # if there are more than 2 edges here, make sure to compare all
        if len(group['geometry']) > 2:
            l = group['geometry'].tolist()
            l.append(l[0])
            geom_pairs = list(zip(l[:-1], l[1:]))
        # otherwise, just compare the first edge to the second edge
        else:
            geom_pairs = [(group['geometry'].iloc[0], group['geometry'].iloc[1])]

        # for each pair of edges to compare
        for geom1, geom2 in geom_pairs:
            # if they don't have the same geometry, flag them as different streets
            if not is_same_geometry(geom1, geom2):
                # add edge uvk, but not edge vuk, otherwise we'll iterate both their keys
                # and they'll still duplicate each other at the end of this process
                different_streets.append((group['u'].iloc[0], group['v'].iloc[0], group['key'].iloc[0]))

    # for each unique different street, iterate its key + 1 so it's unique
    for u, v, k in set(different_streets):
        # filter out key if it appears in data dict as we'll pass it explicitly
        attributes = {k:v for k, v in G[u][v][k].items() if k != 'key'}
        G.add_edge(u, v, key=k+1, **attributes)
        G.remove_edge(u, v, key=k)

    return G



def get_undirected(G):
    """
    Convert a directed graph to an undirected graph that maintains parallel
    edges if geometries differ.

    Parameters
    ----------
    G : networkx multidigraph

    Returns
    -------
    networkx multigraph
    """

    # set from/to nodes before making graph undirected
    G = G.copy()
    for u, v, k, data in G.edges(keys=True, data=True):
        G.edges[u, v, k]['from'] = u
        G.edges[u, v, k]['to'] = v

        # add geometry if it doesn't already exist, to retain parallel
        # edges' distinct geometries
        if 'geometry' not in data:
            point_u = Point((G.nodes[u]['x'], G.nodes[u]['y']))
            point_v = Point((G.nodes[v]['x'], G.nodes[v]['y']))
            data['geometry'] = LineString([point_u, point_v])

    # update edge keys so we don't retain only one edge of sets of parallel edges
    # when we convert from a multidigraph to a multigraph
    G = update_edge_keys(G)

    # now convert multidigraph to a multigraph, retaining all edges in both
    # directions for now, as well as all graph attributes
    H = nx.MultiGraph()
    H.add_nodes_from(G.nodes(data=True))
    H.add_edges_from(G.edges(keys=True, data=True))
    H.graph = G.graph
    H.name = G.name

    # the previous operation added all directed edges from G as undirected
    # edges in H. this means we have duplicate edges for every bi-directional
    # street. so, look through the edges and remove any duplicates
    duplicate_edges = []
    for u, v, key, data in H.edges(keys=True, data=True):

        # if we haven't already flagged this edge as a duplicate
        if not (u, v, key) in duplicate_edges:

            # look at every other edge between u and v, one at a time
            for key_other in H[u][v]:

                # don't compare this edge to itself
                if not key_other == key:

                    # compare the first edge's data to the second's to see if
                    # they are duplicates
                    data_other = H.edges[u, v, key_other]
                    if is_duplicate_edge(data, data_other):

                        # if they match up, flag the duplicate for removal
                        duplicate_edges.append((u, v, key_other))

    H.remove_edges_from(duplicate_edges)
    utils.log('Made undirected graph')

    return H
