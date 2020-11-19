"""Serialize graphs to/from files on disk."""

import ast
import os

import networkx as nx
import pandas as pd
from shapely import wkt

from . import osm_xml
from . import settings
from . import utils
from . import utils_graph


def save_graph_geopackage(G, filepath=None, encoding="utf-8"):
    """
    Save graph nodes and edges to disk as layers in a GeoPackage file.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    filepath : string
        path to the GeoPackage file including extension. if None, use default
        data folder + graph.gpkg
    encoding : string
        the character encoding for the saved file

    Returns
    -------
    None
    """
    # default filepath if none was provided
    if filepath is None:
        filepath = os.path.join(settings.data_folder, "graph.gpkg")

    # if save folder does not already exist, create it
    folder, filename = os.path.split(filepath)
    if not folder == "" and not os.path.exists(folder):
        os.makedirs(folder)

    # convert undirected graph to gdfs and stringify non-numeric columns
    gdf_nodes, gdf_edges = utils_graph.graph_to_gdfs(utils_graph.get_undirected(G))
    gdf_nodes = _stringify_nonnumeric_cols(gdf_nodes)
    gdf_edges = _stringify_nonnumeric_cols(gdf_edges)

    # save the nodes and edges as GeoPackage layers
    gdf_nodes.to_file(filepath, layer="nodes", driver="GPKG", encoding=encoding)
    gdf_edges.to_file(filepath, layer="edges", driver="GPKG", encoding=encoding)
    utils.log(f'Saved graph as GeoPackage at "{filepath}"')


def save_graph_shapefile(G, filepath=None, encoding="utf-8"):
    """
    Save graph nodes and edges to disk as ESRI shapefiles.

    The shapefile format is proprietary and outdated. Whenever possible, you
    should use the superior GeoPackage file format instead, for instance, via
    the save_graph_geopackage function.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    filepath : string
        path to the shapefiles folder (no file extension). if None, use
        default data folder + graph_shapefile
    encoding : string
        the character encoding for the saved files

    Returns
    -------
    None
    """
    # default filepath if none was provided
    if filepath is None:
        filepath = os.path.join(settings.data_folder, "graph_shapefile")

    # if save folder does not already exist, create it (shapefiles
    # get saved as set of files)
    if not filepath == "" and not os.path.exists(filepath):
        os.makedirs(filepath)
    filepath_nodes = os.path.join(filepath, "nodes.shp")
    filepath_edges = os.path.join(filepath, "edges.shp")

    # convert undirected graph to gdfs and stringify non-numeric columns
    gdf_nodes, gdf_edges = utils_graph.graph_to_gdfs(utils_graph.get_undirected(G))
    gdf_nodes = _stringify_nonnumeric_cols(gdf_nodes)
    gdf_edges = _stringify_nonnumeric_cols(gdf_edges)

    # save the nodes and edges as separate ESRI shapefiles
    gdf_nodes.to_file(filepath_nodes, encoding=encoding)
    gdf_edges.to_file(filepath_edges, encoding=encoding)
    utils.log(f'Saved graph as shapefiles at "{filepath}"')


def save_graphml(G, filepath=None, gephi=False, encoding="utf-8"):
    """
    Save graph to disk as GraphML file.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    filepath : string
        path to the GraphML file including extension. if None, use default
        data folder + graph.graphml
    gephi : bool
        if True, give each edge a unique key/id to work around Gephi's
        interpretation of the GraphML specification
    encoding : string
        the character encoding for the saved file

    Returns
    -------
    None
    """
    # default filepath if none was provided
    if filepath is None:
        filepath = os.path.join(settings.data_folder, "graph.graphml")

    # if save folder does not already exist, create it
    folder, filename = os.path.split(filepath)
    if not folder == "" and not os.path.exists(folder):
        os.makedirs(folder)

    if gephi:
        # for gephi compatibility, each edge's key must be unique as an id
        uvdk = zip(G.edges(keys=False, data=True), range(len(G.edges)))
        uvkd = ((u, v, k, d) for (u, v, d), k in uvdk)
        G = nx.MultiDiGraph(uvkd)

    else:
        # make a copy to not mutate original graph object caller passed in
        G = G.copy()

    # stringify all the graph attribute values
    for attr, value in G.graph.items():
        G.graph[attr] = str(value)

    # stringify all the node attribute values
    for _, data in G.nodes(data=True):
        for attr, value in data.items():
            data[attr] = str(value)

    # stringify all the edge attribute values
    for _, _, data in G.edges(keys=False, data=True):
        for attr, value in data.items():
            data[attr] = str(value)

    nx.write_graphml(G, path=filepath, encoding=encoding)
    utils.log(f'Saved graph as GraphML file at "{filepath}"')


def load_graphml(filepath, node_dtypes=None, edge_dtypes=None):
    """
    Load an OSMnx-saved GraphML file from disk.

    Converts the node/edge attributes to appropriate data types, which can be
    customized if needed by passing in node_dtypes or edge_dtypes arguments.

    Parameters
    ----------
    filepath : string
        path to the GraphML file
    node_dtypes : dict
        dict of node attribute names:types to convert values' data types
    edge_dtypes : dict
        dict of edge attribute names:types to convert values' data types

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    # specify default node/edge attribute values' data types
    default_node_dtypes = {
        "elevation": float,
        "elevation_res": float,
        "lat": float,
        "lon": float,
        "osmid": int,
        "x": float,
        "y": float,
    }
    default_edge_dtypes = {
        "bearing": float,
        "grade": float,
        "grade_abs": float,
        "length": float,
        "osmid": int,
        "speed_kph": float,
        "travel_time": float,
    }

    # override default node/edge attr types with user-passed types, if any
    if node_dtypes is not None:
        default_node_dtypes.update(node_dtypes)
    if edge_dtypes is not None:
        default_edge_dtypes.update(edge_dtypes)

    # read the graphml file from disk
    G = nx.read_graphml(filepath, node_type=default_node_dtypes["osmid"], force_multigraph=True)

    # convert node/edge attribute data types
    utils.log("Converting node and edge attribute data types")
    G = _convert_node_attr_types(G, default_node_dtypes)
    G = _convert_edge_attr_types(G, default_edge_dtypes)

    # eval non-string graph attributes to convert to correct types
    for attr in {"simplified", "streets_per_node"}:
        try:
            G.graph[attr] = ast.literal_eval(G.graph[attr])
        except Exception:
            pass

    # remove node_default and edge_default metadata keys if they exist
    if "node_default" in G.graph:
        del G.graph["node_default"]
    if "edge_default" in G.graph:
        del G.graph["edge_default"]

    utils.log(f'Loaded graph with {len(G)} nodes and {len(G.edges)} edges from "{filepath}"')
    return G


def _convert_node_attr_types(G, dtypes=None):
    """
    Convert graph nodes' attributes using a dict of data types.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    dtypes : dict
        dict of node attribute names:types

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    for _, data in G.nodes(data=True):
        for attr in dtypes:
            if attr in data:
                data[attr] = dtypes[attr](data[attr])
    return G


def _convert_edge_attr_types(G, dtypes=None):
    """
    Convert graph edges' attributes using a dict of data types.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    dtypes : dict
        dict of edge attribute names:types

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    # for each edge in the graph, eval attribute value lists and convert types
    for _, _, data in G.edges(data=True, keys=False):

        # edges attributes might have a single value, or a list if simplified
        # first, eval stringified lists to convert them to list objects
        for attr, value in data.items():
            if value.startswith("[") and value.endswith("]"):
                try:
                    data[attr] = ast.literal_eval(value)
                except Exception:
                    pass

        # next, convert attribute value types
        for attr in dtypes:
            if attr in data:
                if isinstance(data[attr], list):
                    # if it's a list, eval it then convert each item
                    data[attr] = [dtypes[attr](item) for item in data[attr]]
                else:
                    # otherwise, just convert the single value
                    data[attr] = dtypes[attr](data[attr])

        # next, eval "oneway" attr to bool
        try:
            data["oneway"] = ast.literal_eval(data["oneway"])
        except (KeyError, ValueError):
            # may lack oneway if settings.all_oneway=True when you created
            # graph, or have values it can't eval if settings.all_oneway=True
            pass

        # if "geometry" attr exists, convert its well-known text to LineString
        if "geometry" in data:
            data["geometry"] = wkt.loads(data["geometry"])

        # delete extraneous "id" attribute added by graphml saving
        if "id" in data:
            del data["id"]

    return G


def _stringify_nonnumeric_cols(gdf):
    """
    Make every non-numeric GeoDataFrame column (besides geometry) a string.

    This allows proper serializing via Fiona of GeoDataFrames with mixed types
    such as strings and ints in the same column.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        gdf to stringify non-numeric columns of

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        gdf with non-numeric columns stringified
    """
    # stringify every non-numeric column other than geometry column
    for col in (c for c in gdf.columns if not c == "geometry"):
        if not pd.api.types.is_numeric_dtype(gdf[col]):
            gdf[col] = gdf[col].fillna("").astype(str)

    return gdf


def save_graph_xml(
    data,
    filepath=None,
    node_tags=settings.osm_xml_node_tags,
    node_attrs=settings.osm_xml_node_attrs,
    edge_tags=settings.osm_xml_way_tags,
    edge_attrs=settings.osm_xml_way_attrs,
    oneway=False,
    merge_edges=True,
    edge_tag_aggs=None,
):
    """
    Save graph to disk as an OSM-formatted XML .osm file.

    This function exists only to allow serialization to the .osm file format
    for applications that require it, and has constraints to conform to that.
    To save/load full-featured OSMnx graphs to/from disk for later use, use
    the save_graphml and load_graphml functions instead.

    Note: for large networks this function can take a long time to run. Before
    using this function, make sure you configured OSMnx as described in the
    example below when you created the graph.

    Example
    -------
    >>> import osmnx as ox
    >>> utn = ox.settings.useful_tags_node
    >>> oxna = ox.settings.osm_xml_node_attrs
    >>> oxnt = ox.settings.osm_xml_node_tags
    >>> utw = ox.settings.useful_tags_way
    >>> oxwa = ox.settings.osm_xml_way_attrs
    >>> oxwt = ox.settings.osm_xml_way_tags
    >>> utn = list(set(utn + oxna + oxnt))
    >>> utw = list(set(utw + oxwa + oxwt))
    >>> ox.config(all_oneway=True, useful_tags_node=utn, useful_tags_way=utw)
    >>> G = ox.graph_from_place('Piedmont, CA, USA', network_type='drive')
    >>> ox.save_graph_xml(G, filepath='./data/graph1.osm')

    Parameters
    ----------
    data : networkx multi(di)graph OR a length 2 iterable of nodes/edges
        geopandas GeoDataFrames
    filepath : string
        path to the .osm file including extension. if None, use default data
        folder + graph.osm
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
    osm_xml._save_graph_xml(
        data,
        filepath,
        node_tags,
        node_attrs,
        edge_tags,
        edge_attrs,
        oneway,
        merge_edges,
        edge_tag_aggs,
    )
