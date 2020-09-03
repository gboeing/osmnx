"""Serialize graphs to/from files on disk."""

import ast
import os
from xml.etree import ElementTree as etree

import networkx as nx
import numpy as np
import pandas as pd
from shapely import wkt

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
        path to the GraphML file including extension. if None, use
        default data folder + graph.graphml
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
        filepath = os.path.join(settings.data_folder, "graph.graphml")

    # if save folder does not already exist, create it
    folder, filename = os.path.split(filepath)
    if not folder == "" and not os.path.exists(folder):
        os.makedirs(folder)

    # make a copy to not edit the original graph object the caller passed in
    G = G.copy()

    if gephi:

        gdf_nodes, gdf_edges = utils_graph.graph_to_gdfs(G)

        # turn each edge's key into a unique ID for Gephi compatibility
        gdf_edges["key"] = range(len(gdf_edges))

        # gephi doesn't handle node attrs named x and y well, so rename
        gdf_nodes["xcoord"] = gdf_nodes["x"]
        gdf_nodes["ycoord"] = gdf_nodes["y"]
        G = utils_graph.graph_from_gdfs(gdf_nodes, gdf_edges)

        # remove graph attributes as Gephi only accepts node and edge attrs
        G.graph = dict()

    else:
        # if not gephi, keep graph attrs and stringify them
        for dict_key in G.graph:
            # convert all the graph attribute values to strings
            G.graph[dict_key] = str(G.graph[dict_key])

    # stringify node and edge attributes
    for _, data in G.nodes(data=True):
        for dict_key in data:
            if gephi and dict_key in {"xcoord", "ycoord"}:
                # don't convert x y values to string if saving for gephi
                continue
            else:
                # convert all the node attribute values to strings
                data[dict_key] = str(data[dict_key])

    for _, _, data in G.edges(keys=False, data=True):
        for dict_key in data:
            # convert all the edge attribute values to strings
            data[dict_key] = str(data[dict_key])

    nx.write_graphml(G, path=filepath, encoding=encoding)
    utils.log(f'Saved graph as GraphML file at "{filepath}"')


def load_graphml(filepath, node_type=int):
    """
    Load an OSMnx-saved GraphML file from disk.

    Converts the node/edge attributes to appropriate data types.

    Parameters
    ----------
    filepath : string
        path to the GraphML file
    node_type : type
        convert node ids to this data type

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    # read the graph from disk
    G = nx.MultiDiGraph(nx.read_graphml(filepath, node_type=node_type))

    # convert node/edge attribute data types
    utils.log("Converting node and edge attribute data types")
    G = _convert_node_attr_types(G, node_type)
    G = _convert_edge_attr_types(G, node_type)

    # convert graph crs attribute from saved string to correct dict data type
    # if it is a stringified dict rather than a proj4 string
    if "crs" in G.graph and G.graph["crs"].startswith("{") and G.graph["crs"].endswith("}"):
        G.graph["crs"] = ast.literal_eval(G.graph["crs"])

    if "streets_per_node" in G.graph:
        G.graph["streets_per_node"] = ast.literal_eval(G.graph["streets_per_node"])

    # remove node_default and edge_default metadata keys if they exist
    if "node_default" in G.graph:
        del G.graph["node_default"]
    if "edge_default" in G.graph:
        del G.graph["edge_default"]

    utils.log(f'Loaded graph with {len(G)} nodes and {len(G.edges)} edges from "{filepath}"')
    return G


def _convert_node_attr_types(G, node_type):
    """
    Convert graph nodes' attributes' types from string to numeric.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    node_type : type
        convert node ID (osmid) to this type

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    for _, data in G.nodes(data=True):

        # convert node ID to user-requested type
        data["osmid"] = node_type(data["osmid"])

        # convert numeric node attributes from string to float
        for attr in {"elevation", "elevation_res", "lat", "lon", "x", "y"}:
            if attr in data:
                data[attr] = float(data[attr])

    return G


def _convert_edge_attr_types(G, node_type):
    """
    Convert graph edges' attributes' types from string to numeric.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    node_type : type
        convert osmid to this type

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    # convert numeric, bool, and list edge attributes from string to correct data types
    for _, _, data in G.edges(data=True, keys=False):

        # parse length to float: should always have only 1 value
        data["length"] = float(data["length"])

        try:
            data["oneway"] = ast.literal_eval(data["oneway"])
            raise ValueError()
        except (KeyError, ValueError):
            # may lack oneway if settings.all_oneway=True when you created
            # graph, or values it can't eval if settings.all_oneway=True
            pass

        # convert to float any possible OSMnx-added edge attributes, which may
        # have multiple values if graph was simplified after they were added
        for attr in {"grade", "grade_abs", "bearing", "speed_kph", "travel_time"}:
            if attr in data:
                if data[attr].startswith("[") and data[attr].endswith("]"):
                    # if it's a list, eval it then convert each item to float
                    data[attr] = [float(a) for a in ast.literal_eval(data[attr])]
                else:
                    data[attr] = float(data[attr])

        # these attributes might have a single value, or a list if edge's
        # topology was simplified
        for attr in {
            "highway",
            "name",
            "bridge",
            "tunnel",
            "lanes",
            "ref",
            "maxspeed",
            "service",
            "access",
            "area",
            "landuse",
            "width",
            "est_width",
        }:
            # if this edge has this attribute, and it starts with '[' and ends
            # with ']', then it's a list to be parsed
            if attr in data and data[attr].startswith("[") and data[attr].endswith("]"):
                # try to convert the string list to a list type, else leave as
                # single-value string (and leave as string if error)
                try:
                    data[attr] = ast.literal_eval(data[attr])
                except Exception:
                    pass

        # osmid might have a single value or a list
        if "osmid" in data:
            if data["osmid"].startswith("[") and data["osmid"].endswith("]"):
                # if it's a list, eval list then convert each element to node_type
                data["osmid"] = [node_type(i) for i in ast.literal_eval(data["osmid"])]
            else:
                # if it's not a list, convert it to the node_type
                data["osmid"] = node_type(data["osmid"])

        # if geometry attribute exists, load the string as well-known text to
        # shapely LineString
        if "geometry" in data:
            data["geometry"] = wkt.loads(data["geometry"])

        # delete the extraneous id attribute added by graphml saving
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
    # default filepath if none was provided
    if filepath is None:
        filepath = os.path.join(settings.data_folder, "graph.osm")

    # if save folder does not already exist, create it
    folder, filename = os.path.split(filepath)
    if not folder == "" and not os.path.exists(folder):
        os.makedirs(folder)

    if not settings.all_oneway:
        raise UserWarning(
            "In order for save_graph_osm to behave properly "
            "the graph must have been created with the "
            "`all_oneway` setting set to True."
        )

    try:
        gdf_nodes, gdf_edges = data
    except ValueError:
        gdf_nodes, gdf_edges = utils_graph.graph_to_gdfs(
            data, node_geometry=False, fill_edge_geometry=False
        )

    # rename columns per osm specification
    gdf_nodes.rename(columns={"osmid": "id", "x": "lon", "y": "lat"}, inplace=True)
    if "id" in gdf_edges.columns:
        gdf_edges = gdf_edges[[col for col in gdf_edges if col != "id"]]
    if "uniqueid" in gdf_edges.columns:
        gdf_edges = gdf_edges.rename(columns={"uniqueid": "id"})
    else:
        gdf_edges = gdf_edges.reset_index().rename(columns={"index": "id"})

    # add default values for required attributes
    for table in (gdf_nodes, gdf_edges):
        table["uid"] = "1"
        table["user"] = "osmnx"
        table["version"] = "1"
        table["changeset"] = "1"
        table["timestamp"] = "2017-01-01T00:00:00Z"

    # convert all datatypes to str
    gdf_nodes = gdf_nodes.applymap(str)
    gdf_edges = gdf_edges.applymap(str)

    # misc. string replacements to meet OSM XML spec
    if "oneway" in gdf_edges.columns:

        # fill blank oneway tags with default (False)
        gdf_edges.loc[pd.isnull(gdf_edges["oneway"]), "oneway"] = oneway
        gdf_edges.loc[:, "oneway"] = gdf_edges["oneway"].astype(str)
        gdf_edges.loc[:, "oneway"] = (
            gdf_edges["oneway"].str.replace("False", "no").replace("True", "yes")
        )

    # initialize XML tree with an OSM root element then append nodes/edges
    root = etree.Element("osm", attrib={"version": "1", "generator": "OSMnx"})
    root = _append_nodes_xml_tree(root, gdf_nodes, node_attrs, node_tags)
    root = _append_edges_xml_tree(
        root, gdf_edges, edge_attrs, edge_tags, edge_tag_aggs, merge_edges
    )

    # write to disk
    etree.ElementTree(root).write(filepath)
    utils.log(f'Saved graph as .osm file at "{filepath}"')


def _append_nodes_xml_tree(root, gdf_nodes, node_attrs, node_tags):
    """
    Append nodes to an XML tree.

    Parameters
    ----------
    root : ElementTree.Element
        xml tree
    gdf_nodes : geopandas.GeoDataFrame
        GeoDataFrame of graph nodes
    node_attrs : list
        osm way attributes to include in output OSM XML
    node_tags : list
        osm way tags to include in output OSM XML

    Returns
    -------
    root : ElementTree.Element
        xml tree with nodes appended
    """
    for _, row in gdf_nodes.iterrows():
        node = etree.SubElement(root, "node", attrib=row[node_attrs].dropna().to_dict())
        for tag in node_tags:
            if tag in gdf_nodes.columns:
                etree.SubElement(node, "tag", attrib={"k": tag, "v": row[tag]})
    return root


def _append_edges_xml_tree(root, gdf_edges, edge_attrs, edge_tags, edge_tag_aggs, merge_edges):
    """
    Append edges to an XML tree.

    Parameters
    ----------
    root : ElementTree.Element
        xml tree
    gdf_edges : geopandas.GeoDataFrame
        GeoDataFrame of graph edges
    edge_attrs : list
        osm way attributes to include in output OSM XML
    edge_tags : list
        osm way tags to include in output OSM XML
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
    merge_edges : bool
        if True merges graph edges such that each OSM way has one entry
        and one entry only in the OSM XML. Otherwise, every OSM way
        will have a separate entry for each node pair it contains.

    Returns
    -------
    root : ElementTree.Element
        xml tree with edges appended
    """
    if merge_edges:

        for e in gdf_edges["id"].unique():
            all_way_edges = gdf_edges[gdf_edges["id"] == e]
            first = all_way_edges.iloc[0]
            edge = etree.SubElement(root, "way", attrib=first[edge_attrs].dropna().to_dict())

            if len(all_way_edges) == 1:
                etree.SubElement(edge, "nd", attrib={"ref": first["u"]})
                etree.SubElement(edge, "nd", attrib={"ref": first["v"]})
            else:
                # topological sort
                ordered_nodes = _get_unique_nodes_ordered_from_way(all_way_edges)
                for node in ordered_nodes:
                    etree.SubElement(edge, "nd", attrib={"ref": node})

            if edge_tag_aggs is None:
                for tag in edge_tags:
                    if tag in all_way_edges.columns:
                        etree.SubElement(edge, "tag", attrib={"k": tag, "v": first[tag]})
            else:
                for tag in edge_tags:
                    if (tag in all_way_edges.columns) and (
                        tag not in (t for t, agg in edge_tag_aggs)
                    ):
                        etree.SubElement(edge, "tag", attrib={"k": tag, "v": first[tag]})

                for tag, agg in edge_tag_aggs:
                    if tag in all_way_edges.columns:
                        etree.SubElement(
                            edge, "tag", attrib={"k": tag, "v": all_way_edges[tag].aggregate(agg)}
                        )
    else:
        # NOTE: this will generate separate OSM ways for each network edge,
        # even if the edges are all part of the same original OSM way. As
        # such, each way will be comprised of two nodes, and there will be
        # many ways with the same OSM id. This does not conform to the
        # OSM XML schema standard, however, the data will still comprise a
        # valid network and will be readable by *most* OSM tools.
        for _, row in gdf_edges.iterrows():
            edge = etree.SubElement(root, "way", attrib=row[edge_attrs].dropna().to_dict())
            etree.SubElement(edge, "nd", attrib={"ref": row["u"]})
            etree.SubElement(edge, "nd", attrib={"ref": row["v"]})
            for tag in edge_tags:
                if tag in gdf_edges.columns:
                    etree.SubElement(edge, "tag", attrib={"k": tag, "v": row[tag]})

    return root


def _get_unique_nodes_ordered_from_way(df_way_edges):
    """
    Recover original node order from df of edges associated w/ single OSM way.

    Parameters
    ----------
    df_way_edges : pandas.DataFrame
        Dataframe containing columns 'u' and 'v' corresponding to
        origin/destination nodes.

    Returns
    -------
    unique_ordered_nodes : list
        An ordered list of unique node IDs.
        Note: If the edges do not all connect (e.g. [(1, 2), (2,3),
        (10, 11), (11, 12), (12, 13)]), then this method will return
        only those nodes associated with the largest component of
        connected edges, even if subsequent connected chunks are contain
        more total nodes. This is done to ensure a proper topological
        representation of nodes in the XML way records because if there
        are unconnected components, the sorting algorithm cannot recover
        their original order. We would not likely ever encounter this
        kind of disconnected structure of nodes within a given way, but
        it is not explicitly forbidden in the OSM XML design schema.
    """
    G = nx.MultiDiGraph()
    all_nodes = list(df_way_edges["u"].values) + list(df_way_edges["v"].values)

    G.add_nodes_from(all_nodes)
    G.add_edges_from(df_way_edges[["u", "v"]].values)

    # copy nodes into new graph
    H = utils_graph.get_largest_component(G, strongly=False)
    unique_ordered_nodes = list(nx.topological_sort(H))
    num_unique_nodes = len(np.unique(all_nodes))

    if len(unique_ordered_nodes) < num_unique_nodes:
        utils.log(f"Recovered order for {len(unique_ordered_nodes)} of {num_unique_nodes} nodes")

    return unique_ordered_nodes
