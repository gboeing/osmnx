"""Read/write .osm formatted XML files."""

import bz2
import xml.sax
from pathlib import Path
from xml.etree import ElementTree as etree

import networkx as nx
import numpy as np
import pandas as pd

from . import settings
from . import utils
from . import utils_graph


class _OSMContentHandler(xml.sax.handler.ContentHandler):
    """
    SAX content handler for OSM XML.

    Used to build an Overpass-like response JSON object in self.object. For
    format notes, see
    http://wiki.openstreetmap.org/wiki/OSM_XML#OSM_XML_file_format_notes and
    http://overpass-api.de/output_formats.html#json
    """

    def __init__(self):
        self._element = None
        self.object = {"elements": []}

    def startElement(self, name, attrs):
        if name == "osm":
            self.object.update({k: v for k, v in attrs.items() if k in {"version", "generator"}})

        elif name in {"node", "way"}:
            self._element = dict(type=name, tags={}, nodes=[], **attrs)
            self._element.update({k: float(v) for k, v in attrs.items() if k in {"lat", "lon"}})
            self._element.update(
                {k: int(v) for k, v in attrs.items() if k in {"id", "uid", "version", "changeset"}}
            )

        elif name == "relation":
            self._element = dict(type=name, tags={}, members=[], **attrs)
            self._element.update(
                {k: int(v) for k, v in attrs.items() if k in {"id", "uid", "version", "changeset"}}
            )

        elif name == "tag":
            self._element["tags"].update({attrs["k"]: attrs["v"]})

        elif name == "nd":
            self._element["nodes"].append(int(attrs["ref"]))

        elif name == "member":
            self._element["members"].append(
                {k: (int(v) if k == "ref" else v) for k, v in attrs.items()}
            )

    def endElement(self, name):
        if name in {"node", "way", "relation"}:
            self.object["elements"].append(self._element)


def _overpass_json_from_file(filepath):
    """
    Read OSM XML from file and return Overpass-like JSON.

    Parameters
    ----------
    filepath : string or pathlib.Path
        path to file containing OSM XML data

    Returns
    -------
    OSMContentHandler object
    """

    def _opener(filepath):
        if filepath.suffix == ".bz2":
            return bz2.BZ2File(filepath)
        else:
            # assume an unrecognized file extension is just XML
            return filepath.open(mode="rb")

    with _opener(Path(filepath)) as f:
        handler = _OSMContentHandler()
        xml.sax.parse(f, handler)
        return handler.object


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
    filepath : string or pathlib.Path
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
        filepath = Path(settings.data_folder) / "graph.osm"
    else:
        filepath = Path(filepath)

    # if save folder does not already exist, create it
    filepath.parent.mkdir(parents=True, exist_ok=True)

    if not settings.all_oneway:  # pragma: no cover
        import warnings

        msg = (
            "In order for save_graph_xml to behave properly the graph must "
            "have been created with the `all_oneway` setting set to True."
        )
        warnings.warn(msg)

    try:
        gdf_nodes, gdf_edges = data
    except ValueError:
        gdf_nodes, gdf_edges = utils_graph.graph_to_gdfs(
            data, node_geometry=False, fill_edge_geometry=False
        )

    # rename columns per osm specification
    gdf_nodes.rename(columns={"x": "lon", "y": "lat"}, inplace=True)
    gdf_nodes = gdf_nodes.reset_index().rename(columns={"osmid": "id"})
    if "id" in gdf_edges.columns:
        gdf_edges = gdf_edges[[col for col in gdf_edges if col != "id"]]
    if "uniqueid" in gdf_edges.columns:
        gdf_edges = gdf_edges.rename(columns={"uniqueid": "id"})
    else:
        gdf_edges = gdf_edges.reset_index().reset_index().rename(columns={"index": "id"})

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
    gdf_edges.reset_index(inplace=True)
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
    df_way_edges.reset_index(inplace=True)
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
