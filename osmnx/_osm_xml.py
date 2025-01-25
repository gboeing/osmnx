"""
Read/write OSM XML files.

For file format information see https://wiki.openstreetmap.org/wiki/OSM_XML
"""

from __future__ import annotations

import bz2
import gzip
import logging as lg
from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING
from typing import Any
from typing import TextIO
from warnings import warn
from xml.etree.ElementTree import Element
from xml.etree.ElementTree import ElementTree
from xml.etree.ElementTree import SubElement
from xml.etree.ElementTree import parse as etree_parse
from xml.sax import parse as sax_parse
from xml.sax.handler import ContentHandler

import networkx as nx
import pandas as pd

from . import convert
from . import projection
from . import settings
from . import truncate
from . import utils
from ._errors import GraphSimplificationError
from ._version import __version__ as osmnx_version

if TYPE_CHECKING:
    from collections.abc import Iterator
    from xml.sax.xmlreader import AttributesImpl

    import geopandas as gpd


# default values for standard "node" and "way" XML subelement attributes
# see: https://wiki.openstreetmap.org/wiki/Elements#Common_attributes
ATTR_DEFAULTS = {
    "changeset": "1",
    "timestamp": utils.ts(style="iso8601"),
    "uid": "1",
    "user": "OSMnx",
    "version": "1",
    "visible": "true",
}

# default values for standard "osm" root XML element attributes
# current OSM editing API version: https://wiki.openstreetmap.org/wiki/API
ROOT_ATTR_DEFAULTS = {
    "attribution": "https://www.openstreetmap.org/copyright",
    "copyright": "OpenStreetMap and contributors",
    "generator": f"OSMnx {osmnx_version}",
    "license": "https://opendatacommons.org/licenses/odbl/1-0/",
    "version": "0.6",
}


class _OSMContentHandler(ContentHandler):
    """
    SAX content handler for OSM XML.

    Builds an Overpass-like response JSON object in self.object. For format
    notes, see https://wiki.openstreetmap.org/wiki/OSM_XML and
    https://overpass-api.de
    """

    def __init__(self) -> None:
        self._element: dict[str, Any] | None = None
        self.object: dict[str, Any] = {"elements": []}

    def startElement(self, name: str, attrs: AttributesImpl) -> None:  # noqa: N802
        # identify node/way/relation attrs to convert from string to numeric
        float_attrs = {"lat", "lon"}
        int_attrs = {"changeset", "id", "uid", "version"}

        if name == "osm":
            self.object.update({k: v for k, v in attrs.items() if k in ROOT_ATTR_DEFAULTS})

        elif name in {"node", "way"}:
            self._element = dict(type=name, tags={}, **attrs)
            if name == "way":
                self._element["nodes"] = []
            self._element.update({k: float(v) for k, v in attrs.items() if k in float_attrs})
            self._element.update({k: int(v) for k, v in attrs.items() if k in int_attrs})

        elif name == "relation":
            self._element = dict(type=name, tags={}, members=[], **attrs)
            self._element.update({k: int(v) for k, v in attrs.items() if k in int_attrs})

        elif name == "tag":
            self._element["tags"].update({attrs["k"]: attrs["v"]})  # type: ignore[index]

        elif name == "nd":
            self._element["nodes"].append(int(attrs["ref"]))  # type: ignore[index]

        elif name == "member":
            self._element["members"].append(  # type: ignore[index]
                {k: (int(v) if k == "ref" else v) for k, v in attrs.items()},
            )

    def endElement(self, name: str) -> None:  # noqa: N802
        if name in {"node", "way", "relation"}:
            self.object["elements"].append(self._element)


@contextmanager
def _open_file(filepath: Path, encoding: str) -> Iterator[TextIO]:
    """
    Open a file and return a file object, optionally handling bz2 or gz files.

    Uses a wrapper context manager to yield the file object to ensure the file
    will always get closed when the caller is finished with it.

    Parameters
    ----------
    filepath
        Path to file.
    encoding
        The file's character encoding.

    Returns
    -------
    file
        The file handle.
    """
    if filepath.suffix == ".bz2":
        with bz2.open(filepath, mode="rt", encoding=encoding) as file:
            yield file
    elif filepath.suffix == ".gz":
        with gzip.open(filepath, mode="rt", encoding=encoding) as file:
            yield file
    else:
        with filepath.open(mode="rt", encoding=encoding) as file:
            yield file


def _overpass_json_from_xml(filepath: Path, encoding: str) -> dict[str, Any]:
    """
    Read OSM XML data from file and return Overpass-like JSON.

    Parameters
    ----------
    filepath
        Path to file containing OSM XML data.
    encoding
        The XML file's character encoding.

    Returns
    -------
    response_json
        A parsed JSON response from the Overpass API.
    """
    with _open_file(filepath, encoding) as file:
        # warn if this XML file was generated by OSMnx itself
        root_attrs = etree_parse(file).getroot().attrib  # noqa: S314
        if "generator" in root_attrs and "OSMnx" in root_attrs["generator"]:
            msg = (
                "The XML file you are loading appears to have been generated "
                "by OSMnx: this use case is not supported and may not behave "
                "as expected. To save/load graphs to/from disk for later use "
                "in OSMnx, use the `io.save_graphml` and `io.load_graphml` "
                "functions instead. Refer to the documentation for details."
            )
            warn(msg, category=UserWarning, stacklevel=2)

        # move back to beginning of file, then parse XML to Overpass-like JSON
        file.seek(0)
        handler = _OSMContentHandler()
        sax_parse(file, handler)  # noqa: S317

    return handler.object


def _save_graph_xml(
    G: nx.MultiDiGraph,
    filepath: str | Path | None,
    way_tag_aggs: dict[str, Any] | None,
    encoding: str = "utf-8",
) -> None:
    """
    Save graph to disk as an OSM XML file.

    Parameters
    ----------
    G
        Unsimplified, unprojected graph to save as an OSM XML file.
    filepath
        Path to the saved file including extension. If None, use default
        `settings.data_folder/graph.osm`.
    way_tag_aggs
        Keys are OSM way tag keys and values are aggregation functions
        (anything accepted as an argument by `pandas.agg`). Allows user to
        aggregate graph edge attribute values into single OSM way values. If
        None, or if some tag's key does not exist in the dict, the way
        attribute will be assigned the value of the first edge of the way.
    encoding
        The character encoding of the saved OSM XML file.

    Returns
    -------
    None
    """
    # default "oneway" value used to fill this tag where missing
    ONEWAY = False

    # round lat/lon coordinates to 7 decimals (approx 5 to 10 mm resolution)
    PRECISION = 7

    # warn user if ox.settings.all_oneway is not currently True (but maybe it
    # was when they created the graph)
    if not settings.all_oneway:
        msg = "Make sure graph was created with `ox.settings.all_oneway=True` to save as OSM XML."
        warn(msg, category=UserWarning, stacklevel=2)

    # warn user if graph is projected
    if projection.is_projected(G.graph["crs"]):
        msg = (
            "Graph should be unprojected to save as OSM XML: the existing "
            "projected x-y coordinates will be saved as lat-lon node attributes. "
            "Project your graph back to lat-lon to avoid this."
        )
        warn(msg, category=UserWarning, stacklevel=2)

    # raise error if graph has been simplified
    if G.graph.get("simplified", False):
        msg = "Graph must be unsimplified to save as OSM XML."
        raise GraphSimplificationError(msg)

    # set default filepath if None was provided
    filepath = Path(settings.data_folder) / "graph.osm" if filepath is None else Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    # convert graph to node/edge gdfs and create dict of spatial bounds
    gdf_nodes, gdf_edges = convert.graph_to_gdfs(G, fill_edge_geometry=False)
    coords = [str(round(c, PRECISION)) for c in gdf_nodes.union_all().bounds]
    bounds = dict(zip(["minlon", "minlat", "maxlon", "maxlat"], coords))

    # add default values (if missing) for standard attrs
    for gdf in (gdf_nodes, gdf_edges):
        for col, value in ATTR_DEFAULTS.items():
            if col not in gdf.columns:
                gdf[col] = value
            else:
                gdf[col] = gdf[col].fillna(value)

    # transform nodes gdf to meet OSM XML spec
    # 1) reset index (osmid) then rename osmid, x, and y columns
    # 2) round lat/lon coordinates
    # 3) drop unnecessary geometry column
    gdf_nodes = gdf_nodes.reset_index().rename(columns={"osmid": "id", "x": "lon", "y": "lat"})
    gdf_nodes[["lon", "lat"]] = gdf_nodes[["lon", "lat"]].round(PRECISION)
    gdf_nodes = gdf_nodes.drop(columns=["geometry"])

    # transform edges gdf to meet OSM XML spec
    # 1) fill and convert oneway bools to strings
    # 2) rename osmid column (but keep (u, v, k) index for processing)
    # 3) drop unnecessary geometry column
    if "oneway" in gdf_edges.columns:
        gdf_edges["oneway"] = gdf_edges["oneway"].fillna(ONEWAY).replace({True: "yes", False: "no"})
    gdf_edges = gdf_edges.rename(columns={"osmid": "id"}).drop(columns=["geometry"])

    # create parent XML element then add bounds, nodes, ways as subelements
    element = Element("osm", attrib=ROOT_ATTR_DEFAULTS)
    _ = SubElement(element, "bounds", attrib=bounds)
    _add_nodes_xml(element, gdf_nodes)
    _add_ways_xml(element, gdf_edges, way_tag_aggs)

    # write to disk
    ElementTree(element).write(filepath, encoding=encoding, xml_declaration=True)
    msg = f"Saved graph as OSM XML file at {str(filepath)!r}"
    utils.log(msg, level=lg.INFO)


def _add_nodes_xml(
    parent: Element,
    gdf_nodes: gpd.GeoDataFrame,
) -> None:
    """
    Add graph nodes as subelements of an XML parent element.

    Parameters
    ----------
    parent
        The XML parent element.
    gdf_nodes
        A GeoDataFrame of graph nodes.

    Returns
    -------
    None
    """
    node_tags = set(settings.useful_tags_node)
    node_attrs = {"id", "lat", "lon"}.union(ATTR_DEFAULTS)

    # add each node attrs dict as a SubElement of parent
    for node in gdf_nodes.to_dict(orient="records"):
        attrs = {k: str(node[k]) for k in node_attrs if pd.notna(node[k])}
        node_element = SubElement(parent, "node", attrib=attrs)

        # add each node tag dict as its own SubElement of the node SubElement
        # for vals that are non-null (or list if node consolidation was done)
        tags = (
            {"k": k, "v": str(node[k])}
            for k in node_tags & node.keys()
            if isinstance(node[k], list) or pd.notna(node[k])
        )
        for tag in tags:
            _ = SubElement(node_element, "tag", attrib=tag)


def _add_ways_xml(
    parent: Element,
    gdf_edges: gpd.GeoDataFrame,
    way_tag_aggs: dict[str, Any] | None,
) -> None:
    """
    Add graph edges (grouped as ways) as subelements of an XML parent element.

    Parameters
    ----------
    parent
        The XML parent element.
    gdf_edges
        A GeoDataFrame of graph edges with OSM way "id" column for grouping
        edges into ways.
    way_tag_aggs
        Keys are OSM way tag keys and values are aggregation functions
        (anything accepted as an argument by `pandas.agg`). Allows user to
        aggregate graph edge attribute values into single OSM way values. If
        None, or if some tag's key does not exist in the dict, the way
        attribute will be assigned the value of the first edge of the way.

    Returns
    -------
    None
    """
    way_tags = set(settings.useful_tags_way)
    way_attrs = list({"id"}.union(ATTR_DEFAULTS))

    for osmid, way in gdf_edges.groupby("id"):
        # STEP 1: add the way and its attrs as a "way" subelement of the
        # parent element
        attrs = way[way_attrs].iloc[0].astype(str).to_dict()
        way_element = SubElement(parent, "way", attrib=attrs)

        # STEP 2: add the way's edges' node IDs as "nd" subelements of the
        # "way" subelement. if way contains more than 1 edge, sort the nodes
        # topologically, otherwise just add node "u" then "v" from index.
        if len(way) == 1:
            nodes = way.index[0][:2]
        else:
            nodes = _sort_nodes(nx.MultiDiGraph(way.index.to_list()), osmid)
        for node in nodes:
            _ = SubElement(way_element, "nd", attrib={"ref": str(node)})

        # STEP 3: add way's edges' tags as "tag" subelements of the "way"
        # subelement. if an agg function was provided for a tag, apply it to
        # the values of the edges in the way. if no agg function was provided
        # for a tag, just use the value from first edge in way.
        for tag in way_tags.intersection(way.columns):
            if way_tag_aggs is not None and tag in way_tag_aggs:
                value = way[tag].agg(way_tag_aggs[tag])
            else:
                value = way[tag].iloc[0]
            if pd.notna(value):
                _ = SubElement(way_element, "tag", attrib={"k": tag, "v": str(value)})


def _sort_nodes(G: nx.MultiDiGraph, osmid: int) -> list[int]:
    """
    Topologically sort the nodes of an OSM way.

    Parameters
    ----------
    G
        The graph representing the OSM way.
    osmid
        The OSM way ID.

    Returns
    -------
    ordered_nodes
        The way's node IDs in topologically sorted order.
    """
    try:
        ordered_nodes = list(nx.topological_sort(G))

    except nx.NetworkXUnfeasible:
        # if it couldn't topologically sort the nodes, the way probably
        # contains a cycle. try removing an edge to break the cycle. first,
        # look for multiple edges emanating from the same source node
        insert_before = True
        edges = [
            edge
            for source in [node for node, degree in G.out_degree() if degree > 1]
            for edge in G.out_edges(source, keys=True)
        ]

        # if none found, then look for multiple edges pointing at the same
        # target node instead
        if len(edges) == 0:
            insert_before = False
            edges = [
                edge
                for target in [node for node, degree in G.in_degree() if degree > 1]
                for edge in G.in_edges(target, keys=True)
            ]

            # if still none, then take the first edge of the way: the entire
            # way could just be a cycle in which each node appears once
            if len(edges) == 0:
                edges = [next(iter(G.edges))]

        # remove one edge at a time and, if the graph remains connected, exit
        # the loop and check if we are able to topologically sort the nodes
        for edge in edges:
            G_ = G.copy()
            G_.remove_edge(*edge)
            if nx.is_weakly_connected(G_):
                break

        try:
            ordered_nodes = list(nx.topological_sort(G_))

            # re-insert (before or after its neighbor as needed) the duplicate
            # source or target node from the edge we removed
            dupe_node = edge[0] if insert_before else edge[1]
            neighbor = edge[1] if insert_before else edge[0]
            position = ordered_nodes.index(neighbor)
            position = position if insert_before else position + 1
            ordered_nodes.insert(position, dupe_node)

        except nx.NetworkXUnfeasible:
            # if it failed again, this way probably contains multiple cycles,
            # so remove a cycle then try to sort the nodes again, recursively.
            # note this is destructive and will be missing in the saved data.
            G_ = G.copy()
            G_.remove_edges_from(nx.find_cycle(G_))
            G_ = truncate.largest_component(G_)
            ordered_nodes = _sort_nodes(G_, osmid)
            msg = f"Had to remove a cycle from way {osmid!r} for topological sort"
            utils.log(msg, level=lg.WARNING)

    return ordered_nodes
