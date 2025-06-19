"""File I/O functions to save/load graphs to/from files on disk."""

from __future__ import annotations

import ast
import contextlib
import logging as lg
from pathlib import Path
from typing import TYPE_CHECKING
from typing import Any

import networkx as nx
import pandas as pd
from shapely import wkt

from . import _osm_xml
from . import convert
from . import settings
from . import utils

if TYPE_CHECKING:
    import geopandas as gpd


def save_graph_geopackage(
    G: nx.MultiDiGraph,
    filepath: str | Path | None = None,
    *,
    directed: bool = False,
    encoding: str = "utf-8",
) -> None:
    """
    Save graph nodes and edges to disk as layers in a GeoPackage file.

    Parameters
    ----------
    G
        The graph to save.
    filepath
        Path to the GeoPackage file including extension. If None, use default
        `settings.data_folder/graph.gpkg`.
    directed
        If False, save one edge for each undirected edge in the graph but
        retain original oneway and to/from information as edge attributes. If
        True, save one edge for each directed edge in the graph.
    encoding
        The character encoding of the saved GeoPackage file.
    """
    # default filepath if none was provided
    filepath = Path(settings.data_folder) / "graph.gpkg" if filepath is None else Path(filepath)

    # if save folder does not already exist, create it
    filepath.parent.mkdir(parents=True, exist_ok=True)

    # convert graph to gdfs and stringify non-numeric columns
    if directed:
        gdf_nodes, gdf_edges = convert.graph_to_gdfs(G)
    else:
        gdf_nodes, gdf_edges = convert.graph_to_gdfs(convert.to_undirected(G))
    gdf_nodes = _stringify_nonnumeric_cols(gdf_nodes)
    gdf_edges = _stringify_nonnumeric_cols(gdf_edges)

    # save the nodes and edges as GeoPackage layers
    gdf_nodes.to_file(filepath, layer="nodes", driver="GPKG", index=True, encoding=encoding)
    gdf_edges.to_file(filepath, layer="edges", driver="GPKG", index=True, encoding=encoding)

    msg = f"Saved graph as GeoPackage at {str(filepath)!r}"
    utils.log(msg, level=lg.INFO)


def save_graphml(
    G: nx.MultiDiGraph,
    filepath: str | Path | None = None,
    *,
    gephi: bool = False,
    encoding: str = "utf-8",
) -> None:
    """
    Save graph to disk as GraphML file.

    Parameters
    ----------
    G
        The graph to save as.
    filepath
        Path to the GraphML file including extension. If None, use default
        `settings.data_folder/graph.graphml`.
    gephi
        If True, give each edge a unique key/id for compatibility with Gephi's
        interpretation of the GraphML specification.
    encoding
        The character encoding of the saved GraphML file.
    """
    # default filepath if none was provided
    filepath = Path(settings.data_folder) / "graph.graphml" if filepath is None else Path(filepath)

    # if save folder does not already exist, create it
    filepath.parent.mkdir(parents=True, exist_ok=True)

    # make a copy to not mutate original graph object caller passed in
    G = G.copy()

    if gephi:
        # for gephi compatibility, each edge's key must be unique as an id
        uvkd = [(u, v, k, d) for k, (u, v, d) in enumerate(G.edges(keys=False, data=True))]
        G.clear_edges()
        G.add_edges_from(uvkd)

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
    msg = f"Saved graph as GraphML file at {str(filepath)!r}"
    utils.log(msg, level=lg.INFO)


def load_graphml(
    filepath: str | Path | None = None,
    *,
    graphml_str: str | None = None,
    node_dtypes: dict[str, Any] | None = None,
    edge_dtypes: dict[str, Any] | None = None,
    graph_dtypes: dict[str, Any] | None = None,
) -> nx.MultiDiGraph:
    """
    Load an OSMnx-saved GraphML file from disk or GraphML string.

    This function converts node, edge, and graph-level attributes (serialized
    as strings) to their appropriate data types. These can be customized as
    needed by passing in dtypes arguments providing types or custom converter
    functions. For example, if you want to convert some attribute's values to
    `bool`, consider using the built-in `ox.io._convert_bool_string` function
    to properly handle "True"/"False" string literals as True/False booleans:
    `ox.load_graphml(fp, node_dtypes={my_attr: ox.io._convert_bool_string})`.

    If you manually configured the `all_oneway=True` setting, you may need to
    manually specify here that edge `oneway` attributes should be type `str`.

    Note that you must pass one and only one of `filepath` or `graphml_str`.
    If passing `graphml_str`, you may need to decode the bytes read from your
    file before converting to string to pass to this function.

    Parameters
    ----------
    filepath
        Path to the GraphML file.
    graphml_str
        Valid and decoded string representation of a GraphML file's contents.
    node_dtypes
        Dict of node attribute names:types to convert values' data types. The
        type can be a type or a custom string converter function.
    edge_dtypes
        Dict of edge attribute names:types to convert values' data types. The
        type can be a type or a custom string converter function.
    graph_dtypes
        Dict of graph-level attribute names:types to convert values' data
        types. The type can be a type or a custom string converter function.

    Returns
    -------
    G
        The loaded MultiDiGraph.
    """
    if (filepath is None and graphml_str is None) or (
        filepath is not None and graphml_str is not None
    ):  # pragma: no cover
        msg = "You must pass one and only one of `filepath` or `graphml_str`."
        raise ValueError(msg)

    # specify default graph/node/edge attribute values' data types
    default_graph_dtypes = {
        "consolidated": _convert_bool_string,
        "simplified": _convert_bool_string,
    }
    default_node_dtypes = {
        "elevation": float,
        "elevation_res": float,
        "osmid": int,
        "street_count": int,
        "x": float,
        "y": float,
    }
    default_edge_dtypes = {
        "bearing": float,
        "grade": float,
        "grade_abs": float,
        "length": float,
        "oneway": _convert_bool_string,
        "osmid": int,
        "reversed": _convert_bool_string,
        "speed_kph": float,
        "travel_time": float,
    }

    # override default graph/node/edge attr types with user-passed types, if any
    if graph_dtypes is not None:
        default_graph_dtypes.update(graph_dtypes)
    if node_dtypes is not None:
        default_node_dtypes.update(node_dtypes)
    if edge_dtypes is not None:
        default_edge_dtypes.update(edge_dtypes)

    if filepath is not None:
        # read the graphml file from disk
        source = filepath
        G = nx.read_graphml(
            Path(filepath),
            node_type=default_node_dtypes["osmid"],
            force_multigraph=True,
        )
    else:
        # parse the graphml string
        source = "string"
        G = nx.parse_graphml(
            graphml_str,
            node_type=default_node_dtypes["osmid"],
            force_multigraph=True,
        )

    # convert graph/node/edge attribute data types
    msg = "Converting node, edge, and graph-level attribute data types"
    utils.log(msg, level=lg.INFO)
    G = _convert_graph_attr_types(G, default_graph_dtypes)
    G = _convert_node_attr_types(G, default_node_dtypes)
    G = _convert_edge_attr_types(G, default_edge_dtypes)

    msg = f"Loaded graph with {len(G)} nodes and {len(G.edges)} edges from {str(source)!r}"
    utils.log(msg, level=lg.INFO)
    return G


def save_graph_xml(
    G: nx.MultiDiGraph,
    filepath: str | Path | None = None,
    *,
    way_tag_aggs: dict[str, Any] | None = None,
    encoding: str = "utf-8",
) -> None:
    """
    Save graph to disk as an OSM XML file.

    This function exists only to allow serialization to the OSM XML format
    for applications that require it, and has constraints to conform to that.
    As such, it has a limited use case which does not include saving/loading
    graphs for subsequent OSMnx analysis. To save/load graphs to/from disk for
    later use in OSMnx, use the `io.save_graphml` and `io.load_graphml`
    functions instead. To load a graph from an OSM XML file that you have
    downloaded or generated elsewhere, use the `graph.graph_from_xml`
    function.

    Use the `settings` module's `useful_tags_node` and `useful_tags_way`
    settings to configure which tags your graph is created and saved with.
    This function merges graph edges such that each OSM way has one entry in
    the XML output, with the way's nodes topologically sorted. `G` must be
    unsimplified to save as OSM XML: otherwise, one edge could comprise
    multiple OSM ways, making it impossible to group and sort edges in way.
    `G` should also have been created with `ox.settings.all_oneway=True` for
    this function to behave properly.

    Parameters
    ----------
    G
        Unsimplified, unprojected graph to save as an OSM XML file.
    filepath
        Path to the saved file including extension. If None, use default
        `settings.data_folder/graph.osm`.
    way_tag_aggs
        Keys are OSM way tag keys and values are aggregation functions
        (anything accepted as an argument by pandas.agg). Allows user to
        aggregate graph edge attribute values into single OSM way values. If
        None, or if some tag's key does not exist in the dict, the way
        attribute will be assigned the value of the first edge of the way.
    encoding
        The character encoding of the saved OSM XML file.
    """
    _osm_xml._save_graph_xml(G, filepath, way_tag_aggs, encoding)


def _convert_graph_attr_types(G: nx.MultiDiGraph, dtypes: dict[str, Any]) -> nx.MultiDiGraph:
    """
    Convert graph-level attributes using a dict of data types.

    Parameters
    ----------
    G
        Graph to convert the graph-level attributes of.
    dtypes
        Dict of graph-level attribute names:types.

    Returns
    -------
    G
        The graph with its graph-level attributes' types converted.
    """
    # remove node_default and edge_default metadata keys if they exist
    G.graph.pop("node_default", None)
    G.graph.pop("edge_default", None)

    for attr in G.graph.keys() & dtypes.keys():
        G.graph[attr] = dtypes[attr](G.graph[attr])

    return G


def _convert_node_attr_types(G: nx.MultiDiGraph, dtypes: dict[str, Any]) -> nx.MultiDiGraph:
    """
    Convert graph nodes' attributes using a dict of data types.

    Parameters
    ----------
    G
        Graph to convert the node attributes of.
    dtypes
        Dict of node attribute names:types.

    Returns
    -------
    G
        The graph with its nodes' attributes' types converted.
    """
    for _, data in G.nodes(data=True):
        # first, eval stringified lists, dicts, or sets to convert them to objects
        # lists, dicts, or sets would be custom attribute types added by a user
        for attr, value in data.items():
            if (value.startswith("[") and value.endswith("]")) or (
                value.startswith("{") and value.endswith("}")
            ):
                with contextlib.suppress(SyntaxError, ValueError):
                    data[attr] = ast.literal_eval(value)

        for attr in data.keys() & dtypes.keys():
            data[attr] = dtypes[attr](data[attr])
    return G


def _convert_edge_attr_types(G: nx.MultiDiGraph, dtypes: dict[str, Any]) -> nx.MultiDiGraph:
    """
    Convert graph edges' attributes using a dict of data types.

    Parameters
    ----------
    G
        Graph to convert the edge attributes of.
    dtypes
        Dict of edge attribute names:types.

    Returns
    -------
    G
        The graph with its edges' attributes' types converted.
    """
    # for each edge in the graph, eval attribute value lists and convert types
    for _, _, data in G.edges(data=True, keys=False):
        # remove extraneous "id" attribute added by graphml saving
        data.pop("id", None)

        # first, eval stringified lists, dicts, or sets to convert them to objects
        # edge attributes might have a single value, or a list if simplified
        # dicts or sets would be custom attribute types added by a user
        for attr, value in data.items():
            if (value.startswith("[") and value.endswith("]")) or (
                value.startswith("{") and value.endswith("}")
            ):
                with contextlib.suppress(SyntaxError, ValueError):
                    data[attr] = ast.literal_eval(value)

        # next, convert attribute value types if attribute appears in dtypes
        for attr in data.keys() & dtypes.keys():
            if isinstance(data[attr], list):
                # if it's a list, eval it then convert each item
                data[attr] = [dtypes[attr](item) for item in data[attr]]
            else:
                # otherwise, just convert the single value
                data[attr] = dtypes[attr](data[attr])

        # if "geometry" attr exists, convert its well-known text to LineString
        if "geometry" in data:
            data["geometry"] = wkt.loads(data["geometry"])

    return G


def _convert_bool_string(value: bool | str) -> bool:  # noqa: FBT001
    """
    Convert a "True" or "False" string literal to corresponding boolean type.

    This is necessary because Python will otherwise parse the string "False"
    to the boolean value True, that is, `bool("False") == True`. This function
    raises a ValueError if a value other than "True" or "False" is passed.

    If the value is already a boolean, this function just returns it, to
    accommodate usage when the value was originally inside a stringified list.

    Parameters
    ----------
    value
        The string to convert to bool.

    Returns
    -------
    bool_value
        The boolean equivalent of the string literal.
    """
    if isinstance(value, bool):
        return value

    if value in {"True", "False"}:
        return value == "True"

    # otherwise the value is not a valid boolean
    msg = f"Invalid literal for boolean: {value!r}."
    raise ValueError(msg)


def _stringify_nonnumeric_cols(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Make every non-numeric GeoDataFrame column (besides geometry) a string.

    This allows proper serializing via Fiona of GeoDataFrames with mixed types
    such as strings and ints in the same column.

    Parameters
    ----------
    gdf
        GeoDataFrame to stringify non-numeric columns of.

    Returns
    -------
    gdf
        GeoDataFrame with non-numeric columns stringified.
    """
    # stringify every non-numeric column other than geometry column
    for col in (c for c in gdf.columns if c != "geometry"):
        if not pd.api.types.is_numeric_dtype(gdf[col]):
            gdf[col] = gdf[col].fillna("").astype(str)

    return gdf
