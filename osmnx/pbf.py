"""
Load OSM PBF files as graphs or features GeoDataFrames.

For file format information see https://wiki.openstreetmap.org/wiki/PBF_Format
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING
from typing import Any

import osmium

from . import graph
from . import simplification
from . import truncate
from . import utils

if TYPE_CHECKING:
    import networkx as nx


def filter_pbf(filepath: str | Path, tags: list[str] | dict[str, list[str]]) -> Path:
    """
    Filter a PBF file based on way tags.

    Parameters
    ----------
    filepath
        Path to PBF file containing OSM data.
    tags
        Tags to filter the elements in the file. If `list`, then only retain
        ways that contain a tag key from the list. If `dict`, then only retain
        ways that contain a tag key from the dict with a matching value from
        that dict key's list of values.

    Returns
    -------
    output_filepath
        Filepath to new filtered PBF file.
    """
    utils.log(f"Filtering PBF file at {str(filepath)!r}...")
    type_err_msg = "`tags` must be a list, or a dict with string keys and list values"

    # make output filepath and delete if it already exists
    filepath = Path(filepath)
    output_filepath = filepath.parent / f"filtered-{filepath.name}"
    if output_filepath.is_file():
        output_filepath.unlink()

    with osmium.SimpleWriter(output_filepath) as f:
        # filter ways based on tags
        fp = osmium.FileProcessor(filepath)
        fp = fp.with_filter(osmium.filter.EntityFilter(osmium.osm.WAY))

        if isinstance(tags, list):
            # if user passed list, only retain ways with those tag keys
            utils.log(f"Filtering ways by tag keys {tags!r}.")
            fp = fp.with_filter(osmium.filter.KeyFilter(*tags))

        elif isinstance(tags, dict):
            # if user passed dict, only retain ways tagged with those key:value pairs
            tag_filters: list[tuple[str, str]] = []
            for tag_key, tag_values in tags.items():
                if not isinstance(tag_values, list):
                    raise TypeError(type_err_msg)
                tag_filters.extend(zip([tag_key] * len(tag_values), tag_values, strict=False))
            utils.log(f"Filtering ways by tag values {tag_filters}.")
            fp = fp.with_filter(osmium.filter.TagFilter(*tag_filters))

        else:
            raise TypeError(type_err_msg)

        way_nodes: list[int] = []
        for way in fp:
            way_nodes.extend(node.ref for node in way.nodes)
            f.add_way(way)

        # filter nodes to retain only those that make up the ways
        way_nodes_unique = set(way_nodes)
        utils.log(f"Filtering way nodes by {len(way_nodes_unique):,} IDs.")
        fp = osmium.FileProcessor(filepath)
        fp = fp.with_filter(osmium.filter.EntityFilter(osmium.osm.NODE))
        fp = fp.with_filter(osmium.filter.IdFilter(way_nodes_unique))

        for node in fp:
            f.add_node(node)

    utils.log(f"Saved filtered PBF file at {str(output_filepath)!r}")
    return output_filepath


def overpass_json_from_pbf(
    filepath: str | Path,
) -> dict[str, Any]:
    """
    Read OSM PBF data from file and return Overpass-like JSON.

    Parameters
    ----------
    filepath
        Path to PBF file containing OSM data.

    Returns
    -------
    response_json
        A parsed JSON response like from the Overpass API.
    """
    elements = []

    # first pass, extract ways and track their node refs
    utils.log(f"Extracting ways from {str(filepath)!r}.")
    fp = osmium.FileProcessor(filepath)
    fp = fp.with_filter(osmium.filter.EntityFilter(osmium.osm.WAY))
    for way in fp:
        way_dict = {
            "type": "way",
            "id": way.id,
            "nodes": [node.ref for node in way.nodes],
            "tags": dict(way.tags),
        }
        elements.append(way_dict)

    # count the ways and extract their constituent nodes' IDs
    way_count = len(elements)
    way_nodes = {n for e in elements for n in e["nodes"]}

    # second pass, filter to nodes that constitute the preceding ways
    utils.log(f"Extracting {len(way_nodes):,} way nodes from {str(filepath)!r}.")
    fp = osmium.FileProcessor(filepath)
    fp = fp.with_filter(osmium.filter.EntityFilter(osmium.osm.NODE))
    fp = fp.with_filter(osmium.filter.IdFilter(way_nodes))
    for node in fp:
        node_dict = {
            "type": "node",
            "id": node.id,
            "lon": node.location.lon,
            "lat": node.location.lat,
            "tags": dict(node.tags),
        }
        elements.append(node_dict)

    node_count = len(elements) - way_count
    utils.log(f"Extracted {node_count:,} nodes and {way_count:,} ways.")
    return {"elements": elements}


def graph_from_pbf(
    filepath: str | Path,
    tags: list[str] | dict[str, list[str]] | None = None,
    *,
    bidirectional: bool = False,
    simplify: bool = True,
    retain_all: bool = False,
) -> nx.MultiDiGraph:
    """
    Create a graph from data in an OSM PBF file.

    Use the `settings` module's `useful_tags_node` and `useful_tags_way`
    settings to configure which OSM node/way tags are added as graph node/edge
    attributes.

    Parameters
    ----------
    filepath
        Path to PBF file containing OSM data.
    tags
        Tags to filter the elements in the file. If `list`, then only retain
        ways that contain a tag key from the list. If `dict`, then only retain
        ways that contain a tag key from the dict with a matching value from
        that dict key's list of values.
    bidirectional
        If True, create bidirectional edges for one-way streets.
    simplify
        If True, simplify graph topology with the `simplify_graph` function.
    retain_all
        If True, return the entire graph even if it is not connected. If
        False, retain only the largest weakly connected component.

    Returns
    -------
    G
        The resulting MultiDiGraph.
    """
    if tags is not None:
        filepath = filter_pbf(filepath, tags)
    response_json = [overpass_json_from_pbf(filepath)]
    G = graph._create_graph(response_json, bidirectional)
    if not retain_all:
        G = truncate.largest_component(G, strongly=False)
    if simplify:
        G = simplification.simplify_graph(G)
    return G
