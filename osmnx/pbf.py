"""Load OSM PBF files as graphs.

For file format information see https://wiki.openstreetmap.org/wiki/PBF_Format
"""

from __future__ import annotations

import logging as lg
from typing import TYPE_CHECKING
from typing import Any
from warnings import warn

# keep osmium optional so importing osmnx does not require the PBF extra
try:
    import osmium
except ImportError:  # pragma: no cover
    osmium = None

from . import _overpass
from . import graph
from . import settings
from . import simplification
from . import truncate
from . import utils

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path

    import networkx as nx


def _remove_incomplete_ways(
    ways: list[dict[str, Any]],
    found_node_ids: set[int],
) -> list[dict[str, Any]]:
    """
    Remove ways with node references that were not found in the PBF.

    Parameters
    ----------
    ways
        Way dictionaries containing node references.
    found_node_ids
        Node IDs found in the PBF.

    Returns
    -------
    complete_ways
        Way dictionaries whose node references were all found.
    """
    complete_ways = []
    incomplete_way_count = 0
    for way in ways:
        if all(node_id in found_node_ids for node_id in way["nodes"]):
            complete_ways.append(way)
        else:
            incomplete_way_count += 1

    if incomplete_way_count > 0:
        msg = (
            "Removed incomplete ways because their node references were missing: "
            f"{incomplete_way_count:,}. This likely resulted from clipped PBF input."
        )
        warn(msg, category=UserWarning, stacklevel=2)
    return complete_ways


def _overpass_json_from_pbf(
    filepath: str | Path,
    network_type: str,
    way_filter: Callable[[dict[str, str]], bool] | None,
) -> dict[str, Any]:
    """
    Read OSM PBF data from file and return Overpass-like JSON.

    Parameters
    ----------
    filepath
        Path to PBF file containing OSM data.
    network_type
        Network type preset used to filter ways when `way_filter` is None.
    way_filter
        Callable used to filter ways, or None to use `network_type`.

    Returns
    -------
    response_json
        A parsed JSON response like from the Overpass API.
    """
    # fail only when the caller actually tries to read a PBF
    if osmium is None:  # pragma: no cover
        msg = "PBF support requires the optional dependency 'osmium'. Install it with 'osmnx[pbf]'."
        raise ImportError(msg)

    ways: list[dict[str, Any]] = []
    way_nodes: set[int] = set()

    # first pass, filter ways before collecting their node refs
    utils.log(f"Extracting ways from {str(filepath)!r}.")
    fp = osmium.FileProcessor(filepath, entities=osmium.osm.WAY)
    for way in fp:
        tags = dict(way.tags)
        if way_filter is not None:
            if not way_filter(tags):
                continue
        elif not _overpass._network_filter_matches(network_type, tags):
            continue

        node_refs = [node.ref for node in way.nodes]
        if len(node_refs) < 2:  # noqa: PLR2004
            continue

        ways.append(
            {
                "type": "way",
                "id": way.id,
                "nodes": node_refs,
                "tags": tags,
            },
        )
        way_nodes.update(node_refs)

    # second pass, filter to nodes that constitute the preceding ways
    utils.log(f"Extracting {len(way_nodes):,} way nodes from {str(filepath)!r}.")
    fp = osmium.FileProcessor(filepath, entities=osmium.osm.NODE)
    fp = fp.with_filter(osmium.filter.IdFilter(way_nodes))
    nodes: list[dict[str, Any]] = []
    found_node_ids: set[int] = set()
    for node in fp:
        nodes.append(
            {
                "type": "node",
                "id": node.id,
                "lon": node.location.lon,
                "lat": node.location.lat,
                "tags": dict(node.tags),
            },
        )
        found_node_ids.add(node.id)

    # discard incomplete ways before handing data to graph construction
    ways = _remove_incomplete_ways(ways, found_node_ids)
    way_nodes = {node_id for way in ways for node_id in way["nodes"]}
    # remove nodes that no retained way uses
    nodes = [node for node in nodes if node["id"] in way_nodes]
    elements = ways + nodes

    utils.log(f"Extracted {len(nodes):,} nodes and {len(ways):,} ways.")
    return {"elements": elements}


def graph_from_pbf(
    filepath: str | Path,
    *,
    network_type: str = "all",
    way_filter: Callable[[dict[str, str]], bool] | None = None,
    bidirectional: bool | None = None,
    simplify: bool = True,
    retain_all: bool = False,
) -> nx.MultiDiGraph:
    """
    Create a graph from data in an OSM PBF file.

    Use the `settings` module's `useful_tags_node` and `useful_tags_way`
    settings to configure which OSM node/way tags are added as graph node/edge
    attributes. When `way_filter` is None, filter ways with the selected
    `network_type` preset. Otherwise, call `way_filter` once for each way with
    a new dictionary containing that way's tags and retain the way when the
    return value is truthy. Arbitrary Overpass `custom_filter` strings are not
    supported for PBF files. PBF presets use the default "private" access
    rule. If `settings.default_access` is customized, this function warns and
    ignores it because PBF filtering cannot evaluate Overpass QL.

    Parameters
    ----------
    filepath
        Path to PBF file containing OSM data.
    network_type
        Network type preset to filter ways. Choose from "all", "all_public",
        "bike", "drive", "drive_service", or "walk". Ignored when
        `way_filter` is provided, but still determines default bidirectionality.
    way_filter
        Callable that receives a new dictionary of each way's tags and returns
        whether to retain the way. If provided, it takes precedence over
        `network_type`.
    bidirectional
        If explicitly True or False, honor that value. If None, derive it
        from `network_type`.
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
    if network_type not in _overpass._NETWORK_REQUIRED_TAGS:
        msg = f"Unrecognized network_type {network_type!r}."
        raise ValueError(msg)

    excluded = _overpass._NETWORK_EXCLUDED_PATTERNS.get(network_type, {})
    access = excluded.get("access")
    if (
        way_filter is None
        and access == _overpass._DEFAULT_PBF_ACCESS
        and settings.default_access != _overpass._DEFAULT_OVERPASS_ACCESS
    ):
        # warn once here because the matcher runs once per PBF way
        msg = (
            "PBF preset filtering ignores a customized "
            "settings.default_access and uses the default 'private' rule. "
            "Provide way_filter to customize access filtering."
        )
        warn(msg, category=UserWarning, stacklevel=2)

    if bidirectional is None:
        bidirectional = network_type in settings.bidirectional_network_types

    # reuse the normal graph builder so PBF and Overpass graphs stay consistent
    response_json = [_overpass_json_from_pbf(filepath, network_type, way_filter)]
    G = graph._create_graph(response_json, bidirectional)
    if not retain_all:
        G = truncate.largest_component(G, strongly=False)
    if simplify:
        G = simplification.simplify_graph(G)

    msg = f"graph_from_pbf returned graph with {len(G):,} nodes and {len(G.edges):,} edges"
    utils.log(msg, level=lg.INFO)
    return G
