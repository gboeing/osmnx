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
    network_type: str | None,
) -> dict[str, Any]:
    """
    Read OSM PBF data from file and return Overpass-like JSON.

    Parameters
    ----------
    filepath
        Path to PBF file containing OSM data.
    network_type
        Network type preset used to filter ways, or None to retain all ways.

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
    way_node_ids: set[int] = set()

    # first pass, collect ways and their node IDs
    utils.log(f"Extracting ways from {str(filepath)!r}.")
    ways_fp = osmium.FileProcessor(filepath, entities=osmium.osm.WAY)
    if network_type is not None:
        ways_fp = ways_fp.with_filter(osmium.filter.KeyFilter("highway"))
    for way in ways_fp:
        tags = dict(way.tags)
        if network_type is not None and not _overpass._network_filter_matches(network_type, tags):
            continue

        node_ids = [node.ref for node in way.nodes]
        if len(node_ids) < 2:  # noqa: PLR2004
            continue

        ways.append(
            {
                "type": "way",
                "id": way.id,
                "nodes": node_ids,
                "tags": tags,
            },
        )
        way_node_ids.update(node_ids)

    # second pass, filter to nodes that constitute the preceding ways
    utils.log(f"Extracting {len(way_node_ids):,} way nodes from {str(filepath)!r}.")
    nodes_fp = osmium.FileProcessor(filepath, entities=osmium.osm.NODE)
    nodes_fp = nodes_fp.with_filter(osmium.filter.IdFilter(way_node_ids))
    nodes: list[dict[str, Any]] = []
    found_node_ids: set[int] = set()
    for node in nodes_fp:
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
    way_node_ids = {node_id for way in ways for node_id in way["nodes"]}
    # remove nodes that no retained way uses
    nodes = [node for node in nodes if node["id"] in way_node_ids]

    utils.log(f"Extracted {len(nodes):,} nodes and {len(ways):,} ways.")
    return {"elements": ways + nodes}


def graph_from_pbf(
    filepath: str | Path,
    *,
    network_type: str | None = "all",
    simplify: bool = True,
    retain_all: bool = False,
) -> nx.MultiDiGraph:
    """
    Create a graph from data in an OSM PBF file.

    Use the `settings` module's `useful_tags_node` and `useful_tags_way`
    settings to configure which OSM node/way tags are added as graph node/edge
    attributes. Filter ways with the selected `network_type` preset, or pass
    None to retain all ways in a PBF that was pre-filtered with external OSM
    tooling. Add any desired non-highway tags, such as "railway", to
    `settings.useful_tags_way` before loading. PBF presets that filter access
    use the default "private" rule. If `settings.default_access` is
    customized, this function warns and ignores it because PBF filtering
    cannot evaluate Overpass QL.

    Parameters
    ----------
    filepath
        Path to PBF file containing OSM data.
    network_type
        Network type preset to filter ways. Choose from "all", "all_public",
        "bike", "drive", "drive_service", or "walk". If None, retain all
        ways. This unfiltered path is intended for pre-filtered PBF files and
        can be slow for general regional extracts.
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
    if network_type is not None:
        if network_type not in _overpass._NETWORK_EXCLUDED_PATTERNS:
            msg = f"Unrecognized network_type {network_type!r}."
            raise ValueError(msg)

        excluded = _overpass._NETWORK_EXCLUDED_PATTERNS[network_type]
        if "access" in excluded and settings.default_access != _overpass._DEFAULT_OVERPASS_ACCESS:
            msg = (
                "PBF preset filtering ignores a customized "
                "settings.default_access and uses the default 'private' rule. "
                "To customize access filtering, pre-filter the PBF with external "
                "OSM tooling and pass network_type=None."
            )
            warn(msg, category=UserWarning, stacklevel=2)

    bidirectional = (
        network_type is not None and network_type in settings.bidirectional_network_types
    )

    # reuse the normal graph builder so PBF and Overpass graphs stay consistent
    response_json = [_overpass_json_from_pbf(filepath, network_type)]
    G = graph._create_graph(response_json, bidirectional)
    if not retain_all:
        G = truncate.largest_component(G, strongly=False)
    if simplify:
        G = simplification.simplify_graph(G)

    msg = f"graph_from_pbf returned graph with {len(G):,} nodes and {len(G.edges):,} edges"
    utils.log(msg, level=lg.INFO)
    return G
