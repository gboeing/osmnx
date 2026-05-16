"""Convert local OSM PBF data to Overpass-like JSON."""

from __future__ import annotations

import logging as lg
from typing import TYPE_CHECKING
from typing import Any

from . import utils

if TYPE_CHECKING:
    from pathlib import Path

    from ._osm_filters import _WayFilter

# osmium is an optional dependency for PBF querying
try:
    import osmium
except ImportError:  # pragma: no cover
    osmium = None


def _overpass_json_from_pbf(filepath: Path, way_filter: _WayFilter) -> dict[str, Any]:
    """
    Read OSM PBF network data from file and return Overpass-like JSON.

    Parameters
    ----------
    filepath
        Path to PBF file containing OSM data.
    way_filter
        Filter specifying which OSM ways to include.

    Returns
    -------
    response_json
        A parsed JSON response like from the Overpass API.
    """
    way_elements, node_refs = _extract_ways(filepath, way_filter)
    node_elements, found_node_refs = _extract_nodes(filepath, node_refs)
    way_elements = _discard_ways_missing_nodes(way_elements, found_node_refs)

    used_node_refs = {node_id for way in way_elements for node_id in way["nodes"]}
    node_elements = [node for node in node_elements if node["id"] in used_node_refs]

    utils.log(f"Extracted {len(node_elements):,} nodes and {len(way_elements):,} ways.")
    return {"elements": [*node_elements, *way_elements]}


def _extract_ways(filepath: Path, way_filter: _WayFilter) -> tuple[list[dict[str, Any]], set[int]]:
    """Extract matching ways from a PBF file and collect their node references."""
    utils.log(f"Extracting ways from {str(filepath)!r}.")
    way_elements: list[dict[str, Any]] = []
    node_refs: set[int] = set()

    if osmium is None:  # pragma: no cover
        msg = "PBF support requires the optional dependency `osmium`: install `osmnx[pbf]`."
        raise ImportError(msg)

    fp = osmium.FileProcessor(str(filepath))
    fp = fp.with_filter(osmium.filter.EntityFilter(osmium.osm.WAY))
    for way in fp:
        tags = dict(way.tags)
        if way_filter.matches(tags):
            way_node_refs = [node.ref for node in way.nodes]
            if len(way_node_refs) > 1:
                node_refs.update(way_node_refs)
                way_elements.append(
                    {"type": "way", "id": way.id, "nodes": way_node_refs, "tags": tags},
                )

    return way_elements, node_refs


def _extract_nodes(filepath: Path, node_refs: set[int]) -> tuple[list[dict[str, Any]], set[int]]:
    """Extract the requested nodes from a PBF file."""
    utils.log(f"Extracting {len(node_refs):,} way nodes from {str(filepath)!r}.")
    node_elements: list[dict[str, Any]] = []
    found_node_refs: set[int] = set()

    if osmium is None:  # pragma: no cover
        msg = "PBF support requires the optional dependency `osmium`: install `osmnx[pbf]`."
        raise ImportError(msg)

    fp = osmium.FileProcessor(str(filepath))
    fp = fp.with_filter(osmium.filter.EntityFilter(osmium.osm.NODE))
    fp = fp.with_filter(osmium.filter.IdFilter(node_refs))
    for node in fp:
        if node.location.valid():
            found_node_refs.add(node.id)
            node_elements.append(
                {
                    "type": "node",
                    "id": node.id,
                    "lon": node.location.lon,
                    "lat": node.location.lat,
                    "tags": dict(node.tags),
                },
            )

    return node_elements, found_node_refs


def _discard_ways_missing_nodes(
    way_elements: list[dict[str, Any]],
    found_node_refs: set[int],
) -> list[dict[str, Any]]:
    """Discard ways with missing node references to avoid invalid graph edges."""
    valid_way_elements = [
        way for way in way_elements if all(node_id in found_node_refs for node_id in way["nodes"])
    ]
    discarded_count = len(way_elements) - len(valid_way_elements)

    if discarded_count > 0:
        msg = (
            f"Discarded {discarded_count:,} way(s) with missing node references, likely due to "
            "clipped PBF input data."
        )
        utils.log(msg, level=lg.WARNING)

    return valid_way_elements
