"""
Download and create graphs from OpenStreetMap data.

Refer to the Getting Started guide for usage limitations.
"""

from __future__ import annotations

import logging as lg
from collections.abc import Iterable
from itertools import groupby
from pathlib import Path
from typing import TYPE_CHECKING
from typing import Any

import networkx as nx
from shapely import MultiPolygon
from shapely import Polygon

from . import _osm_xml
from . import _overpass
from . import distance
from . import geocoder
from . import projection
from . import settings
from . import simplification
from . import stats
from . import truncate
from . import utils
from . import utils_geo
from ._errors import CacheOnlyInterruptError
from ._errors import InsufficientResponseError
from ._version import __version__

if TYPE_CHECKING:
    from collections.abc import Iterable
    from pathlib import Path


def graph_from_bbox(
    bbox: tuple[float, float, float, float],
    *,
    network_type: str = "all",
    simplify: bool = True,
    retain_all: bool = False,
    truncate_by_edge: bool = False,
    custom_filter: str | list[str] | None = None,
) -> nx.MultiDiGraph:
    """
    Download and create a graph within a lat-lon bounding box.

    This function uses filters to query the Overpass API: you can either
    specify a pre-defined `network_type` or provide your own `custom_filter`
    with Overpass QL.

    Use the `settings` module's `useful_tags_node` and `useful_tags_way`
    settings to configure which OSM node/way tags are added as graph node/edge
    attributes. You can also use the `settings` module to retrieve a snapshot
    of historical OSM data as of a certain date, or to configure the Overpass
    server timeout, memory allocation, and other custom settings.

    Parameters
    ----------
    bbox
        Bounding box as `(left, bottom, right, top)`. Coordinates should be in
        unprojected latitude-longitude degrees (EPSG:4326).
    network_type
        {"all", "all_public", "bike", "drive", "drive_service", "walk"}
        What type of street network to retrieve if `custom_filter` is None.
    simplify
        If True, simplify graph topology via the `simplify_graph` function.
    retain_all
        If True, return the entire graph even if it is not connected. If
        False, retain only the largest weakly connected component.
    truncate_by_edge
        If True, retain nodes outside bounding box if at least one of node's
        neighbors is within the bounding box.
    custom_filter
        A custom ways filter to be used instead of the `network_type` presets,
        e.g. `'["power"~"line"]' or '["highway"~"motorway|trunk"]'`. If `str`,
        the intersection of keys/values will be used, e.g., `'[maxspeed=50][lanes=2]'`
        will return all ways having both maxspeed of 50 and two lanes. If
        `list`, the union of the `list` items will be used, e.g.,
        `['[maxspeed=50]', '[lanes=2]']` will return all ways having either
        maximum speed of 50 or two lanes. Also pass in a `network_type` that
        is in `settings.bidirectional_network_types` if you want the graph to
        be fully bidirectional.

    Returns
    -------
    G

    Notes
    -----
    Very large query areas use the `utils_geo._consolidate_subdivide_geometry`
    function to automatically make multiple requests: see that function's
    documentation for caveats.
    """
    # convert bounding box to a polygon
    polygon = utils_geo.bbox_to_poly(bbox)

    # create graph using this polygon geometry
    G = graph_from_polygon(
        polygon,
        network_type=network_type,
        simplify=simplify,
        retain_all=retain_all,
        truncate_by_edge=truncate_by_edge,
        custom_filter=custom_filter,
    )

    msg = f"graph_from_bbox returned graph with {len(G):,} nodes and {len(G.edges):,} edges"
    utils.log(msg, level=lg.INFO)
    return G


def graph_from_point(
    center_point: tuple[float, float],
    dist: float,
    *,
    dist_type: str = "bbox",
    network_type: str = "all",
    simplify: bool = True,
    retain_all: bool = False,
    truncate_by_edge: bool = False,
    custom_filter: str | list[str] | None = None,
) -> nx.MultiDiGraph:
    """
    Download and create a graph within some distance of a lat-lon point.

    This function uses filters to query the Overpass API: you can either
    specify a pre-defined `network_type` or provide your own `custom_filter`
    with Overpass QL.

    Use the `settings` module's `useful_tags_node` and `useful_tags_way`
    settings to configure which OSM node/way tags are added as graph node/edge
    attributes. You can also use the `settings` module to retrieve a snapshot
    of historical OSM data as of a certain date, or to configure the Overpass
    server timeout, memory allocation, and other custom settings.

    Parameters
    ----------
    center_point
        The `(lat, lon)` center point around which to construct the graph.
        Coordinates should be in unprojected latitude-longitude degrees
        (EPSG:4326).
    dist
        Retain only those nodes within this many meters of `center_point`,
        measuring distance according to `dist_type`.
    dist_type
        {"bbox", "network"}
        If "bbox", retain only those nodes within a bounding box of `dist`
        length/width. If "network", retain only those nodes within `dist`
        network distance of the nearest node to `center_point`.
    network_type
        {"all", "all_public", "bike", "drive", "drive_service", "walk"}
        What type of street network to retrieve if `custom_filter` is None.
    simplify
        If True, simplify graph topology with the `simplify_graph` function.
    retain_all
        If True, return the entire graph even if it is not connected. If
        False, retain only the largest weakly connected component.
    truncate_by_edge
        If True, retain nodes outside bounding box if at least one of node's
        neighbors is within the bounding box.
    custom_filter
        A custom ways filter to be used instead of the `network_type` presets,
        e.g. `'["power"~"line"]' or '["highway"~"motorway|trunk"]'`. If `str`,
        the intersection of keys/values will be used, e.g., `'[maxspeed=50][lanes=2]'`
        will return all ways having both maxspeed of 50 and two lanes. If
        `list`, the union of the `list` items will be used, e.g.,
        `['[maxspeed=50]', '[lanes=2]']` will return all ways having either
        maximum speed of 50 or two lanes. Also pass in a `network_type` that
        is in `settings.bidirectional_network_types` if you want the graph to
        be fully bidirectional.

    Returns
    -------
    G

    Notes
    -----
    Very large query areas use the `utils_geo._consolidate_subdivide_geometry`
    function to automatically make multiple requests: see that function's
    documentation for caveats.
    """
    if dist_type not in {"bbox", "network"}:  # pragma: no cover
        msg = "`dist_type` must be 'bbox' or 'network'."
        raise ValueError(msg)

    # create bounding box from center point and distance in each direction
    bbox = utils_geo.bbox_from_point(center_point, dist)

    # create a graph from the bounding box
    G = graph_from_bbox(
        bbox,
        network_type=network_type,
        simplify=simplify,
        retain_all=retain_all,
        truncate_by_edge=truncate_by_edge,
        custom_filter=custom_filter,
    )

    if dist_type == "network":
        # find node nearest to center then truncate graph by dist from it
        node = distance.nearest_nodes(G, X=center_point[1], Y=center_point[0])
        G = truncate.truncate_graph_dist(G, node, dist)

    msg = f"graph_from_point returned graph with {len(G):,} nodes and {len(G.edges):,} edges"
    utils.log(msg, level=lg.INFO)
    return G


def graph_from_address(
    address: str,
    dist: float,
    *,
    dist_type: str = "bbox",
    network_type: str = "all",
    simplify: bool = True,
    retain_all: bool = False,
    truncate_by_edge: bool = False,
    custom_filter: str | list[str] | None = None,
) -> nx.MultiDiGraph | tuple[nx.MultiDiGraph, tuple[float, float]]:
    """
    Download and create a graph within some distance of an address.

    This function uses filters to query the Overpass API: you can either
    specify a pre-defined `network_type` or provide your own `custom_filter`
    with Overpass QL.

    Use the `settings` module's `useful_tags_node` and `useful_tags_way`
    settings to configure which OSM node/way tags are added as graph node/edge
    attributes. You can also use the `settings` module to retrieve a snapshot
    of historical OSM data as of a certain date, or to configure the Overpass
    server timeout, memory allocation, and other custom settings.

    Parameters
    ----------
    address
        The address to geocode and use as the central point around which to
        construct the graph.
    dist
        Retain only those nodes within this many meters of `center_point`,
        measuring distance according to `dist_type`.
    dist_type
        {"network", "bbox"}
        If "bbox", retain only those nodes within a bounding box of `dist`. If
        "network", retain only those nodes within `dist` network distance from
        the centermost node.
    network_type
        {"all", "all_public", "bike", "drive", "drive_service", "walk"}
        What type of street network to retrieve if `custom_filter` is None.
    simplify
        If True, simplify graph topology with the `simplify_graph` function.
    retain_all
        If True, return the entire graph even if it is not connected. If
        False, retain only the largest weakly connected component.
    truncate_by_edge
        If True, retain nodes outside bounding box if at least one of node's
        neighbors is within the bounding box.
    custom_filter
        A custom ways filter to be used instead of the `network_type` presets,
        e.g. `'["power"~"line"]' or '["highway"~"motorway|trunk"]'`. If `str`,
        the intersection of keys/values will be used, e.g., `'[maxspeed=50][lanes=2]'`
        will return all ways having both maxspeed of 50 and two lanes. If
        `list`, the union of the `list` items will be used, e.g.,
        `['[maxspeed=50]', '[lanes=2]']` will return all ways having either
        maximum speed of 50 or two lanes. Also pass in a `network_type` that
        is in `settings.bidirectional_network_types` if you want the graph to
        be fully bidirectional.

    Returns
    -------
    G or (G, (lat, lon))

    Notes
    -----
    Very large query areas use the `utils_geo._consolidate_subdivide_geometry`
    function to automatically make multiple requests: see that function's
    documentation for caveats.
    """
    # geocode the address string to a (lat, lon) point
    point = geocoder.geocode(address)

    # then create a graph from this point
    G = graph_from_point(
        point,
        dist,
        dist_type=dist_type,
        network_type=network_type,
        simplify=simplify,
        retain_all=retain_all,
        truncate_by_edge=truncate_by_edge,
        custom_filter=custom_filter,
    )

    msg = f"graph_from_address returned graph with {len(G):,} nodes and {len(G.edges):,} edges"
    utils.log(msg, level=lg.INFO)
    return G


def graph_from_place(
    query: str | dict[str, str] | list[str | dict[str, str]],
    *,
    network_type: str = "all",
    simplify: bool = True,
    retain_all: bool = False,
    truncate_by_edge: bool = False,
    which_result: int | None | list[int | None] = None,
    custom_filter: str | list[str] | None = None,
) -> nx.MultiDiGraph:
    """
    Download and create a graph within the boundaries of some place(s).

    The query must be geocodable and OSM must have polygon boundaries for the
    geocode result. If OSM does not have a polygon for this place, you can
    instead get its street network using the `graph_from_address` function,
    which geocodes the place name to a point and gets the network within some
    distance of that point.

    If OSM does have polygon boundaries for this place but you're not finding
    it, try to vary the query string, pass in a structured query dict, or vary
    the `which_result` argument to use a different geocode result. If you know
    the OSM ID of the place, you can retrieve its boundary polygon using the
    `geocode_to_gdf` function, then pass it to the `features_from_polygon`
    function.

    This function uses filters to query the Overpass API: you can either
    specify a pre-defined `network_type` or provide your own `custom_filter`
    with Overpass QL.

    Use the `settings` module's `useful_tags_node` and `useful_tags_way`
    settings to configure which OSM node/way tags are added as graph node/edge
    attributes. You can also use the `settings` module to retrieve a snapshot
    of historical OSM data as of a certain date, or to configure the Overpass
    server timeout, memory allocation, and other custom settings.

    Parameters
    ----------
    query
        The query or queries to geocode to retrieve place boundary polygon(s).
    network_type
        {"all", "all_public", "bike", "drive", "drive_service", "walk"}
        What type of street network to retrieve if `custom_filter` is None.
    simplify
        If True, simplify graph topology with the `simplify_graph` function.
    retain_all
        If True, return the entire graph even if it is not connected. If
        False, retain only the largest weakly connected component.
    truncate_by_edge
        If True, retain nodes outside bounding box if at least one of node's
        neighbors is within the bounding box.
    which_result
        which geocoding result to use. if None, auto-select the first
        (Multi)Polygon or raise an error if OSM doesn't return one.
    custom_filter
        A custom ways filter to be used instead of the `network_type` presets,
        e.g. `'["power"~"line"]' or '["highway"~"motorway|trunk"]'`. If `str`,
        the intersection of keys/values will be used, e.g., `'[maxspeed=50][lanes=2]'`
        will return all ways having both maxspeed of 50 and two lanes. If
        `list`, the union of the `list` items will be used, e.g.,
        `['[maxspeed=50]', '[lanes=2]']` will return all ways having either
        maximum speed of 50 or two lanes. Also pass in a `network_type` that
        is in `settings.bidirectional_network_types` if you want the graph to
        be fully bidirectional.

    Returns
    -------
    G

    Notes
    -----
    Very large query areas use the `utils_geo._consolidate_subdivide_geometry`
    function to automatically make multiple requests: see that function's
    documentation for caveats.
    """
    # extract the geometry from the GeoDataFrame to use in query
    polygon = geocoder.geocode_to_gdf(query, which_result=which_result).union_all()
    msg = "Constructed place geometry polygon(s) to query Overpass"
    utils.log(msg, level=lg.INFO)

    # create graph using this polygon(s) geometry
    G = graph_from_polygon(
        polygon,
        network_type=network_type,
        simplify=simplify,
        retain_all=retain_all,
        truncate_by_edge=truncate_by_edge,
        custom_filter=custom_filter,
    )

    msg = f"graph_from_place returned graph with {len(G):,} nodes and {len(G.edges):,} edges"
    utils.log(msg, level=lg.INFO)
    return G


def graph_from_polygon(
    polygon: Polygon | MultiPolygon,
    *,
    network_type: str = "all",
    simplify: bool = True,
    retain_all: bool = False,
    truncate_by_edge: bool = False,
    custom_filter: str | list[str] | None = None,
) -> nx.MultiDiGraph:
    """
    Download and create a graph within the boundaries of a (Multi)Polygon.

    This function uses filters to query the Overpass API: you can either
    specify a pre-defined `network_type` or provide your own `custom_filter`
    with Overpass QL.

    Use the `settings` module's `useful_tags_node` and `useful_tags_way`
    settings to configure which OSM node/way tags are added as graph node/edge
    attributes. You can also use the `settings` module to retrieve a snapshot
    of historical OSM data as of a certain date, or to configure the Overpass
    server timeout, memory allocation, and other custom settings.

    Parameters
    ----------
    polygon
        The geometry within which to construct the graph. Coordinates should
        be in unprojected latitude-longitude degrees (EPSG:4326).
    network_type
        {"all", "all_public", "bike", "drive", "drive_service", "walk"}
        What type of street network to retrieve if `custom_filter` is None.
    simplify
        If True, simplify graph topology with the `simplify_graph` function.
    retain_all
        If True, return the entire graph even if it is not connected. If
        False, retain only the largest weakly connected component.
    truncate_by_edge
        If True, retain nodes outside bounding box if at least one of node's
        neighbors is within the bounding box.
    custom_filter
        A custom ways filter to be used instead of the `network_type` presets,
        e.g. `'["power"~"line"]' or '["highway"~"motorway|trunk"]'`. If `str`,
        the intersection of keys/values will be used, e.g., `'[maxspeed=50][lanes=2]'`
        will return all ways having both maxspeed of 50 and two lanes. If
        `list`, the union of the `list` items will be used, e.g.,
        `['[maxspeed=50]', '[lanes=2]']` will return all ways having either
        maximum speed of 50 or two lanes. Also pass in a `network_type` that
        is in `settings.bidirectional_network_types` if you want the graph to
        be fully bidirectional.

    Returns
    -------
    G

    Notes
    -----
    Very large query areas use the `utils_geo._consolidate_subdivide_geometry`
    function to automatically make multiple requests: see that function's
    documentation for caveats.
    """
    # verify that the geometry is valid and is a shapely Polygon/MultiPolygon
    # before proceeding
    if not polygon.is_valid:  # pragma: no cover
        msg = "The geometry of `polygon` is invalid."
        raise ValueError(msg)
    if not isinstance(polygon, (Polygon, MultiPolygon)):  # pragma: no cover
        msg = (
            "Geometry must be a shapely Polygon or MultiPolygon. If you "
            "requested graph from place name, make sure your query resolves "
            "to a Polygon or MultiPolygon, and not some other geometry, like "
            "a Point. See OSMnx documentation for details."
        )
        raise TypeError(msg)

    # create a new buffered polygon 0.5km around the desired one
    poly_proj, crs_utm = projection.project_geometry(polygon)
    poly_proj_buff = poly_proj.buffer(500)
    poly_buff, _ = projection.project_geometry(poly_proj_buff, crs=crs_utm, to_latlong=True)

    # download the network data from OSM within buffered polygon
    response_jsons = _overpass._download_overpass_network(poly_buff, network_type, custom_filter)

    # create buffered graph from the downloaded data
    bidirectional = network_type in settings.bidirectional_network_types
    G_buff = _create_graph(response_jsons, bidirectional)

    # truncate buffered graph to the buffered polygon and retain_all for
    # now. needed because overpass returns entire ways that also include
    # nodes outside the poly if the way (that is, a way with a single OSM
    # ID) has a node inside the poly at some point.
    G_buff = truncate.truncate_graph_polygon(G_buff, poly_buff, truncate_by_edge=truncate_by_edge)

    # keep only the largest weakly connected component if retain_all is False
    if not retain_all:
        G_buff = truncate.largest_component(G_buff, strongly=False)

    # simplify the graph topology
    if simplify:
        G_buff = simplification.simplify_graph(G_buff)

    # truncate graph by original polygon to return graph within polygon
    # caller wants. don't simplify again: this allows us to retain
    # intersections along the street that may now only connect 2 street
    # segments in the network, but in reality also connect to an
    # intersection just outside the polygon
    G = truncate.truncate_graph_polygon(G_buff, polygon, truncate_by_edge=truncate_by_edge)

    # keep only the largest weakly connected component if retain_all is False
    # we're doing this again in case the last truncate disconnected anything
    # on the periphery
    if not retain_all:
        G = truncate.largest_component(G, strongly=False)

    # count how many physical streets in buffered graph connect to each
    # intersection in un-buffered graph, to retain true counts for each
    # intersection, even if some of its neighbors are outside the polygon
    spn = stats.count_streets_per_node(G_buff, nodes=G.nodes)
    nx.set_node_attributes(G, values=spn, name="street_count")

    msg = f"graph_from_polygon returned graph with {len(G):,} nodes and {len(G.edges):,} edges"
    utils.log(msg, level=lg.INFO)
    return G


def graph_from_xml(
    filepath: str | Path,
    *,
    bidirectional: bool = False,
    simplify: bool = True,
    retain_all: bool = False,
    encoding: str = "utf-8",
) -> nx.MultiDiGraph:
    """
    Create a graph from data in an OSM XML file.

    Do not load an XML file previously generated by OSMnx: this use case is
    not supported and may not behave as expected. To save/load graphs to/from
    disk for later use in OSMnx, use the `io.save_graphml` and
    `io.load_graphml` functions instead.

    Use the `settings` module's `useful_tags_node` and `useful_tags_way`
    settings to configure which OSM node/way tags are added as graph node/edge
    attributes.

    Parameters
    ----------
    filepath
        Path to file containing OSM XML data.
    bidirectional
        If True, create bidirectional edges for one-way streets.
    simplify
        If True, simplify graph topology with the `simplify_graph` function.
    retain_all
        If True, return the entire graph even if it is not connected. If
        False, retain only the largest weakly connected component.
    encoding
        The OSM XML file's character encoding.

    Returns
    -------
    G
    """
    # transmogrify file of OSM XML data into JSON
    response_jsons = [_osm_xml._overpass_json_from_xml(filepath, encoding)]

    # create graph using this response JSON
    G = _create_graph(response_jsons, bidirectional)

    # keep only the largest weakly connected component if retain_all is False
    if not retain_all:
        G = truncate.largest_component(G, strongly=False)

    # simplify the graph topology as the last step
    if simplify:
        G = simplification.simplify_graph(G)

    msg = f"graph_from_xml returned graph with {len(G):,} nodes and {len(G.edges):,} edges"
    utils.log(msg, level=lg.INFO)
    return G


def _create_graph(
    response_jsons: Iterable[dict[str, Any]],
    bidirectional: bool,  # noqa: FBT001
) -> nx.MultiDiGraph:
    """
    Create a NetworkX MultiDiGraph from Overpass API responses.

    Adds length attributes in meters (great-circle distance between endpoints)
    to all of the graph's (pre-simplified, straight-line) edges via the
    `distance.add_edge_lengths` function.

    Parameters
    ----------
    response_jsons
        Iterable of JSON responses from the Overpass API.
    retain_all
        If True, return the entire graph even if it is not connected.
        Otherwise, retain only the largest weakly connected component.
    bidirectional
        If True, create bidirectional edges for one-way streets.

    Returns
    -------
    G
    """
    # each dict's keys are OSM IDs and values are dicts of attributes
    nodes: dict[int, dict[str, Any]] = {}
    paths: dict[int, dict[str, Any]] = {}

    # consume response_jsons generator to download data from server
    response_count = 0
    for response_json in response_jsons:
        response_count += 1

        # if cache_only_mode, consume response_jsons then continue next loop
        if settings.cache_only_mode:  # pragma: no cover
            continue

        # otherwise, extract nodes and paths from the downloaded OSM data
        nodes_temp, paths_temp = _parse_nodes_paths(response_json)
        nodes.update(nodes_temp)
        paths.update(paths_temp)

    msg = f"Retrieved all data from API in {response_count} request(s)"
    utils.log(msg, level=lg.INFO)
    if settings.cache_only_mode:  # pragma: no cover
        # after consuming all response_jsons in loop, raise exception to catch
        msg = "Interrupted because `settings.cache_only_mode=True`."
        raise CacheOnlyInterruptError(msg)

    # ensure we got some node/way data back from the server request(s)
    if (len(nodes) == 0) and (len(paths) == 0):  # pragma: no cover
        msg = "No data elements in server response. Check query location/filters and log."
        raise InsufficientResponseError(msg)

    # create the MultiDiGraph and set its graph-level attributes
    metadata = {
        "created_date": utils.ts(),
        "created_with": f"OSMnx {__version__}",
        "crs": settings.default_crs,
    }
    G = nx.MultiDiGraph(**metadata)

    # add each OSM node and way (a path of edges) to the graph
    msg = f"Creating graph from {len(nodes):,} OSM nodes and {len(paths):,} OSM ways..."
    utils.log(msg, level=lg.INFO)
    G.add_nodes_from(nodes.items())
    _add_paths(G, paths.values(), bidirectional)

    msg = f"Created graph with {len(G):,} nodes and {len(G.edges):,} edges"
    utils.log(msg, level=lg.INFO)

    # add length (great-circle distance between nodes) attribute to each edge
    if len(G.edges) > 0:
        G = distance.add_edge_lengths(G)

    return G


def _convert_node(element: dict[str, Any]) -> dict[str, Any]:
    """
    Convert an OSM node element into the format for a NetworkX node.

    Parameters
    ----------
    element
        OSM element of type "node".

    Returns
    -------
    node
    """
    node = {"y": element["lat"], "x": element["lon"]}
    if "tags" in element:
        for useful_tag in settings.useful_tags_node:
            if useful_tag in element["tags"]:
                node[useful_tag] = element["tags"][useful_tag]
    return node


def _convert_path(element: dict[str, Any]) -> dict[str, Any]:
    """
    Convert an OSM way element into the format for a NetworkX path.

    Parameters
    ----------
    element
        OSM element of type "way".

    Returns
    -------
    path
    """
    path = {"osmid": element["id"]}

    # remove any consecutive duplicate elements in the list of nodes
    path["nodes"] = [group[0] for group in groupby(element["nodes"])]

    if "tags" in element:
        for useful_tag in settings.useful_tags_way:
            if useful_tag in element["tags"]:
                path[useful_tag] = element["tags"][useful_tag]
    return path


def _parse_nodes_paths(
    response_json: dict[str, Any],
) -> tuple[dict[int, dict[str, Any]], dict[int, dict[str, Any]]]:
    """
    Construct dicts of nodes and paths from an Overpass response.

    Parameters
    ----------
    response_json
        JSON response from the Overpass API.

    Returns
    -------
    nodes, paths
        Each dict's keys are OSM IDs and values are dicts of attributes.
    """
    nodes = {}
    paths = {}
    for element in response_json["elements"]:
        if element["type"] == "node":
            nodes[element["id"]] = _convert_node(element)
        elif element["type"] == "way":
            paths[element["id"]] = _convert_path(element)

    return nodes, paths


def _is_path_one_way(attrs: dict[str, Any], bidirectional: bool, oneway_values: set[str]) -> bool:  # noqa: FBT001
    """
    Determine if a path of nodes allows travel in only one direction.

    Parameters
    ----------
    attrs
        A path's `tag:value` attribute data.
    bidirectional
        Whether this is a bidirectional network type.
    oneway_values
        The values OSM uses in its "oneway" tag to denote True.

    Returns
    -------
    is_one_way
    """
    # rule 1
    if settings.all_oneway:
        # if globally configured to set every edge one-way, then it's one-way
        return True

    # rule 2
    if bidirectional:
        # if this is a bidirectional network type, then nothing in it is
        # considered one-way. eg, if this is a walking network, this may very
        # well be a one-way street (as cars/bikes go), but in a walking-only
        # network it is a bidirectional edge (you can walk both directions on
        # a one-way street). so we will add this path (in both directions) to
        # the graph and set its oneway attribute to False.
        return False

    # rule 3
    if "oneway" in attrs and attrs["oneway"] in oneway_values:
        # if this path is tagged as one-way and if it is not a bidirectional
        # network type then we'll add the path in one direction only
        return True

    # rule 4
    if "junction" in attrs and attrs["junction"] == "roundabout":  # noqa: SIM103
        # roundabouts are also one-way but are not explicitly tagged as such
        return True

    # otherwise, if no rule passed then this path is not tagged as a one-way
    return False


def _is_path_reversed(attrs: dict[str, Any], reversed_values: set[str]) -> bool:
    """
    Determine if the order of nodes in a path should be reversed.

    Parameters
    ----------
    attrs
        A path's `tag:value` attribute data.
    reversed_values
        The values OSM uses in its 'oneway' tag to denote travel can only
        occur in the opposite direction of the node order.

    Returns
    -------
    is_reversed
    """
    return "oneway" in attrs and attrs["oneway"] in reversed_values


def _add_paths(
    G: nx.MultiDiGraph,
    paths: Iterable[dict[str, Any]],
    bidirectional: bool,  # noqa: FBT001
) -> None:
    """
    Add OSM paths to the graph as edges.

    Parameters
    ----------
    G
        The graph to add paths to.
    paths
        Iterable of paths' `tag:value` attribute data dicts.
    bidirectional
        If True, create bidirectional edges for one-way streets.

    Returns
    -------
    None
    """
    # the values OSM uses in its 'oneway' tag to denote True, and to denote
    # travel can only occur in the opposite direction of the node order. see:
    # https://wiki.openstreetmap.org/wiki/Key:oneway
    # https://www.geofabrik.de/de/data/geofabrik-osm-gis-standard-0.7.pdf
    oneway_values = {"yes", "true", "1", "-1", "reverse", "T", "F"}
    reversed_values = {"-1", "reverse", "T"}

    for path in paths:
        # extract/remove the ordered list of nodes from this path element so
        # we don't add it as a superfluous attribute to the edge later
        nodes = path.pop("nodes")

        # reverse the order of nodes in the path if this path is both one-way
        # and only allows travel in the opposite direction of nodes' order
        is_one_way = _is_path_one_way(path, bidirectional, oneway_values)
        if is_one_way and _is_path_reversed(path, reversed_values):
            nodes.reverse()

        # set the oneway attribute, but only if when not forcing all edges to
        # oneway with the all_oneway setting. With the all_oneway setting, you
        # want to preserve the original OSM oneway attribute for later clarity
        if not settings.all_oneway:
            path["oneway"] = is_one_way

        # zip path nodes to get (u, v) tuples like [(0,1), (1,2), (2,3)].
        edges = list(zip(nodes[:-1], nodes[1:]))

        # add all the edge tuples and give them the path's tag:value attrs
        path["reversed"] = False
        G.add_edges_from(edges, **path)

        # if the path is NOT one-way, reverse direction of each edge and add
        # this path going the opposite direction too
        if not is_one_way:
            path["reversed"] = True
            G.add_edges_from([(v, u) for u, v in edges], **path)
