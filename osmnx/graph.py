"""
Download and create graphs from OpenStreetMap data.

This module uses filters to query the Overpass API: you can either specify a
built-in network type or provide your own custom filter with Overpass QL.

Refer to the Getting Started guide for usage limitations.
"""

import itertools
from warnings import warn

import networkx as nx
from shapely.geometry import MultiPolygon
from shapely.geometry import Polygon

from . import _overpass
from . import distance
from . import geocoder
from . import osm_xml
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


def graph_from_bbox(
    north=None,
    south=None,
    east=None,
    west=None,
    bbox=None,
    network_type="all",
    simplify=True,
    retain_all=False,
    truncate_by_edge=False,
    clean_periphery=None,
    custom_filter=None,
):
    """
    Download and create a graph within some bounding box.

    You can use the `settings` module to retrieve a snapshot of historical OSM
    data as of a certain date, or to configure the Overpass server timeout,
    memory allocation, and other custom settings.

    Parameters
    ----------
    north : float
        deprecated, do not use
    south : float
        deprecated, do not use
    east : float
        deprecated, do not use
    west : float
        deprecated, do not use
    bbox : tuple of floats
        bounding box as (north, south, east, west)
    network_type : string {"all", "all_public", "bike", "drive", "drive_service", "walk"}
        what type of street network to get if custom_filter is None
    simplify : bool
        if True, simplify graph topology with the `simplify_graph` function
    retain_all : bool
        if True, return the entire graph even if it is not connected.
        otherwise, retain only the largest weakly connected component.
    truncate_by_edge : bool
        if True, retain nodes outside bounding box if at least one of node's
        neighbors is within the bounding box
    clean_periphery : bool
        deprecated, do not use
    custom_filter : string
        a custom ways filter to be used instead of the network_type presets
        e.g., '["power"~"line"]' or '["highway"~"motorway|trunk"]'. Also pass
        in a network_type that is in settings.bidirectional_network_types if
        you want graph to be fully bi-directional.

    Returns
    -------
    G : networkx.MultiDiGraph

    Notes
    -----
    Very large query areas use the `utils_geo._consolidate_subdivide_geometry`
    function to automatically make multiple requests: see that function's
    documentation for caveats.
    """
    if not (north is None and south is None and east is None and west is None):
        msg = (
            "The `north`, `south`, `east`, and `west` parameters are deprecated and "
            "will be removed in the v2.0.0 release. Use the `bbox` parameter instead. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
        )
        warn(msg, FutureWarning, stacklevel=2)
        bbox = (north, south, east, west)

    msg = (
        "The expected order of coordinates in `bbox` will change in the "
        "v2.0.0 release to `(left, bottom, right, top)`."
    )
    warn(msg, FutureWarning, stacklevel=2)

    # convert bounding box to a polygon
    polygon = utils_geo.bbox_to_poly(bbox=bbox)

    # create graph using this polygon geometry
    G = graph_from_polygon(
        polygon,
        network_type=network_type,
        simplify=simplify,
        retain_all=retain_all,
        truncate_by_edge=truncate_by_edge,
        clean_periphery=clean_periphery,
        custom_filter=custom_filter,
    )

    utils.log(f"graph_from_bbox returned graph with {len(G):,} nodes and {len(G.edges):,} edges")
    return G


def graph_from_point(
    center_point,
    dist=1000,
    dist_type="bbox",
    network_type="all",
    simplify=True,
    retain_all=False,
    truncate_by_edge=False,
    clean_periphery=None,
    custom_filter=None,
):
    """
    Download and create a graph within some distance of a (lat, lon) point.

    You can use the `settings` module to retrieve a snapshot of historical OSM
    data as of a certain date, or to configure the Overpass server timeout,
    memory allocation, and other custom settings.

    Parameters
    ----------
    center_point : tuple
        the (lat, lon) center point around which to construct the graph
    dist : int
        retain only those nodes within this many meters of the center of the
        graph, with distance determined according to dist_type argument
    dist_type : string {"network", "bbox"}
        if "bbox", retain only those nodes within a bounding box of the
        distance parameter. if "network", retain only those nodes within some
        network distance from the center-most node.
    network_type : string, {"all", "all_public", "bike", "drive", "drive_service", "walk"}
        what type of street network to get if custom_filter is None
    simplify : bool
        if True, simplify graph topology with the `simplify_graph` function
    retain_all : bool
        if True, return the entire graph even if it is not connected.
        otherwise, retain only the largest weakly connected component.
    truncate_by_edge : bool
        if True, retain nodes outside bounding box if at least one of node's
        neighbors is within the bounding box
    clean_periphery : bool,
        deprecated, do not use
    custom_filter : string
        a custom ways filter to be used instead of the network_type presets
        e.g., '["power"~"line"]' or '["highway"~"motorway|trunk"]'. Also pass
        in a network_type that is in settings.bidirectional_network_types if
        you want graph to be fully bi-directional.

    Returns
    -------
    G : networkx.MultiDiGraph

    Notes
    -----
    Very large query areas use the `utils_geo._consolidate_subdivide_geometry`
    function to automatically make multiple requests: see that function's
    documentation for caveats.
    """
    if dist_type not in {"bbox", "network"}:  # pragma: no cover
        msg = 'dist_type must be "bbox" or "network"'
        raise ValueError(msg)

    # create bounding box from center point and distance in each direction
    bbox = utils_geo.bbox_from_point(center_point, dist)

    # create a graph from the bounding box
    G = graph_from_bbox(
        bbox=bbox,
        network_type=network_type,
        simplify=simplify,
        retain_all=retain_all,
        truncate_by_edge=truncate_by_edge,
        clean_periphery=clean_periphery,
        custom_filter=custom_filter,
    )

    if dist_type == "network":
        # if dist_type is network, find node in graph nearest to center point
        # then truncate graph by network dist from it
        node = distance.nearest_nodes(G, X=[center_point[1]], Y=[center_point[0]])[0]
        G = truncate.truncate_graph_dist(G, node, max_dist=dist)

    utils.log(f"graph_from_point returned graph with {len(G):,} nodes and {len(G.edges):,} edges")
    return G


def graph_from_address(
    address,
    dist=1000,
    dist_type="bbox",
    network_type="all",
    simplify=True,
    retain_all=False,
    truncate_by_edge=False,
    return_coords=None,
    clean_periphery=None,
    custom_filter=None,
):
    """
    Download and create a graph within some distance of an address.

    You can use the `settings` module to retrieve a snapshot of historical OSM
    data as of a certain date, or to configure the Overpass server timeout,
    memory allocation, and other custom settings.

    Parameters
    ----------
    address : string
        the address to geocode and use as the central point around which to
        construct the graph
    dist : int
        retain only those nodes within this many meters of the center of the
        graph
    dist_type : string {"network", "bbox"}
        if "bbox", retain only those nodes within a bounding box of the
        distance parameter. if "network", retain only those nodes within some
        network distance from the center-most node.
    network_type : string {"all", "all_public", "bike", "drive", "drive_service", "walk"}
        what type of street network to get if custom_filter is None
    simplify : bool
        if True, simplify graph topology with the `simplify_graph` function
    retain_all : bool
        if True, return the entire graph even if it is not connected.
        otherwise, retain only the largest weakly connected component.
    truncate_by_edge : bool
        if True, retain nodes outside bounding box if at least one of node's
        neighbors is within the bounding box
    return_coords : bool
        deprecated, do not use
    clean_periphery : bool
        deprecated, do not use
    custom_filter : string
        a custom ways filter to be used instead of the network_type presets
        e.g., '["power"~"line"]' or '["highway"~"motorway|trunk"]'. Also pass
        in a network_type that is in settings.bidirectional_network_types if
        you want graph to be fully bi-directional.

    Returns
    -------
    networkx.MultiDiGraph or optionally (networkx.MultiDiGraph, (lat, lon))

    Notes
    -----
    Very large query areas use the `utils_geo._consolidate_subdivide_geometry`
    function to automatically make multiple requests: see that function's
    documentation for caveats.
    """
    if return_coords is None:
        return_coords = False
    else:
        warn(
            "The `return_coords` argument has been deprecated and will be removed in "
            "the v2.0.0 release. Future behavior will be as though `return_coords=False`. "
            "If you want the address's geocoded coordinates, use the `geocode` function. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123",
            FutureWarning,
            stacklevel=2,
        )
    # geocode the address string to a (lat, lon) point
    point = geocoder.geocode(query=address)

    # then create a graph from this point
    G = graph_from_point(
        point,
        dist,
        dist_type,
        network_type=network_type,
        simplify=simplify,
        retain_all=retain_all,
        truncate_by_edge=truncate_by_edge,
        clean_periphery=clean_periphery,
        custom_filter=custom_filter,
    )
    utils.log(f"graph_from_address returned graph with {len(G):,} nodes and {len(G.edges):,} edges")

    if return_coords:
        return G, point

    # otherwise
    return G


def graph_from_place(
    query,
    network_type="all",
    simplify=True,
    retain_all=False,
    truncate_by_edge=False,
    which_result=None,
    buffer_dist=None,
    clean_periphery=None,
    custom_filter=None,
):
    """
    Download and create a graph within the boundaries of some place(s).

    The query must be geocodable and OSM must have polygon boundaries for the
    geocode result. If OSM does not have a polygon for this place, you can
    instead get its street network using the graph_from_address function,
    which geocodes the place name to a point and gets the network within some
    distance of that point.

    If OSM does have polygon boundaries for this place but you're not finding
    it, try to vary the query string, pass in a structured query dict, or vary
    the which_result argument to use a different geocode result. If you know
    the OSM ID of the place, you can retrieve its boundary polygon using the
    geocode_to_gdf function, then pass it to the graph_from_polygon function.

    You can use the `settings` module to retrieve a snapshot of historical OSM
    data as of a certain date, or to configure the Overpass server timeout,
    memory allocation, and other custom settings.

    Parameters
    ----------
    query : string or dict or list
        the query or queries to geocode to get place boundary polygon(s)
    network_type : string {"all", "all_public", "bike", "drive", "drive_service", "walk"}
        what type of street network to get if custom_filter is None
    simplify : bool
        if True, simplify graph topology with the `simplify_graph` function
    retain_all : bool
        if True, return the entire graph even if it is not connected.
        otherwise, retain only the largest weakly connected component.
    truncate_by_edge : bool
        if True, retain nodes outside boundary polygon if at least one of
        node's neighbors is within the polygon
    which_result : int
        which geocoding result to use. if None, auto-select the first
        (Multi)Polygon or raise an error if OSM doesn't return one.
    buffer_dist : float
        deprecated, do not use
    clean_periphery : bool
        deprecated, do not use
    custom_filter : string
        a custom ways filter to be used instead of the network_type presets
        e.g., '["power"~"line"]' or '["highway"~"motorway|trunk"]'. Also pass
        in a network_type that is in settings.bidirectional_network_types if
        you want graph to be fully bi-directional.

    Returns
    -------
    G : networkx.MultiDiGraph

    Notes
    -----
    Very large query areas use the `utils_geo._consolidate_subdivide_geometry`
    function to automatically make multiple requests: see that function's
    documentation for caveats.
    """
    if buffer_dist is not None:
        warn(
            "The buffer_dist argument has been deprecated and will be removed "
            "in the v2.0.0 release. Buffer your query area directly, if desired. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123",
            FutureWarning,
            stacklevel=2,
        )

    # create a GeoDataFrame with the spatial boundaries of the place(s)
    if isinstance(query, (str, dict)):
        # if it is a string (place name) or dict (structured place query),
        # then it is a single place
        gdf_place = geocoder.geocode_to_gdf(
            query, which_result=which_result, buffer_dist=buffer_dist
        )
    elif isinstance(query, list):
        # if it is a list, it contains multiple places to get
        gdf_place = geocoder.geocode_to_gdf(query, buffer_dist=buffer_dist)
    else:  # pragma: no cover
        msg = "query must be dict, string, or list of strings"
        raise TypeError(msg)

    # extract the geometry from the GeoDataFrame to use in API query
    polygon = gdf_place["geometry"].unary_union
    utils.log("Constructed place geometry polygon(s) to query API")

    # create graph using this polygon(s) geometry
    G = graph_from_polygon(
        polygon,
        network_type=network_type,
        simplify=simplify,
        retain_all=retain_all,
        truncate_by_edge=truncate_by_edge,
        clean_periphery=clean_periphery,
        custom_filter=custom_filter,
    )

    utils.log(f"graph_from_place returned graph with {len(G):,} nodes and {len(G.edges):,} edges")
    return G


def graph_from_polygon(
    polygon,
    network_type="all",
    simplify=True,
    retain_all=False,
    truncate_by_edge=False,
    clean_periphery=None,
    custom_filter=None,
):
    """
    Download and create a graph within the boundaries of a (multi)polygon.

    You can use the `settings` module to retrieve a snapshot of historical OSM
    data as of a certain date, or to configure the Overpass server timeout,
    memory allocation, and other custom settings.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        the shape to get network data within. coordinates should be in
        unprojected latitude-longitude degrees (EPSG:4326).
    network_type : string {"all", "all_public", "bike", "drive", "drive_service", "walk"}
        what type of street network to get if custom_filter is None
    simplify : bool
        if True, simplify graph topology with the `simplify_graph` function
    retain_all : bool
        if True, return the entire graph even if it is not connected.
        otherwise, retain only the largest weakly connected component.
    truncate_by_edge : bool
        if True, retain nodes outside boundary polygon if at least one of
        node's neighbors is within the polygon
    clean_periphery : bool
        deprecated, do not use
    custom_filter : string
        a custom ways filter to be used instead of the network_type presets
        e.g., '["power"~"line"]' or '["highway"~"motorway|trunk"]'. Also pass
        in a network_type that is in settings.bidirectional_network_types if
        you want graph to be fully bi-directional.

    Returns
    -------
    G : networkx.MultiDiGraph

    Notes
    -----
    Very large query areas use the `utils_geo._consolidate_subdivide_geometry`
    function to automatically make multiple requests: see that function's
    documentation for caveats.
    """
    if clean_periphery is None:
        clean_periphery = True
    else:
        warn(
            "The clean_periphery argument has been deprecated and will be removed in "
            "the v2.0.0 release. Future behavior will be as though clean_periphery=True. "
            "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123",
            FutureWarning,
            stacklevel=2,
        )

    # verify that the geometry is valid and is a shapely Polygon/MultiPolygon
    # before proceeding
    if not polygon.is_valid:  # pragma: no cover
        msg = "The geometry to query within is invalid"
        raise ValueError(msg)
    if not isinstance(polygon, (Polygon, MultiPolygon)):  # pragma: no cover
        msg = (
            "Geometry must be a shapely Polygon or MultiPolygon. If you "
            "requested graph from place name, make sure your query resolves "
            "to a Polygon or MultiPolygon, and not some other geometry, like "
            "a Point. See OSMnx documentation for details."
        )
        raise TypeError(msg)

    if clean_periphery:
        # create a new buffered polygon 0.5km around the desired one
        buffer_dist = 500
        poly_proj, crs_utm = projection.project_geometry(polygon)
        poly_proj_buff = poly_proj.buffer(buffer_dist)
        poly_buff, _ = projection.project_geometry(poly_proj_buff, crs=crs_utm, to_latlong=True)

        # download the network data from OSM within buffered polygon
        response_jsons = _overpass._download_overpass_network(
            poly_buff, network_type, custom_filter
        )

        # create buffered graph from the downloaded data
        bidirectional = network_type in settings.bidirectional_network_types
        G_buff = _create_graph(response_jsons, retain_all=True, bidirectional=bidirectional)

        # truncate buffered graph to the buffered polygon and retain_all for
        # now. needed because overpass returns entire ways that also include
        # nodes outside the poly if the way (that is, a way with a single OSM
        # ID) has a node inside the poly at some point.
        G_buff = truncate.truncate_graph_polygon(G_buff, poly_buff, True, truncate_by_edge)

        # simplify the graph topology
        if simplify:
            G_buff = simplification.simplify_graph(G_buff)

        # truncate graph by original polygon to return graph within polygon
        # caller wants. don't simplify again: this allows us to retain
        # intersections along the street that may now only connect 2 street
        # segments in the network, but in reality also connect to an
        # intersection just outside the polygon
        G = truncate.truncate_graph_polygon(G_buff, polygon, retain_all, truncate_by_edge)

        # count how many physical streets in buffered graph connect to each
        # intersection in un-buffered graph, to retain true counts for each
        # intersection, even if some of its neighbors are outside the polygon
        spn = stats.count_streets_per_node(G_buff, nodes=G.nodes)
        nx.set_node_attributes(G, values=spn, name="street_count")

    # if clean_periphery=False, just use the polygon as provided
    else:
        # download the network data from OSM
        response_jsons = _overpass._download_overpass_network(polygon, network_type, custom_filter)

        # create graph from the downloaded data
        bidirectional = network_type in settings.bidirectional_network_types
        G = _create_graph(response_jsons, retain_all=True, bidirectional=bidirectional)

        # truncate the graph to the extent of the polygon
        G = truncate.truncate_graph_polygon(G, polygon, retain_all, truncate_by_edge)

        # simplify the graph topology after truncation. don't truncate after
        # simplifying or you may have simplified out to an endpoint beyond the
        # truncation distance, which would strip out the entire edge
        if simplify:
            G = simplification.simplify_graph(G)

        # count how many physical streets connect to each intersection/deadend
        # note this will be somewhat inaccurate due to periphery effects, so
        # it's best to parameterize function with clean_periphery=True
        spn = stats.count_streets_per_node(G)
        nx.set_node_attributes(G, values=spn, name="street_count")
        warn(
            "the graph-level street_count attribute will likely be inaccurate "
            "when you set clean_periphery=False",
            stacklevel=2,
        )

    utils.log(f"graph_from_polygon returned graph with {len(G):,} nodes and {len(G.edges):,} edges")
    return G


def graph_from_xml(
    filepath, bidirectional=False, simplify=True, retain_all=False, encoding="utf-8"
):
    """
    Create a graph from data in a .osm formatted XML file.

    Do not load an XML file generated by OSMnx: this use case is not supported
    and may not behave as expected. To save/load graphs to/from disk for later
    use in OSMnx, use the `io.save_graphml` and `io.load_graphml` functions
    instead.

    Parameters
    ----------
    filepath : string or pathlib.Path
        path to file containing OSM XML data
    bidirectional : bool
        if True, create bi-directional edges for one-way streets
    simplify : bool
        if True, simplify graph topology with the `simplify_graph` function
    retain_all : bool
        if True, return the entire graph even if it is not connected.
        otherwise, retain only the largest weakly connected component.
    encoding : string
        the XML file's character encoding

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    # transmogrify file of OSM XML data into JSON
    response_jsons = [osm_xml._overpass_json_from_file(filepath, encoding)]

    # create graph using this response JSON
    G = _create_graph(response_jsons, bidirectional=bidirectional, retain_all=retain_all)

    # simplify the graph topology as the last step
    if simplify:
        G = simplification.simplify_graph(G)

    utils.log(f"graph_from_xml returned graph with {len(G):,} nodes and {len(G.edges):,} edges")
    return G


def _create_graph(response_jsons, retain_all=False, bidirectional=False):
    """
    Create a networkx MultiDiGraph from Overpass API responses.

    Adds length attributes in meters (great-circle distance between endpoints)
    to all of the graph's (pre-simplified, straight-line) edges via the
    `distance.add_edge_lengths` function.

    Parameters
    ----------
    response_jsons : iterable
        iterable of dicts of JSON responses from from the Overpass API
    retain_all : bool
        if True, return the entire graph even if it is not connected.
        otherwise, retain only the largest weakly connected component.
    bidirectional : bool
        if True, create bi-directional edges for one-way streets

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    response_count = 0
    nodes = {}
    paths = {}

    # consume response_jsons generator to download data from server
    for response_json in response_jsons:
        response_count += 1

        # if cache_only_mode, consume response_jsons then continue next loop
        if settings.cache_only_mode:  # pragma: no cover
            continue

        # otherwise, extract nodes and paths from the downloaded OSM data
        nodes_temp, paths_temp = _parse_nodes_paths(response_json)
        nodes.update(nodes_temp)
        paths.update(paths_temp)

    utils.log(f"Retrieved all data from API in {response_count} request(s)")
    if settings.cache_only_mode:  # pragma: no cover
        # after consuming all response_jsons in loop, raise exception to catch
        msg = "Interrupted because `settings.cache_only_mode=True`"
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
    utils.log(f"Creating graph from {len(nodes):,} OSM nodes and {len(paths):,} OSM ways...")
    G.add_nodes_from(nodes.items())
    _add_paths(G, paths.values(), bidirectional)

    # retain only the largest connected component if retain_all=False
    if not retain_all:
        G = truncate.largest_component(G)

    utils.log(f"Created graph with {len(G):,} nodes and {len(G.edges):,} edges")

    # add length (great-circle distance between nodes) attribute to each edge
    if len(G.edges) > 0:
        G = distance.add_edge_lengths(G)

    return G


def _convert_node(element):
    """
    Convert an OSM node element into the format for a networkx node.

    Parameters
    ----------
    element : dict
        an OSM node element

    Returns
    -------
    node : dict
    """
    node = {"y": element["lat"], "x": element["lon"]}
    if "tags" in element:
        for useful_tag in settings.useful_tags_node:
            if useful_tag in element["tags"]:
                node[useful_tag] = element["tags"][useful_tag]
    return node


def _convert_path(element):
    """
    Convert an OSM way element into the format for a networkx path.

    Parameters
    ----------
    element : dict
        an OSM way element

    Returns
    -------
    path : dict
    """
    path = {"osmid": element["id"]}

    # remove any consecutive duplicate elements in the list of nodes
    path["nodes"] = [group[0] for group in itertools.groupby(element["nodes"])]

    if "tags" in element:
        for useful_tag in settings.useful_tags_way:
            if useful_tag in element["tags"]:
                path[useful_tag] = element["tags"][useful_tag]
    return path


def _parse_nodes_paths(response_json):
    """
    Construct dicts of nodes and paths from an Overpass response.

    Parameters
    ----------
    response_json : dict
        JSON response from the Overpass API

    Returns
    -------
    nodes, paths : tuple of dicts
        dicts' keys = osmid and values = dict of attributes
    """
    nodes = {}
    paths = {}
    for element in response_json["elements"]:
        if element["type"] == "node":
            nodes[element["id"]] = _convert_node(element)
        elif element["type"] == "way":
            paths[element["id"]] = _convert_path(element)

    return nodes, paths


def _is_path_one_way(path, bidirectional, oneway_values):
    """
    Determine if a path of nodes allows travel in only one direction.

    Parameters
    ----------
    path : dict
        a path's `tag:value` attribute data
    bidirectional : bool
        whether this is a bi-directional network type
    oneway_values : set
        the values OSM uses in its 'oneway' tag to denote True

    Returns
    -------
    bool
    """
    # rule 1
    if settings.all_oneway:
        # if globally configured to set every edge one-way, then it's one-way
        return True

    # rule 2
    if bidirectional:
        # if this is a bi-directional network type, then nothing in it is
        # considered one-way. eg, if this is a walking network, this may very
        # well be a one-way street (as cars/bikes go), but in a walking-only
        # network it is a bi-directional edge (you can walk both directions on
        # a one-way street). so we will add this path (in both directions) to
        # the graph and set its oneway attribute to False.
        return False

    # rule 3
    if "oneway" in path and path["oneway"] in oneway_values:
        # if this path is tagged as one-way and if it is not a bi-directional
        # network type then we'll add the path in one direction only
        return True

    # rule 4
    if "junction" in path and path["junction"] == "roundabout":
        # roundabouts are also one-way but are not explicitly tagged as such
        return True

    # otherwise, if no rule passed then this path is not tagged as a one-way
    return False


def _is_path_reversed(path, reversed_values):
    """
    Determine if the order of nodes in a path should be reversed.

    Parameters
    ----------
    path : dict
        a path's `tag:value` attribute data
    reversed_values : set
        the values OSM uses in its 'oneway' tag to denote travel can only
        occur in the opposite direction of the node order

    Returns
    -------
    bool
    """
    return "oneway" in path and path["oneway"] in reversed_values


def _add_paths(G, paths, bidirectional=False):
    """
    Add a list of paths to the graph as edges.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        graph to add paths to
    paths : list
        list of paths' `tag:value` attribute data dicts
    bidirectional : bool
        if True, create bi-directional edges for one-way streets

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
