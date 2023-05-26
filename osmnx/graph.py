"""Graph creation functions."""

import itertools
import warnings

import networkx as nx
from shapely.geometry import MultiPolygon
from shapely.geometry import Polygon

from . import distance
from . import downloader
from . import geocoder
from . import osm_xml
from . import projection
from . import settings
from . import simplification
from . import stats
from . import truncate
from . import utils
from . import utils_geo
from . import utils_graph
from ._errors import EmptyOverpassResponse
from ._version import __version__


def graph_from_bbox(
    north,
    south,
    east,
    west,
    network_type="all_private",
    simplify=True,
    retain_all=False,
    truncate_by_edge=False,
    clean_periphery=True,
    custom_filter=None,
):
    """
    Create a graph from OSM within some bounding box.

    Parameters
    ----------
    north : float
        northern latitude of bounding box
    south : float
        southern latitude of bounding box
    east : float
        eastern longitude of bounding box
    west : float
        western longitude of bounding box
    network_type : string {"all_private", "all", "bike", "drive", "drive_service", "walk"}
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
        if True, buffer 500m to get a graph larger than requested, then
        simplify, then truncate it to requested spatial boundaries
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
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via the `settings` module. Very large query areas
    will use the utils_geo._consolidate_subdivide_geometry function to perform
    multiple queries: see that function's documentation for caveats.
    """
    # convert bounding box to a polygon
    polygon = utils_geo.bbox_to_poly(north, south, east, west)

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

    utils.log(f"graph_from_bbox returned graph with {len(G)} nodes and {len(G.edges)} edges")
    return G


def graph_from_point(
    center_point,
    dist=1000,
    dist_type="bbox",
    network_type="all_private",
    simplify=True,
    retain_all=False,
    truncate_by_edge=False,
    clean_periphery=True,
    custom_filter=None,
):
    """
    Create a graph from OSM within some distance of some (lat, lng) point.

    Parameters
    ----------
    center_point : tuple
        the (lat, lng) center point around which to construct the graph
    dist : int
        retain only those nodes within this many meters of the center of the
        graph, with distance determined according to dist_type argument
    dist_type : string {"network", "bbox"}
        if "bbox", retain only those nodes within a bounding box of the
        distance parameter. if "network", retain only those nodes within some
        network distance from the center-most node (requires that scikit-learn
        is installed as an optional dependency).
    network_type : string, {"all_private", "all", "bike", "drive", "drive_service", "walk"}
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
        if True, buffer 500m to get a graph larger than requested, then
        simplify, then truncate it to requested spatial boundaries
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
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via the `settings` module. Very large query areas
    will use the utils_geo._consolidate_subdivide_geometry function to perform
    multiple queries: see that function's documentation for caveats.
    """
    if dist_type not in {"bbox", "network"}:  # pragma: no cover
        raise ValueError('dist_type must be "bbox" or "network"')

    # create bounding box from center point and distance in each direction
    north, south, east, west = utils_geo.bbox_from_point(center_point, dist)

    # create a graph from the bounding box
    G = graph_from_bbox(
        north,
        south,
        east,
        west,
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

    utils.log(f"graph_from_point returned graph with {len(G)} nodes and {len(G.edges)} edges")
    return G


def graph_from_address(
    address,
    dist=1000,
    dist_type="bbox",
    network_type="all_private",
    simplify=True,
    retain_all=False,
    truncate_by_edge=False,
    return_coords=False,
    clean_periphery=True,
    custom_filter=None,
):
    """
    Create a graph from OSM within some distance of some address.

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
        network distance from the center-most node (requires that scikit-learn
        is installed as an optional dependency).
    network_type : string {"all_private", "all", "bike", "drive", "drive_service", "walk"}
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
        optionally also return the geocoded coordinates of the address
    clean_periphery : bool,
        if True, buffer 500m to get a graph larger than requested, then
        simplify, then truncate it to requested spatial boundaries
    custom_filter : string
        a custom ways filter to be used instead of the network_type presets
        e.g., '["power"~"line"]' or '["highway"~"motorway|trunk"]'. Also pass
        in a network_type that is in settings.bidirectional_network_types if
        you want graph to be fully bi-directional.

    Returns
    -------
    networkx.MultiDiGraph or optionally (networkx.MultiDiGraph, (lat, lng))

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via the `settings` module. Very large query areas
    will use the utils_geo._consolidate_subdivide_geometry function to perform
    multiple queries: see that function's documentation for caveats.
    """
    # geocode the address string to a (lat, lng) point
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
    utils.log(f"graph_from_address returned graph with {len(G)} nodes and {len(G.edges)} edges")

    if return_coords:
        return G, point
    else:
        return G


def graph_from_place(
    query,
    network_type="all_private",
    simplify=True,
    retain_all=False,
    truncate_by_edge=False,
    which_result=None,
    buffer_dist=None,
    clean_periphery=True,
    custom_filter=None,
):
    """
    Create graph from OSM within the boundaries of some geocodable place(s).

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

    Parameters
    ----------
    query : string or dict or list
        the query or queries to geocode to get place boundary polygon(s)
    network_type : string {"all_private", "all", "bike", "drive", "drive_service", "walk"}
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
        distance to buffer around the place geometry, in meters
    clean_periphery : bool
        if True, buffer 500m to get a graph larger than requested, then
        simplify, then truncate it to requested spatial boundaries
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
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via the `settings` module. Very large query areas
    will use the utils_geo._consolidate_subdivide_geometry function to perform
    multiple queries: see that function's documentation for caveats.
    """
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
        raise TypeError("query must be dict, string, or list of strings")

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

    utils.log(f"graph_from_place returned graph with {len(G)} nodes and {len(G.edges)} edges")
    return G


def graph_from_polygon(
    polygon,
    network_type="all_private",
    simplify=True,
    retain_all=False,
    truncate_by_edge=False,
    clean_periphery=True,
    custom_filter=None,
):
    """
    Create a graph from OSM within the boundaries of some shapely polygon.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        the shape to get network data within. coordinates should be in
        unprojected latitude-longitude degrees (EPSG:4326).
    network_type : string {"all_private", "all", "bike", "drive", "drive_service", "walk"}
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
        if True, buffer 500m to get a graph larger than requested, then
        simplify, then truncate it to requested spatial boundaries
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
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via the `settings` module. Very large query areas
    will use the utils_geo._consolidate_subdivide_geometry function to perform
    multiple queries: see that function's documentation for caveats.
    """
    # verify that the geometry is valid and is a shapely Polygon/MultiPolygon
    # before proceeding
    if not polygon.is_valid:  # pragma: no cover
        raise ValueError("The geometry to query within is invalid")
    if not isinstance(polygon, (Polygon, MultiPolygon)):  # pragma: no cover
        raise TypeError(
            "Geometry must be a shapely Polygon or MultiPolygon. If you requested "
            "graph from place name, make sure your query resolves to a Polygon or "
            "MultiPolygon, and not some other geometry, like a Point. See OSMnx "
            "documentation for details."
        )

    if clean_periphery:
        # create a new buffered polygon 0.5km around the desired one
        buffer_dist = 500
        poly_proj, crs_utm = projection.project_geometry(polygon)
        poly_proj_buff = poly_proj.buffer(buffer_dist)
        poly_buff, _ = projection.project_geometry(poly_proj_buff, crs=crs_utm, to_latlong=True)

        # download the network data from OSM within buffered polygon
        response_jsons = downloader._osm_network_download(poly_buff, network_type, custom_filter)

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
        response_jsons = downloader._osm_network_download(polygon, network_type, custom_filter)

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
        msg = (
            "the graph-level street_count attribute will likely be inaccurate "
            "when you set clean_periphery=False"
        )
        warnings.warn(msg, stacklevel=1)

    utils.log(f"graph_from_polygon returned graph with {len(G)} nodes and {len(G.edges)} edges")
    return G


def graph_from_xml(filepath, bidirectional=False, simplify=True, retain_all=False):
    """
    Create a graph from data in a .osm formatted XML file.

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

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    # transmogrify file of OSM XML data into JSON
    response_jsons = [osm_xml._overpass_json_from_file(filepath)]

    # create graph using this response JSON
    G = _create_graph(response_jsons, bidirectional=bidirectional, retain_all=retain_all)

    # simplify the graph topology as the last step
    if simplify:
        G = simplification.simplify_graph(G)

    utils.log(f"graph_from_xml returned graph with {len(G)} nodes and {len(G.edges)} edges")
    return G


def _create_graph(response_jsons, retain_all=False, bidirectional=False):
    """
    Create a networkx MultiDiGraph from Overpass API responses.

    Adds length attributes in meters (great-circle distance between endpoints)
    to all of the graph's (pre-simplified, straight-line) edges via the
    `distance.add_edge_lengths` function.

    Parameters
    ----------
    response_jsons : list
        list of dicts of JSON responses from from the Overpass API
    retain_all : bool
        if True, return the entire graph even if it is not connected.
        otherwise, retain only the largest weakly connected component.
    bidirectional : bool
        if True, create bi-directional edges for one-way streets

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    utils.log("Creating graph from downloaded OSM data...")

    # make sure we got data back from the server request(s)
    if not any(rj["elements"] for rj in response_jsons):  # pragma: no cover
        msg = "There are no data elements in the server response. Check log and query location/filters."
        raise EmptyOverpassResponse(msg)

    # create the graph as a MultiDiGraph and set its meta-attributes
    metadata = {
        "created_date": utils.ts(),
        "created_with": f"OSMnx {__version__}",
        "crs": settings.default_crs,
    }
    G = nx.MultiDiGraph(**metadata)

    # extract nodes and paths from the downloaded osm data
    nodes = {}
    paths = {}
    for response_json in response_jsons:
        nodes_temp, paths_temp = _parse_nodes_paths(response_json)
        nodes.update(nodes_temp)
        paths.update(paths_temp)

    # add each osm node to the graph
    for node, data in nodes.items():
        G.add_node(node, **data)

    # add each osm way (ie, a path of edges) to the graph
    _add_paths(G, paths.values(), bidirectional)

    # retain only the largest connected component if retain_all is False
    if not retain_all:
        G = utils_graph.get_largest_component(G)

    utils.log(f"Created graph with {len(G)} nodes and {len(G.edges)} edges")

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
        a path's tag:value attribute data
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
    elif bidirectional:
        # if this is a bi-directional network type, then nothing in it is
        # considered one-way. eg, if this is a walking network, this may very
        # well be a one-way street (as cars/bikes go), but in a walking-only
        # network it is a bi-directional edge (you can walk both directions on
        # a one-way street). so we will add this path (in both directions) to
        # the graph and set its oneway attribute to False.
        return False

    # rule 3
    elif "oneway" in path and path["oneway"] in oneway_values:
        # if this path is tagged as one-way and if it is not a bi-directional
        # network type then we'll add the path in one direction only
        return True

    # rule 4
    elif "junction" in path and path["junction"] == "roundabout":
        # roundabouts are also one-way but are not explicitly tagged as such
        return True

    else:
        # otherwise this path is not tagged as a one-way
        return False


def _is_path_reversed(path, reversed_values):
    """
    Determine if the order of nodes in a path should be reversed.

    Parameters
    ----------
    path : dict
        a path's tag:value attribute data
    reversed_values : set
        the values OSM uses in its 'oneway' tag to denote travel can only
        occur in the opposite direction of the node order

    Returns
    -------
    bool
    """
    if "oneway" in path and path["oneway"] in reversed_values:
        return True
    else:
        return False


def _add_paths(G, paths, bidirectional=False):
    """
    Add a list of paths to the graph as edges.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        graph to add paths to
    paths : list
        list of paths' tag:value attribute data dicts
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
