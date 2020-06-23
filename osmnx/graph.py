"""Graph creation functions."""

import bz2
import os
import xml.sax
from itertools import groupby

import networkx as nx
from shapely.geometry import MultiPolygon
from shapely.geometry import Polygon

from . import distance
from . import downloader
from . import geocoder
from . import projection
from . import settings
from . import simplification
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
    network_type : string
        what type of street network to get if custom_filter is None. One of
        'walk', 'bike', 'drive', 'drive_service', 'all', or 'all_private'.
    simplify : bool
        if True, simplify the graph topology
    retain_all : bool
        if True, return the entire graph even if it is not connected
    truncate_by_edge : bool
        if True, retain node if it's outside bounding box but at least one of
        node's neighbors are within the bounding box
    clean_periphery : bool
        if True, buffer 500m to get a graph larger than requested, then
        simplify, then truncate it to requested spatial boundaries
    custom_filter : string
        a custom network filter to be used instead of the network_type presets,
        e.g., '["power"~"line"]' or '["highway"~"motorway|trunk"]'. Also pass
        in a network_type that is in settings.bidirectional_network_types if
        you want graph to be fully bidirectional.

    Returns
    -------
    G : networkx.MultiDiGraph

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
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

    utils.log(f"graph_from_bbox returned graph with {len(G)} nodes and {len(G.edges())} edges")
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
    dist_type : string
        {'network', 'bbox'} if 'bbox', retain only those nodes within a bounding
        box of the distance parameter. if 'network', retain only those nodes
        within some network distance from the center-most node.
    network_type : string
        what type of street network to get if custom_filter is None. One of
        'walk', 'bike', 'drive', 'drive_service', 'all', or 'all_private'.
    simplify : bool
        if True, simplify the graph topology
    retain_all : bool
        if True, return the entire graph even if it is not connected
    truncate_by_edge : bool
        if True, retain node if it's outside bounding box but at least one of
        node's neighbors are within bounding box
    clean_periphery : bool,
        if True, buffer 500m to get a graph larger than requested, then
        simplify, then truncate it to requested spatial boundaries
    custom_filter : string
        a custom network filter to be used instead of the network_type presets,
        e.g., '["power"~"line"]' or '["highway"~"motorway|trunk"]'. Also pass
        in a network_type that is in settings.bidirectional_network_types if
        you want graph to be fully bidirectional.

    Returns
    -------
    G : networkx.MultiDiGraph

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """
    if dist_type not in {"bbox", "network"}:
        raise ValueError('dist_type must be "bbox" or "network"')

    # create a bounding box from the center point and the distance in each
    # direction
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

    # if the network dist_type is network, find the node in the graph
    # nearest to the center point, and truncate the graph by network distance
    # from this node
    if dist_type == "network":
        centermost_node = distance.get_nearest_node(G, center_point)
        G = truncate.truncate_graph_dist(G, centermost_node, max_dist=dist)

    utils.log(f"graph_from_point returned graph with {len(G)} nodes and {len(G.edges())} edges")
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
    dist_type : string
        {'network', 'bbox'} if 'bbox', retain only those nodes within a bounding
        box of the distance parameter.
        if 'network', retain only those nodes within some network distance from
        the center-most node.
    network_type : string
        what type of street network to get if custom_filter is None. One of
        'walk', 'bike', 'drive', 'drive_service', 'all', or 'all_private'.
    simplify : bool
        if True, simplify the graph topology
    retain_all : bool
        if True, return the entire graph even if it is not connected
    truncate_by_edge : bool
        if True, retain node if it's outside bounding box but at least one of
        node's neighbors are within bounding box
    return_coords : bool
        optionally also return the geocoded coordinates of the address
    clean_periphery : bool,
        if True, buffer 500m to get a graph larger than requested, then
        simplify, then truncate it to requested spatial boundaries
    custom_filter : string
        a custom network filter to be used instead of the network_type presets,
        e.g., '["power"~"line"]' or '["highway"~"motorway|trunk"]'. Also pass
        in a network_type that is in settings.bidirectional_network_types if
        you want graph to be fully bidirectional.

    Returns
    -------
    networkx.MultiDiGraph or optionally (networkx.MultiDiGraph, (lat, lng))

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
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
    utils.log(f"graph_from_address returned graph with {len(G)} nodes and {len(G.edges())} edges")

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
    which_result=1,
    buffer_dist=None,
    clean_periphery=True,
    custom_filter=None,
):
    """
    Create graph from OSM within the boundaries of some geocodable place(s).

    The query must be geocodable and OSM must have polygon boundaries for the
    geocode result. If OSM does not have a polygon for this place, you can
    instead get its street network using the graph_from_address function, which
    geocodes the place name to a point and gets the network within some distance
    of that point. Alternatively, you might try to vary the which_result
    parameter to use a different geocode result. For example, the first geocode
    result (ie, the default) might resolve to a point geometry, but the second
    geocode result for this query might resolve to a polygon, in which case you
    can use graph_from_place with which_result=2.

    Parameters
    ----------
    query : string or dict or list
        the place(s) to geocode/download data for
    network_type : string
        what type of street network to get if custom_filter is None. One of
        'walk', 'bike', 'drive', 'drive_service', 'all', or 'all_private'.
    simplify : bool
        if True, simplify the graph topology
    retain_all : bool
        if True, return the entire graph even if it is not connected
    truncate_by_edge : bool
        if True, retain node if it's outside polygon but at least one of
        node's neighbors are within bbox
    which_result : int
        max number of results to return and which to process upon receipt
    buffer_dist : float
        distance to buffer around the place geometry, in meters
    clean_periphery : bool
        if True, buffer 500m to get a graph larger than requested, then
        simplify, then truncate it to requested spatial boundaries
    custom_filter : string
        a custom network filter to be used instead of the network_type presets,
        e.g., '["power"~"line"]' or '["highway"~"motorway|trunk"]'. Also pass
        in a network_type that is in settings.bidirectional_network_types if
        you want graph to be fully bidirectional.

    Returns
    -------
    G : networkx.MultiDiGraph

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """
    # create a GeoDataFrame with the spatial boundaries of the place(s)
    if isinstance(query, (str, dict)):
        # if it is a string (place name) or dict (structured place query), then
        # it is a single place
        gdf_place = geocoder.geocode_to_gdf(
            query, which_result=which_result, buffer_dist=buffer_dist
        )
    elif isinstance(query, list):
        # if it is a list, it contains multiple places to get
        gdf_place = geocoder.geocode_to_gdf(query, buffer_dist=buffer_dist)
    else:
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

    utils.log(f"graph_from_place returned graph with {len(G)} nodes and {len(G.edges())} edges")
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
        the shape to get network data within. coordinates should be in units of
        latitude-longitude degrees.
    network_type : string
        what type of street network to get if custom_filter is None. One of
        'walk', 'bike', 'drive', 'drive_service', 'all', or 'all_private'.
    simplify : bool
        if True, simplify the graph topology
    retain_all : bool
        if True, return the entire graph even if it is not connected
    truncate_by_edge : bool
        if True, retain node if it's outside polygon but at least one of
        node's neighbors are within polygon
    clean_periphery : bool
        if True, buffer 500m to get a graph larger than requested, then
        simplify, then truncate it to requested spatial boundaries
    custom_filter : string
        a custom network filter to be used instead of the network_type presets,
        e.g., '["power"~"line"]' or '["highway"~"motorway|trunk"]'. Also pass
        in a network_type that is in settings.bidirectional_network_types if
        you want graph to be fully bidirectional.

    Returns
    -------
    G : networkx.MultiDiGraph

    Notes
    -----
    You can configure the Overpass server timeout, memory allocation, and
    other custom settings via ox.config().
    """
    # verify that the geometry is valid and is a shapely Polygon/MultiPolygon
    # before proceeding
    if not polygon.is_valid:
        raise ValueError("The geometry to query within is invalid")
    if not isinstance(polygon, (Polygon, MultiPolygon)):
        raise TypeError(
            "Geometry must be a shapely Polygon or MultiPolygon. If you requested "
            "graph from place name or address, make sure your query resolves to a "
            "Polygon or MultiPolygon, and not some other geometry, like a Point. "
            "See OSMnx documentation for details."
        )

    if clean_periphery:
        # create a new buffered polygon 0.5km around the desired one
        buffer_dist = 500
        poly_proj, crs_utm = projection.project_geometry(polygon)
        poly_proj_buff = poly_proj.buffer(buffer_dist)
        poly_buff, _ = projection.project_geometry(poly_proj_buff, crs=crs_utm, to_latlong=True)

        # download the network data from OSM within buffered polygon
        response_jsons = downloader._osm_net_download(poly_buff, network_type, custom_filter)

        # create buffered graph from the downloaded data
        G_buff = _create_graph(
            response_jsons,
            retain_all=True,
            bidirectional=network_type in settings.bidirectional_network_types,
        )

        # truncate buffered graph to the buffered polygon and retain_all for
        # now. needed because overpass returns entire ways that also include
        # nodes outside the poly if the way (that is, a way with a single OSM
        # ID) has a node inside the poly at some point.
        G_buff = truncate.truncate_graph_polygon(
            G_buff, poly_buff, retain_all=True, truncate_by_edge=truncate_by_edge
        )

        # simplify the graph topology
        if simplify:
            G_buff = simplification.simplify_graph(G_buff)

        # truncate graph by original polygon to return graph within polygon
        # caller wants. don't simplify again: this allows us to retain
        # intersections along the street that may now only connect 2 street
        # segments in the network, but in reality also connect to an
        # intersection just outside the polygon
        G = truncate.truncate_graph_polygon(
            G_buff, polygon, retain_all=retain_all, truncate_by_edge=truncate_by_edge
        )

        # count how many street segments in buffered graph emanate from each
        # intersection in un-buffered graph, to retain true counts for each
        # intersection, even if some of its neighbors are outside the polygon
        G.graph["streets_per_node"] = utils_graph.count_streets_per_node(G_buff, nodes=G.nodes())

    # if clean_periphery=False, just use the polygon as provided
    else:
        # download the network data from OSM
        response_jsons = downloader._osm_net_download(polygon, network_type, custom_filter)

        # create graph from the downloaded data
        G = _create_graph(
            response_jsons,
            retain_all=True,
            bidirectional=network_type in settings.bidirectional_network_types,
        )

        # truncate the graph to the extent of the polygon
        G = truncate.truncate_graph_polygon(
            G, polygon, retain_all=retain_all, truncate_by_edge=truncate_by_edge
        )

        # simplify the graph topology as the last step. don't truncate after
        # simplifying or you may have simplified out to an endpoint beyond the
        # truncation distance, in which case you will then strip out your entire
        # edge
        if simplify:
            G = simplification.simplify_graph(G)

    utils.log(f"graph_from_polygon returned graph with {len(G)} nodes and {len(G.edges())} edges")
    return G


def graph_from_xml(filepath, bidirectional=False, simplify=True, retain_all=False):
    """
    Create a graph from data in an OSM-formatted XML file.

    Parameters
    ----------
    filepath : string
        path to file containing OSM XML data
    bidirectional : bool
        if True, create bidirectional edges for one-way streets
    simplify : bool
        if True, simplify the graph topology
    retain_all : bool
        if True, return the entire graph even if it is not connected

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    # transmogrify file of OSM XML data into JSON
    response_jsons = [_overpass_json_from_file(filepath)]

    # create graph using this response JSON
    G = _create_graph(response_jsons, bidirectional=bidirectional, retain_all=retain_all)

    # simplify the graph topology as the last step.
    if simplify:
        G = simplification.simplify_graph(G)

    utils.log(f"graph_from_xml returned graph with {len(G)} nodes and {len(G.edges())} edges")
    return G


def _overpass_json_from_file(filepath):
    """
    Read OSM XML from file and return Overpass-like JSON.

    Parameters
    ----------
    filepath : string
        path to file containing OSM XML data

    Returns
    -------
    OSMContentHandler object
    """

    def _opener(filepath):
        _, ext = os.path.splitext(filepath)
        if ext == ".bz2":
            # Use Python 2/3 compatible BZ2File()
            return bz2.BZ2File(filepath)
        else:
            # Assume an unrecognized file extension is just XML
            return open(filepath, mode="rb")

    with _opener(filepath) as file:
        handler = _OSMContentHandler()
        xml.sax.parse(file, handler)
        return handler.object


def _create_graph(response_jsons, retain_all=False, bidirectional=False):
    """
    Create a networkx MultiDiGraph from Overpass API responses.

    Parameters
    ----------
    response_jsons : list
        list of dicts of JSON responses from from the Overpass API
    retain_all : bool
        if True, return the entire graph even if it is not connected
    bidirectional : bool
        if True, create bidirectional edges for one-way streets

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    utils.log("Creating graph from downloaded OSM data...")

    # make sure we got data back from the server requests
    elements = []
    for response_json in response_jsons:
        elements.extend(response_json["elements"])
    if len(elements) < 1:
        raise EmptyOverpassResponse("There are no data elements in the response JSON")

    # create the graph as a MultiDiGraph and set its meta-attributes
    G = nx.MultiDiGraph(
        created_date=utils.ts(), created_with=f"OSMnx {__version__}", crs=settings.default_crs
    )

    # extract nodes and paths from the downloaded osm data
    nodes = {}
    paths = {}
    for osm_data in response_jsons:
        nodes_temp, paths_temp = _parse_osm_nodes_paths(osm_data)
        for key, value in nodes_temp.items():
            nodes[key] = value
        for key, value in paths_temp.items():
            paths[key] = value

    # add each osm node to the graph
    for node, data in nodes.items():
        G.add_node(node, **data)

    # add each osm way (ie, a path of edges) to the graph
    G = _add_paths(G, paths, bidirectional=bidirectional)

    # retain only the largest connected component, if caller did not
    # set retain_all=True
    if not retain_all:
        G = utils_graph.get_largest_component(G)

    utils.log(f"Created graph with {len(G)} nodes and {len(G.edges())} edges")

    # add length (great circle distance between nodes) attribute to each edge to
    # use as weight
    if len(G.edges) > 0:
        G = utils_graph.add_edge_lengths(G)

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
    node = {}
    node["y"] = element["lat"]
    node["x"] = element["lon"]
    node["osmid"] = element["id"]
    if "tags" in element:
        for useful_tag in settings.useful_tags_node:
            if useful_tag in element["tags"]:
                node[useful_tag] = element["tags"][useful_tag]
    return node


def _convert_path(element):
    """
    Convert an OSM way element into the format for a networkx graph path.

    Parameters
    ----------
    element : dict
        an OSM way element

    Returns
    -------
    path : dict
    """
    path = {}
    path["osmid"] = element["id"]

    # remove any consecutive duplicate elements in the list of nodes
    grouped_list = groupby(element["nodes"])
    path["nodes"] = [group[0] for group in grouped_list]

    if "tags" in element:
        for useful_tag in settings.useful_tags_way:
            if useful_tag in element["tags"]:
                path[useful_tag] = element["tags"][useful_tag]
    return path


def _parse_osm_nodes_paths(osm_data):
    """
    Construct dicts of nodes and paths.

    Dict key=osmid and value=dict of attributes.

    Parameters
    ----------
    osm_data : dict
        JSON response from from the Overpass API

    Returns
    -------
    nodes, paths : tuple of dicts
    """
    nodes = {}
    paths = {}
    for element in osm_data["elements"]:
        if element["type"] == "node":
            key = element["id"]
            nodes[key] = _convert_node(element)
        elif element["type"] == "way":
            key = element["id"]
            paths[key] = _convert_path(element)

    return nodes, paths


def _add_path(G, data, one_way):
    """
    Add a path to the graph.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    data : dict
        the attributes of the path
    one_way : bool
        if this path is one-way or if it is bi-directional

    Returns
    -------
    None
    """
    # extract the ordered list of nodes from this path element, then delete it
    # so we don't add it as an attribute to the edge later
    path_nodes = data["nodes"]
    del data["nodes"]

    # set the oneway attribute to the passed-in value, to make it consistent
    # True/False values, but only do this if you aren't forcing all edges to
    # oneway with the all_oneway setting. With the all_oneway setting, you
    # likely still want to preserve the original OSM oneway attribute.
    if not settings.all_oneway:
        data["oneway"] = one_way

    # zip together the path nodes so you get tuples like (0,1), (1,2), (2,3)
    # and so on
    path_edges = list(zip(path_nodes[:-1], path_nodes[1:]))
    G.add_edges_from(path_edges, **data)

    # if the path is NOT one-way
    if not one_way:
        # reverse the direction of each edge and add this path going the
        # opposite direction
        path_edges_opposite_direction = [(v, u) for u, v in path_edges]
        G.add_edges_from(path_edges_opposite_direction, **data)


def _add_paths(G, paths, bidirectional=False):
    """
    Add a collection of paths to the graph.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    paths : dict
        the paths from OSM
    bidirectional : bool
        if True, create bidirectional edges for one-way streets

    Returns
    -------
    None
    """
    # the list of values OSM uses in its 'oneway' tag to denote True
    # https://www.geofabrik.de/de/data/geofabrik-osm-gis-standard-0.7.pdf
    osm_oneway_values = ["yes", "true", "1", "-1", "T", "F"]

    for data in paths.values():

        if settings.all_oneway is True:
            _add_path(G, data, one_way=True)
        # if this path is tagged as one-way and if it is not a walking network,
        # then we'll add the path in one direction only
        elif ("oneway" in data and data["oneway"] in osm_oneway_values) and not bidirectional:
            if data["oneway"] == "-1" or data["oneway"] == "T":
                # paths with a one-way value of -1 or T are one-way, but in the
                # reverse direction of the nodes' order, see osm documentation
                data["nodes"] = list(reversed(data["nodes"]))
            # add this path (in only one direction) to the graph
            _add_path(G, data, one_way=True)

        elif ("junction" in data and data["junction"] == "roundabout") and not bidirectional:
            # roundabout are also oneway but not tagged as is
            _add_path(G, data, one_way=True)

        # else, this path is not tagged as one-way or it is a walking network
        # (you can walk both directions on a one-way street)
        else:
            # add this path (in both directions) to the graph and set its
            # 'oneway' attribute to False. if this is a walking network, this
            # may very well be a one-way street (as cars/bikes go), but in a
            # walking-only network it is a bi-directional edge
            _add_path(G, data, one_way=False)

    return G


class _OSMContentHandler(xml.sax.handler.ContentHandler):
    """
    SAX content handler for OSM XML.

    Used to build an Overpass-like response JSON object in self.object. For format
    notes, see http://wiki.openstreetmap.org/wiki/OSM_XML#OSM_XML_file_format_notes
    and http://overpass-api.de/output_formats.html#json
    """

    def __init__(self):
        self._element = None
        self.object = {"elements": []}

    def startElement(self, name, attrs):
        if name == "osm":
            self.object.update({k: attrs[k] for k in attrs.keys() if k in {"version", "generator"}})

        elif name in {"node", "way"}:
            self._element = dict(type=name, tags={}, nodes=[], **attrs)
            self._element.update({k: float(attrs[k]) for k in attrs.keys() if k in {"lat", "lon"}})
            self._element.update(
                {
                    k: int(attrs[k])
                    for k in attrs.keys()
                    if k in {"id", "uid", "version", "changeset"}
                }
            )

        elif name == "tag":
            self._element["tags"].update({attrs["k"]: attrs["v"]})

        elif name == "nd":
            self._element["nodes"].append(int(attrs["ref"]))

        elif name == "relation":
            # Placeholder for future relation support.
            # Look for nested members and tags.
            pass

    def endElement(self, name):
        if name in {"node", "way"}:
            self.object["elements"].append(self._element)
