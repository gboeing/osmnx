"""Graph utility functions."""

from warnings import warn

from . import convert
from . import routing
from . import truncate
from . import utils


def graph_to_gdfs(G, nodes=True, edges=True, node_geometry=True, fill_edge_geometry=True):
    """
    Do not use: deprecated.

    Use the `convert.graph_to_gdfs` function instead.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        deprecated, do not use
    nodes : bool
        deprecated, do not use
    edges : bool
        deprecated, do not use
    node_geometry : bool
        deprecated, do not use
    fill_edge_geometry : bool
        deprecated, do not use

    Returns
    -------
    geopandas.GeoDataFrame or tuple
    """
    msg = (
        "The `graph_to_gdfs` function has moved to the `convert` module. Calling "
        "`utils_graph.graph_to_gdfs` is deprecated and will be removed in the "
        "v2.0.0 release. Call it via `convert.graph_to_gdfs` instead. "
        "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
    )
    warn(msg, FutureWarning, stacklevel=2)
    return convert.graph_to_gdfs(G, nodes, edges, node_geometry, fill_edge_geometry)


def graph_from_gdfs(gdf_nodes, gdf_edges, graph_attrs=None):
    """
    Do not use: deprecated.

    Use the `convert.graph_from_gdfs` function instead.

    Parameters
    ----------
    gdf_nodes : geopandas.GeoDataFrame
        deprecated, do not use
    gdf_edges : geopandas.GeoDataFrame
        deprecated, do not use
    graph_attrs : dict
        deprecated, do not use

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    msg = (
        "The `graph_from_gdfs` function has moved to the `convert` module. Calling "
        "`utils_graph.graph_from_gdfs` is deprecated and will be removed in the "
        "v2.0.0 release. Call it via `convert.graph_from_gdfs` instead. "
        "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
    )
    warn(msg, FutureWarning, stacklevel=2)
    return convert.graph_from_gdfs(gdf_nodes, gdf_edges, graph_attrs)


def route_to_gdf(G, route, weight="length"):
    """
    Do not use: deprecated.

    Use the `routing.route_to_gdf` function instead.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        deprecated, do not use
    route : list
        deprecated, do not use
    weight : string
        deprecated, do not use

    Returns
    -------
    gdf_edges : geopandas.GeoDataFrame
    """
    msg = (
        "The `route_to_gdf` function has moved to the `routing` module. Calling "
        "`utils_graph.route_to_gdf` is deprecated and will be removed in the "
        "v2.0.0 release. Call it via `routing.route_to_gdf` instead. "
        "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
    )
    warn(msg, FutureWarning, stacklevel=2)
    return routing.route_to_gdf(G, route, weight)


def get_route_edge_attributes(
    G, route, attribute=None, minimize_key="length", retrieve_default=None
):
    """
    Do not use: deprecated.

    Use the `routing.route_to_gdf` function instead.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        deprecated, do not use
    route : list
        deprecated, do not use
    attribute : string
        deprecated, do not use
    minimize_key : string
        deprecated, do not use
    retrieve_default : Callable[Tuple[Any, Any], Any]
        deprecated, do not use

    Returns
    -------
    attribute_values : list
    """
    warn(
        "The `get_route_edge_attributes` function has been deprecated and will "
        "be removed in the v2.0.0 release. Use the `routing.route_to_gdf` function instead. "
        "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123",
        FutureWarning,
        stacklevel=2,
    )
    attribute_values = []
    for u, v in zip(route[:-1], route[1:]):
        # if there are parallel edges between two nodes, select the one with the
        # lowest value of minimize_key
        data = min(G.get_edge_data(u, v).values(), key=lambda x: x[minimize_key])
        if attribute is None:
            attribute_value = data
        elif retrieve_default is not None:
            attribute_value = data.get(attribute, retrieve_default(u, v))
        else:
            attribute_value = data[attribute]
        attribute_values.append(attribute_value)
    return attribute_values


def remove_isolated_nodes(G, warn=True):
    """
    Do not use: deprecated.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        deprecated, do not use
    warn : bool
        deprecated, do not use

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    if warn:
        msg = "The `remove_isolated_nodes` function is deprecated and will be removed in the v2.0.0 release."
        warn(msg, FutureWarning, stacklevel=2)

    # make a copy to not mutate original graph object caller passed in
    G = G.copy()

    # get the set of all isolated nodes, then remove them
    isolated_nodes = {node for node, degree in G.degree() if degree < 1}
    G.remove_nodes_from(isolated_nodes)
    utils.log(f"Removed {len(isolated_nodes):,} isolated nodes")
    return G


def get_largest_component(G, strongly=False):
    """
    Do not use: deprecated.

    Use the `truncate.largest_component` function instead.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        deprecated, do not use
    strongly : bool
        deprecated, do not use

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    msg = (
        "The `get_largest_component` function is deprecated and will be removed in the "
        "v2.0.0 release. Replace it with `truncate.largest_component` instead. "
        "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
    )
    warn(msg, FutureWarning, stacklevel=2)
    return truncate.largest_component(G, strongly)


def get_digraph(G, weight="length"):
    """
    Do not use: deprecated.

    Use the `convert.to_digraph` function instead.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        deprecated, do not use
    weight : string
        deprecated, do not use

    Returns
    -------
    networkx.DiGraph
    """
    msg = (
        "The `get_digraph` function is deprecated and will be removed in the "
        "v2.0.0 release. Replace it with `convert.to_digraph` instead. "
        "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
    )
    warn(msg, FutureWarning, stacklevel=2)
    return convert.to_digraph(G, weight)


def get_undirected(G):
    """
    Do not use: deprecated.

    Use the `convert.to_undirected` function instead.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        deprecated, do not use

    Returns
    -------
    networkx.MultiGraph
    """
    msg = (
        "The `get_undirected` function is deprecated and will be removed in the "
        "v2.0.0 release. Replace it with `convert.to_undirected` instead. "
        "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
    )
    warn(msg, FutureWarning, stacklevel=2)
    return convert.to_undirected(G)
