"""Create leaflet web maps via folium."""

import json

from . import utils_graph

# folium is an optional dependency for the folium plotting functions
try:
    import folium
except ImportError:
    folium = None


def _make_folium_polyline(
    edge, edge_color, edge_width, edge_opacity, popup_attribute=None, **kwargs
):
    """
    Turn row GeoDataFrame into a folium PolyLine with attributes.

    Parameters
    ----------
    edge : GeoSeries
        a row from the gdf_edges GeoDataFrame
    edge_color : string
        color of the edge lines
    edge_width : numeric
        width of the edge lines
    edge_opacity : numeric
        opacity of the edge lines
    popup_attribute : string
        edge attribute to display in a pop-up when an edge is clicked, if
        None, no popup
    kwargs : dict
        Extra parameters passed through to folium

    Returns
    -------
    pl : folium.PolyLine
    """
    # check if we were able to import folium successfully
    if not folium:
        raise ImportError("The folium package must be installed to use this optional feature.")

    # locations is a list of points for the polyline
    # folium takes coords in lat,lon but geopandas provides them in lon,lat
    # so we have to flip them around
    locations = list([(lat, lng) for lng, lat in edge["geometry"].coords])

    # if popup_attribute is None, then create no pop-up
    if popup_attribute is None:
        popup = None
    else:
        # folium doesn't interpret html in the html argument (weird), so can't
        # do newlines without an iframe
        popup_text = json.dumps(edge[popup_attribute])
        popup = folium.Popup(html=popup_text)

    # create a folium polyline with attributes
    pl = folium.PolyLine(
        locations=locations,
        popup=popup,
        color=edge_color,
        weight=edge_width,
        opacity=edge_opacity,
        **kwargs,
    )
    return pl


def plot_graph_folium(
    G,
    graph_map=None,
    popup_attribute=None,
    tiles="cartodbpositron",
    zoom=1,
    fit_bounds=True,
    edge_color="#333333",
    edge_width=5,
    edge_opacity=1,
    **kwargs,
):
    """
    Plot a graph on an interactive folium web map.

    Note that anything larger than a small city can take a long time to plot
    and create a large web map file that is very slow to load as JavaScript.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    graph_map : folium.folium.Map or folium.FeatureGroup
        if not None, plot the graph on this preexisting folium map object
    popup_attribute : string
        edge attribute to display in a pop-up when an edge is clicked
    tiles : string
        name of a folium tileset
    zoom : int
        initial zoom level for the map
    fit_bounds : bool
        if True, fit the map to the boundaries of the route's edges
    edge_color : string
        color of the edge lines
    edge_width : numeric
        width of the edge lines
    edge_opacity : numeric
        opacity of the edge lines
    kwargs : dict
        Extra keyword arguments passed through to folium

    Returns
    -------
    graph_map : folium.folium.Map
    """
    # check if we were able to import folium successfully
    if not folium:
        raise ImportError("The folium package must be installed to use this optional feature.")

    # create gdf of the graph edges
    gdf_edges = utils_graph.graph_to_gdfs(G, nodes=False, fill_edge_geometry=True)

    # get graph centroid
    x, y = gdf_edges.unary_union.centroid.xy
    graph_centroid = (y[0], x[0])

    # create the folium web map if one wasn't passed-in
    if graph_map is None:
        graph_map = folium.Map(location=graph_centroid, zoom_start=zoom, tiles=tiles)

    # add each graph edge to the map
    for _, row in gdf_edges.iterrows():
        pl = _make_folium_polyline(
            edge=row,
            edge_color=edge_color,
            edge_width=edge_width,
            edge_opacity=edge_opacity,
            popup_attribute=popup_attribute,
            **kwargs,
        )
        pl.add_to(graph_map)

    # if fit_bounds is True, fit the map to the bounds of the route by passing
    # list of lat-lng points as [southwest, northeast]
    if fit_bounds and isinstance(graph_map, folium.Map):
        tb = gdf_edges.total_bounds
        bounds = [[tb[1], tb[0]], [tb[3], tb[2]]]
        graph_map.fit_bounds(bounds)

    return graph_map


def plot_route_folium(
    G,
    route,
    route_map=None,
    popup_attribute=None,
    tiles="cartodbpositron",
    zoom=1,
    fit_bounds=True,
    route_color="#cc0000",
    route_width=5,
    route_opacity=1,
    **kwargs,
):
    """
    Plot a route on an interactive folium web map.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    route : list
        the route as a list of nodes
    route_map : folium.folium.Map
        if not None, plot the route on this preexisting folium map object
    popup_attribute : string
        edge attribute to display in a pop-up when an edge is clicked
    tiles : string
        name of a folium tileset
    zoom : int
        initial zoom level for the map
    fit_bounds : bool
        if True, fit the map to the boundaries of the route's edges
    route_color : string
        color of the route's line
    route_width : numeric
        width of the route's line
    route_opacity : numeric
        opacity of the route lines
    kwargs : dict
        Extra parameters passed through to folium

    Returns
    -------
    route_map : folium.folium.Map
    """
    # check if we were able to import folium successfully
    if not folium:
        raise ImportError("The folium package must be installed to use this optional feature.")

    # create gdf of the route edges
    gdf_edges = utils_graph.graph_to_gdfs(G, nodes=False, fill_edge_geometry=True)
    route_nodes = list(zip(route[:-1], route[1:]))
    index = [
        gdf_edges[(gdf_edges["u"] == u) & (gdf_edges["v"] == v)].index[0] for u, v in route_nodes
    ]
    gdf_route_edges = gdf_edges.loc[index]

    # get route centroid
    x, y = gdf_route_edges.unary_union.centroid.xy
    route_centroid = (y[0], x[0])

    # create the folium web map if one wasn't passed-in
    if route_map is None:
        route_map = folium.Map(location=route_centroid, zoom_start=zoom, tiles=tiles)

    # add each route edge to the map
    for _, row in gdf_route_edges.iterrows():
        pl = _make_folium_polyline(
            edge=row,
            edge_color=route_color,
            edge_width=route_width,
            edge_opacity=route_opacity,
            popup_attribute=popup_attribute,
            **kwargs,
        )
        pl.add_to(route_map)

    # if fit_bounds is True, fit the map to the bounds of the route by passing
    # list of lat-lng points as [southwest, northeast]
    if fit_bounds and isinstance(route_map, folium.Map):
        tb = gdf_route_edges.total_bounds
        bounds = [(tb[1], tb[0]), (tb[3], tb[2])]
        route_map.fit_bounds(bounds)

    return route_map
