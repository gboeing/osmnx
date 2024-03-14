"""Create interactive Leaflet web maps of graphs and routes via folium.

This module is deprecated. Do not use. It will be removed in the v2.0.0 release.
You can generate and explore interactive web maps of graph nodes, edges,
and/or routes automatically using GeoPandas.GeoDataFrame.explore instead, for
example like: `ox.graph_to_gdfs(G, nodes=False).explore()`. See the OSMnx
examples gallery for complete details and demonstrations.
"""

import json
from warnings import warn

from . import convert

# folium is an optional dependency for the folium plotting functions
try:
    import folium
except ImportError:  # pragma: no cover
    folium = None


def plot_graph_folium(
    G,
    graph_map=None,
    popup_attribute=None,
    tiles="cartodbpositron",
    zoom=1,
    fit_bounds=True,
    **kwargs,
):
    """
    Do not use: deprecated.

    You can generate and explore interactive web maps of graph nodes, edges,
    and/or routes automatically using GeoPandas.GeoDataFrame.explore instead,
    for example like: `ox.graph_to_gdfs(G, nodes=False).explore()`. See the
    OSMnx examples gallery for complete details and demonstrations.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        deprecated
    graph_map : folium.folium.Map
        deprecated
    popup_attribute : string
        deprecated
    tiles : string
        deprecated
    zoom : int
        deprecated
    fit_bounds : bool
        deprecated
    kwargs
        deprecated

    Returns
    -------
    folium.folium.Map
    """
    warn(
        "The `folium` module has been deprecated and will be removed in the v2.0.0 release. "
        "You can generate and explore interactive web maps of graph nodes, edges, "
        "and/or routes automatically using GeoPandas.GeoDataFrame.explore instead, "
        "for example like: `ox.graph_to_gdfs(G, nodes=False).explore()`. See the "
        "OSMnx examples gallery for complete details and demonstrations.",
        FutureWarning,
        stacklevel=2,
    )
    # create gdf of all graph edges
    gdf_edges = convert.graph_to_gdfs(G, nodes=False)
    return _plot_folium(gdf_edges, graph_map, popup_attribute, tiles, zoom, fit_bounds, **kwargs)


def plot_route_folium(
    G,
    route,
    route_map=None,
    popup_attribute=None,
    tiles="cartodbpositron",
    zoom=1,
    fit_bounds=True,
    **kwargs,
):
    """
    Do not use: deprecated.

    You can generate and explore interactive web maps of graph nodes, edges,
    and/or routes automatically using GeoPandas.GeoDataFrame.explore instead,
    for example like: `ox.graph_to_gdfs(G, nodes=False).explore()`. See the
    OSMnx examples gallery for complete details and demonstrations.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        deprecated
    route : list
        deprecated
    route_map : folium.folium.Map
        deprecated
    popup_attribute : string
        deprecated
    tiles : string
        deprecated
    zoom : int
        deprecated
    fit_bounds : bool
        deprecated
    kwargs
        deprecated

    Returns
    -------
    folium.folium.Map
    """
    warn(
        "The `folium` module has been deprecated and will be removed in the v2.0.0 release. "
        "You can generate and explore interactive web maps of graph nodes, edges, "
        "and/or routes automatically using GeoPandas.GeoDataFrame.explore instead, "
        "for example like: `ox.graph_to_gdfs(G, nodes=False).explore()`. See the "
        "OSMnx examples gallery for complete details and demonstrations.",
        FutureWarning,
        stacklevel=2,
    )
    # create gdf of the route edges in order
    node_pairs = zip(route[:-1], route[1:])
    uvk = ((u, v, min(G[u][v].items(), key=lambda k: k[1]["length"])[0]) for u, v in node_pairs)
    gdf_edges = convert.graph_to_gdfs(G.subgraph(route), nodes=False).loc[uvk]
    return _plot_folium(gdf_edges, route_map, popup_attribute, tiles, zoom, fit_bounds, **kwargs)


def _plot_folium(gdf, m, popup_attribute, tiles, zoom, fit_bounds, **kwargs):
    """
    Plot a GeoDataFrame of LineStrings on a folium map object.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        a GeoDataFrame of LineString geometries and attributes
    m : folium.folium.Map or folium.FeatureGroup
        if not None, plot on this preexisting folium map object
    popup_attribute : string
        attribute to display in pop-up on-click, if None, no popup
    tiles : string
        name of a folium tileset
    zoom : int
        initial zoom level for the map
    fit_bounds : bool
        if True, fit the map to gdf's boundaries
    kwargs
        keyword arguments to pass to folium.PolyLine()

    Returns
    -------
    m : folium.folium.Map
    """
    # check if we were able to import folium successfully
    if folium is None:  # pragma: no cover
        msg = "folium must be installed to use this optional feature"
        raise ImportError(msg)

    # get centroid
    x, y = gdf.unary_union.centroid.xy
    centroid = (y[0], x[0])

    # create the folium web map if one wasn't passed-in
    if m is None:
        m = folium.Map(location=centroid, zoom_start=zoom, tiles=tiles)

    # identify the geometry and popup columns
    attrs = ["geometry"] if popup_attribute is None else ["geometry", popup_attribute]

    # add each edge to the map
    for vals in gdf[attrs].to_numpy():
        params = dict(zip(["geom", "popup_val"], vals))
        pl = _make_folium_polyline(**params, **kwargs)
        pl.add_to(m)

    # if fit_bounds is True, fit the map to the bounds of the route by passing
    # list of lat-lon points as [southwest, northeast]
    if fit_bounds and isinstance(m, folium.Map):
        tb = gdf.total_bounds
        m.fit_bounds([(tb[1], tb[0]), (tb[3], tb[2])])

    return m


def _make_folium_polyline(geom, popup_val=None, **kwargs):
    """
    Turn LineString geometry into a folium PolyLine with attributes.

    Parameters
    ----------
    geom : shapely LineString
        geometry of the line
    popup_val : string
        text to display in pop-up when a line is clicked, if None, no popup
    kwargs
        keyword arguments to pass to folium.PolyLine()

    Returns
    -------
    pl : folium.PolyLine
    """
    # locations is a list of points for the polyline folium takes coords in
    # lat,lon but geopandas provides them in lon,lat so we must reverse them
    locations = [(lat, lon) for lon, lat in geom.coords]

    # create popup if popup_val is not None
    # folium doesn't interpret html, so can't do newlines without iframe
    popup = None if popup_val is None else folium.Popup(html=json.dumps(popup_val))

    # create a folium polyline with attributes
    return folium.PolyLine(locations=locations, popup=popup, **kwargs)
