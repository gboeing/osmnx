"""Visualize street networks, routes, orientations, and geospatial features."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from typing import Literal
from typing import overload

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd

from . import bearing
from . import graph
from . import projection
from . import settings
from . import simplification
from . import utils
from . import utils_geo
from . import utils_graph

# matplotlib is an optional dependency needed for visualization
try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib import colormaps
    from matplotlib import colors
    from matplotlib.axes._axes import Axes
    from matplotlib.figure import Figure
    from matplotlib.projections.polar import PolarAxes
except ImportError:  # pragma: no cover
    plt = cm = colormaps = colors = None  # type: ignore[assignment]


def get_colors(
    n: int,
    cmap: str = "viridis",
    start: float = 0,
    stop: float = 1,
    alpha: float = 1,
    return_hex: bool = False,
) -> list:
    """
    Get `n` evenly-spaced colors from a matplotlib colormap.

    Parameters
    ----------
    n : int
        number of colors
    cmap : string
        name of a matplotlib colormap
    start : float
        where to start in the colorspace
    stop : float
        where to end in the colorspace
    alpha : float
        opacity, the alpha channel for the RGBa colors
    return_hex : bool
        if True, convert RGBa colors to HTML-like hexadecimal RGB strings. if
        False, return colors as (R, G, B, alpha) tuples.

    Returns
    -------
    color_list : list
    """
    _verify_mpl()
    color_list = [colormaps[cmap](x) for x in np.linspace(start, stop, n)]
    if return_hex:
        color_list = [colors.to_hex(c) for c in color_list]
    else:
        color_list = [(r, g, b, alpha) for r, g, b, _ in color_list]
    return color_list


def get_node_colors_by_attr(
    G: nx.MultiDiGraph,
    attr: str,
    num_bins: int = None,
    cmap: str = "viridis",
    start: float = 0,
    stop: float = 1,
    na_color: str = "none",
    equal_size: bool = False,
) -> pd.Series:
    """
    Get colors based on node attribute values.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    attr : string
        name of a numerical node attribute
    num_bins : int
        if None, linearly map a color to each value. otherwise, assign values
        to this many bins then assign a color to each bin.
    cmap : string
        name of a matplotlib colormap
    start : float
        where to start in the colorspace
    stop : float
        where to end in the colorspace
    na_color : string
        what color to assign nodes with missing attr values
    equal_size : bool
        ignored if num_bins is None. if True, bin into equal-sized quantiles
        (requires unique bin edges). if False, bin into equal-spaced bins.

    Returns
    -------
    node_colors : pandas.Series
        series labels are node IDs and values are colors
    """
    _verify_mpl()
    vals = pd.Series(nx.get_node_attributes(G, attr))
    return _get_colors_by_value(vals, num_bins, cmap, start, stop, na_color, equal_size)


def get_edge_colors_by_attr(
    G: nx.MultiDiGraph,
    attr: str,
    num_bins: int = None,
    cmap: str = "viridis",
    start: float = 0,
    stop: float = 1,
    na_color: str = "none",
    equal_size: bool = False,
) -> pd.Series:
    """
    Get colors based on edge attribute values.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    attr : string
        name of a numerical edge attribute
    num_bins : int
        if None, linearly map a color to each value. otherwise, assign values
        to this many bins then assign a color to each bin.
    cmap : string
        name of a matplotlib colormap
    start : float
        where to start in the colorspace
    stop : float
        where to end in the colorspace
    na_color : string
        what color to assign edges with missing attr values
    equal_size : bool
        ignored if num_bins is None. if True, bin into equal-sized quantiles
        (requires unique bin edges). if False, bin into equal-spaced bins.

    Returns
    -------
    edge_colors : pandas.Series
        series labels are edge IDs (u, v, key) and values are colors
    """
    _verify_mpl()
    vals = pd.Series(nx.get_edge_attributes(G, attr))
    return _get_colors_by_value(vals, num_bins, cmap, start, stop, na_color, equal_size)


def plot_graph(
    G: nx.MultiDiGraph,
    ax: Axes = None,
    figsize: tuple = (8, 8),
    bgcolor: str = "#111111",
    node_color: str | list = "w",
    node_size: int | list = 15,
    node_alpha: float = None,
    node_edgecolor: str | list = "none",
    node_zorder: int = 1,
    edge_color: str | list = "#999999",
    edge_linewidth: int | list = 1,
    edge_alpha: float = None,
    show: bool = True,
    close: bool = False,
    save: bool = False,
    filepath: str = None,
    dpi: int = 300,
    bbox: tuple = None,
) -> tuple:
    """
    Visualize a graph.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    ax : matplotlib axes instance
        if not None, plot on this pre-existing axes instance
    figsize : tuple
        if ax is None, create new figure with size (width, height)
    bgcolor : string
        background color of plot
    node_color : string or list
        color(s) of the nodes
    node_size : int or list
        size(s) of the nodes: if 0, then skip plotting the nodes
    node_alpha : float
        opacity of the nodes, note: if you passed RGBA values to node_color,
        set node_alpha=None to use the alpha channel in node_color
    node_edgecolor : string or list
        color(s) of the nodes' markers' borders
    node_zorder : int
        zorder to plot nodes: edges are always 1, so set node_zorder=0 to plot
        nodes below edges
    edge_color : string or list
        color(s) of the edges' lines
    edge_linewidth : int or list
        width(s) of the edges' lines: if 0, then skip plotting the edges
    edge_alpha : float
        opacity of the edges, note: if you passed RGBA values to edge_color,
        set edge_alpha=None to use the alpha channel in edge_color
    show : bool
        if True, call pyplot.show() to show the figure
    close : bool
        if True, call pyplot.close() to close the figure
    save : bool
        if True, save the figure to disk at filepath
    filepath : string
        if save is True, the path to the file. file format determined from
        extension. if None, use settings.imgs_folder/image.png
    dpi : int
        if save is True, the resolution of saved file
    bbox : tuple
        bounding box as (north, south, east, west). if None, will calculate
        from spatial extents of plotted geometries.

    Returns
    -------
    fig, ax : tuple
        tuple of matplotlib (Figure, Axes)
    """
    _verify_mpl()
    max_node_size = max(node_size) if hasattr(node_size, "__iter__") else node_size
    max_edge_lw = max(edge_linewidth) if hasattr(edge_linewidth, "__iter__") else edge_linewidth
    if max_node_size <= 0 and max_edge_lw <= 0:  # pragma: no cover
        msg = "Either `node_size` or `edge_linewidth` must be > 0 to plot something."
        raise ValueError(msg)

    # create fig, ax as needed
    utils.log("Begin plotting the graph...")
    fig, ax = _get_fig_ax(ax=ax, figsize=figsize, bgcolor=bgcolor, polar=False)

    if max_edge_lw > 0:
        # plot the edges' geometries
        gdf_edges = utils_graph.graph_to_gdfs(G, nodes=False)["geometry"]
        ax = gdf_edges.plot(ax=ax, color=edge_color, lw=edge_linewidth, alpha=edge_alpha, zorder=1)

    if max_node_size > 0:
        # scatter plot the nodes' x/y coordinates
        gdf_nodes = utils_graph.graph_to_gdfs(G, edges=False, node_geometry=False)[["x", "y"]]
        ax.scatter(  # type: ignore[union-attr]
            x=gdf_nodes["x"],
            y=gdf_nodes["y"],
            s=node_size,
            c=node_color,
            alpha=node_alpha,
            edgecolor=node_edgecolor,
            zorder=node_zorder,
        )

    # get spatial extents from bbox parameter or the edges' geometries
    padding = 0.0
    if bbox is None:
        try:
            west, south, east, north = gdf_edges.total_bounds
        except NameError:
            west, south = gdf_nodes.min()
            east, north = gdf_nodes.max()
        bbox = north, south, east, west
        padding = 0.02  # pad 2% to not cut off peripheral nodes' circles

    # configure axes appearance, save/show figure as specified, and return
    ax = _config_ax(ax, G.graph["crs"], bbox, padding)  # type: ignore[arg-type]
    fig, ax = _save_and_show(fig, ax, save, show, close, filepath, dpi)
    utils.log("Finished plotting the graph")
    return fig, ax


def plot_graph_route(
    G: nx.MultiDiGraph,
    route: list,
    route_color: str = "r",
    route_linewidth: int = 4,
    route_alpha: float = 0.5,
    orig_dest_size: int = 100,
    ax: Axes = None,
    **pg_kwargs: Any,
) -> tuple:
    """
    Visualize a route along a graph.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    route : list
        route as a list of node IDs
    route_color : string
        color of the route
    route_linewidth : int
        width of the route line
    route_alpha : float
        opacity of the route line
    orig_dest_size : int
        size of the origin and destination nodes
    ax : matplotlib axes instance
        if not None, plot on this pre-existing axes instance
    pg_kwargs
        keyword arguments to pass to plot_graph

    Returns
    -------
    fig, ax : tuple
        tuple of matplotlib (Figure, Axes)
    """
    _verify_mpl()
    if ax is None:
        # plot the graph but not the route, and override any user show/close
        # args for now: we'll do that later
        overrides = {"show", "save", "close"}
        kwargs = {k: v for k, v in pg_kwargs.items() if k not in overrides}
        fig, ax = plot_graph(G, show=False, save=False, close=False, **kwargs)
    else:
        fig = ax.figure

    # scatterplot origin and destination points (first/last nodes in route)
    od_x = (G.nodes[route[0]]["x"], G.nodes[route[-1]]["x"])
    od_y = (G.nodes[route[0]]["y"], G.nodes[route[-1]]["y"])
    ax.scatter(od_x, od_y, s=orig_dest_size, c=route_color, alpha=route_alpha, edgecolor="none")

    # assemble the route edge geometries' x and y coords then plot the line
    x = []
    y = []
    for u, v in zip(route[:-1], route[1:]):
        # if there are parallel edges, select the shortest in length
        data = min(G.get_edge_data(u, v).values(), key=lambda d: d["length"])
        if "geometry" in data:
            # if geometry attribute exists, add all its coords to list
            xs, ys = data["geometry"].xy
            x.extend(xs)
            y.extend(ys)
        else:
            # otherwise, the edge is a straight line from node to node
            x.extend((G.nodes[u]["x"], G.nodes[v]["x"]))
            y.extend((G.nodes[u]["y"], G.nodes[v]["y"]))
    ax.plot(x, y, c=route_color, lw=route_linewidth, alpha=route_alpha)

    # save and show the figure as specified, passing relevant kwargs
    sas_kwargs = {"save", "show", "close", "filepath", "file_format", "dpi"}
    kwargs = {k: v for k, v in pg_kwargs.items() if k in sas_kwargs}
    fig, ax = _save_and_show(fig, ax, **kwargs)
    return fig, ax


def plot_graph_routes(
    G: nx.MultiDiGraph,
    routes: list[list],
    route_colors: str | list = "r",
    route_linewidths: int | list = 4,
    **pgr_kwargs: Any,
) -> tuple:
    """
    Visualize several routes along a graph.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    routes : list
        a list of lists of node IDs
    route_colors : string or list
        if string, 1 color for all routes. if list, the colors for each
        route.
    route_linewidths : int or list
        if int, 1 linewidth for all routes. if list, the linewidth for
        each route.
    pgr_kwargs
        keyword arguments to pass to plot_graph_route

    Returns
    -------
    fig, ax : tuple
        tuple of matplotlib (Figure, Axes)
    """
    _verify_mpl()

    # check for valid arguments
    if not all(isinstance(r, list) for r in routes):  # pragma: no cover
        msg = "routes must be a list of route lists"
        raise ValueError(msg)
    if len(routes) <= 1:  # pragma: no cover
        msg = "You must pass more than 1 route"
        raise ValueError(msg)
    if isinstance(route_colors, str):
        route_colors = [route_colors] * len(routes)
    if len(routes) != len(route_colors):  # pragma: no cover
        msg = "route_colors list must have same length as routes"
        raise ValueError(msg)
    if isinstance(route_linewidths, int):
        route_linewidths = [route_linewidths] * len(routes)
    if len(routes) != len(route_linewidths):  # pragma: no cover
        msg = "route_linewidths list must have same length as routes"
        raise ValueError(msg)

    # plot the graph and the first route
    overrides = {"route", "route_color", "route_linewidth", "show", "save", "close"}
    kwargs = {k: v for k, v in pgr_kwargs.items() if k not in overrides}
    fig, ax = plot_graph_route(
        G,
        route=routes[0],
        route_color=route_colors[0],
        route_linewidth=route_linewidths[0],
        show=False,
        save=False,
        close=False,
        **kwargs,
    )

    # plot the subsequent routes on top of existing ax
    overrides.update({"ax"})
    kwargs = {k: v for k, v in pgr_kwargs.items() if k not in overrides}
    r_rc_rlw = zip(routes[1:], route_colors[1:], route_linewidths[1:])
    for route, route_color, route_linewidth in r_rc_rlw:
        fig, ax = plot_graph_route(
            G,
            route=route,
            route_color=route_color,
            route_linewidth=route_linewidth,
            show=False,
            save=False,
            close=False,
            ax=ax,
            **kwargs,
        )

    # save and show the figure as specified, passing relevant kwargs
    sas_kwargs = {"save", "show", "close", "filepath", "file_format", "dpi"}
    kwargs = {k: v for k, v in pgr_kwargs.items() if k in sas_kwargs}
    fig, ax = _save_and_show(fig, ax, **kwargs)
    return fig, ax


def plot_figure_ground(
    G: nx.MultiDiGraph = None,
    address: str = None,
    point: tuple = None,
    dist: float = 805,
    network_type: str = "drive_service",
    street_widths: dict = None,
    default_width: int = 4,
    figsize: tuple = (8, 8),
    edge_color: str = "w",
    smooth_joints: bool = True,
    **pg_kwargs: Any,
) -> tuple:
    """
    Plot a figure-ground diagram of a street network.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph, must be unprojected
    address : string
        address to geocode as the center point if G is not passed in
    point : tuple
        the (lat, lon) center point if address and G are not passed in
    dist : float
        how many meters to extend north, south, east, west from center point
    network_type : string
        what type of street network to get
    street_widths : dict
        dict keys are street types and values are widths to plot in pixels
    default_width : int
        fallback width for any street type not in street_widths
    figsize : tuple
        tuple of (width, height) of figure (values should be equal)
    edge_color : string
        color of the edges' lines
    smooth_joints : bool
        if True, plot nodes same width as streets to smooth line joints and
        prevent cracks between them from showing
    pg_kwargs
        keyword arguments to pass to plot_graph

    Returns
    -------
    fig, ax : tuple
        tuple of matplotlib (Figure, Axes)
    """
    _verify_mpl()
    multiplier = 1.2

    # if user did not pass in custom street widths, create a dict of defaults
    if street_widths is None:
        street_widths = {
            "footway": 1.5,
            "steps": 1.5,
            "pedestrian": 1.5,
            "service": 1.5,
            "path": 1.5,
            "track": 1.5,
            "motorway": 6,
        }

    # if G was passed in, plot it centered on its node centroid
    if G is not None:
        gdf_nodes = utils_graph.graph_to_gdfs(G, edges=False, node_geometry=True)
        lonlat_point = gdf_nodes.unary_union.centroid.coords[0]
        point = tuple(reversed(lonlat_point))

    # otherwise get network by address or point (whichever was passed) using a
    # dist multiplier to ensure we get more than enough network. simplify in
    # non-strict mode to not combine multiple street types into single edge
    elif address is not None:
        G, point = graph.graph_from_address(
            address,
            dist=dist * multiplier,
            dist_type="bbox",
            network_type=network_type,
            simplify=False,
            truncate_by_edge=True,
            return_coords=True,
        )
        G = simplification.simplify_graph(G, strict=False)
    elif point is not None:
        G = graph.graph_from_point(
            point,
            dist=dist * multiplier,
            dist_type="bbox",
            network_type=network_type,
            simplify=False,
            truncate_by_edge=True,
        )
        G = simplification.simplify_graph(G, strict=False)
    else:  # pragma: no cover
        msg = "You must pass an address or lat-lon point or graph."
        raise ValueError(msg)

    # we need an undirected graph to find every edge incident on a node
    Gu = utils_graph.get_undirected(G)

    # for each edge, get a linewidth according to street type
    edge_linewidths = []
    for _, _, d in Gu.edges(keys=False, data=True):
        street_type = d["highway"][0] if isinstance(d["highway"], list) else d["highway"]
        if street_type in street_widths:
            edge_linewidths.append(street_widths[street_type])
        else:
            edge_linewidths.append(default_width)

    if smooth_joints:
        # for each node, get a nodesize according to the narrowest incident edge
        node_widths = {}
        for node in Gu.nodes:
            # first, identify all the highway types of this node's incident edges
            ie_data = [Gu.get_edge_data(node, nbr) for nbr in Gu.neighbors(node)]
            edge_types = [d[min(d)]["highway"] for d in ie_data]
            if not edge_types:
                # if node has no incident edges, make size zero
                node_widths[node] = 0
            else:
                # flatten the list of edge types
                et_flat = []
                for et in edge_types:
                    if isinstance(et, list):
                        et_flat.extend(et)
                    else:
                        et_flat.append(et)

                # lookup corresponding width for each edge type in flat list
                edge_widths = [street_widths.get(et, default_width) for et in et_flat]

                # node diameter should equal largest edge width to make joints
                # perfectly smooth. alternatively use min(?) to prevent
                # anything larger from extending past smallest street's line.
                # circle marker sizes are in area, so use diameter squared.
                circle_diameter = max(edge_widths)
                circle_area = circle_diameter**2
                node_widths[node] = circle_area

        # assign the node size to each node in the graph
        node_sizes: list | int = [node_widths[node] for node in Gu.nodes]
    else:
        node_sizes = 0

    # define the view extents of the plotting figure
    bbox = utils_geo.bbox_from_point(point, dist, project_utm=False)  # type: ignore[arg-type]

    # plot the figure
    overrides = {"bbox", "node_size", "node_color", "edge_linewidth"}
    kwargs = {k: v for k, v in pg_kwargs.items() if k not in overrides}
    fig, ax = plot_graph(
        G=Gu,
        bbox=bbox,
        figsize=figsize,
        node_size=node_sizes,
        node_color=edge_color,
        edge_color=edge_color,
        edge_linewidth=edge_linewidths,
        **kwargs,
    )
    return fig, ax


def plot_footprints(
    gdf: gpd.GeoDataFrame,
    ax: Axes = None,
    figsize: tuple = (8, 8),
    color: str = "orange",
    edge_color: str = "none",
    edge_linewidth: float = 0,
    alpha: float = None,
    bgcolor: str = "#111111",
    bbox: tuple = None,
    save: bool = False,
    show: bool = True,
    close: bool = False,
    filepath: str = None,
    dpi: int = 600,
) -> tuple:
    """
    Visualize a GeoDataFrame of geospatial features' footprints.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame of footprints (shapely Polygons and MultiPolygons)
    ax : matplotlib axes instance
        if not None, plot on this pre-existing axes instance
    figsize : tuple
        if ax is None, create new figure with size (width, height)
    color : string
        color of the footprints
    edge_color : string
        color of the edge of the footprints
    edge_linewidth : float
        width of the edge of the footprints
    alpha : float
        opacity of the footprints
    bgcolor : string
        background color of the plot
    bbox : tuple
        bounding box as (north, south, east, west). if None, will calculate
        from the spatial extents of the geometries in gdf
    save : bool
        if True, save the figure to disk at filepath
    show : bool
        if True, call pyplot.show() to show the figure
    close : bool
        if True, call pyplot.close() to close the figure
    filepath : string
        if save is True, the path to the file. file format determined from
        extension. if None, use settings.imgs_folder/image.png
    dpi : int
        if save is True, the resolution of saved file

    Returns
    -------
    fig, ax : tuple
        tuple of matplotlib (Figure, Axes)
    """
    _verify_mpl()
    fig, ax = _get_fig_ax(ax=ax, figsize=figsize, bgcolor=bgcolor, polar=False)

    # retain only Polygons and MultiPolygons, then plot
    gdf = gdf[gdf["geometry"].type.isin({"Polygon", "MultiPolygon"})]
    ax = gdf.plot(
        ax=ax, facecolor=color, edgecolor=edge_color, linewidth=edge_linewidth, alpha=alpha
    )

    # determine figure extents
    if bbox is None:
        west, south, east, north = gdf.total_bounds
    else:
        north, south, east, west = bbox

    # configure axes appearance, save/show figure as specified, and return
    ax = _config_ax(ax, gdf.crs, (north, south, east, west), 0)  # type: ignore[arg-type]
    fig, ax = _save_and_show(fig, ax, save, show, close, filepath, dpi)
    return fig, ax


def plot_orientation(
    Gu: nx.MultiGraph,
    num_bins: int = 36,
    min_length: float = 0,
    weight: str = None,
    ax: PolarAxes = None,
    figsize: tuple = (5, 5),
    area: bool = True,
    color: str = "#003366",
    edgecolor: str = "k",
    linewidth: float = 0.5,
    alpha: float = 0.7,
    title: str = None,
    title_y: float = 1.05,
    title_font: dict = None,
    xtick_font: dict = None,
) -> tuple:
    """
    Plot a polar histogram of a spatial network's bidirectional edge bearings.

    Ignores self-loop edges as their bearings are undefined.

    For more info see: Boeing, G. 2019. "Urban Spatial Order: Street Network
    Orientation, Configuration, and Entropy." Applied Network Science, 4 (1),
    67. https://doi.org/10.1007/s41109-019-0189-1

    Parameters
    ----------
    Gu : networkx.MultiGraph
        undirected, unprojected graph with `bearing` attributes on each edge
    num_bins : int
        number of bins; for example, if `num_bins=36` is provided, then each
        bin will represent 10 degrees around the compass
    min_length : float
        ignore edges with `length` attributes less than `min_length`
    weight : string
        if not None, weight edges' bearings by this (non-null) edge attribute
    ax : matplotlib.projections.polar.PolarAxes
        if not None, plot on this pre-existing axes instance (must have
        projection=polar)
    figsize : tuple
        if ax is None, create new figure with size (width, height)
    area : bool
        if True, set bar length so area is proportional to frequency,
        otherwise set bar length so height is proportional to frequency
    color : string
        color of histogram bars
    edgecolor : string
        color of histogram bar edges
    linewidth : float
        width of histogram bar edges
    alpha : float
        opacity of histogram bars
    title : string
        title for plot
    title_y : float
        y position to place title
    title_font : dict
        the title's fontdict to pass to matplotlib
    xtick_font : dict
        the xtick labels' fontdict to pass to matplotlib

    Returns
    -------
    fig, ax : tuple
        tuple of matplotlib (Figure, PolarAxes)
    """
    _verify_mpl()

    if title_font is None:
        title_font = {"family": "DejaVu Sans", "size": 24, "weight": "bold"}
    if xtick_font is None:
        xtick_font = {
            "family": "DejaVu Sans",
            "size": 10,
            "weight": "bold",
            "alpha": 1.0,
            "zorder": 3,
        }

    # get the bearings' distribution's bin counts and edges
    bin_counts, bin_edges = bearing._bearings_distribution(Gu, num_bins, min_length, weight)

    # positions: where to center each bar. ignore the last bin edge, because
    # it's the same as the first (i.e., 0 degrees = 360 degrees)
    positions = np.radians(bin_edges[:-1])

    # width: make bars fill the circumference without gaps or overlaps
    width = 2 * np.pi / num_bins

    # radius: how long to make each bar. set bar length so either the bar area
    # (ie, via sqrt) or the bar height is proportional to the bin's frequency
    bin_frequency = bin_counts / bin_counts.sum()
    radius = np.sqrt(bin_frequency) if area else bin_frequency

    # create PolarAxes (if not passed-in) then set N at top and go clockwise
    fig, ax = _get_fig_ax(ax=ax, figsize=figsize, bgcolor=None, polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction("clockwise")
    ax.set_ylim(top=radius.max())

    # configure the y-ticks and remove their labels
    ax.set_yticks(np.linspace(0, radius.max(), 5))
    ax.set_yticklabels(labels="")

    # configure the x-ticks and their labels
    xticklabels = ["N", "", "E", "", "S", "", "W", ""]
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(labels=xticklabels, fontdict=xtick_font)
    ax.tick_params(axis="x", which="major", pad=-2)

    # draw the bars
    ax.bar(
        positions,
        height=radius,
        width=width,
        align="center",
        bottom=0,
        zorder=2,
        color=color,
        edgecolor=edgecolor,
        linewidth=linewidth,
        alpha=alpha,
    )

    if title:
        ax.set_title(title, y=title_y, fontdict=title_font)
    fig.tight_layout()
    return fig, ax


def _get_colors_by_value(
    vals: pd.Series,
    num_bins: int | None,
    cmap: str,
    start: float,
    stop: float,
    na_color: str,
    equal_size: bool,
) -> pd.Series:
    """
    Map colors to the values in a series.

    Parameters
    ----------
    vals : pandas.Series
        series labels are node/edge IDs and values are attribute values
    num_bins : int
        if None, linearly map a color to each value. otherwise, assign values
        to this many bins then assign a color to each bin.
    cmap : string
        name of a matplotlib colormap
    start : float
        where to start in the colorspace
    stop : float
        where to end in the colorspace
    na_color : string
        what color to assign to missing values
    equal_size : bool
        ignored if num_bins is None. if True, bin into equal-sized quantiles
        (requires unique bin edges). if False, bin into equal-spaced bins.

    Returns
    -------
    color_series : pandas.Series
        series labels are node/edge IDs and values are colors
    """
    if len(vals) == 0:
        msg = "There are no attribute values."
        raise ValueError(msg)

    if num_bins is None:
        # calculate min/max values based on start/stop and data range
        vals_min = vals.dropna().min()
        vals_max = vals.dropna().max()
        full_range = (vals_max - vals_min) / (stop - start)
        full_min = vals_min - full_range * start
        full_max = full_min + full_range

        # linearly map a color to each attribute value
        normalizer = colors.Normalize(full_min, full_max)
        scalar_mapper = cm.ScalarMappable(normalizer, colormaps[cmap])
        color_series = vals.map(scalar_mapper.to_rgba)
        color_series.loc[pd.isna(vals)] = na_color

    else:
        # otherwise, bin values then assign colors to bins
        cut_func = pd.qcut if equal_size else pd.cut
        bins = cut_func(vals, num_bins, labels=range(num_bins))  # type: ignore[operator]
        bin_colors = get_colors(num_bins, cmap, start, stop)
        color_list = [bin_colors[b] if pd.notna(b) else na_color for b in bins]
        color_series = pd.Series(color_list, index=bins.index)

    return color_series


def _save_and_show(
    fig: Figure,
    ax: Axes,
    save: bool = False,
    show: bool = True,
    close: bool = True,
    filepath: str = None,
    dpi: int = 300,
) -> tuple:
    """
    Save a figure to disk and/or show it, as specified by args.

    Parameters
    ----------
    fig : matplotlib figure
        matplotlib figure
    ax : matplotlib axes instance
        matplotlib axes
    save : bool
        if True, save the figure to disk at filepath
    show : bool
        if True, call pyplot.show() to show the figure
    close : bool
        if True, call pyplot.close() to close the figure
    filepath : string
        if save is True, the path to the file. file format determined from
        extension. if None, use settings.imgs_folder/image.png
    dpi : int
        if save is True, the resolution of saved file

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axes
    """
    fig.canvas.draw()
    fig.canvas.flush_events()

    if save:
        # default filepath, if none provided
        fp = Path(settings.imgs_folder) / "image.png" if filepath is None else Path(filepath)

        # if save folder does not already exist, create it
        fp.parent.mkdir(parents=True, exist_ok=True)

        # get the file extension and figure facecolor
        ext = fp.suffix.strip(".")
        fc = fig.get_facecolor()

        if ext == "svg":
            # if the file format is svg, prep the fig/ax for saving
            ax.axis("off")
            ax.set_position((0, 0, 1, 1))
            ax.patch.set_alpha(0.0)
            fig.patch.set_alpha(0.0)
            fig.savefig(fp, bbox_inches=0, format=ext, facecolor=fc, transparent=True)
        else:
            # constrain saved figure's extent to interior of the axes
            extent = ax.bbox.transformed(fig.dpi_scale_trans.inverted())

            # temporarily turn figure frame on to save with facecolor
            fig.set_frameon(True)
            fig.savefig(fp, dpi=dpi, bbox_inches=extent, format=ext, facecolor=fc, transparent=True)
            fig.set_frameon(False)  # and turn it back off again
        utils.log(f"Saved figure to disk at {fp}")

    if show:
        plt.show()

    if close:
        plt.close()

    return fig, ax


def _config_ax(ax: Axes, crs: Any, bbox: tuple, padding: float) -> Axes:
    """
    Configure matplotlib axes instance for display.

    Parameters
    ----------
    ax : matplotlib axes instance
        the axes instance containing the plot
    crs : dict or string or pyproj.CRS
        the CRS of the plotted geometries
    bbox : tuple
        bounding box as (north, south, east, west)
    padding : float
        relative padding to add around the plot's bbox

    Returns
    -------
    ax : matplotlib axes
        the configured/styled axes instance
    """
    # set the axes view limits to bbox + relative padding
    north, south, east, west = bbox
    padding_ns = (north - south) * padding
    padding_ew = (east - west) * padding
    ax.set_ylim((south - padding_ns, north + padding_ns))
    ax.set_xlim((west - padding_ew, east + padding_ew))

    # set margins to zero, point ticks inward, turn off ax border and x/y axis
    # so there is no space around the plot
    ax.margins(0)
    ax.tick_params(which="both", direction="in")
    _ = [s.set_visible(False) for s in ax.spines.values()]  # type: ignore[func-returns-value]
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    # set aspect ratio
    if projection.is_projected(crs):
        # if projected, make equal aspect ratio
        ax.set_aspect("equal")
    else:
        # if not projected, conform aspect ratio to not stretch plot
        cos_lat = np.cos((south + north) / 2 / 180 * np.pi)
        ax.set_aspect(1 / cos_lat)

    return ax


# if polar = False, return Axes
@overload
def _get_fig_ax(
    ax: Axes | None, figsize: tuple, bgcolor: str | None, polar: Literal[False]
) -> tuple[Figure, Axes]:
    ...


# if polar = True, return PolarAxes
@overload
def _get_fig_ax(
    ax: Axes | None, figsize: tuple, bgcolor: str | None, polar: Literal[True]
) -> tuple[Figure, PolarAxes]:
    ...


def _get_fig_ax(
    ax: Axes | None, figsize: tuple, bgcolor: str | None, polar: bool
) -> tuple[Figure, Axes | PolarAxes]:
    """
    Generate a matplotlib Figure and (Polar)Axes or return existing ones.

    Parameters
    ----------
    ax : matplotlib axes instance
        if not None, plot on this pre-existing axes instance
    figsize : tuple
        if ax is None, create new figure with size (width, height)
    bgcolor : string
        background color of plot
    polar : bool
        if True, generate PolarAxes instead of Axes

    Returns
    -------
    fig, ax : tuple
        tuple of matplotlib (Figure, Axes) or (Figure, PolarAxes)
    """
    if ax is None:
        if polar:
            # make PolarAxes
            fig, ax = plt.subplots(figsize=figsize, subplot_kw={"projection": "polar"})
        else:
            # make regular Axes
            fig, ax = plt.subplots(figsize=figsize, facecolor=bgcolor, frameon=False)
            ax.set_facecolor(bgcolor)
    else:
        fig = ax.figure  # type: ignore[assignment]

    return fig, ax


def _verify_mpl() -> None:
    """
    Verify that matplotlib is installed and successfully imported.

    Returns
    -------
    None
    """
    if cm is None or colors is None or plt is None or colormaps is None:  # pragma: no cover
        msg = "matplotlib must be installed as an optional dependency for visualization"
        raise ImportError(msg)
