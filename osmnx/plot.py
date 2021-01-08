"""Plot spatial geometries, street networks, and routes."""

from pathlib import Path

import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from . import graph
from . import settings
from . import simplification
from . import utils
from . import utils_geo
from . import utils_graph


def get_colors(n, cmap="viridis", start=0.0, stop=1.0, alpha=1.0, return_hex=False):
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
    color_list = [cm.get_cmap(cmap)(x) for x in np.linspace(start, stop, n)]
    if return_hex:
        color_list = [colors.to_hex(c) for c in color_list]
    else:
        color_list = [(r, g, b, alpha) for r, g, b, _ in color_list]
    return color_list


def get_node_colors_by_attr(
    G, attr, num_bins=None, cmap="viridis", start=0, stop=1, na_color="none", equal_size=False
):
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
    vals = pd.Series(nx.get_node_attributes(G, attr))
    return _get_colors_by_value(vals, num_bins, cmap, start, stop, na_color, equal_size)


def get_edge_colors_by_attr(
    G, attr, num_bins=None, cmap="viridis", start=0, stop=1, na_color="none", equal_size=False
):
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
    vals = pd.Series(nx.get_edge_attributes(G, attr))
    return _get_colors_by_value(vals, num_bins, cmap, start, stop, na_color, equal_size)


def plot_graph(
    G,
    ax=None,
    figsize=(8, 8),
    bgcolor="#111111",
    node_color="w",
    node_size=15,
    node_alpha=None,
    node_edgecolor="none",
    node_zorder=1,
    edge_color="#999999",
    edge_linewidth=1,
    edge_alpha=None,
    show=True,
    close=False,
    save=False,
    filepath=None,
    dpi=300,
    bbox=None,
):
    """
    Plot a graph.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    ax : matplotlib axis
        if not None, plot on this preexisting axis
    figsize : tuple
        if ax is None, create new figure with size (width, height)
    bgcolor : string
        background color of plot
    node_color : string or list
        color(s) of the nodes
    node_size : int
        size of the nodes: if 0, then skip plotting the nodes
    node_alpha : float
        opacity of the nodes, note: if you passed RGBA values to node_color,
        set node_alpha=None to use the alpha channel in node_color
    node_edgecolor : string
        color of the nodes' markers' borders
    node_zorder : int
        zorder to plot nodes: edges are always 1, so set node_zorder=0 to plot
        nodes below edges
    edge_color : string or list
        color(s) of the edges' lines
    edge_linewidth : float
        width of the edges' lines: if 0, then skip plotting the edges
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
        matplotlib figure, axis
    """
    max_node_size = max(node_size) if hasattr(node_size, "__iter__") else node_size
    max_edge_lw = max(edge_linewidth) if hasattr(edge_linewidth, "__iter__") else edge_linewidth
    if max_node_size <= 0 and max_edge_lw <= 0:
        raise ValueError("Either node_size or edge_linewidth must be > 0 to plot something.")

    # create fig, ax as needed
    utils.log("Begin plotting the graph...")
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, facecolor=bgcolor, frameon=False)
        ax.set_facecolor(bgcolor)
    else:
        fig = ax.figure

    if max_edge_lw > 0:
        # plot the edges' geometries
        gdf_edges = utils_graph.graph_to_gdfs(G, nodes=False)["geometry"]
        ax = gdf_edges.plot(ax=ax, color=edge_color, lw=edge_linewidth, alpha=edge_alpha, zorder=1)

    if max_node_size > 0:
        # scatter plot the nodes' x/y coordinates
        gdf_nodes = utils_graph.graph_to_gdfs(G, edges=False, node_geometry=False)[["x", "y"]]
        ax.scatter(
            x=gdf_nodes["x"],
            y=gdf_nodes["y"],
            s=node_size,
            c=node_color,
            alpha=node_alpha,
            edgecolor=node_edgecolor,
            zorder=node_zorder,
        )

    # get spatial extents from bbox parameter or the edges' geometries
    padding = 0
    if bbox is None:
        try:
            west, south, east, north = gdf_edges.total_bounds
        except NameError:
            west, south = gdf_nodes.min()
            east, north = gdf_nodes.max()
        bbox = north, south, east, west
        padding = 0.02  # pad 2% to not cut off peripheral nodes' circles

    # configure axis appearance, save/show figure as specified, and return
    ax = _config_ax(ax, G.graph["crs"], bbox, padding)
    fig, ax = _save_and_show(fig, ax, save, show, close, filepath, dpi)
    utils.log("Finished plotting the graph")
    return fig, ax


def plot_graph_route(
    G,
    route,
    route_color="r",
    route_linewidth=4,
    route_alpha=0.5,
    orig_dest_size=100,
    ax=None,
    **pg_kwargs,
):
    """
    Plot a route along a graph.

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
    ax : matplotlib axis
        if not None, plot route on this preexisting axis instead of creating a
        new fig, ax and drawing the underlying graph
    pg_kwargs
        keyword arguments to pass to plot_graph

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axis
    """
    if ax is None:
        # plot the graph but not the route, and override any user show/close
        # args for now: we'll do that later
        override = {"show", "save", "close"}
        kwargs = {k: v for k, v in pg_kwargs.items() if k not in override}
        fig, ax = plot_graph(G, show=False, save=False, close=False, **kwargs)
    else:
        fig = ax.figure

    # scatterplot origin and destination points (first/last nodes in route)
    x = (G.nodes[route[0]]["x"], G.nodes[route[-1]]["x"])
    y = (G.nodes[route[0]]["y"], G.nodes[route[-1]]["y"])
    ax.scatter(x, y, s=orig_dest_size, c=route_color, alpha=route_alpha, edgecolor="none")

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


def plot_graph_routes(G, routes, route_colors="r", **pgr_kwargs):
    """
    Plot several routes along a graph.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    routes : list
        routes as a list of lists of node IDs
    route_colors : string or list
        if string, 1 color for all routes. if list, the colors for each route.
    pgr_kwargs
        keyword arguments to pass to plot_graph_route

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axis
    """
    # check for valid arguments
    if not all(isinstance(r, list) for r in routes):
        raise ValueError("routes must be a list of route lists")
    if len(routes) < 2:
        raise ValueError("You must pass more than 1 route")
    if isinstance(route_colors, str):
        route_colors = [route_colors] * len(routes)
    if len(routes) != len(route_colors):
        raise ValueError("route_colors list must have same length as routes")

    # plot the graph and the first route
    override = {"route", "route_color", "show", "save", "close"}
    kwargs = {k: v for k, v in pgr_kwargs.items() if k not in override}
    fig, ax = plot_graph_route(
        G,
        route=routes[0],
        route_color=route_colors[0],
        show=False,
        save=False,
        close=False,
        **kwargs,
    )

    # plot the subsequent routes on top of existing ax
    override.update({"ax"})
    kwargs = {k: v for k, v in pgr_kwargs.items() if k not in override}
    for route, route_color in zip(routes[1:], route_colors[1:]):
        fig, ax = plot_graph_route(
            G,
            route=route,
            route_color=route_color,
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
    G=None,
    address=None,
    point=None,
    dist=805,
    network_type="drive_service",
    street_widths=None,
    default_width=4,
    figsize=(8, 8),
    edge_color="w",
    smooth_joints=True,
    **pg_kwargs,
):
    """
    Plot a figure-ground diagram of a street network.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph, must be unprojected
    address : string
        address to geocode as the center point if G is not passed in
    point : tuple
        center point if address and G are not passed in
    dist : numeric
        how many meters to extend north, south, east, west from center point
    network_type : string
        what type of street network to get
    street_widths : dict
        dict keys are street types and values are widths to plot in pixels
    default_width : numeric
        fallback width in pixels for any street type not in street_widths
    figsize : numeric
        (width, height) of figure, should be equal
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
        matplotlib figure, axis
    """
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
        lnglat_point = gdf_nodes.unary_union.centroid.coords[0]
        point = tuple(reversed(lnglat_point))

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
    else:
        raise ValueError("You must pass an address or lat-lng point or graph.")

    # we need an undirected graph to find every edge incident to a node
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
        node_widths = dict()
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
                edge_widths = [
                    street_widths[et] if et in street_widths else default_width for et in et_flat
                ]

                # node diameter should equal largest edge width to make joints
                # perfectly smooth. alternatively use min(?) to prevent
                # anything larger from extending past smallest street's line.
                # circle marker sizes are in area, so use diameter squared.
                circle_diameter = max(edge_widths)
                circle_area = circle_diameter ** 2
                node_widths[node] = circle_area

        # assign the node size to each node in the graph
        node_sizes = [node_widths[node] for node in Gu.nodes]
    else:
        node_sizes = 0

    # define the view extents of the plotting figure
    bbox = utils_geo.bbox_from_point(point, dist, project_utm=False)

    # plot the figure
    override = {"bbox", "node_size", "node_color", "edge_linewidth"}
    kwargs = {k: v for k, v in pg_kwargs.items() if k not in override}
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
    gdf,
    ax=None,
    figsize=(8, 8),
    color="orange",
    alpha=None,
    bgcolor="#111111",
    bbox=None,
    save=False,
    show=True,
    close=False,
    filepath=None,
    dpi=600,
):
    """
    Plot a GeoDataFrame of geospatial entities' footprints.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame of footprints (shapely Polygons and MultiPolygons)
    ax : axis
        if not None, plot on this preexisting axis
    figsize : tuple
        if ax is None, create new figure with size (width, height)
    color : string
        color of the footprints
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
        matplotlib figure, axis
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, facecolor=bgcolor, frameon=False)
        ax.set_facecolor(bgcolor)
    else:
        fig = ax.figure

    # retain only Polygons and MultiPolygons, then plot
    gdf = gdf[gdf["geometry"].type.isin({"Polygon", "MultiPolygon"})]
    ax = gdf.plot(ax=ax, facecolor=color, edgecolor="none", linewidth=0, alpha=alpha)

    # determine figure extents
    if bbox is None:
        west, south, east, north = gdf.total_bounds
    else:
        north, south, east, west = bbox

    # configure axis appearance, save/show figure as specified, and return
    ax = _config_ax(ax, gdf.crs, (north, south, east, west), 0)
    fig, ax = _save_and_show(fig, ax, save, show, close, filepath, dpi)
    return fig, ax


def _get_colors_by_value(vals, num_bins, cmap, start, stop, na_color, equal_size):
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
        raise ValueError("There are no attribute values.")

    if num_bins is None:
        # calculate min/max values based on start/stop and data range
        vals_min = vals.dropna().min()
        vals_max = vals.dropna().max()
        full_range = (vals_max - vals_min) / (stop - start)
        full_min = vals_min - full_range * start
        full_max = full_min + full_range

        # linearly map a color to each attribute value
        normalizer = colors.Normalize(full_min, full_max)
        scalar_mapper = cm.ScalarMappable(normalizer, cm.get_cmap(cmap))
        color_series = vals.map(scalar_mapper.to_rgba)
        color_series.loc[pd.isnull(vals)] = na_color

    else:
        # otherwise, bin values then assign colors to bins
        cut_func = pd.qcut if equal_size else pd.cut
        bins = cut_func(vals, num_bins, labels=range(num_bins))
        bin_colors = get_colors(num_bins, cmap, start, stop)
        color_list = [bin_colors[b] if pd.notnull(b) else na_color for b in bins]
        color_series = pd.Series(color_list, index=bins.index)

    return color_series


def _save_and_show(fig, ax, save=False, show=True, close=True, filepath=None, dpi=300):
    """
    Save a figure to disk and/or show it, as specified by args.

    Parameters
    ----------
    fig : figure
        matplotlib figure
    ax : axis
        matplotlib axis
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
        matplotlib figure, axis
    """
    fig.canvas.draw()
    fig.canvas.flush_events()

    if save:

        # default filepath, if none provided
        if filepath is None:
            filepath = Path(settings.imgs_folder) / "image.png"
        else:
            filepath = Path(filepath)

        # if save folder does not already exist, create it
        filepath.parent.mkdir(parents=True, exist_ok=True)

        # get the file extension and figure facecolor
        ext = filepath.suffix.strip(".")
        fc = fig.get_facecolor()

        if ext == "svg":
            # if the file format is svg, prep the fig/ax for saving
            ax.axis("off")
            ax.set_position([0, 0, 1, 1])
            ax.patch.set_alpha(0.0)
            fig.patch.set_alpha(0.0)
            fig.savefig(filepath, bbox_inches=0, format=ext, facecolor=fc, transparent=True)
        else:
            # constrain saved figure's extent to interior of the axis
            extent = ax.bbox.transformed(fig.dpi_scale_trans.inverted())

            # temporarily turn figure frame on to save with facecolor
            fig.set_frameon(True)
            fig.savefig(
                filepath, dpi=dpi, bbox_inches=extent, format=ext, facecolor=fc, transparent=True
            )
            fig.set_frameon(False)  # and turn it back off again
        utils.log(f"Saved figure to disk at {filepath}")

    if show:
        plt.show()

    if close:
        plt.close()

    return fig, ax


def _config_ax(ax, crs, bbox, padding):
    """
    Configure axis for display.

    Parameters
    ----------
    ax : matplotlib axis
        the axis containing the plot
    crs : dict or string or pyproj.CRS
        the CRS of the plotted geometries
    bbox : tuple
        bounding box as (north, south, east, west)
    padding : float
        relative padding to add around the plot's bbox

    Returns
    -------
    ax : matplotlib axis
        the configured/styled axis
    """
    # set the axis view limits to bbox + relative padding
    north, south, east, west = bbox
    padding_ns = (north - south) * padding
    padding_ew = (east - west) * padding
    ax.set_ylim((south - padding_ns, north + padding_ns))
    ax.set_xlim((west - padding_ew, east + padding_ew))

    # set margins to zero, point ticks inward, turn off ax border and x/y axis
    # so there is no space around the plot
    ax.margins(0)
    ax.tick_params(which="both", direction="in")
    _ = [s.set_visible(False) for s in ax.spines.values()]
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    # set aspect ratio
    if crs == settings.default_crs:
        # if data are not projected, conform aspect ratio to not stretch plot
        coslat = np.cos((south + north) / 2 / 180 * np.pi)
        ax.set_aspect(1 / coslat)
    else:
        # if projected, make everything square
        ax.set_aspect("equal")

    return ax
