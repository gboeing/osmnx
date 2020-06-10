"""Plot spatial geometries, street networks, and routes."""

import os

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from descartes import PolygonPatch
from matplotlib.collections import LineCollection
from matplotlib.collections import PatchCollection
from shapely.geometry import MultiPolygon
from shapely.geometry import Polygon

from . import graph
from . import settings
from . import simplification
from . import utils
from . import utils_geo
from . import utils_graph


def _rgb_color_list_to_hex(color_list):
    """
    Convert a list of RGBa colors to a list of hexadecimal color codes.

    Parameters
    ----------
    color_list : list
        list of RGBa colors

    Returns
    -------
    color_list_hex : list
    """
    color_list_rgb = [[int(x * 255) for x in c[0:3]] for c in color_list]
    color_list_hex = [f"#{rgb[0]:02X}{rgb[1]:02X}{rgb[2]:02X}" for rgb in color_list_rgb]
    return color_list_hex


def get_colors(n, cmap="viridis", start=0.0, stop=1.0, alpha=1.0, return_hex=False):
    """
    Return n-length list of RGBa colors from passed colormap name and alpha.

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
        if True, convert RGBa colors to a hexadecimal string

    Returns
    -------
    colors : list
    """
    colors = [cm.get_cmap(cmap)(x) for x in np.linspace(start, stop, n)]
    colors = [(r, g, b, alpha) for r, g, b, _ in colors]
    if return_hex:
        colors = _rgb_color_list_to_hex(colors)
    return colors


def get_node_colors_by_attr(
    G, attr, num_bins=None, cmap="viridis", start=0, stop=1, na_color="none"
):
    """
    Get a list of node colors by binning continuous attribute into quantiles.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    attr : string
        name of the attribute
    num_bins : int
        how many quantiles (default None assigns each node to its own bin)
    cmap : string
        name of a matplotlib colormap
    start : float
        where to start in the colorspace
    stop : float
        where to end in the colorspace
    na_color : string
        what color to assign nodes with null attribute values

    Returns
    -------
    node_colors : list
    """
    if num_bins is None:
        num_bins = len(G)
    bin_labels = range(num_bins)
    attr_values = pd.Series([data[attr] for node, data in G.nodes(data=True)])
    cats = pd.qcut(x=attr_values, q=num_bins, labels=bin_labels)
    colors = get_colors(num_bins, cmap, start, stop)
    node_colors = [colors[int(cat)] if pd.notnull(cat) else na_color for cat in cats]
    return node_colors


def get_edge_colors_by_attr(
    G, attr, num_bins=None, cmap="viridis", start=0, stop=1, na_color="none"
):
    """
    Get a list of edge colors by binning continuous attribute into quantiles.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    attr : string
        name of the continuous-variable attribute
    num_bins : int
        how many quantiles
    cmap : string
        name of a matplotlib colormap
    start : float
        where to start in the colorspace
    stop : float
        where to end in the colorspace
    na_color : string
        what color to assign nodes with null attribute values

    Returns
    -------
    edge_colors : list
    """
    if num_bins is None:
        num_bins = len(G.edges())
    bin_labels = range(num_bins)
    attr_values = pd.Series([data[attr] for u, v, key, data in G.edges(keys=True, data=True)])
    cats = pd.qcut(x=attr_values, q=num_bins, labels=bin_labels)
    colors = get_colors(num_bins, cmap, start, stop)
    edge_colors = [colors[int(cat)] if pd.notnull(cat) else na_color for cat in cats]
    return edge_colors


def _node_list_to_coordinate_lines(G, route):
    """
    Create list of lines that follow the route defined by the list of nodes.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    route : list
        list of node IDs composing a valid path in G

    Returns
    -------
    lines : list
        list of lines as pairs ((x_start, y_start), (x_stop, y_stop))
    """
    edge_nodes = list(zip(route[:-1], route[1:]))
    lines = []
    for u, v in edge_nodes:
        # if there are parallel edges, select the shortest in length
        data = min(G.get_edge_data(u, v).values(), key=lambda x: x["length"])

        # if it has a geometry attribute (ie, a list of line segments)
        if "geometry" in data:
            # add them to the list of lines to plot
            xs, ys = data["geometry"].xy
            lines.append(list(zip(xs, ys)))
        else:
            # if it doesn't have a geometry attribute, the edge is a straight
            # line from node to node
            x1 = G.nodes[u]["x"]
            y1 = G.nodes[u]["y"]
            x2 = G.nodes[v]["x"]
            y2 = G.nodes[v]["y"]
            line = [(x1, y1), (x2, y2)]
            lines.append(line)
    return lines


def _save_and_show(fig, ax, save=False, show=True, close=True, filepath=None, dpi=300):
    """
    Save a figure to disk and/or show it, as specified by args.

    Parameters
    ----------
    fig : figure
        matplotlib figure
    ax : axes
        matplotlib axes
    save : bool
        if True, save the figure to disk at filepath
    show : bool
        if True, call pyplot.show() to show the figure
    close : bool
        if True, call pyplot.close() to close the figure
    filepath : string
        if save is True, the path to the file. file format determined from
        extension. if None, use default image folder + image.png
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
        if filepath is None:
            filepath = os.path.join(settings.imgs_folder, "image.png")

        # if save folder does not already exist, create it
        folder, _ = os.path.split(filepath)
        if not folder == "" and not os.path.exists(folder):
            os.makedirs(folder)

        fc = fig.get_facecolor()

        _, extension = os.path.splitext(filepath)
        extension = extension.strip(".")

        if extension == "svg":
            # if the file format is svg, prep the fig/ax a bit for saving
            ax.axis("off")
            ax.set_position([0, 0, 1, 1])
            ax.patch.set_alpha(0.0)
            fig.patch.set_alpha(0.0)
            fig.savefig(
                filepath, bbox_inches=0, format=extension, facecolor=fc, transparent=True,
            )
        else:
            # constrain saved figure's extent to interior of the axes
            extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            fig.savefig(
                filepath,
                dpi=dpi,
                bbox_inches=extent,
                format=extension,
                facecolor=fc,
                transparent=True,
            )
        utils.log(f"Saved figure to disk at {filepath}")

    if show:
        plt.show()

    if close:
        plt.close()

    return fig, ax


def _config_ax(ax, crs, bbox, padding=0):
    """
    Configure axes for display.

    Parameters
    ----------
    ax : matplotlib axes
        the axes containing the plot
    crs : dict or string or pyproj.CRS
        the CRS of the plotted geometries
    bbox : tuple
        bounding box as (north, south, east, west)
    padding : float
        relative padding to add around the plot's bbox

    Returns
    -------
    ax : matplotlib axes
        the configured/styled axes
    """
    # set the axes view limits
    north, south, east, west = bbox
    padding_ns = (north - south) * padding
    padding_ew = (east - west) * padding
    ax.set_ylim((south - padding_ns, north + padding_ns))
    ax.set_xlim((west - padding_ew, east + padding_ew))

    # turn off the axes display, set margins to zero, and point ticks inward
    # so there's no space around the plot
    ax.axis("off")
    ax.margins(0)
    ax.tick_params(which="both", direction="in")

    xaxis = ax.get_xaxis()
    yaxis = ax.get_yaxis()
    xaxis.get_major_formatter().set_useOffset(False)
    yaxis.get_major_formatter().set_useOffset(False)
    xaxis.set_visible(False)
    yaxis.set_visible(False)

    # set aspect ratio
    if crs == settings.default_crs:
        # if data are not projected, conform aspect ratio to not stretch plot
        coslat = np.cos((south + north) / 2.0 / 180.0 * np.pi)
        ax.set_aspect(1.0 / coslat)
    else:
        # if projected, make everything square
        ax.set_aspect("equal")

    return ax


def plot_graph(
    G,
    ax=None,
    figsize=(8, 8),
    bgcolor="#222222",
    node_color="w",
    node_size=15,
    node_alpha=1,
    node_edgecolor="none",
    node_zorder=1,
    edge_color="#999999",
    edge_linewidth=1,
    edge_alpha=1,
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
    ax : matplotlib axes
        if not None, plot on this preexisting axes
    figsize : tuple
        if ax is None, sets figure size as (width, height). if None, default
        to 6 inch height with proportional width
    bgcolor : string
        background color of plot
    node_color : string
        color of the nodes
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
    edge_color : string
        color of the edges' lines
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
        extension. if None, use default image folder + image.png
    dpi : int
        if save is True, the resolution of saved file
    bbox : tuple
        bounding box as (north, south, east, west). if None, will calculate
        from spatial extents of edges.

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axes
    """
    if node_size <= 0 and edge_linewidth <= 0:
        raise ValueError("Either node_size or edge_linewidth must be > 0 to plot something.")

    # create fig, ax as needed
    utils.log("Begin plotting the graph...")
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, facecolor=bgcolor)
        ax.set_facecolor(bgcolor)
    else:
        fig = ax.figure

    if edge_linewidth > 0:
        # plot the edges' geometries
        gdf_edges = utils_graph.graph_to_gdfs(G, nodes=False)
        ax = gdf_edges.plot(ax=ax, color=edge_color, lw=edge_linewidth, alpha=edge_alpha, zorder=1)

    if node_size > 0:
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
        padding = 0.02  # pad 2% to not trim nodes' circles

    # configure axes appearance, save/show figure as specified, and return
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
    ax : matplotlib axes
        if not None, plot route on this preexisting axes instead of creating a
        new fig, ax and drawing the underlying graph
    pg_kwargs
        keyword arguments to pass to plot_graph

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axes
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

    # add the routes to the axes as a LineCollection
    lines = _node_list_to_coordinate_lines(G, route)
    lc = LineCollection(lines, colors=route_color, linewidths=route_linewidth, alpha=route_alpha)
    ax.add_collection(lc)

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
    route_colors : list
        if string, the color for all routes; if list, the color for each route
    pgr_kwargs
        keyword arguments to pass to plot_graph_route

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axes
    """
    # check for valid arguments
    if not all([isinstance(r, list) for r in routes]):
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
        how many meters to extend north, south, east, and west from the center
        point
    network_type : string
        what type of network to get
    street_widths : dict
        dict keys are street types and values are widths to plot in pixels
    default_width : numeric
        fallback street width in pixels for any street type not found in
        street_widths dict
    figsize : numeric
        (width, height) of figure, should be equal
    edge_color : string
        color of the edges' lines
    smooth_joints : bool
        if True, plot nodes same width as streets to smooth line joints and
        prevent cracks between them from showing

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axes
    """
    multiplier = 1.2

    # if G was passed-in, use this graph in the plot, centered on the centroid
    # of its nodes
    if G is not None:
        gdf_nodes = utils_graph.graph_to_gdfs(G, edges=False, node_geometry=True)
        lnglat_point = gdf_nodes.unary_union.centroid.coords[0]
        point = tuple(reversed(lnglat_point))

    # otherwise, get the network by either address or point, whichever was
    # passed-in, using a distance multiplier to make sure we get more than
    # enough network. simplify in non-strict mode to not combine multiple street
    # types into single edge
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

    # if user did not pass in custom street widths, create a dict of default
    # values
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

    # we need an undirected graph to find every edge incident to a node
    G_undir = utils_graph.get_undirected(G)

    # for each network edge, get a linewidth according to street type (the OSM
    # 'highway' value)
    edge_linewidths = []
    for _, _, data in G_undir.edges(keys=False, data=True):
        street_type = data["highway"][0] if isinstance(data["highway"], list) else data["highway"]
        if street_type in street_widths:
            edge_linewidths.append(street_widths[street_type])
        else:
            edge_linewidths.append(default_width)

    if smooth_joints:
        # for each node, get a nodesize according to the narrowest incident edge
        node_widths = {}
        for node in G_undir.nodes():
            # first, identify all the highway types of this node's incident edges
            incident_edges_data = [
                G_undir.get_edge_data(node, neighbor) for neighbor in G_undir.neighbors(node)
            ]
            edge_types = [data[min(data)]["highway"] for data in incident_edges_data]
            if len(edge_types) < 1:
                # if node has no incident edges, make size zero
                node_widths[node] = 0
            else:
                # flatten the list of edge types
                edge_types_flat = []
                for et in edge_types:
                    if isinstance(et, list):
                        edge_types_flat.extend(et)
                    else:
                        edge_types_flat.append(et)

                # for each edge type in the flattened list, lookup the
                # corresponding width
                edge_widths = [
                    street_widths[edge_type] if edge_type in street_widths else default_width
                    for edge_type in edge_types_flat
                ]

                # the node diameter will be the biggest of the edge widths, to
                # make joints perfectly smooth alternatively, use min (?) to
                # prevent anything larger from extending past smallest
                # street's line
                circle_diameter = max(edge_widths)

                # mpl circle marker sizes are in area, so it is the diameter
                # squared
                circle_area = circle_diameter ** 2
                node_widths[node] = circle_area

        # assign the node size to each node in the graph
        node_sizes = [node_widths[node] for node in G_undir.nodes()]
    else:
        node_sizes = 0

    # define the spatial extents of the plotting figure to make it square, in
    # projected units, and cropped to the desired area
    bbox = utils_geo.bbox_from_point(point, dist, project_utm=False)

    # plot the figure
    override = {"bbox", "node_size", "node_color", "edge_linewidth"}
    kwargs = {k: v for k, v in pg_kwargs.items() if k not in override}
    fig, ax = plot_graph(
        G=G_undir,
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
    figsize=None,
    color="orange",
    bgcolor="#222222",
    bbox=None,
    save=False,
    show=True,
    close=False,
    filepath=None,
    dpi=600,
):
    """
    Plot a GeoDataFrame of footprints.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame of footprints (shapely Polygons and MultiPolygons)
    ax : axes
        if not None, plot on this preexisting axes
    figsize : tuple
        if ax is None, use figsize (width, height) to determine figure size
    color : string
        color of the footprints
    bgcolor : string
        background color of the plot
    set_bounds : bool
        if True, set bounds from either passed-in bbox or the spatial extent
        of gdf
    bbox : tuple
        if True, set the figure's extent to this bounding box, otherwise use
        the spatial extents of the geometries in gdf
    save : bool
        if True, save the figure to disk at filepath
    show : bool
        if True, call pyplot.show() to show the figure
    close : bool
        if True, call pyplot.close() to close the figure
    filepath : string
        if save is True, the path to the file. file format determined from
        extension. if None, use default image folder + image.png
    dpi : int
        if save is True, the resolution of saved file

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, facecolor=bgcolor)
        ax.set_facecolor(bgcolor)
    else:
        fig = ax.figure

    # extract each polygon as a descartes patch, and add to a matplotlib patch
    # collection
    patches = []
    for geometry in gdf["geometry"]:
        if isinstance(geometry, Polygon):
            patches.append(PolygonPatch(geometry))
        elif isinstance(geometry, MultiPolygon):
            for (
                subpolygon
            ) in geometry:  # if geometry is multipolygon, go through each constituent subpolygon
                patches.append(PolygonPatch(subpolygon))
    pc = PatchCollection(patches, facecolor=color, edgecolor=color, linewidth=0, alpha=1)
    ax.add_collection(pc)

    # determine figure extents
    if bbox is None:
        # set the figure bounds to the polygons' bounds
        west, south, east, north = gdf.total_bounds
    else:
        north, south, east, west = bbox

    # configure axes appearance, save/show figure as specified, and return
    ax = _config_ax(ax, gdf.crs, (north, south, east, west))
    fig, ax = _save_and_show(fig, ax, save, show, close, filepath, dpi)
    return fig, ax
