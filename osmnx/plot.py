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


def plot_shape(
    gdf,
    fc="#cbe0f0",
    ec="#999999",
    linewidth=1,
    alpha=1,
    figsize=(6, 6),
    margin=0.02,
    axis_off=True,
):
    """
    Plot a GeoDataFrame of place boundary geometries.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        the gdf containing the geometries to plot
    fc : string or list
        the facecolor (or list of facecolors) for the polygons
    ec : string or list
        the edgecolor (or list of edgecolors) for the polygons
    linewidth : numeric
        the width of the polygon edge lines
    alpha : numeric
        the opacity
    figsize : tuple
        width, height of plotting figure
    margin : numeric
        the size of the figure margins
    axis_off : bool
        if True, disable the matplotlib axes display

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axes
    """
    # if facecolor or edgecolor is a string instead of a list, make sure we have
    # as many colors as gdf elements
    if isinstance(fc, str):
        fc = [fc] * len(gdf)
    if isinstance(ec, str):
        ec = [ec] * len(gdf)

    # plot the geometries one at a time
    fig, ax = plt.subplots(figsize=figsize)
    for geometry, facecolor, edgecolor in zip(gdf["geometry"], fc, ec):
        if isinstance(geometry, (Polygon, MultiPolygon)):
            if isinstance(geometry, Polygon):
                geometry = MultiPolygon([geometry])
            for polygon in geometry:
                patch = PolygonPatch(
                    polygon, fc=facecolor, ec=edgecolor, linewidth=linewidth, alpha=alpha
                )
                ax.add_patch(patch)
        else:
            raise ValueError(
                "All geometries in GeoDataFrame must be shapely Polygons or MultiPolygons"
            )

    # adjust the axes margins and limits around the image and make axes
    # equal-aspect
    west, south, east, north = gdf.unary_union.bounds
    margin_ns = (north - south) * margin
    margin_ew = (east - west) * margin
    ax.set_ylim((south - margin_ns, north + margin_ns))
    ax.set_xlim((west - margin_ew, east + margin_ew))
    ax.set_aspect(aspect="equal", adjustable="box")
    if axis_off:
        ax.axis("off")

    plt.show()
    return fig, ax


def _rgb_color_list_to_hex(color_list):
    """
    Convert a list of RGBa colors to a list of hexadecimal color codes.

    Parameters
    ----------
    color_list : list
        the list of RGBa colors

    Returns
    -------
    color_list_hex : list
    """
    color_list_rgb = [[int(x * 255) for x in c[0:3]] for c in color_list]
    color_list_hex = [f"#{rgb[0]:02X}{rgb[1]:02X}{rgb[2]:02X}" for rgb in color_list_rgb]
    return color_list_hex


def get_colors(n, cmap="viridis", start=0.0, stop=1.0, alpha=1.0, return_hex=False):
    """
    Return n-length list of RGBa colors from the passed colormap name and alpha.

    Parameters
    ----------
    n : int
        number of colors
    cmap : string
        name of a colormap
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
        the name of the attribute
    num_bins : int
        how many quantiles (default None assigns each node to its own bin)
    cmap : string
        name of a colormap
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
        the name of the continuous-variable attribute
    num_bins : int
        how many quantiles
    cmap : string
        name of a colormap
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


def _node_list_to_coordinate_lines(G, route, use_geom=True):
    """
    Make list of lines that follow the route defined by the list of nodes.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    route : list
        a valid graph path as a list of node IDs
    use_geom : bool
        if True, use the spatial geometry attribute of the edges to draw
        geographically accurate edges, rather than just lines straight from
        node to node

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
        if "geometry" in data and use_geom:
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


def _save_and_show(
    fig, ax, save=False, show=True, close=True, filepath=None, dpi=300, axis_off=True,
):
    """
    Save a figure to disk and show it, as specified.

    Parameters
    ----------
    fig : figure
        matplotlib figure
    ax : axes
        matplotlib axes
    save : bool
        if True, save the figure as an image file to disk
    show : bool
        if True, call plt.show() to show the figure
    close : bool
        if True, call plt.close() to close the figure
    filepath : string
        if save is True, the path to the file. file format determined from
        extension. if None, use default image folder + image.png
    dpi : int
        if save is True, the resolution of the image file
    axis_off : bool
        if True, matplotlib axes was turned off by plot_graph so constrain the
        saved figure's extent to the interior of the axes

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axes
    """
    fig.canvas.draw()
    fig.canvas.flush_events()

    if save:

        # default filepath, if none was provided
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
            if axis_off:
                # if axis is turned off, constrain the saved figure's extent to
                # the interior of the axes
                extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            else:
                extent = "tight"
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


def _config_ax(ax, axis_off, equal_aspect, top, bottom, crs):

    xaxis = ax.get_xaxis()
    yaxis = ax.get_yaxis()

    xaxis.get_major_formatter().set_useOffset(False)
    yaxis.get_major_formatter().set_useOffset(False)

    # if axis_off, turn off the axes display set the margins to zero and point
    # the ticks in so there's no space around the plot
    if axis_off:
        ax.axis("off")
        ax.margins(0)
        ax.tick_params(which="both", direction="in")
        xaxis.set_visible(False)
        yaxis.set_visible(False)

    if equal_aspect:
        # make everything square
        ax.set_aspect("equal")

    elif crs == settings.default_crs:
        # if data are not projected, conform aspect ratio to not stretch plot
        coslat = np.cos((bottom + top) / 2.0 / 180.0 * np.pi)
        ax.set_aspect(1.0 / coslat)

    return ax


def plot_graph(
    G,
    ax=None,
    figsize=None,
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
    margin=0.02,
    equal_aspect=False,
    axis_off=True,
    annotate=False,
    draw_graph=True,
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
        if ax is None, use figsize (width, height) to determine figure size;
        if None, default to 6 inch height with proportional width
    bgcolor : string
        the background color of the plot
    node_color : string
        the color of the nodes. color is passed to matplotlib
    node_size : int
        the size of the nodes
    node_alpha : float
        the opacity of the nodes. if you passed RGBA values to node_color,
        then set this to None to use the alpha channel in node_color. if 0,
        then skip plotting the nodes.
    node_edgecolor : string
        the color of the node's marker's border
    node_zorder : int
        zorder to plot nodes, edges are always 1, so make node_zorder 0 to
        plot nodes beneath them or 2 to plot nodes atop them
    edge_color : string
        the color of the edges' lines. color is passed to matplotlib.
    edge_linewidth : float
        the width of the edges' lines
    edge_alpha : float
        the opacity of the edges' lines. if you passed RGBA values to
        edge_color, then set this to None to use the alpha channel in
        edge_color. if 0, then skip plotting the edges.
    show : bool
        if True, call plt.show() to show the figure
    close : bool
        if True, call plt.close() to close the figure
    save : bool
        if True, save the figure as an image file to disk
    filepath : string
        if save is True, the path to the file. file format determined from
        extension. if None, use default image folder + image.png
    dpi : int
        if save is True, the resolution of the image file
    bbox : tuple
        bounding box as (north, south, east, west); if None will calculate
        from spatial extents of edges; if passing bbox, you probably also want
        to pass margin=0 to constrain it
    margin : float
        relative margin around the figure
    equal_aspect : bool
        if True, set the axes aspect ratio equal
    axis_off : bool
        if True, turn off the matplotlib axes
    annotate : bool
        if True, annotate the nodes with their IDs
    draw_graph : bool
        if False, create and style fig and ax, but do not draw graph on them

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axes
    """
    utils.log("Begin plotting the graph...")
    node_Xs = [float(x) for _, x in G.nodes(data="x")]
    node_Ys = [float(y) for _, y in G.nodes(data="y")]

    # get north, south, east, west values either from bbox parameter or from the
    # spatial extent of the edges' geometries
    if bbox is None:
        edges = utils_graph.graph_to_gdfs(G, nodes=False, fill_edge_geometry=True)
        west, south, east, north = edges.total_bounds
    else:
        north, south, east, west = bbox

    if ax is not None:
        fig = ax.figure
    else:
        if figsize is None:
            # if figsize not passed in, make height=6in and proportional width
            bbox_aspect_ratio = (north - south) / (east - west)
            fig_height = 6
            fig_width = fig_height / bbox_aspect_ratio
            figsize = (fig_width, fig_height)

        # create the figure and axes
        fig, ax = plt.subplots(figsize=figsize, facecolor=bgcolor)
        ax.set_facecolor(bgcolor)

    if draw_graph:

        if edge_alpha > 0:
            # draw the edges as lines from node to node
            lines = []
            for u, v, data in G.edges(keys=False, data=True):
                if "geometry" in data:
                    # if edge has geometry attr, add it to list of lines to plot
                    xs, ys = data["geometry"].xy
                    lines.append(list(zip(xs, ys)))
                else:
                    # if it doesn't have geometry attr, draw straight line from
                    # node to node
                    x1 = G.nodes[u]["x"]
                    y1 = G.nodes[u]["y"]
                    x2 = G.nodes[v]["x"]
                    y2 = G.nodes[v]["y"]
                    line = [(x1, y1), (x2, y2)]
                    lines.append(line)

            # add the lines to the axes as a LineCollection
            lc = LineCollection(
                lines, colors=edge_color, linewidths=edge_linewidth, alpha=edge_alpha, zorder=1
            )
            ax.add_collection(lc)
            utils.log("Drew the graph edges")

        if node_alpha > 0:
            # scatter plot the nodes
            ax.scatter(
                node_Xs,
                node_Ys,
                s=node_size,
                c=node_color,
                alpha=node_alpha,
                edgecolor=node_edgecolor,
                zorder=node_zorder,
            )

    # set the extent of the figure
    margin_ns = (north - south) * margin
    margin_ew = (east - west) * margin
    ax.set_ylim((south - margin_ns, north + margin_ns))
    ax.set_xlim((west - margin_ew, east + margin_ew))

    # configure axes appearance
    ax = _config_ax(ax, axis_off, equal_aspect, north, south, crs=G.graph["crs"])

    # annotate the nodes with their IDs if annotate=True
    if annotate:
        for node, data in G.nodes(data=True):
            ax.annotate(node, xy=(data["x"], data["y"]))

    # save and show the figure as specified
    fig, ax = _save_and_show(fig, ax, save, show, close, filepath, dpi, axis_off)
    return fig, ax


def plot_graph_route(
    G, route, route_color="r", route_linewidth=4, route_alpha=0.5, orig_dest_size=100, **pg_kwargs,
):
    """
    Plot a route along a graph.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    route : list
        the route as a list of nodes
    route_color : string
        the color of the route
    route_linewidth : int
        the width of the route line
    route_alpha : float
        the opacity of the route line
    orig_dest_size : int
        the size of the origin and destination nodes
    **pg_kwargs
        keyword arguments to pass to plot_graph

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axes
    """
    # plot the graph but not the route, and override any user show/close
    # args for now: we'll do that later
    override = {"G", "show", "save", "close"}
    kwargs = {k: v for k, v in pg_kwargs.items() if k not in override}
    fig, ax = plot_graph(G, show=False, save=False, close=False, **kwargs)

    # scatterplot origin and destination points (first/last nodes in route)
    orig_dest_Xs = (G.nodes[route[0]]["x"], G.nodes[route[-1]]["x"])
    orig_dest_Ys = (G.nodes[route[0]]["y"], G.nodes[route[-1]]["y"])
    ax.scatter(
        orig_dest_Xs,
        orig_dest_Ys,
        s=orig_dest_size,
        c=route_color,
        alpha=route_alpha,
        edgecolor="none",
        zorder=3,
    )

    # add the routes to the axes as a LineCollection
    lines = _node_list_to_coordinate_lines(G, route)
    lc = LineCollection(
        lines, colors=route_color, linewidths=route_linewidth, alpha=route_alpha, zorder=2
    )
    ax.add_collection(lc)

    # save and show the figure as specified, passing relevant kwargs
    sas_kwargs = {"save", "show", "close", "filepath", "file_format", "dpi", "axis_off"}
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
        the routes as a list of lists of nodes
    route_colors : list
        if string, the color for all routes; if list, the color for each route
    pgr_kwargs
        keyword arguments to pass to plot_graph_route

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axes
    """
    # initial checks
    if not all([isinstance(r, list) for r in routes]):
        raise ValueError("routes must be a list of route lists")
    if len(routes) < 2:
        raise ValueError("You must pass more than 1 route")
    if isinstance(route_colors, str):
        route_colors = [route_colors] * len(routes)
    if len(routes) != len(route_colors):
        raise ValueError("route_colors list must have same length as routes")

    # plot the first route on top of the graph
    override = {"G", "route", "route_color", "show", "save", "close"}
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

    # plot the subsequent routes on top of ax, without replotting the graph
    override.update({"fig", "ax", "draw_graph"})
    kwargs = {k: v for k, v in pgr_kwargs.items() if k not in override}
    for route, route_color in zip(routes[1:], route_colors[1:]):
        fig, ax = plot_graph_route(
            G,
            route=route,
            route_color=route_color,
            ax=ax,
            draw_graph=False,
            show=False,
            save=False,
            close=False,
            **kwargs,
        )

    # save and show the figure as specified, passing relevant kwargs
    sas_kwargs = {"save", "show", "close", "filepath", "file_format", "dpi", "axis_off"}
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
    Plot figure-ground diagram of a street network.

    Defaults to one square mile.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph, must be unprojected
    address : string
        the address to geocode as the center point if G is not passed in
    point : tuple
        the center point if address and G are not passed in
    dist : numeric
        how many meters to extend north, south, east, and west from the center
        point
    network_type : string
        what type of network to get
    street_widths : dict
        where keys are street types and values are widths to plot in pixels
    default_width : numeric
        the default street width in pixels for any street type not found in
        street_widths dict
    figsize : numeric
        width, height of figure (should be equal)
    edge_color : string
        the color of the edges' lines
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
    G_undir = G.to_undirected()

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
            edge_types = [data[0]["highway"] for data in incident_edges_data]
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
    fig, ax = plot_graph(
        G_undir,
        bbox=bbox,
        figsize=figsize,
        margin=0,
        axis_off=True,
        equal_aspect=False,
        node_size=node_sizes,
        node_color=edge_color,
        edge_color=edge_color,
        edge_linewidth=edge_linewidths,
        **pg_kwargs,
    )

    return fig, ax


def plot_footprints(
    gdf,
    ax=None,
    figsize=None,
    color="orange",
    bgcolor="#222222",
    set_bounds=True,
    bbox=None,
    save=False,
    show=True,
    close=False,
    filepath=None,
    dpi=600,
    equal_aspect=False,
    axis_off=True,
):
    """
    Plot a GeoDataFrame of footprints.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame of footprints
    ax : axes
        matplotlib axes
    figsize : tuple
        if ax is None, use figsize (width, height) to determine figure size
    color : string
        the color of the footprints
    bgcolor : string
        the background color of the plot
    set_bounds : bool
        if True, set bounds from either passed-in bbox or the spatial extent
        of gdf
    bbox : tuple
        if True and if set_bounds is True, set the display bounds to this bbox
    save : bool
        if True, save the figure as an image file to disk
    show : bool
        if True, call plt.show() to show the figure
    close : bool
        if True, call plt.close() to close the figure
    filepath : string
        if save is True, the path to the file. file format determined from
        extension. if None, use default image folder + image.png
    dpi : int
        if save is True, the resolution of the image file
    equal_aspect : bool
        if True, set the axes aspect ratio equal
    axis_off : bool
        if True, turn off the matplotlib axes

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

    if set_bounds:
        if bbox is None:
            # set the figure bounds to the polygons' bounds
            left, bottom, right, top = gdf.total_bounds
        else:
            top, bottom, right, left = bbox
        ax.set_xlim((left, right))
        ax.set_ylim((bottom, top))

    # configure axes appearance
    ax = _config_ax(ax, axis_off, equal_aspect, top, bottom, crs=gdf.crs)

    fig, ax = _save_and_show(
        fig=fig,
        ax=ax,
        save=save,
        show=show,
        close=close,
        filepath=filepath,
        dpi=dpi,
        axis_off=axis_off,
    )

    return fig, ax
