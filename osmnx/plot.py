"""Plot spatial geometries, street networks, and routes."""

import os
import warnings

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


def _warn_deprecated_params(**kwargs):
    """
    Warn about deprecated params.

    Parameters
    ----------
    kwargs

    Returns
    -------
    None
    """
    params = [k for k, v in kwargs.items() if v is not None]
    if len(params) > 0:

        param_str = ", ".join(params)
        msg = f"The {param_str} parameter(s) have been deprecated and will be removed in the next release. "

        cond1 = "fig_height" in kwargs and kwargs["fig_height"] is not None
        cond2 = "fig_width" in kwargs and kwargs["fig_width"] is not None
        if cond1 or cond2:
            msg += "Note, fig_height and fig_width are replaced by the figsize parameter, use that instead. "

        cond1 = "filename" in kwargs and kwargs["filename"] is not None
        cond2 = "file_format" in kwargs and kwargs["file_format"] is not None
        if cond1 or cond2:
            msg += "Note, filename and file_format are replaced by the filepath parameter, use that instead. "

        if "orig_dest_node_size" in kwargs and kwargs["orig_dest_node_size"] is not None:
            msg += "Note, orig_dest_node_size is replaced by the orig_dest_size parameter, use that instead. "

        warnings.warn(msg)


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
        the size of the plotting figure
    margin : numeric
        the size of the figure margins
    axis_off : bool
        if True, disable the matplotlib axes display

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axis
    """
    warnings.warn(
        "the plot_shape function has been deprecated and will be removed in the next release, use gdf.plot() instead"
    )

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

    # adjust the axis margins and limits around the image and make axes
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


def get_edge_colors_by_attr(G, attr, num_bins=5, cmap="viridis", start=0, stop=1, na_color="none"):
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


def _save_and_show(fig, ax, save, show, close, filename, file_format, dpi, axis_off):
    """
    Save a figure to disk and show it, as specified.

    Parameters
    ----------
    fig : figure
        matplotlib figure
    ax : axis
        matplotlib axis
    save : bool
        whether to save the figure to disk or not
    show : bool
        whether to display the figure or not
    close : bool
        close the figure (only if show equals False) to prevent display
    filename : string
        the name of the file to save
    file_format : string
        the format of the file to save (e.g., 'jpg', 'png', 'svg')
    dpi : int
        the resolution of the image file if saving
    axis_off : bool
        if True matplotlib axis was turned off by plot_graph so constrain the
        saved figure's extent to the interior of the axis

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axis
    """
    # save the figure if specified
    if save:
        # create the save folder if it doesn't already exist
        if not os.path.exists(settings.imgs_folder):
            os.makedirs(settings.imgs_folder)
        path_filename = os.path.join(settings.imgs_folder, os.extsep.join([filename, file_format]))

        if file_format == "svg":
            # if the file_format is svg, prep the fig/ax a bit for saving
            ax.axis("off")
            ax.set_position([0, 0, 1, 1])
            ax.patch.set_alpha(0.0)
            fig.patch.set_alpha(0.0)
            fig.savefig(
                path_filename,
                bbox_inches=0,
                format=file_format,
                facecolor=fig.get_facecolor(),
                transparent=True,
            )
        else:
            if axis_off:
                # if axis is turned off, constrain the saved figure's extent to
                # the interior of the axis
                extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            else:
                extent = "tight"
            fig.savefig(
                path_filename,
                dpi=dpi,
                bbox_inches=extent,
                format=file_format,
                facecolor=fig.get_facecolor(),
                transparent=True,
            )
        utils.log("Saved the figure to disk")

    # show the figure if specified
    if show:
        plt.show()
        utils.log("Showed the plot")
    # if show=False, close the figure if close=True to prevent display
    elif close:
        plt.close()

    return fig, ax


def plot_graph(
    G,
    bbox=None,
    fig_height=None,
    fig_width=None,
    margin=0.02,
    axis_off=None,
    equal_aspect=None,
    bgcolor="w",
    show=True,
    save=False,
    close=True,
    file_format=None,
    filename=None,
    dpi=300,
    annotate=None,
    node_color="#66ccff",
    node_size=15,
    node_alpha=1,
    node_edgecolor="none",
    node_zorder=1,
    edge_color="#999999",
    edge_linewidth=1,
    edge_alpha=1,
    use_geom=None,
    figsize=None,
    filepath=None,
):
    """
    Plot a networkx spatial graph.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    bbox : tuple
        bounding box as north,south,east,west - if None will calculate from
        spatial extents of data. if passing a bbox, you probably also want to
        pass margin=0 to constrain it.
    fig_height : int
        deprecated, do not use
    fig_width : int
        deprecated, do not use
    margin : float
        relative margin around the figure
    axis_off : bool
        deprecated, do not use
    equal_aspect : bool
        deprecated, do not use
    bgcolor : string
        the background color of the figure and axis
    show : bool
        if True, show the figure
    save : bool
        if True, save the figure as an image file to disk
    close : bool
        close the figure (only if show equals False) to prevent display
    file_format : string
        deprecated, do not use
    filename : string
        deprecated, do not use
    dpi : int
        the resolution of the image file if saving
    annotate : bool
        deprecated, do not use
    node_color : string
        the color of the nodes. color is passed to matplotlib
    node_size : int
        the size of the nodes
    node_alpha : float
        the opacity of the nodes. if you passed RGBA values to node_color,
        then set this to None to use the alpha channel in node_color
    node_edgecolor : string
        the color of the node's marker's border
    node_zorder : int
        zorder to plot nodes, edges are always 2, so make node_zorder 1 to plot
        nodes beneath them or 3 to plot nodes atop them
    edge_color : string
        the color of the edges' lines. color is passed to matplotlib.
    edge_linewidth : float
        the width of the edges' lines
    edge_alpha : float
        the opacity of the edges' lines. if you passed RGBA values to edge_color,
        then set this to None to use the alpha channel in edge_color
    use_geom : bool
        deprecated, do not use
    figsize : tuple
        figure (width, height)
    filepath : string
        filename.ext to save image in settings.imgs_folder

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axis
    """
    _warn_deprecated_params(
        fig_height=fig_height,
        fig_width=fig_width,
        file_format=file_format,
        filename=filename,
        annotate=annotate,
        use_geom=use_geom,
        axis_off=axis_off,
        equal_aspect=equal_aspect,
    )
    if axis_off is None:
        axis_off = True
    if equal_aspect is None:
        equal_aspect = False
    if fig_height is None:
        fig_height = 6
    if file_format is None:
        file_format = "png"
    if filename is None:
        filename = "temp"
    if annotate is None:
        annotate = False
    if use_geom is None:
        use_geom = True
    if figsize is not None:
        fig_width, fig_height = figsize
    if filepath is not None:
        folders, filename_ext = os.path.split(filepath)
        filename, file_format = os.path.splitext(filename_ext)
        file_format = file_format.strip(".")

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

    # if caller did not pass in a fig_width, calculate it proportionately from
    # the fig_height and bounding box aspect ratio
    bbox_aspect_ratio = (north - south) / (east - west)
    if fig_width is None:
        fig_width = fig_height / bbox_aspect_ratio

    # create the figure and axis
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), facecolor=bgcolor)
    ax.set_facecolor(bgcolor)

    # draw the edges as lines from node to node
    lines = []
    for u, v, data in G.edges(keys=False, data=True):
        if "geometry" in data and use_geom:
            # if it has a geometry attribute (a list of line segments), add them
            # to the list of lines to plot
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

    # add the lines to the axis as a linecollection
    lc = LineCollection(
        lines, colors=edge_color, linewidths=edge_linewidth, alpha=edge_alpha, zorder=2
    )
    ax.add_collection(lc)
    utils.log("Drew the graph edges")

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

    # configure axis appearance
    xaxis = ax.get_xaxis()
    yaxis = ax.get_yaxis()

    xaxis.get_major_formatter().set_useOffset(False)
    yaxis.get_major_formatter().set_useOffset(False)

    # if axis_off, turn off the axis display set the margins to zero and point
    # the ticks in so there's no space around the plot
    if axis_off:
        ax.axis("off")
        ax.margins(0)
        ax.tick_params(which="both", direction="in")
        xaxis.set_visible(False)
        yaxis.set_visible(False)
        fig.canvas.draw()

    if equal_aspect:
        # make everything square
        ax.set_aspect("equal")
        fig.canvas.draw()
    else:
        # if the graph is not projected, conform the aspect ratio to not stretch the plot
        if G.graph["crs"] == settings.default_crs:
            coslat = np.cos((min(node_Ys) + max(node_Ys)) / 2.0 / 180.0 * np.pi)
            ax.set_aspect(1.0 / coslat)
            fig.canvas.draw()

    # annotate the axis with node IDs if annotate=True
    if annotate:
        for node, data in G.nodes(data=True):
            ax.annotate(node, xy=(data["x"], data["y"]))

    # save and show the figure as specified
    fig, ax = _save_and_show(fig, ax, save, show, close, filename, file_format, dpi, axis_off)
    return fig, ax


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


def plot_graph_route(
    G,
    route,
    bbox=None,
    fig_height=None,
    fig_width=None,
    margin=0.02,
    bgcolor="w",
    axis_off=None,
    show=True,
    save=False,
    close=True,
    file_format=None,
    filename=None,
    dpi=300,
    annotate=None,
    node_color="#999999",
    node_size=15,
    node_alpha=1,
    node_edgecolor="none",
    node_zorder=1,
    edge_color="#999999",
    edge_linewidth=1,
    edge_alpha=1,
    use_geom=None,
    origin_point=None,
    destination_point=None,
    route_color="r",
    route_linewidth=4,
    route_alpha=0.5,
    orig_dest_node_alpha=None,
    orig_dest_node_size=None,
    orig_dest_node_color=None,
    orig_dest_point_color=None,
    figsize=None,
    filepath=None,
    orig_dest_size=None,
):
    """
    Plot a route along a networkx spatial graph.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    route : list
        the route as a list of nodes
    bbox : tuple
        bounding box as north,south,east,west - if None will calculate from
        spatial extents of data. if passing a bbox, you probably also want to
        pass margin=0 to constrain it.
    fig_height : int
        deprecated, do not use
    fig_width : int
        deprecated, do not use
    margin : float
        relative margin around the figure
    axis_off : bool
        deprecated, do not use
    bgcolor : string
        the background color of the figure and axis
    show : bool
        if True, show the figure
    save : bool
        if True, save the figure as an image file to disk
    close : bool
        close the figure (only if show equals False) to prevent display
    file_format : string
        deprecated, do not use
    filename : string
        deprecated, do not use
    dpi : int
        the resolution of the image file if saving
    annotate : bool
        deprecated, do not use
    node_color : string
        the color of the nodes
    node_size : int
        the size of the nodes
    node_alpha : float
        the opacity of the nodes
    node_edgecolor : string
        the color of the node's marker's border
    node_zorder : int
        zorder to plot nodes, edges are always 2, so make node_zorder 1 to plot
        nodes beneath them or 3 to plot nodes atop them
    edge_color : string
        the color of the edges' lines
    edge_linewidth : float
        the width of the edges' lines
    edge_alpha : float
        the opacity of the edges' lines
    use_geom : bool
        deprecated, do not use
    origin_point : tuple
        deprecated, do not use
    destination_point : tuple
        deprecated, do not use
    route_color : string
        the color of the route
    route_linewidth : int
        the width of the route line
    route_alpha : float
        the opacity of the route line
    orig_dest_node_alpha : float
        deprecated, do not use
    orig_dest_node_size : int
        deprecated, do not use
    orig_dest_node_color : string
        deprecated, do not use
    orig_dest_point_color : string
        deprecated, do not use
    figsize : tuple
        figure (width, height)
    filepath : string
        filename.ext to save image in settings.imgs_folder
    orig_dest_size : int
        the size of the origin and destination nodes

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axis
    """
    _warn_deprecated_params(
        fig_height=fig_height,
        fig_width=fig_width,
        file_format=file_format,
        filename=filename,
        annotate=annotate,
        use_geom=use_geom,
        origin_point=origin_point,
        destination_point=destination_point,
        orig_dest_node_alpha=orig_dest_node_alpha,
        orig_dest_node_color=orig_dest_node_color,
        orig_dest_point_color=orig_dest_point_color,
        orig_dest_node_size=orig_dest_node_size,
        axis_off=axis_off,
    )
    if axis_off is None:
        axis_off = True
    if fig_height is None:
        fig_height = 6
    if file_format is None:
        file_format = "png"
    if filename is None:
        filename = "temp"
    if annotate is None:
        annotate = False
    if use_geom is None:
        use_geom = True
    if orig_dest_node_alpha is None:
        orig_dest_node_alpha = 0.5
    if orig_dest_node_size is None:
        orig_dest_node_size = 100
    if orig_dest_node_color is None:
        orig_dest_node_color = "r"
    if orig_dest_point_color is None:
        orig_dest_point_color = "b"
    if figsize is not None:
        fig_width, fig_height = figsize
    if filepath is not None:
        folders, filename_ext = os.path.split(filepath)
        filename, file_format = os.path.splitext(filename_ext)
        file_format = file_format.strip(".")
    if orig_dest_size is not None:
        orig_dest_node_size = orig_dest_size

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # plot the graph but not the route
        fig, ax = plot_graph(
            G,
            bbox=bbox,
            fig_height=fig_height,
            fig_width=fig_width,
            margin=margin,
            axis_off=axis_off,
            bgcolor=bgcolor,
            show=False,
            save=False,
            close=False,
            filename=filename,
            dpi=dpi,
            annotate=annotate,
            node_color=node_color,
            node_size=node_size,
            node_alpha=node_alpha,
            node_edgecolor=node_edgecolor,
            node_zorder=node_zorder,
            edge_color=edge_color,
            edge_linewidth=edge_linewidth,
            edge_alpha=edge_alpha,
            use_geom=use_geom,
        )

    # the origin and destination nodes are the first and last nodes in the route
    origin_node = route[0]
    destination_node = route[-1]

    if origin_point is None or destination_point is None:
        # if caller didn't pass points, use the first and last node in route as
        # origin/destination
        origin_destination_lats = (G.nodes[origin_node]["y"], G.nodes[destination_node]["y"])
        origin_destination_lons = (G.nodes[origin_node]["x"], G.nodes[destination_node]["x"])
    else:
        # otherwise, use the passed points as origin/destination
        origin_destination_lats = (origin_point[0], destination_point[0])
        origin_destination_lons = (origin_point[1], destination_point[1])
        orig_dest_node_color = orig_dest_point_color

    # scatter the origin and destination points
    ax.scatter(
        origin_destination_lons,
        origin_destination_lats,
        s=orig_dest_node_size,
        c=orig_dest_node_color,
        alpha=orig_dest_node_alpha,
        edgecolor=node_edgecolor,
        zorder=4,
    )

    # plot the route lines
    lines = _node_list_to_coordinate_lines(G, route, use_geom)

    # add the lines to the axis as a linecollection
    lc = LineCollection(
        lines, colors=route_color, linewidths=route_linewidth, alpha=route_alpha, zorder=3
    )
    ax.add_collection(lc)

    # save and show the figure as specified
    fig, ax = _save_and_show(fig, ax, save, show, close, filename, file_format, dpi, axis_off)
    return fig, ax


def plot_graph_routes(
    G,
    routes,
    bbox=None,
    fig_height=None,
    fig_width=None,
    margin=0.02,
    bgcolor="w",
    axis_off=None,
    show=True,
    save=False,
    close=True,
    file_format=None,
    filename=None,
    dpi=300,
    annotate=None,
    node_color="#999999",
    node_size=15,
    node_alpha=1,
    node_edgecolor="none",
    node_zorder=1,
    edge_color="#999999",
    edge_linewidth=1,
    edge_alpha=1,
    use_geom=None,
    orig_dest_points=None,
    route_color="r",
    route_linewidth=4,
    route_alpha=0.5,
    orig_dest_node_alpha=None,
    orig_dest_node_size=None,
    orig_dest_node_color=None,
    orig_dest_point_color=None,
    figsize=None,
    filepath=None,
    orig_dest_size=None,
):
    """
    Plot several routes along a networkx spatial graph.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    routes : list
        the routes as a list of lists of nodes
    bbox : tuple
        bounding box as north,south,east,west - if None will calculate from
        spatial extents of data. if passing a bbox, you probably also want to
        pass margin=0 to constrain it.
    fig_height : int
        deprecated, do not use
    fig_width : int
        deprecated, do not use
    margin : float
        relative margin around the figure
    axis_off : bool
        if True turn off the matplotlib axis
    bgcolor : string
        the background color of the figure and axis
    show : bool
        if True, show the figure
    save : bool
        if True, save the figure as an image file to disk
    close : bool
        close the figure (only if show equals False) to prevent display
    file_format : string
        deprecated, do not use
    filename : string
        deprecated, do not use
    dpi : int
        the resolution of the image file if saving
    annotate : bool
        deprecated, do not use
    node_color : string
        the color of the nodes
    node_size : int
        the size of the nodes
    node_alpha : float
        the opacity of the nodes
    node_edgecolor : string
        the color of the node's marker's border
    node_zorder : int
        zorder to plot nodes, edges are always 2, so make node_zorder 1 to plot
        nodes beneath them or 3 to plot nodes atop them
    edge_color : string
        the color of the edges' lines
    edge_linewidth : float
        the width of the edges' lines
    edge_alpha : float
        the opacity of the edges' lines
    use_geom : bool
        deprecated, do not use
    orig_dest_points : list of tuples
        deprecated, do not use
    route_color : string
        route color (note: will be renamed `route_colors` and take list in
        next release)
    route_linewidth : int
        the width of the route line
    route_alpha : float
        the opacity of the route line
    orig_dest_node_alpha : float
        deprecated, do not use
    orig_dest_node_size : int
        deprecated, do not use
    orig_dest_node_color : string
        deprecated, do not use
    orig_dest_point_color : string
        deprecated, do not use
    figsize : tuple
        figure (width, height)
    filepath : string
        filename.ext to save image in settings.imgs_folder
    orig_dest_size : int
        the size of the origin and destination nodes

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axis
    """
    _warn_deprecated_params(
        fig_height=fig_height,
        fig_width=fig_width,
        file_format=file_format,
        filename=filename,
        annotate=annotate,
        use_geom=use_geom,
        orig_dest_points=orig_dest_points,
        orig_dest_node_alpha=orig_dest_node_alpha,
        orig_dest_node_color=orig_dest_node_color,
        orig_dest_point_color=orig_dest_point_color,
        orig_dest_node_size=orig_dest_node_size,
        axis_off=axis_off,
    )
    if axis_off is None:
        axis_off = True
    if fig_height is None:
        fig_height = 6
    if file_format is None:
        file_format = "png"
    if filename is None:
        filename = "temp"
    if annotate is None:
        annotate = False
    if use_geom is None:
        use_geom = True
    if orig_dest_node_alpha is None:
        orig_dest_node_alpha = 0.5
    if orig_dest_node_size is None:
        orig_dest_node_size = 100
    if orig_dest_node_color is None:
        orig_dest_node_color = "r"
    if orig_dest_point_color is None:
        orig_dest_point_color = "b"
    if figsize is not None:
        fig_width, fig_height = figsize
    if filepath is not None:
        folders, filename_ext = os.path.split(filepath)
        filename, file_format = os.path.splitext(filename_ext)
        file_format = file_format.strip(".")
    if orig_dest_size is not None:
        orig_dest_node_size = orig_dest_size

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # plot the graph but not the routes
        fig, ax = plot_graph(
            G,
            bbox=bbox,
            fig_height=fig_height,
            fig_width=fig_width,
            margin=margin,
            axis_off=axis_off,
            bgcolor=bgcolor,
            show=False,
            save=False,
            close=False,
            filename=filename,
            dpi=dpi,
            annotate=annotate,
            node_color=node_color,
            node_size=node_size,
            node_alpha=node_alpha,
            node_edgecolor=node_edgecolor,
            node_zorder=node_zorder,
            edge_color=edge_color,
            edge_linewidth=edge_linewidth,
            edge_alpha=edge_alpha,
            use_geom=use_geom,
        )

    # save coordinates of the given reference points
    orig_dest_points_lats = []
    orig_dest_points_lons = []

    if orig_dest_points is None:
        # if caller didn't pass points, use the first and last node in each route as
        # origin/destination points
        for route in routes:
            origin_node = route[0]
            destination_node = route[-1]
            orig_dest_points_lats.append(G.nodes[origin_node]["y"])
            orig_dest_points_lats.append(G.nodes[destination_node]["y"])
            orig_dest_points_lons.append(G.nodes[origin_node]["x"])
            orig_dest_points_lons.append(G.nodes[destination_node]["x"])

    else:
        # otherwise, use the passed points as origin/destination points
        for point in orig_dest_points:
            orig_dest_points_lats.append(point[0])
            orig_dest_points_lons.append(point[1])
        orig_dest_node_color = orig_dest_point_color

    # scatter the origin and destination points
    ax.scatter(
        orig_dest_points_lons,
        orig_dest_points_lats,
        s=orig_dest_node_size,
        c=orig_dest_node_color,
        alpha=orig_dest_node_alpha,
        edgecolor=node_edgecolor,
        zorder=4,
    )

    # plot the routes lines
    lines = []
    for route in routes:
        lines.extend(_node_list_to_coordinate_lines(G, route, use_geom))

    # add the lines to the axis as a linecollection
    lc = LineCollection(
        lines, colors=route_color, linewidths=route_linewidth, alpha=route_alpha, zorder=3
    )
    ax.add_collection(lc)

    # save and show the figure as specified
    fig, ax = _save_and_show(fig, ax, save, show, close, filename, file_format, dpi, axis_off)
    return fig, ax


def plot_figure_ground(
    G=None,
    address=None,
    point=None,
    dist=805,
    network_type="drive_service",
    street_widths=None,
    default_width=4,
    fig_length=None,
    edge_color="w",
    bgcolor="#333333",
    smooth_joints=True,
    filename=None,
    file_format=None,
    show=False,
    save=True,
    close=True,
    dpi=300,
    figsize=None,
    filepath=None,
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
    fig_length : numeric
        deprecated, do not use
    edge_color : string
        the color of the streets
    bgcolor : string
        the color of the background
    smooth_joints : bool
        if True, plot nodes same width as streets to smooth line joints and
        prevent cracks between them from showing
    filename : string
        deprecated, do not use
    file_format : string
        deprecated, do not use
    show : bool
        if True, show the figure
    save : bool
        if True, save the figure as an image file to disk
    close : bool
        close the figure (only if show equals False) to prevent display
    dpi : int
        the resolution of the image file if saving
    figsize : tuple
        figure width, height (should be equal to each other)
    filepath : string
        filename.ext to save image in settings.imgs_folder

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axis
    """
    _warn_deprecated_params(fig_length=fig_length, file_format=file_format, filename=filename)
    if fig_length is None:
        fig_length = 8
    if file_format is None:
        file_format = "png"
    if filename is None:
        filename = "temp"
    if figsize is not None:
        fig_length = figsize[0]
    if filepath is not None:
        folders, filename_ext = os.path.split(filepath)
        filename, file_format = os.path.splitext(filename_ext)
        file_format = file_format.strip(".")

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

    # create a filename if one was not passed
    if filename is None and save:
        filename = "figure_ground"

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # plot the figure
        fig, ax = plot_graph(
            G_undir,
            bbox=bbox,
            fig_height=fig_length,
            margin=0,
            axis_off=True,
            equal_aspect=False,
            bgcolor=bgcolor,
            node_size=node_sizes,
            node_color=edge_color,
            edge_linewidth=edge_linewidths,
            edge_color=edge_color,
            show=show,
            save=save,
            close=close,
            filename=filename,
            file_format=file_format,
            dpi=dpi,
        )

    return fig, ax


def plot_footprints(
    gdf,
    fig=None,
    ax=None,
    figsize=None,
    color="#333333",
    bgcolor="w",
    set_bounds=True,
    bbox=None,
    save=False,
    show=True,
    close=False,
    filename=None,
    file_format=None,
    dpi=600,
    filepath=None,
):
    """
    Plot a GeoDataFrame of footprints.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame of footprints
    fig : figure
        matplotlib figure
    ax : axis
        matplotlib axis
    figsize : tuple
        (width, height) size of matplotlib figure
    color : string
        the color of the footprints
    bgcolor : string
        the background color of the plot
    set_bounds : bool
        if True, set bounds from either passed-in bbox or the spatial extent of the gdf
    bbox : tuple
        if True and if set_bounds is True, set the display bounds to this bbox
    save : bool
        whether to save the figure to disk or not
    show : bool
        whether to display the figure or not
    close : bool
        close the figure (only if show equals False) to prevent display
    filename : string
        deprecated, do not use
    file_format : string
        deprecated, do not use
    dpi : int
        the resolution of the image file if saving
    filepath : string
        filename.ext to save image in settings.imgs_folder

    Returns
    -------
    fig, ax : tuple
        matplotlib figure, axis
    """
    _warn_deprecated_params(file_format=file_format, filename=filename)
    if file_format is None:
        file_format = "png"
    if filename is None:
        filename = "temp"
    if filepath is not None:
        folders, filename_ext = os.path.split(filepath)
        filename, file_format = os.path.splitext(filename_ext)
        file_format = file_format.strip(".")

    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=figsize, facecolor=bgcolor)
        ax.set_facecolor(bgcolor)

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

    # turn off the axis display set the margins to zero and point the ticks in
    # so there's no space around the plot
    ax.axis("off")
    ax.margins(0)
    ax.tick_params(which="both", direction="in")
    fig.canvas.draw()

    # make everything square
    ax.set_aspect("equal")
    fig.canvas.draw()

    fig, ax = _save_and_show(
        fig=fig,
        ax=ax,
        save=save,
        show=show,
        close=close,
        filename=filename,
        file_format=file_format,
        dpi=dpi,
        axis_off=True,
    )

    return fig, ax
