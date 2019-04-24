################################################################################
# Module: plot.py
# Description: Plot spatial geometries, street networks, and routes
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

import time
import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
from descartes import PolygonPatch
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon

from . import settings
from .core import graph_from_address
from .core import graph_from_point
from .core import bbox_from_point
from .projection import project_graph
from .save_load import graph_to_gdfs
from .simplify import simplify_graph
from .utils import log


# folium is an optional dependency for the folium plotting functions
try:
    import folium
except ImportError as e:
    folium = None


def plot_shape(gdf, fc='#cbe0f0', ec='#999999', linewidth=1, alpha=1,
               figsize=(6,6), margin=0.02, axis_off=True):
    """
    Plot a GeoDataFrame of place boundary geometries.

    Parameters
    ----------
    gdf : GeoDataFrame
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
    """

    # if facecolor or edgecolor is a string instead of a list, make sure we have
    # as many colors as gdf elements
    if isinstance(fc, str):
        fc = [fc] * len(gdf)
    if isinstance(ec, str):
        ec = [ec] * len(gdf)

    # plot the geometries one at a time
    fig, ax = plt.subplots(figsize=figsize)
    for geometry, facecolor, edgecolor in zip(gdf['geometry'], fc, ec):
        if isinstance(geometry, (Polygon, MultiPolygon)):
            if isinstance(geometry, Polygon):
                geometry = MultiPolygon([geometry])
            for polygon in geometry:
                patch = PolygonPatch(polygon, fc=facecolor, ec=edgecolor, linewidth=linewidth, alpha=alpha)
                ax.add_patch(patch)
        else:
            raise ValueError('All geometries in GeoDataFrame must be shapely Polygons or MultiPolygons')

    # adjust the axis margins and limits around the image and make axes
    # equal-aspect
    west, south, east, north = gdf.unary_union.bounds
    margin_ns = (north - south) * margin
    margin_ew = (east - west) * margin
    ax.set_ylim((south - margin_ns, north + margin_ns))
    ax.set_xlim((west - margin_ew, east + margin_ew))
    ax.set_aspect(aspect='equal', adjustable='box')
    if axis_off:
        ax.axis('off')

    plt.show()
    return fig, ax


def rgb_color_list_to_hex(color_list):
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
    color_list_rgb = [[int(x*255) for x in c[0:3]] for c in color_list]
    color_list_hex = ['#{:02X}{:02X}{:02X}'.format(rgb[0], rgb[1], rgb[2]) for rgb in color_list_rgb]
    return color_list_hex


def get_colors(n, cmap='viridis', start=0., stop=1., alpha=1., return_hex=False):
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
        colors = rgb_color_list_to_hex(colors)
    return colors


def get_node_colors_by_attr(G, attr, num_bins=None, cmap='viridis', start=0, stop=1, na_color='none'):
    """
    Get a list of node colors by binning some continuous-variable attribute into
    quantiles.

    Parameters
    ----------
    G : networkx multidigraph
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
    list
    """
    if num_bins is None:
        num_bins=len(G.nodes())
    bin_labels = range(num_bins)
    attr_values = pd.Series([data[attr] for node, data in G.nodes(data=True)])
    cats = pd.qcut(x=attr_values, q=num_bins, labels=bin_labels)
    colors = get_colors(num_bins, cmap, start, stop)
    node_colors = [colors[int(cat)] if pd.notnull(cat) else na_color for cat in cats]
    return node_colors


def get_edge_colors_by_attr(G, attr, num_bins=5, cmap='viridis', start=0, stop=1, na_color='none'):
    """
    Get a list of edge colors by binning some continuous-variable attribute into
    quantiles.

    Parameters
    ----------
    G : networkx multidigraph
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
    list
    """
    if num_bins is None:
        num_bins=len(G.edges())
    bin_labels = range(num_bins)
    attr_values = pd.Series([data[attr] for u, v, key, data in G.edges(keys=True, data=True)])
    cats = pd.qcut(x=attr_values, q=num_bins, labels=bin_labels)
    colors = get_colors(num_bins, cmap, start, stop)
    edge_colors = [colors[int(cat)] if pd.notnull(cat) else na_color for cat in cats]
    return edge_colors


def save_and_show(fig, ax, save, show, close, filename, file_format, dpi, axis_off):
    """
    Save a figure to disk and show it, as specified.

    Parameters
    ----------
    fig : figure
    ax : axis
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
    """
    # save the figure if specified
    if save:
        start_time = time.time()

        # create the save folder if it doesn't already exist
        if not os.path.exists(settings.imgs_folder):
            os.makedirs(settings.imgs_folder)
        path_filename = os.path.join(settings.imgs_folder, os.extsep.join([filename, file_format]))

        if file_format == 'svg':
            # if the file_format is svg, prep the fig/ax a bit for saving
            ax.axis('off')
            ax.set_position([0, 0, 1, 1])
            ax.patch.set_alpha(0.)
            fig.patch.set_alpha(0.)
            fig.savefig(path_filename, bbox_inches=0, format=file_format, facecolor=fig.get_facecolor(), transparent=True)
        else:
            if axis_off:
                # if axis is turned off, constrain the saved figure's extent to
                # the interior of the axis
                extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            else:
                extent = 'tight'
            fig.savefig(path_filename, dpi=dpi, bbox_inches=extent, format=file_format, facecolor=fig.get_facecolor(), transparent=True)
        log('Saved the figure to disk in {:,.2f} seconds'.format(time.time()-start_time))

    # show the figure if specified
    if show:
        start_time = time.time()
        plt.show()
        log('Showed the plot in {:,.2f} seconds'.format(time.time()-start_time))
    # if show=False, close the figure if close=True to prevent display
    elif close:
        plt.close()

    return fig, ax


def plot_graph(G, bbox=None, fig_height=6, fig_width=None, margin=0.02,
               axis_off=True, equal_aspect=False, bgcolor='w', show=True,
               save=False, close=True, file_format='png', filename='temp',
               dpi=300, annotate=False, node_color='#66ccff', node_size=15,
               node_alpha=1, node_edgecolor='none', node_zorder=1,
               edge_color='#999999', edge_linewidth=1, edge_alpha=1,
               use_geom=True):
    """
    Plot a networkx spatial graph.

    Parameters
    ----------
    G : networkx multidigraph
    bbox : tuple
        bounding box as north,south,east,west - if None will calculate from
        spatial extents of data. if passing a bbox, you probably also want to
        pass margin=0 to constrain it.
    fig_height : int
        matplotlib figure height in inches
    fig_width : int
        matplotlib figure width in inches
    margin : float
        relative margin around the figure
    axis_off : bool
        if True turn off the matplotlib axis
    equal_aspect : bool
        if True set the axis aspect ratio equal
    bgcolor : string
        the background color of the figure and axis
    show : bool
        if True, show the figure
    save : bool
        if True, save the figure as an image file to disk
    close : bool
        close the figure (only if show equals False) to prevent display
    file_format : string
        the format of the file to save (e.g., 'jpg', 'png', 'svg')
    filename : string
        the name of the file if saving
    dpi : int
        the resolution of the image file if saving
    annotate : bool
        if True, annotate the nodes in the figure
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
        if True, use the spatial geometry attribute of the edges to draw
        geographically accurate edges, rather than just lines straight from node
        to node

    Returns
    -------
    fig, ax : tuple
    """

    log('Begin plotting the graph...')
    node_Xs = [float(x) for _, x in G.nodes(data='x')]
    node_Ys = [float(y) for _, y in G.nodes(data='y')]

    # get north, south, east, west values either from bbox parameter or from the
    # spatial extent of the edges' geometries
    if bbox is None:
        edges = graph_to_gdfs(G, nodes=False, fill_edge_geometry=True)
        west, south, east, north = edges.total_bounds
    else:
        north, south, east, west = bbox

    # if caller did not pass in a fig_width, calculate it proportionately from
    # the fig_height and bounding box aspect ratio
    bbox_aspect_ratio = (north-south)/(east-west)
    if fig_width is None:
        fig_width = fig_height / bbox_aspect_ratio

    # create the figure and axis
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), facecolor=bgcolor)
    ax.set_facecolor(bgcolor)

    # draw the edges as lines from node to node
    start_time = time.time()
    lines = []
    for u, v, data in G.edges(keys=False, data=True):
        if 'geometry' in data and use_geom:
            # if it has a geometry attribute (a list of line segments), add them
            # to the list of lines to plot
            xs, ys = data['geometry'].xy
            lines.append(list(zip(xs, ys)))
        else:
            # if it doesn't have a geometry attribute, the edge is a straight
            # line from node to node
            x1 = G.nodes[u]['x']
            y1 = G.nodes[u]['y']
            x2 = G.nodes[v]['x']
            y2 = G.nodes[v]['y']
            line = [(x1, y1), (x2, y2)]
            lines.append(line)

    # add the lines to the axis as a linecollection
    lc = LineCollection(lines, colors=edge_color, linewidths=edge_linewidth, alpha=edge_alpha, zorder=2)
    ax.add_collection(lc)
    log('Drew the graph edges in {:,.2f} seconds'.format(time.time()-start_time))

    # scatter plot the nodes
    ax.scatter(node_Xs, node_Ys, s=node_size, c=node_color, alpha=node_alpha, edgecolor=node_edgecolor, zorder=node_zorder)

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
        ax.axis('off')
        ax.margins(0)
        ax.tick_params(which='both', direction='in')
        xaxis.set_visible(False)
        yaxis.set_visible(False)
        fig.canvas.draw()

    if equal_aspect:
        # make everything square
        ax.set_aspect('equal')
        fig.canvas.draw()
    else:
        # if the graph is not projected, conform the aspect ratio to not stretch the plot
        if G.graph['crs'] == settings.default_crs:
            coslat = np.cos((min(node_Ys) + max(node_Ys)) / 2. / 180. * np.pi)
            ax.set_aspect(1. / coslat)
            fig.canvas.draw()

    # annotate the axis with node IDs if annotate=True
    if annotate:
        for node, data in G.nodes(data=True):
            ax.annotate(node, xy=(data['x'], data['y']))

    # save and show the figure as specified
    fig, ax = save_and_show(fig, ax, save, show, close, filename, file_format, dpi, axis_off)
    return fig, ax


def node_list_to_coordinate_lines(G, node_list, use_geom=True):
    """
    Given a list of nodes, return a list of lines that together follow the path
    defined by the list of nodes.

    Parameters
    ----------
    G : networkx multidigraph
    route : list
        the route as a list of nodes
    use_geom : bool
        if True, use the spatial geometry attribute of the edges to draw
        geographically accurate edges, rather than just lines straight from node
        to node

    Returns
    -------
    lines : list of lines given as pairs ( (x_start, y_start), (x_stop, y_stop) )
    """
    edge_nodes = list(zip(node_list[:-1], node_list[1:]))
    lines = []
    for u, v in edge_nodes:
        # if there are parallel edges, select the shortest in length
        data = min(G.get_edge_data(u, v).values(), key=lambda x: x['length'])

        # if it has a geometry attribute (ie, a list of line segments)
        if 'geometry' in data and use_geom:
            # add them to the list of lines to plot
            xs, ys = data['geometry'].xy
            lines.append(list(zip(xs, ys)))
        else:
            # if it doesn't have a geometry attribute, the edge is a straight
            # line from node to node
            x1 = G.nodes[u]['x']
            y1 = G.nodes[u]['y']
            x2 = G.nodes[v]['x']
            y2 = G.nodes[v]['y']
            line = [(x1, y1), (x2, y2)]
            lines.append(line)
    return lines

def plot_graph_route(G, route, bbox=None, fig_height=6, fig_width=None,
                     margin=0.02, bgcolor='w', axis_off=True, show=True,
                     save=False, close=True, file_format='png', filename='temp',
                     dpi=300, annotate=False, node_color='#999999',
                     node_size=15, node_alpha=1, node_edgecolor='none',
                     node_zorder=1, edge_color='#999999', edge_linewidth=1,
                     edge_alpha=1, use_geom=True, origin_point=None,
                     destination_point=None, route_color='r', route_linewidth=4,
                     route_alpha=0.5, orig_dest_node_alpha=0.5,
                     orig_dest_node_size=100, orig_dest_node_color='r',
                     orig_dest_point_color='b'):
    """
    Plot a route along a networkx spatial graph.

    Parameters
    ----------
    G : networkx multidigraph
    route : list
        the route as a list of nodes
    bbox : tuple
        bounding box as north,south,east,west - if None will calculate from
        spatial extents of data. if passing a bbox, you probably also want to
        pass margin=0 to constrain it.
    fig_height : int
        matplotlib figure height in inches
    fig_width : int
        matplotlib figure width in inches
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
        the format of the file to save (e.g., 'jpg', 'png', 'svg')
    filename : string
        the name of the file if saving
    dpi : int
        the resolution of the image file if saving
    annotate : bool
        if True, annotate the nodes in the figure
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
        if True, use the spatial geometry attribute of the edges to draw
        geographically accurate edges, rather than just lines straight from node
        to node
    origin_point : tuple
        optional, an origin (lat, lon) point to plot instead of the origin node
    destination_point : tuple
        optional, a destination (lat, lon) point to plot instead of the
        destination node
    route_color : string
        the color of the route
    route_linewidth : int
        the width of the route line
    route_alpha : float
        the opacity of the route line
    orig_dest_node_alpha : float
        the opacity of the origin and destination nodes
    orig_dest_node_size : int
        the size of the origin and destination nodes
    orig_dest_node_color : string
        the color of the origin and destination nodes
    orig_dest_point_color : string
        the color of the origin and destination points if being plotted instead
        of nodes

    Returns
    -------
    fig, ax : tuple
    """

    # plot the graph but not the route
    fig, ax = plot_graph(G, bbox=bbox, fig_height=fig_height, fig_width=fig_width,
                         margin=margin, axis_off=axis_off, bgcolor=bgcolor,
                         show=False, save=False, close=False, filename=filename,
                         dpi=dpi, annotate=annotate, node_color=node_color,
                         node_size=node_size, node_alpha=node_alpha,
                         node_edgecolor=node_edgecolor, node_zorder=node_zorder,
                         edge_color=edge_color, edge_linewidth=edge_linewidth,
                         edge_alpha=edge_alpha, use_geom=use_geom)

    # the origin and destination nodes are the first and last nodes in the route
    origin_node = route[0]
    destination_node = route[-1]

    if origin_point is None or destination_point is None:
        # if caller didn't pass points, use the first and last node in route as
        # origin/destination
        origin_destination_lats = (G.nodes[origin_node]['y'], G.nodes[destination_node]['y'])
        origin_destination_lons = (G.nodes[origin_node]['x'], G.nodes[destination_node]['x'])
    else:
        # otherwise, use the passed points as origin/destination
        origin_destination_lats = (origin_point[0], destination_point[0])
        origin_destination_lons = (origin_point[1], destination_point[1])
        orig_dest_node_color = orig_dest_point_color

    # scatter the origin and destination points
    ax.scatter(origin_destination_lons, origin_destination_lats, s=orig_dest_node_size,
               c=orig_dest_node_color, alpha=orig_dest_node_alpha, edgecolor=node_edgecolor, zorder=4)

    # plot the route lines
    lines = node_list_to_coordinate_lines(G, route, use_geom)

    # add the lines to the axis as a linecollection
    lc = LineCollection(lines, colors=route_color, linewidths=route_linewidth, alpha=route_alpha, zorder=3)
    ax.add_collection(lc)

    # save and show the figure as specified
    fig, ax = save_and_show(fig, ax, save, show, close, filename, file_format, dpi, axis_off)
    return fig, ax

def plot_graph_routes(G, routes, bbox=None, fig_height=6, fig_width=None,
                      margin=0.02, bgcolor='w', axis_off=True, show=True,
                      save=False, close=True, file_format='png', filename='temp',
                      dpi=300, annotate=False, node_color='#999999',
                      node_size=15, node_alpha=1, node_edgecolor='none',
                      node_zorder=1, edge_color='#999999', edge_linewidth=1,
                      edge_alpha=1, use_geom=True, orig_dest_points=None,
                      route_color='r', route_linewidth=4,
                      route_alpha=0.5, orig_dest_node_alpha=0.5,
                      orig_dest_node_size=100, orig_dest_node_color='r',
                      orig_dest_point_color='b'):
    """
    Plot several routes along a networkx spatial graph.

    Parameters
    ----------
    G : networkx multidigraph
    routes : list
        the routes as a list of lists of nodes
    bbox : tuple
        bounding box as north,south,east,west - if None will calculate from
        spatial extents of data. if passing a bbox, you probably also want to
        pass margin=0 to constrain it.
    fig_height : int
        matplotlib figure height in inches
    fig_width : int
        matplotlib figure width in inches
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
        the format of the file to save (e.g., 'jpg', 'png', 'svg')
    filename : string
        the name of the file if saving
    dpi : int
        the resolution of the image file if saving
    annotate : bool
        if True, annotate the nodes in the figure
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
        if True, use the spatial geometry attribute of the edges to draw
        geographically accurate edges, rather than just lines straight from node
        to node
    orig_dest_points : list of tuples
        optional, a group of (lat, lon) points to plot instead of the
        origins and destinations of each route nodes
    route_color : string
        the color of the route
    route_linewidth : int
        the width of the route line
    route_alpha : float
        the opacity of the route line
    orig_dest_node_alpha : float
        the opacity of the origin and destination nodes
    orig_dest_node_size : int
        the size of the origin and destination nodes
    orig_dest_node_color : string
        the color of the origin and destination nodes
    orig_dest_point_color : string
        the color of the origin and destination points if being plotted instead
        of nodes

    Returns
    -------
    fig, ax : tuple
    """

    # plot the graph but not the routes
    fig, ax = plot_graph(G, bbox=bbox, fig_height=fig_height, fig_width=fig_width,
                         margin=margin, axis_off=axis_off, bgcolor=bgcolor,
                         show=False, save=False, close=False, filename=filename,
                         dpi=dpi, annotate=annotate, node_color=node_color,
                         node_size=node_size, node_alpha=node_alpha,
                         node_edgecolor=node_edgecolor, node_zorder=node_zorder,
                         edge_color=edge_color, edge_linewidth=edge_linewidth,
                         edge_alpha=edge_alpha, use_geom=use_geom)

    # save coordinates of the given reference points
    orig_dest_points_lats = []
    orig_dest_points_lons = []

    if orig_dest_points is None:
        # if caller didn't pass points, use the first and last node in each route as
        # origin/destination points
        for route in routes:
            origin_node = route[0]
            destination_node = route[-1]
            orig_dest_points_lats.append(G.nodes[origin_node]['y'])
            orig_dest_points_lats.append(G.nodes[destination_node]['y'])
            orig_dest_points_lons.append(G.nodes[origin_node]['x'])
            orig_dest_points_lons.append(G.nodes[destination_node]['x'])

    else:
        # otherwise, use the passed points as origin/destination points
        for point in orig_dest_points:
            orig_dest_points_lats.append(point[0])
            orig_dest_points_lons.append(point[1])
        orig_dest_node_color = orig_dest_point_color

    # scatter the origin and destination points
    ax.scatter(orig_dest_points_lons, orig_dest_points_lats, s=orig_dest_node_size,
               c=orig_dest_node_color, alpha=orig_dest_node_alpha, edgecolor=node_edgecolor, zorder=4)

    # plot the routes lines
    lines = []
    for route in routes:
        lines.extend(node_list_to_coordinate_lines(G, route, use_geom))

    # add the lines to the axis as a linecollection
    lc = LineCollection(lines, colors=route_color, linewidths=route_linewidth, alpha=route_alpha, zorder=3)
    ax.add_collection(lc)

    # save and show the figure as specified
    fig, ax = save_and_show(fig, ax, save, show, close, filename, file_format, dpi, axis_off)
    return fig, ax

def make_folium_polyline(edge, edge_color, edge_width, edge_opacity, popup_attribute=None):

    """
    Turn a row from the gdf_edges GeoDataFrame into a folium PolyLine with
    attributes.

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
        edge attribute to display in a pop-up when an edge is clicked, if None,
        no popup

    Returns
    -------
    pl : folium.PolyLine
    """

    # check if we were able to import folium successfully
    if not folium:
        raise ImportError('The folium package must be installed to use this optional feature.')

    # locations is a list of points for the polyline
    # folium takes coords in lat,lon but geopandas provides them in lon,lat
    # so we have to flip them around
    locations = list([(lat, lon) for lon, lat in edge['geometry'].coords])

    # if popup_attribute is None, then create no pop-up
    if popup_attribute is None:
        popup = None
    else:
        # folium doesn't interpret html in the html argument (weird), so can't
        # do newlines without an iframe
        popup_text = json.dumps(edge[popup_attribute])
        popup = folium.Popup(html=popup_text)

    # create a folium polyline with attributes
    pl = folium.PolyLine(locations=locations, popup=popup,
                         color=edge_color, weight=edge_width, opacity=edge_opacity)
    return pl


def plot_graph_folium(G, graph_map=None, popup_attribute=None,
                      tiles='cartodbpositron', zoom=1, fit_bounds=True,
                      edge_color='#333333', edge_width=5, edge_opacity=1):
    """
    Plot a graph on an interactive folium web map.

    Note that anything larger than a small city can take a long time to plot and
    create a large web map file that is very slow to load as JavaScript.

    Parameters
    ----------
    G : networkx multidigraph
    graph_map : folium.folium.Map
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

    Returns
    -------
    graph_map : folium.folium.Map
    """

    # check if we were able to import folium successfully
    if not folium:
        raise ImportError('The folium package must be installed to use this optional feature.')

    # create gdf of the graph edges
    gdf_edges = graph_to_gdfs(G, nodes=False, fill_edge_geometry=True)

    # get graph centroid
    x, y = gdf_edges.unary_union.centroid.xy
    graph_centroid = (y[0], x[0])

    # create the folium web map if one wasn't passed-in
    if graph_map is None:
        graph_map = folium.Map(location=graph_centroid, zoom_start=zoom, tiles=tiles)

    # add each graph edge to the map
    for _, row in gdf_edges.iterrows():
        pl = make_folium_polyline(edge=row, edge_color=edge_color, edge_width=edge_width,
                                  edge_opacity=edge_opacity, popup_attribute=popup_attribute)
        pl.add_to(graph_map)

    # if fit_bounds is True, fit the map to the bounds of the route by passing
    # list of lat-lng points as [southwest, northeast]
    if fit_bounds:
        tb = gdf_edges.total_bounds
        bounds = [(tb[1], tb[0]), (tb[3], tb[2])]
        graph_map.fit_bounds(bounds)

    return graph_map


def plot_route_folium(G, route, route_map=None, popup_attribute=None,
                      tiles='cartodbpositron', zoom=1, fit_bounds=True,
                      route_color='#cc0000', route_width=5, route_opacity=1):
    """
    Plot a route on an interactive folium web map.

    Parameters
    ----------
    G : networkx multidigraph
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

    Returns
    -------
    route_map : folium.folium.Map
    """

    # check if we were able to import folium successfully
    if not folium:
        raise ImportError('The folium package must be installed to use this optional feature.')

    # create gdf of the route edges
    gdf_edges = graph_to_gdfs(G, nodes=False, fill_edge_geometry=True)
    route_nodes = list(zip(route[:-1], route[1:]))
    index = [gdf_edges[(gdf_edges['u']==u) & (gdf_edges['v']==v)].index[0] for u, v in route_nodes]
    gdf_route_edges = gdf_edges.loc[index]

    # get route centroid
    x, y = gdf_route_edges.unary_union.centroid.xy
    route_centroid = (y[0], x[0])

    # create the folium web map if one wasn't passed-in
    if route_map is None:
        route_map = folium.Map(location=route_centroid, zoom_start=zoom, tiles=tiles)

    # add each route edge to the map
    for _, row in gdf_route_edges.iterrows():
        pl = make_folium_polyline(edge=row, edge_color=route_color, edge_width=route_width,
                                  edge_opacity=route_opacity, popup_attribute=popup_attribute)
        pl.add_to(route_map)

    # if fit_bounds is True, fit the map to the bounds of the route by passing
    # list of lat-lng points as [southwest, northeast]
    if fit_bounds:
        tb = gdf_route_edges.total_bounds
        bounds = [(tb[1], tb[0]), (tb[3], tb[2])]
        route_map.fit_bounds(bounds)

    return route_map


def plot_figure_ground(G=None, address=None, point=None, dist=805,
                       network_type='drive_service', street_widths=None,
                       default_width=4, fig_length=8, edge_color='w',
                       bgcolor='#333333', smooth_joints=True, filename=None,
                       file_format='png', show=False, save=True, close=True,
                       dpi=300):
    """
    Plot a figure-ground diagram of a street network, defaulting to one square
    mile.

    Parameters
    ----------
    G : networkx multidigraph
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
        the height and width of this square diagram
    edge_color : string
        the color of the streets
    bgcolor : string
        the color of the background
    smooth_joints : bool
        if True, plot nodes same width as streets to smooth line joints and
        prevent cracks between them from showing
    filename : string
        filename to save the image as
    file_format : string
        the format of the file to save (e.g., 'jpg', 'png', 'svg')
    show : bool
        if True, show the figure
    save : bool
        if True, save the figure as an image file to disk
    close : bool
        close the figure (only if show equals False) to prevent display
    dpi : int
        the resolution of the image file if saving

    Returns
    -------
    fig, ax : tuple
    """

    multiplier = 1.2

    # if G was passed-in, use this graph in the plot, centered on the centroid
    # of its nodes
    if G is not None:
        gdf_nodes = graph_to_gdfs(G, edges=False, node_geometry=True)
        lnglat_point = gdf_nodes.unary_union.centroid.coords[0]
        point = tuple(reversed(lnglat_point))

    # otherwise, get the network by either address or point, whichever was
    # passed-in, using a distance multiplier to make sure we get more than
    # enough network. simplify in non-strict mode to not combine multiple street
    # types into single edge
    elif address is not None:
        G, point = graph_from_address(address, distance=dist*multiplier, distance_type='bbox', network_type=network_type,
                                      simplify=False, truncate_by_edge=True, return_coords=True)
        G = simplify_graph(G, strict=False)
    elif point is not None:
        G = graph_from_point(point, distance=dist*multiplier, distance_type='bbox', network_type=network_type,
                             simplify=False, truncate_by_edge=True)
        G = simplify_graph(G, strict=False)
    else:
        raise ValueError('You must pass an address or lat-long point or graph.')

    # if user did not pass in custom street widths, create a dict of default
    # values
    if street_widths is None:
        street_widths = {'footway' : 1.5,
                         'steps' : 1.5,
                         'pedestrian' : 1.5,
                         'service' : 1.5,
                         'path' : 1.5,
                         'track' : 1.5,
                         'motorway' : 6}

    # we need an undirected graph to find every edge incident to a node
    G_undir = G.to_undirected()

    # for each network edge, get a linewidth according to street type (the OSM
    # 'highway' value)
    edge_linewidths = []
    for _, _, data in G_undir.edges(keys=False, data=True):
        street_type = data['highway'][0] if isinstance(data['highway'], list) else data['highway']
        if street_type in street_widths:
            edge_linewidths.append(street_widths[street_type])
        else:
            edge_linewidths.append(default_width)

    if smooth_joints:
        # for each node, get a nodesize according to the narrowest incident edge
        node_widths = {}
        for node in G_undir.nodes():
            # first, identify all the highway types of this node's incident edges
            incident_edges_data = [G_undir.get_edge_data(node, neighbor) for neighbor in G_undir.neighbors(node)]
            edge_types = [data[0]['highway'] for data in incident_edges_data]
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
                edge_widths = [street_widths[edge_type] if edge_type in street_widths else default_width for edge_type in edge_types_flat]

                # the node diameter will be the biggest of the edge widths, to make joints perfectly smooth
                # alternatively, use min (?) to pervent anything larger from extending past smallest street's line
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
    bbox = bbox_from_point(point, dist, project_utm=False)

    # create a filename if one was not passed
    if filename is None and save:
        filename = 'figure_ground_{}_{}'.format(point, network_type)

    # plot the figure
    fig, ax = plot_graph(G_undir, bbox=bbox, fig_height=fig_length,
                         margin=0, axis_off=True, equal_aspect=False,
                         bgcolor=bgcolor, node_size=node_sizes,
                         node_color=edge_color, edge_linewidth=edge_linewidths,
                         edge_color=edge_color, show=show, save=save,
                         close=close, filename=filename, file_format=file_format,
                         dpi=dpi)

    return fig, ax
