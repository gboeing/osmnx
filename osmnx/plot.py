###################################################################################################
# Module: plot.py
# Description: Plot spatial geometries, street networks, and routes
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
###################################################################################################

import time
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
from descartes import PolygonPatch
from shapely.geometry import Polygon, MultiPolygon

from . import globals
from .utils import log


def plot_shape(gdf, fc='#cbe0f0', ec='#999999', linewidth=1, alpha=1, figsize=(6,6), margin=0.02, axis_off=True):
    """
    Plot a GeoDataFrame of place boundary geometries.
    
    Parameters
    ----------
    gdf : GeoDataFrame, the gdf containing the geometries to plot
    fc : string, the facecolor for the polygons
    ec : string, the edgecolor for the polygons
    linewidth : numeric, the width of the polygon edge lines
    alpha : numeric, the opacity
    figsize : tuple, the size of the plotting figure
    margin : numeric, the size of the figure margins
    axis_off : bool, if True, disable the matplotlib axes display
    
    Returns
    -------
    fig, ax : tuple
    """
    # plot the geometries one at a time
    fig, ax = plt.subplots(figsize=figsize)
    for geometry in gdf['geometry'].tolist():
        if isinstance(geometry, (Polygon, MultiPolygon)):
            if isinstance(geometry, Polygon):
                geometry = MultiPolygon([geometry])
            for polygon in geometry:
                patch = PolygonPatch(polygon, fc=fc, ec=ec, linewidth=linewidth, alpha=alpha)
                ax.add_patch(patch)
        else:
            raise ValueError('All geometries in GeoDataFrame must be shapely Polygons or MultiPolygons')

    # adjust the axis margins and limits around the image and make axes equal-aspect
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


def get_edge_colors_by_attr(G, attr, num_bins=5, cmap='spectral', start=0.1, stop=0.9):
    """
    Get a list of edge colors by binning some continuous-variable attribute into quantiles
    
    Parameters
    ----------
    G : graph
    attr : string, the name of the continuous-variable attribute
    num_bins : int, how many quantiles
    cmap : string, name of a colormap
    start : float, where to start in the colorspace
    stop : float, where to end in the colorspace
    
    Returns
    -------
    colors : list
    """
    bin_labels = range(num_bins)
    attr_values = pd.Series([data[attr] for u, v, key, data in G.edges(keys=True, data=True)])
    cats = pd.qcut(x=attr_values, q=num_bins, labels=bin_labels)
    color_list = [cm.get_cmap(cmap)(x) for x in np.linspace(start, stop, num_bins)]
    colors = [color_list[cat] for cat in cats]
    return colors
    
    
def save_and_show(fig, ax, save, show, close, filename, file_format, dpi, axis_off):
    """
    Save a figure to disk and show it, as specified.
    
    Parameters
    ----------
    fig : figure
    ax : axis
    save : bool, whether to save the figure to disk or not
    show : bool, whether to display the figure or not
    close : close the figure (only if show equals False) to prevent display
    filename : string, the name of the file to save
    file_format : string, the format of the file to save (e.g., 'jpg', 'png', 'svg')
    dpi : int, the resolution of the image file if saving
    axis_off : bool, if True matplotlib axis was turned off by plot_graph so constrain the saved figure's extent to the interior of the axis
    
    Returns
    -------
    fig, ax : figure, axis
    """
    # save the figure if specified
    if save:
        start_time = time.time()
        
        # create the save folder if it doesn't already exist
        if not os.path.exists(globals.imgs_folder):
            os.makedirs(globals.imgs_folder)
        path_filename = '{}/{}.{}'.format(globals.imgs_folder, filename, file_format)
        
        if file_format == 'svg':
            # if the file_format is svg, prep the fig/ax a bit for saving
            ax.axis('off')
            ax.set_position([0, 0, 1, 1])
            ax.patch.set_alpha(0.)
            fig.patch.set_alpha(0.)
            fig.savefig(path_filename, bbox_inches=0, format=file_format, facecolor=fig.get_facecolor(), transparent=True)
        else:
            if axis_off:
                # if axis is turned off, constrain the saved figure's extent to the interior of the axis
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
    
    
def plot_graph(G, bbox=None, fig_height=6, fig_width=None, margin=0.02, axis_off=True, bgcolor='w',
               show=True, save=False, close=True, file_format='jpg', filename='temp', dpi=300, annotate=False,
               node_color='#66ccff', node_size=15, node_alpha=1, node_edgecolor='none', node_zorder=1,
               edge_color='#999999', edge_linewidth=1, edge_alpha=1, use_geom=True):
    """
    Plot a networkx spatial graph.
    
    Parameters
    ----------
    G : graph
    bbox : tuple, bounding box as north,south,east,west - if None will calculate from spatial extents of data
    fig_height : int, matplotlib figure height in inches
    fig_width : int, matplotlib figure width in inches
    margin : float, relative margin around the figure
    axis_off : bool, if True turn off the matplotlib axis
    bgcolor : string, the background color of the figure and axis
    show : bool, if True, show the figure
    save : bool, if True, save the figure as an image file to disk
    close : close the figure (only if show equals False) to prevent display
    file_format : string, the format of the file to save (e.g., 'jpg', 'png', 'svg')
    filename : string, the name of the file if saving
    dpi : int, the resolution of the image file if saving
    annotate : bool, if True, annotate the nodes in the figure
    node_color : string, the color of the nodes
    node_size : int, the size of the nodes
    node_alpha : float, the opacity of the nodes
    node_edgecolor : string, the color of the node's marker's border
    node_zorder : int, zorder to plot nodes, edges are always 2, so make node_zorder 1 to plot nodes beneath them or 3 to plot nodes atop them
    edge_color : string, the color of the edges' lines
    edge_linewidth : float, the width of the edges' lines
    edge_alpha : float, the opacity of the edges' lines
    use_geom : bool, if True, use the spatial geometry attribute of the edges to draw geographically accurate edges, rather than just lines straight from node to node
    
    Returns
    -------
    fig, ax : figure, axis    
    """
    
    log('Begin plotting the graph...')
    node_Xs = [float(node['x']) for node in G.node.values()]
    node_Ys = [float(node['y']) for node in G.node.values()]
    
    # get north, south, east, west values either from bbox parameter or from min/max node coordinate values
    if bbox is None:
        north = max(node_Ys)
        south = min(node_Ys)
        east = max(node_Xs)
        west = min(node_Xs)
    else:
        north, south, east, west = bbox
    
    # if caller did not pass in a fig_width, calculate it proportionately from the fig_height and bounding box aspect ratio
    bbox_aspect_ratio = (north-south)/(east-west)
    if fig_width is None:
        fig_width = fig_height / bbox_aspect_ratio
    
    # create the figure and axis
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), facecolor=bgcolor)
    ax.set_axis_bgcolor(bgcolor)
    
    # draw the edges as lines from node to node
    start_time = time.time()
    lines = []
    for u, v, key, data in G.edges(keys=True, data=True):
        if 'geometry' in data and use_geom:
            # if it has a geometry attribute (a list of line segments), add them to the list of lines to plot
            xs, ys = data['geometry'].xy
            lines.append(list(zip(xs, ys)))
        else:
            # if it doesn't have a geometry attribute, the edge is a straight line from node to node
            x1 = G.node[u]['x']
            y1 = G.node[u]['y']
            x2 = G.node[v]['x']
            y2 = G.node[v]['y']
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
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    if axis_off:
        ax.axis('off')
    
    # annotate the axis with node IDs if annotate=True
    if annotate:
        for node, data in G.nodes(data=True):
            ax.annotate(node, xy=(data['x'], data['y']))
            
    # save and show the figure as specified
    fig, ax = save_and_show(fig, ax, save, show, close, filename, file_format, dpi, axis_off)
    return fig, ax


def plot_graph_route(G, route, bbox=None, fig_height=6, fig_width=None, margin=0.02, bgcolor='w',
                     axis_off=True, show=True, save=False, close=True, file_format='jpg', filename='temp', dpi=300, annotate=False,
                     node_color='#999999', node_size=15, node_alpha=1, node_edgecolor='none', node_zorder=1,
                     edge_color='#999999', edge_linewidth=1, edge_alpha=1, use_geom=True,
                     origin_point=None, destination_point=None,
                     route_color='r', route_linewidth=4, route_alpha=0.5, orig_dest_node_alpha=0.5,
                     orig_dest_node_size=100, orig_dest_node_color='r', orig_dest_point_color='b'):
    """
    Plot a route along a networkx spatial graph.
    
    Parameters
    ----------
    G : graph
    route : list, the route as a list of nodes
    bbox : tuple, bounding box as north,south,east,west - if None will calculate from spatial extents of data
    fig_height : int, matplotlib figure height in inches
    fig_width : int, matplotlib figure width in inches
    margin : float, relative margin around the figure
    axis_off : bool, if True turn off the matplotlib axis
    bgcolor : string, the background color of the figure and axis
    show : bool, if True, show the figure
    save : bool, if True, save the figure as an image file to disk
    close : close the figure (only if show equals False) to prevent display
    file_format : string, the format of the file to save (e.g., 'jpg', 'png', 'svg')
    filename : string, the name of the file if saving
    dpi : int, the resolution of the image file if saving
    annotate : bool, if True, annotate the nodes in the figure
    node_color : string, the color of the nodes
    node_size : int, the size of the nodes
    node_alpha : float, the opacity of the nodes
    node_edgecolor : string, the color of the node's marker's border
    node_zorder : int, zorder to plot nodes, edges are always 2, so make node_zorder 1 to plot nodes beneath them or 3 to plot nodes atop them
    edge_color : string, the color of the edges' lines
    edge_linewidth : float, the width of the edges' lines
    edge_alpha : float, the opacity of the edges' lines
    use_geom : bool, if True, use the spatial geometry attribute of the edges to draw geographically accurate edges, rather than just lines straight from node to node
    origin_point : tuple, optional, an origin (lat, lon) point to plot instead of the origin node
    destination_point : tuple, optional, a destination (lat, lon) point to plot instead of the destination node
    route_color : string, the color of the route
    route_linewidth : int, the width of the route line
    route_alpha : float, the opacity of the route line
    orig_dest_node_alpha : float, the opacity of the origin and destination nodes
    orig_dest_node_size : int, the size of the origin and destination nodes
    orig_dest_node_color : string, the color of the origin and destination nodes
    orig_dest_point_color : string, the color of the origin and destination points if being plotted instead of nodes
    
    Returns
    -------
    fig, ax : figure, axis    
    """
    
    # plot the graph but not the route
    fig, ax = plot_graph(G, bbox=bbox, fig_height=fig_height, fig_width=fig_width, margin=margin, axis_off=axis_off, bgcolor=bgcolor,
                         show=False, save=False, close=False, filename=filename, dpi=dpi, annotate=annotate,
                         node_color=node_color, node_size=node_size, node_alpha=node_alpha, node_edgecolor=node_edgecolor, node_zorder=node_zorder,
                         edge_color=edge_color, edge_linewidth=edge_linewidth, edge_alpha=edge_alpha, use_geom=use_geom)
    
    # the origin and destination nodes are the first and last nodes in the route
    origin_node = route[0]
    destination_node = route[-1]
        
    if origin_point is None or destination_point is None:
        # if caller didn't pass points, use the first and last node in route as origin/destination    
        origin_destination_lats = (G.node[origin_node]['y'], G.node[destination_node]['y'])
        origin_destination_lons = (G.node[origin_node]['x'], G.node[destination_node]['x'])
    else:
        # otherwise, use the passed points as origin/destination
        origin_destination_lats = (origin_point[0], destination_point[0])
        origin_destination_lons = (origin_point[1], destination_point[1])
        orig_dest_node_color = orig_dest_point_color
    
    # scatter the origin and destination points
    ax.scatter(origin_destination_lons, origin_destination_lats, s=orig_dest_node_size, 
               c=orig_dest_node_color, alpha=orig_dest_node_alpha, edgecolor=node_edgecolor, zorder=4)
    
    # plot the route lines
    edge_nodes = list(zip(route[:-1], route[1:]))
    lines = []
    for u, v in edge_nodes:
        # if there are parallel edges, select the shortest in length
        data = min([data for data in G.edge[u][v].values()], key=lambda x: x['length'])
        
        # if it has a geometry attribute (ie, a list of line segments)
        if 'geometry' in data and use_geom:
            # add them to the list of lines to plot
            xs, ys = data['geometry'].xy
            lines.append(list(zip(xs, ys)))
        else:
            # if it doesn't have a geometry attribute, the edge is a straight line from node to node
            x1 = G.node[u]['x']
            y1 = G.node[u]['y']
            x2 = G.node[v]['x']
            y2 = G.node[v]['y']
            line = [(x1, y1), (x2, y2)]
            lines.append(line)
    
    # add the lines to the axis as a linecollection    
    lc = LineCollection(lines, colors=route_color, linewidths=route_linewidth, alpha=route_alpha, zorder=3)
    ax.add_collection(lc)
    
    # save and show the figure as specified
    fig, ax = save_and_show(fig, ax, save, show, close, filename, file_format, dpi, axis_off)
    return fig, ax
    
    