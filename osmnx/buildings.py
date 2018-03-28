################################################################################
# Module: buildings.py
# Description: Download and plot building footprints from OpenStreetMap
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

import time
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
from descartes import PolygonPatch

from . import settings
from .core import consolidate_subdivide_geometry
from .core import get_polygons_coordinates
from .core import overpass_request
from .core import bbox_from_point
from .core import gdf_from_place
from .plot import save_and_show
from .projection import project_geometry
from .utils import log
from .utils import geocode


def osm_bldg_download(polygon=None, north=None, south=None, east=None, west=None,
                      timeout=180, memory=None, max_query_area_size=50*1000*50*1000):
    """
    Download OpenStreetMap building footprint data.

    Parameters
    ----------
    polygon : shapely Polygon or MultiPolygon
        geographic shape to fetch the building footprints within
    north : float
        northern latitude of bounding box
    south : float
        southern latitude of bounding box
    east : float
        eastern longitude of bounding box
    west : float
        western longitude of bounding box
    timeout : int
        the timeout interval for requests and to pass to API
    memory : int
        server memory allocation size for the query, in bytes. If none, server
        will use its default allocation size
    max_query_area_size : float
        max area for any part of the geometry, in the units the geometry is in:
        any polygon bigger will get divided up for multiple queries to API
        (default is 50,000 * 50,000 units (ie, 50km x 50km in area, if units are
        meters))

    Returns
    -------
    list
        list of response_json dicts
    """

    # check if we're querying by polygon or by bounding box based on which
    # argument(s) where passed into this function
    by_poly = polygon is not None
    by_bbox = not (north is None or south is None or east is None or west is None)
    if not (by_poly or by_bbox):
        raise ValueError('You must pass a polygon or north, south, east, and west')

    response_jsons = []

    # pass server memory allocation in bytes for the query to the API
    # if None, pass nothing so the server will use its default allocation size
    # otherwise, define the query's maxsize parameter value as whatever the
    # caller passed in
    if memory is None:
        maxsize = ''
    else:
        maxsize = '[maxsize:{}]'.format(memory)

    # define the query to send the API
    if by_bbox:
        # turn bbox into a polygon and project to local UTM
        polygon = Polygon([(west, south), (east, south), (east, north), (west, north)])
        geometry_proj, crs_proj = project_geometry(polygon)

        # subdivide it if it exceeds the max area size (in meters), then project
        # back to lat-long
        geometry_proj_consolidated_subdivided = consolidate_subdivide_geometry(geometry_proj, max_query_area_size=max_query_area_size)
        geometry, _ = project_geometry(geometry_proj_consolidated_subdivided, crs=crs_proj, to_latlong=True)
        log('Requesting building footprints data within bounding box from API in {:,} request(s)'.format(len(geometry)))
        start_time = time.time()

        # loop through each polygon rectangle in the geometry (there will only
        # be one if original bbox didn't exceed max area size)
        for poly in geometry:
            # represent bbox as south,west,north,east and round lat-longs to 8
            # decimal places (ie, within 1 mm) so URL strings aren't different
            # due to float rounding issues (for consistent caching)
            west, south, east, north = poly.bounds
            query_template = ('[out:json][timeout:{timeout}]{maxsize};((way["building"]({south:.8f},'
                              '{west:.8f},{north:.8f},{east:.8f});(._;>;););(relation["building"]'
                              '({south:.8f},{west:.8f},{north:.8f},{east:.8f});(._;>;);));out;')
            query_str = query_template.format(north=north, south=south, east=east, west=west, timeout=timeout, maxsize=maxsize)
            response_json = overpass_request(data={'data':query_str}, timeout=timeout)
            response_jsons.append(response_json)
        msg = ('Got all building footprints data within bounding box from '
               'API in {:,} request(s) and {:,.2f} seconds')
        log(msg.format(len(geometry), time.time()-start_time))

    elif by_poly:
        # project to utm, divide polygon up into sub-polygons if area exceeds a
        # max size (in meters), project back to lat-long, then get a list of polygon(s) exterior coordinates
        geometry_proj, crs_proj = project_geometry(polygon)
        geometry_proj_consolidated_subdivided = consolidate_subdivide_geometry(geometry_proj, max_query_area_size=max_query_area_size)
        geometry, _ = project_geometry(geometry_proj_consolidated_subdivided, crs=crs_proj, to_latlong=True)
        polygon_coord_strs = get_polygons_coordinates(geometry)
        log('Requesting building footprints data within polygon from API in {:,} request(s)'.format(len(polygon_coord_strs)))
        start_time = time.time()

        # pass each polygon exterior coordinates in the list to the API, one at
        # a time
        for polygon_coord_str in polygon_coord_strs:
            query_template = ('[out:json][timeout:{timeout}]{maxsize};(way'
                              '(poly:"{polygon}")["building"];(._;>;);relation'
                              '(poly:"{polygon}")["building"];(._;>;));out;')
            query_str = query_template.format(polygon=polygon_coord_str, timeout=timeout, maxsize=maxsize)
            response_json = overpass_request(data={'data':query_str}, timeout=timeout)
            response_jsons.append(response_json)
        msg = ('Got all building footprints data within polygon from API in '
               '{:,} request(s) and {:,.2f} seconds')
        log(msg.format(len(polygon_coord_strs), time.time()-start_time))

    return response_jsons


def create_buildings_gdf(polygon=None, north=None, south=None, east=None,
                         west=None, retain_invalid=False):
    """
    Get building footprint data from OSM then assemble it into a GeoDataFrame.

    Parameters
    ----------
    polygon : shapely Polygon or MultiPolygon
        geographic shape to fetch the building footprints within
    north : float
        northern latitude of bounding box
    south : float
        southern latitude of bounding box
    east : float
        eastern longitude of bounding box
    west : float
        western longitude of bounding box
    retain_invalid : bool
        if False discard any building footprints with an invalid geometry

    Returns
    -------
    GeoDataFrame
    """

    responses = osm_bldg_download(polygon, north, south, east, west)

    vertices = {}
    for response in responses:
        for result in response['elements']:
            if 'type' in result and result['type']=='node':
                vertices[result['id']] = {'lat' : result['lat'],
                                          'lon' : result['lon']}

    buildings = {}
    for response in responses:
        for result in response['elements']:
            if 'type' in result and result['type']=='way':
                nodes = result['nodes']
                try:
                    polygon = Polygon([(vertices[node]['lon'], vertices[node]['lat']) for node in nodes])
                except Exception:
                    log('Polygon has invalid geometry: {}'.format(nodes))
                building = {'nodes' : nodes,
                            'geometry' : polygon}

                if 'tags' in result:
                    for tag in result['tags']:
                        building[tag] = result['tags'][tag]

                buildings[result['id']] = building

    gdf = gpd.GeoDataFrame(buildings).T
    gdf.crs = settings.default_crs

    if not retain_invalid:
        # drop all invalid geometries
        gdf = gdf[gdf['geometry'].is_valid]

    return gdf


def buildings_from_point(point, distance, retain_invalid=False):
    """
    Get building footprints within some distance north, south, east, and west of
    a lat-long point.

    Parameters
    ----------
    point : tuple
        a lat-long point
    distance : numeric
        distance in meters
    retain_invalid : bool
        if False discard any building footprints with an invalid geometry

    Returns
    -------
    GeoDataFrame
    """

    bbox = bbox_from_point(point=point, distance=distance)
    north, south, east, west = bbox
    return create_buildings_gdf(north=north, south=south, east=east, west=west, retain_invalid=retain_invalid)


def buildings_from_address(address, distance, retain_invalid=False):
    """
    Get building footprints within some distance north, south, east, and west of
    an address.

    Parameters
    ----------
    address : string
        the address to geocode to a lat-long point
    distance : numeric
        distance in meters
    retain_invalid : bool
        if False discard any building footprints with an invalid geometry

    Returns
    -------
    GeoDataFrame
    """

    # geocode the address string to a (lat, lon) point
    point = geocode(query=address)

    # get buildings within distance of this point
    return buildings_from_point(point, distance, retain_invalid=retain_invalid)


def buildings_from_polygon(polygon, retain_invalid=False):
    """
    Get building footprints within some polygon.

    Parameters
    ----------
    polygon : Polygon

    retain_invalid : bool
        if False discard any building footprints with an invalid geometry

    Returns
    -------
    GeoDataFrame
    """

    return create_buildings_gdf(polygon=polygon, retain_invalid=retain_invalid)


def buildings_from_place(place, retain_invalid=False):
    """
    Get building footprints within the boundaries of some place.

    Parameters
    ----------
    place : string
        the query to geocode to get geojson boundary polygon
    retain_invalid : bool
        if False discard any building footprints with an invalid geometry

    Returns
    -------
    GeoDataFrame
    """

    city = gdf_from_place(place)
    polygon = city['geometry'].iloc[0]
    return create_buildings_gdf(polygon, retain_invalid=retain_invalid)


def plot_buildings(gdf, fig=None, ax=None, figsize=None, color='#333333', bgcolor='w', set_bounds=True, bbox=None,
                   save=False, show=True, close=False, filename='image', file_format='png', dpi=600):
    """
    Plot a GeoDataFrame of building footprints.

    Parameters
    ----------
    gdf : GeoDataFrame
        building footprints
    fig : figure
    ax : axis
    figsize : tuple
    color : string
        the color of the building footprints
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
        the name of the file to save
    file_format : string
        the format of the file to save (e.g., 'jpg', 'png', 'svg')
    dpi : int
        the resolution of the image file if saving

    Returns
    -------
    fig, ax : tuple

    """

    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=figsize, facecolor=bgcolor)
        ax.set_facecolor(bgcolor)

    # extract each polygon as a descartes patch, and add to a matplotlib patch
    # collection
    patches = []
    for geometry in gdf['geometry']:
        if isinstance(geometry, Polygon):
            patches.append(PolygonPatch(geometry))
        elif isinstance(geometry, MultiPolygon):
            for subpolygon in geometry: #if geometry is multipolygon, go through each constituent subpolygon
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
    ax.axis('off')
    ax.margins(0)
    ax.tick_params(which='both', direction='in')
    fig.canvas.draw()

    # make everything square
    ax.set_aspect('equal')
    fig.canvas.draw()

    fig, ax = save_and_show(fig=fig, ax=ax, save=save, show=show, close=close,
                            filename=filename, file_format=file_format, dpi=dpi, axis_off=True)

    return fig, ax
