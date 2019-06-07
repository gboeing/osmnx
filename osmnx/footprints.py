################################################################################
# Module: footprints.py
# Description: Download and plot footprints from OpenStreetMap
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

import geopandas as gpd
import matplotlib.pyplot as plt
import time
from descartes import PolygonPatch
from matplotlib.collections import PatchCollection
from shapely.geometry import LineString
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
from shapely.ops import polygonize

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


def osm_footprints_download(polygon=None, north=None, south=None, east=None, west=None,
                            footprint_type='building', timeout=180, memory=None,
                            max_query_area_size=50*1000*50*1000):
    """
    Download OpenStreetMap footprint data as a list of json responses.

    Parameters
    ----------
    polygon : shapely Polygon or MultiPolygon
        geographic shape to fetch the footprints within
    north : float
        northern latitude of bounding box
    south : float
        southern latitude of bounding box
    east : float
        eastern longitude of bounding box
    west : float
        western longitude of bounding box
    footprint_type : string
        type of footprint to be downloaded. OSM tag key e.g. 'building', 'landuse', 'place', etc.
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
        log('Requesting footprints data within bounding box from API in {:,} request(s)'.format(len(geometry)))
        start_time = time.time()

        # loop through each polygon rectangle in the geometry (there will only
        # be one if original bbox didn't exceed max area size)
        for poly in geometry:
            # represent bbox as south,west,north,east and round lat-longs to 8
            # decimal places (ie, within 1 mm) so URL strings aren't different
            # due to float rounding issues (for consistent caching)
            west, south, east, north = poly.bounds
            query_template = ('[out:json][timeout:{timeout}]{maxsize};'
                              '((way["{footprint_type}"]({south:.8f},{west:.8f},{north:.8f},{east:.8f});'
                              '(._;>;););'
                              '(relation["{footprint_type}"]({south:.8f},{west:.8f},{north:.8f},{east:.8f});'
                              '(._;>;);););out;')
            query_str = query_template.format(north=north, south=south, east=east, west=west, timeout=timeout,
                                              maxsize=maxsize, footprint_type=footprint_type)
            response_json = overpass_request(data={'data':query_str}, timeout=timeout)
            response_jsons.append(response_json)
        msg = ('Got all footprint data within bounding box from '
               'API in {:,} request(s) and {:,.2f} seconds')
        log(msg.format(len(geometry), time.time()-start_time))

    elif by_poly:
        # project to utm, divide polygon up into sub-polygons if area exceeds a
        # max size (in meters), project back to lat-long, then get a list of polygon(s) exterior coordinates
        geometry_proj, crs_proj = project_geometry(polygon)
        geometry_proj_consolidated_subdivided = consolidate_subdivide_geometry(geometry_proj, max_query_area_size=max_query_area_size)
        geometry, _ = project_geometry(geometry_proj_consolidated_subdivided, crs=crs_proj, to_latlong=True)
        polygon_coord_strs = get_polygons_coordinates(geometry)
        log('Requesting footprint data within polygon from API in {:,} request(s)'.format(len(polygon_coord_strs)))
        start_time = time.time()

        # pass each polygon exterior coordinates in the list to the API, one at
        # a time
        for polygon_coord_str in polygon_coord_strs:
            query_template = ('[out:json][timeout:{timeout}]{maxsize};('
                              'way(poly:"{polygon}")["{footprint_type}"];(._;>;);'
                              'relation(poly:"{polygon}")["{footprint_type}"];(._;>;););out;')
            query_str = query_template.format(polygon=polygon_coord_str, timeout=timeout, maxsize=maxsize,
                                              footprint_type=footprint_type)
            response_json = overpass_request(data={'data':query_str}, timeout=timeout)
            response_jsons.append(response_json)
        msg = ('Got all footprint data within polygon from API in '
               '{:,} request(s) and {:,.2f} seconds')
        log(msg.format(len(polygon_coord_strs), time.time()-start_time))

    return response_jsons


def create_footprints_gdf(polygon=None, north=None, south=None, east=None, west=None,
                          footprint_type='building', retain_invalid=False, responses=None):
    """
    Get footprint (polygon) data from OSM and convert it into a GeoDataFrame.

    Parameters
    ----------
    polygon : shapely Polygon or MultiPolygon
        geographic shape to fetch the footprints within
    north : float
        northern latitude of bounding box
    south : float
        southern latitude of bounding box
    east : float
        eastern longitude of bounding box
    west : float
        western longitude of bounding box
    footprint_type : string
        type of footprint to be downloaded. OSM tag key e.g. 'building', 'landuse', 'place', etc.
    retain_invalid : bool
        if False discard any footprints with an invalid geometry
    responses : list
        list of response jsons

    Returns
    -------
    GeoDataFrame
    """
    # allow pickling between downloading footprints and converting them to a GeoDataFrame
    if responses is None:
        responses = osm_footprints_download(polygon, north, south, east, west, footprint_type)

    # parse the list of responses into separate dicts of vertices, footprints and relations
    # create a set of ways not directly tagged with footprint_type
    vertices, footprints, relations, untagged_ways = responses_to_dicts(responses, footprint_type)

    # create simple Shapely geometries (Polygon or LineString) for all of the ways in footprints
    for footprint_key, footprint_val in footprints.items():
        footprint_val['geometry'] = create_footprint_geometry(footprint_key, footprint_val, vertices)

    # create a complex Shapely Polygon or MultiPolygon for each relation
    for relation_key, relation_val in relations.items():
        relation_val['geometry'] = create_relation_geometry(relation_key, relation_val, footprints)
    
    # merge relations into the footprints dictionary
    footprints.update(relations)

    # delete supporting geometry not directly tagged with footprint_type from the footprints dictionary
    for untagged_way in untagged_ways:
        try:
            del footprints[untagged_way]
        except KeyError:
            log('untagged_way {} not found in footprints dict'.format(untagged_way))

    # Convert footprints dictionary to a GeoDataFrame
    gdf = gpd.GeoDataFrame.from_dict(footprints, orient='index')
    gdf.crs = settings.default_crs

    # filter the gdf to only include valid Polygons or MultiPolygons
    if not retain_invalid:    
        filter1 = gdf['geometry'].is_valid
        filter2 = (gdf['geometry'].geom_type == 'Polygon') | (gdf['geometry'].geom_type == 'MultiPolygon')
        filter_combined = filter1 & filter2
        gdf = gdf[filter_combined]
    
    return gdf


def responses_to_dicts(responses, footprint_type):
    """
    Parse a list of json responses into dictionaries of vertices, footprints, and relations.

    Note: OSM's data model and the Overpass API will return open ways (lines) as part of
    a 'polygon' query. These may be fragments of the inner and outer rings of relations or
    they may be open ways mistakenly tagged with 'polygon' type tags.

    Ways not directly tagged with the footprint type are added to the untagged_ways set for
    removal from the footprints dictionary at the end of the process.

    Some inner ways of relations may be tagged with the footprint type in their own right e.g.
    landuse=meadow as an inner way in a landuse=forest relation and need to be kept. These are
    created here.

    Parameters
    ----------
    responses : list
        list of json responses
    footprint_type : string
        type of footprint downloaded. OSM tag key e.g. 'building', 'landuse', 'place', etc.

    Returns
    -------
    vertices
        dictionary of OSM nodes including their lat, lon coordinates
    footprints
        dictionary of OSM ways including their nodes and tags
    relations
        dictionary of OSM relations including member ids and tags
    untagged_footprints
        set of ids for ways or relations not directly tagged with footprint_type
    """
    # create dictionaries to hold vertices, footprints and relations
    vertices = {}
    footprints = {}
    relations = {}
    # create a set to hold the ids of ways not directly tagged as footprint_type
    untagged_footprints = set()

    # loop through each response once adding each element to one of the dicts
    for response in responses:
        for element in response['elements']:
            # NODES - only keep coordinates
            if 'type' in element and element['type']=='node':
                vertices[element['id']] = {'lat' : element['lat'],
                                           'lon' : element['lon']}
            # WAYS - both open and closed
            elif 'type' in element and element['type']=='way':
                footprint = {'nodes' : element['nodes']}
                if 'tags' in element:
                    for tag in element['tags']:
                        footprint[tag] = element['tags'][tag]
                footprints[element['id']] = footprint
                # add ways not individually tagged with footprint_type to the untagged_footprints set
                if ('tags' not in element) or (footprint_type not in element['tags']):
                    untagged_footprints.add(element['id'])
            # RELATIONS
            elif 'type' in element and element['type']=='relation':
                relation = {'members' : {}}
                for member in element['members']:
                    if 'type' in member and member['type']=='way':
                        relation['members'].update({member['ref']:member.get('role')})
                if 'tags' in element:
                    for tag in element['tags']:
                        relation[tag] = element['tags'][tag]
                relations[element['id']] = relation
                # add relations not individually tagged with footprint_type to the untagged_footprints set
                if ('tags' not in element) or (footprint_type not in element['tags']):
                    untagged_footprints.add(element['id'])
            # Log any other Elements found in the response
            else:
                log('Element {} is not a node, way or relation'.format(element['id']))

    return vertices, footprints, relations, untagged_footprints


def create_footprint_geometry(footprint_key, footprint_val, vertices):
    """
    Create Shapely geometry for open or closed ways in the initial footprints dictionary.

    Closed ways are converted directly to Shapely Polygons, open ways (fragments that will
    form the outer and inner rings of relations) are converted to LineStrings.

    Parameters
    ----------
    footprint_key : int
        the id of the way/footprint to process
    footprint_val : dict
        the nodes and tags of the footprint
    vertices : dict
        the dictionary of OSM nodes with their coordinates

    Returns
    -------
    Shapely Polygon or LineString
    """
    # CLOSED WAYS
    if footprint_val['nodes'][0] == footprint_val['nodes'][-1]:
        try:
            footprint_geometry = Polygon([(vertices[node]['lon'], vertices[node]['lat']) for node in footprint_val['nodes']])
        except Exception:
            log('Polygon has invalid geometry: {}'.format(footprint_key))
    # OPEN WAYS    
    else:
        try:
            footprint_geometry = LineString([(vertices[node]['lon'], vertices[node]['lat']) for node in footprint_val['nodes']])
        except Exception:
            log('LineString has invalid geometry: {}'.format(footprint_key))

    return footprint_geometry


def create_relation_geometry(relation_key, relation_val, footprints):
    """
    Create Shapely geometry for relations - Polygons with holes or MultiPolygons

    OSM relations are used to define complex polygons - polygons with holes or
    multi-polygons. The polygons' outer and inner rings may be made up of chains
    of LineStrings. https://wiki.openstreetmap.org/wiki/Relation:multipolygon 
    requires that multipolygon rings have an outer or inner 'role'.
    
    OSM's data model allows a polygon type tag e.g. 'building' to be added to 
    any OSM element. This can include non-polygon relations e.g. bus routes.
    Relations that do not have at least one closed ring with an outer role 
    are filtered out.

    Inner rings that are tagged with the footprint type in their own right e.g.
    landuse=meadow as an inner ring of landuse=forest will have been included in
    the footprints dictionary as part of the original parsing and are not dealt
    with here.

    Parameters
    ----------
    relation_key : int
        the id of the relation to process
    relation_val : dict
        members and tags of the relation
    footprints : dictionary
        dictionary of all footprints (including open and closed ways)

    Returns
    -------
    Shapely Polygon or MultiPolygon
    """

    # create empty lists to hold member geometries
    multipoly = []
    outer_polys = []
    outer_lines = []
    inner_polys = []
    inner_lines = []

    # add each members geometry to a list according to its role and geometry type
    for member_id, member_role in relation_val['members'].items():
        if member_role == 'outer':
            if footprints[member_id]['geometry'].geom_type == 'Polygon':
                outer_polys.append(footprints[member_id]['geometry'])
            elif footprints[member_id]['geometry'].geom_type == 'LineString':
                outer_lines.append(footprints[member_id]['geometry'])
        elif member_role == 'inner':
            if footprints[member_id]['geometry'].geom_type == 'Polygon':
                inner_polys.append(footprints[member_id]['geometry'])
            elif footprints[member_id]['geometry'].geom_type == 'LineString':
                inner_lines.append(footprints[member_id]['geometry'])

    # try to polygonize open outer ways and concatenate them to outer_polys
    if len(outer_lines) > 0:
        try:
            result = list(polygonize(outer_lines))
        except Exception:
            log("polygonize failed for 'outer' ways in relation: {}".format(relation_key))
        else:
            outer_polys += result

    # try to polygonize open inner ways and concatenate them to inner_polys
    if len(inner_lines) > 0:
        try:
            result = list(polygonize(inner_lines))
        except Exception:
            log("polygonize failed for 'inner' ways in relation: {}".format(relation_key))
        else:
            inner_polys += result

    # filter out relations missing both 'outer' and 'inner' polygons or just 'outer'
    if len(outer_polys + inner_polys) == 0:
        log("Relation {} missing 'outer' and 'inner' closed ways".format(relation_key))
    elif len(outer_polys) == 0:
        log("Relation {} missing 'outer' closed ways".format(relation_key))
    # process the others to multipolygons
    else:
        for outer_poly in outer_polys:
            temp_poly = outer_poly
            for inner_poly in inner_polys:
                if inner_poly.within(outer_poly):
                    temp_poly=temp_poly.difference(inner_poly)
            multipoly.append(temp_poly)

    # return relations with one outer way as Polygons, multiple outer ways as MultiPolygons
    if len(multipoly) == 1:
        return Polygon(multipoly[0])
    elif len(multipoly) > 1:    
        return MultiPolygon(multipoly)
    else:
        log('relation {} could not be converted to a complex footprint'.format(relation_key))


def footprints_from_point(point, distance, footprint_type='building', retain_invalid=False):
    """
    Get footprints within some distance north, south, east, and west of
    a lat-long point.

    Parameters
    ----------
    point : tuple
        a lat-long point
    distance : numeric
        distance in meters
    footprint_type : string
        type of footprint to be downloaded. OSM tag key e.g. 'building', 'landuse', 'place', etc.
    retain_invalid : bool
        if False discard any footprints with an invalid geometry

    Returns
    -------
    GeoDataFrame
    """

    bbox = bbox_from_point(point=point, distance=distance)
    north, south, east, west = bbox
    return create_footprints_gdf(north=north, south=south, east=east, west=west,
                                 footprint_type=footprint_type, retain_invalid=retain_invalid)


def footprints_from_address(address, distance, footprint_type='building', retain_invalid=False):
    """
    Get footprints within some distance north, south, east, and west of
    an address.

    Parameters
    ----------
    address : string
        the address to geocode to a lat-long point
    distance : numeric
        distance in meters
    footprint_type : string
        type of footprint to be downloaded. OSM tag key e.g. 'building', 'landuse', 'place', etc.
    retain_invalid : bool
        if False discard any footprints with an invalid geometry

    Returns
    -------
    GeoDataFrame
    """

    # geocode the address string to a (lat, lon) point
    point = geocode(query=address)

    # get footprints within distance of this point
    return footprints_from_point(point, distance, footprint_type=footprint_type,
                                 retain_invalid=retain_invalid)


def footprints_from_polygon(polygon, footprint_type='building', retain_invalid=False):
    """
    Get footprints within some polygon.

    Parameters
    ----------
    polygon : shapely Polygon or MultiPolygon
        the shape to get data within. coordinates should be in units of
        latitude-longitude degrees.
    footprint_type : string
        type of footprint to be downloaded. OSM tag key e.g. 'building', 'landuse', 'place', etc.
    retain_invalid : bool
        if False discard any footprints with an invalid geometry

    Returns
    -------
    GeoDataFrame
    """

    return create_footprints_gdf(polygon=polygon, footprint_type=footprint_type,
                                 retain_invalid=retain_invalid)


def footprints_from_place(place, footprint_type='building', retain_invalid=False, which_result=1):
    """
    Get footprints within the boundaries of some place.

    The query must be geocodable and OSM must have polygon boundaries for the
    geocode result. If OSM does not have a polygon for this place, you can
    instead get its footprints using the footprints_from_address function, which
    geocodes the place name to a point and gets the footprints within some distance
    of that point.

    Parameters
    ----------
    place : string
        the query to geocode to get geojson boundary polygon
    footprint_type : string
        type of footprint to be downloaded. OSM tag key e.g. 'building', 'landuse', 'place', etc.
    retain_invalid : bool
        if False discard any footprints with an invalid geometry
    which_result : int
        max number of results to return and which to process upon receipt

    Returns
    -------
    GeoDataFrame
    """

    city = gdf_from_place(place, which_result=which_result)
    polygon = city['geometry'].iloc[0]
    return create_footprints_gdf(polygon, retain_invalid=retain_invalid,
                                 footprint_type=footprint_type)


def plot_footprints(gdf, fig=None, ax=None, figsize=None, color='#333333', bgcolor='w',
                   set_bounds=True, bbox=None, save=False, show=True, close=False,
                   filename='image', file_format='png', dpi=600):
    """
    Plot a GeoDataFrame of footprints.

    Parameters
    ----------
    gdf : GeoDataFrame
        footprints
    fig : figure
    ax : axis
    figsize : tuple
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
