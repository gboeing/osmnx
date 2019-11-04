################################################################################
# Module: elevation.py
# Description: Get node elevations and edge grades from the Google Maps
#              Elevation API
# License: MIT, see full license in LICENSE.txt
# Web: https://github.com/gboeing/osmnx
################################################################################

import math
import time

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
import rasterio
import requests
from geopandas.tools import sjoin
from scipy.interpolate import interp2d
from shapely.geometry import box

from .core import get_from_cache, save_to_cache
from .utils import log


def add_node_elevations_from_dem(G,
                                 dem_paths,
                                 num_buffer=1,
                                 interp_kind='linear'):
    """
    Add elevation for each node in network using a Digital Elevation Model (DEM)
    and add to node as an attribute. Measurement unit is unchanged from the
    source DEM; most DEM files are in meters.

    DEM sources include SRTM, which has worldwide 1 arc-second (30m) coverage,
    and USGS 3DEP, which has 1/3 arc-second (10m) and 1 arc-second (30m)
    coverage for the continental United States.

    Parameters
    ----------
    G : networkx multidigraph
    dem_paths : str, pathlib.Path or list of str or pathlib.Path
        path to one or more DEM files with elevation data. Can be any file
        format readable by GDAL.
    num_buffer : int
        number of bordering cells around (lon, lat) to use when interpolating
    interp_kind : str
        kind of interpolation. Passed to scipy.interpolate.interp2d. Can be
        ['linear', 'cubic', 'quintic']. Note that 'cubic' requires 'num_buffer'
        of at least 3 and 'quintic' requires 'num_buffer' of at least 5.

    Returns
    -------
    G : networkx multidigraph
    """

    if type(dem_paths) == str:
        dem_paths = [dem_paths]

    # make a DataFrame of all the nodes' coordinates, with node id as index
    node_points = pd.DataFrame([(node, data['x'], data['y'])
                                for node, data in G.nodes(data=True)],
                               columns=['node', 'lon', 'lat'])
    node_points = node_points.set_index('node')

    # Get bounding box of each DEM file in dem_paths
    extents = get_dem_extents(dem_paths)

    # For each node, find the files within 6 arc-seconds (this should be
    # suitable for most elevation DEMs, since they are generally a resolution of
    # max 1 or 3 arc-second.)
    # This is additionally scaled by `num_buffer`, to find enough boxes around
    # the point
    # So first create a box of this distance around each node, to then intersect
    # with the DEM extents
    # TODO use datset.res to make this exact
    # For example: (9.259259252584599e-05, 9.259259252584599e-05) for 1/3 arcsec
    dist = 1 / 60 / 60 * 6 * num_buffer
    node_boxes = gpd.GeoDataFrame(
        node_points,
        geometry=node_points.apply(lambda r: box(r['lon'] - dist, r[
            'lat'] - dist, r['lon'] + dist, r['lat'] + dist),
                                   axis=1))

    # Now find the DEMs that intersect each point's bounding box (according to
    # `dist`)
    intersection = sjoin(node_boxes, extents, how='left')

    # After the spatial join, there is no longer only one row per index value
    # (node key). Nodes that are near the edge of two DEMs may have two rows. So
    # create a grouped DataFrame and then iterate over the groups, passing each
    # to `query_dem`
    grouped = intersection.groupby(intersection.index)
    results = []
    for name, group in grouped:
        # Within each group, all lon/lat values should be the same
        msg = 'lat/lon values changing within group'
        assert all(group[['lon', 'lat']].nunique() == 1), msg

        value = query_dem(dem_paths=group['path'].values,
                          lon=group['lon'].iloc[0],
                          lat=group['lat'].iloc[0])
        results.append([name, value])

    results_df = pd.DataFrame(results, columns=['node',
                                                'elevation']).set_index('node')

    # add elevation as an attribute to the nodes
    node_points = node_points.join(results_df, how='left')
    msg = 'missing elevation for some points'
    assert sum(node_points['elevation'].isna()) == 0, msg

    # round elevation to millimeter
    node_points['elevation'] = node_points['elevation'].round(3)
    nx.set_node_attributes(G,
                           name='elevation',
                           values=node_points['elevation'].to_dict())
    log('Added elevation data to all nodes.')

    return G


def get_dem_extents(dem_paths):
    """Find bounding boxes for each DEM file

    Parameters
    ----------
    dem_paths : str, pathlib.Path or list of str or pathlib.Path
        path to one or more DEM files with elevation data. Can be any file
        format readable by GDAL.

    Returns
    -------
    extents : GeoDataFrame with DEM file path and its bounding box as geometry
    """
    # Create dict with bounding boxes of each file
    extents = []

    for dem_path in dem_paths:
        with rasterio.open(dem_path) as dataset:
            extents.append([dem_path, box(*dataset.bounds)])

    extents = gpd.GeoDataFrame(extents,
                               columns=['path', 'geometry'],
                               geometry='geometry')
    return extents


def query_dem(dem_paths, lon, lat, num_buffer=1, interp_kind='linear'):
    """Query elevation data for given point

    Parameters
    ----------
    dem_paths : str, pathlib.Path or list of str or pathlib.Path
        path to one or more DEM files with elevation data. Can be any file
        format readable by GDAL.
    lon : float
        longitude
    lat : float
        latitude
    num_buffer : int
        number of bordering cells around (lon, lat) to use when interpolating
    interp_kind : str
        kind of interpolation. Passed to scipy.interpolate.interp2d. Can be
        ['linear', 'cubic', 'quintic']. Note that 'cubic' requires 'num_buffer'
        of at least 3 and 'quintic' requires 'num_buffer' of at least 5.

    Returns
    -------
    value : float
        elevation in terms of the unit of the DEM (usually meters)
    """
    values_list = []
    for dem_path in dem_paths:
        values = get_data_from_file(dem_path=dem_path,
                                    lon=lon,
                                    lat=lat,
                                    num_buffer=num_buffer)
        if values is not None:
            values_list.append(values)

    # Take responses and create lists of lat/lons/values to interpolate over
    x = []
    y = []
    z = []
    for values in values_list:
        x.extend(values[0].flatten())
        y.extend(values[1].flatten())
        z.extend(values[2].flatten())

    # Interpolate over the values
    fun = interp2d(x=x, y=y, z=z, kind=interp_kind)
    return fun(lon, lat)[0]


def get_data_from_file(dem_path, lon, lat, num_buffer):
    """Get array of longitude, latitude, and elevation values from DEM file

    Parameters
    ----------
    dem_path : str or pathlib.Path
        path to DEM file with elevation data. Can be any fileformat readable by
        GDAL.
    lon : float
        longitude
    lat : float
        latitude
    num_buffer : int
        number of bordering cells around (lon, lat) to use when interpolating

    Returns
    -------
    array : None or 3D Numpy array
        returns None if requested cells are out of bounds of the file
        otherwise, returns 3D array containing array of longitude values, array
        of latitude values, and array of elevation values.
    """
    with rasterio.open(dem_path) as dataset:
        # Find row, column of elevation square inside raster
        row, col = dataset.index(lon, lat)

        # Make window include cells around it
        # The number of additional cells depends on the value of num_buffer
        # When num_buffer==1, an additional 8 cells will be loaded and
        # interpolated on;
        # When num_buffer==2, an additional 24 cells will be loaded and
        # interpolated on, etc.
        # When using kind='linear' interpolation, I'm not sure if having the
        # extra cells makes a difference; ie if it creates the plane based only
        # on the closest cells or from all. When using kind='cubic', it's
        # probably more accurate with more cells.
        minrow = row - num_buffer if row >= num_buffer else row
        maxrow = row + num_buffer if row + num_buffer <= dataset.width else row
        mincol = col - num_buffer if col >= num_buffer else col
        maxcol = col + num_buffer if col + num_buffer <= dataset.width else col

        # Restrict to within dataset height/width
        maxrow = dataset.height - 1 if maxrow >= dataset.height else maxrow
        maxcol = dataset.width - 1 if maxcol >= dataset.width else maxcol
        minrow = 0 if minrow < 0 else minrow
        mincol = 0 if mincol < 0 else mincol

        if (minrow >= dataset.height) or (mincol >= dataset.width):
            return None
        if (maxrow <= minrow) or (maxcol <= mincol):
            return None

        # Add +1 to deal with range() not including end
        maxrow += 1
        maxcol += 1

        # Retrieve just this window of values from the DEM
        # NOTE: Assumes that the elevation data is in band 1 of the file
        window = ([minrow, maxrow], [mincol, maxcol])
        val_arr = dataset.read(1, window=window)

        msg = 'array has too few or too many values'
        max_num = 2 * num_buffer + 1
        assert (1 <= val_arr.shape[0] <= max_num) and (1 <= val_arr.shape[1] <=
                                                       max_num), msg

        # Create array of latitude/longitude pairs for each cell center
        lons = []
        lats = []
        for row in range(minrow, maxrow):
            lon_cols = []
            lat_cols = []
            for col in range(mincol, maxcol):
                lon, lat = dataset.xy(row, col)
                lon_cols.append(lon)
                lat_cols.append(lat)

            lons.append(lon_cols)
            lats.append(lat_cols)

        # Array with longitudes, latitudes, values
        # I.e. x, y, z
        return np.array([np.array(lons), np.array(lats), val_arr])


def add_node_elevations(G,
                        api_key,
                        max_locations_per_batch=350,
                        pause_duration=0.02):  # pragma: no cover
    """
    Get the elevation (meters) of each node in the network and add it to the
    node as an attribute.

    Parameters
    ----------
    G : networkx multidigraph
    api_key : string
        your google maps elevation API key
    max_locations_per_batch : int
        max number of coordinate pairs to submit in each API call (if this is
        too high, the server will reject the request because its character
        limit exceeds the max)
    pause_duration : float
        time to pause between API calls

    Returns
    -------
    G : networkx multidigraph
    """

    # google maps elevation API endpoint
    url_template = 'https://maps.googleapis.com/maps/api/elevation/json?locations={}&key={}'

    # make a pandas series of all the nodes' coordinates as 'lat,lng'
    # round coorindates to 5 decimal places (approx 1 meter) to be able to fit
    # in more locations per API call
    node_points = pd.Series({node:'{:.5f},{:.5f}'.format(data['y'], data['x']) for node, data in G.nodes(data=True)})
    log('Requesting node elevations from the API in {} calls.'.format(math.ceil(len(node_points) / max_locations_per_batch)))

    # break the series of coordinates into chunks of size max_locations_per_batch
    # API format is locations=lat,lng|lat,lng|lat,lng|lat,lng...
    results = []
    for i in range(0, len(node_points), max_locations_per_batch):
        chunk = node_points.iloc[i : i + max_locations_per_batch]
        locations = '|'.join(chunk)
        url = url_template.format(locations, api_key)

        # check if this request is already in the cache (if global use_cache=True)
        cached_response_json = get_from_cache(url)
        if cached_response_json is not None:
            response_json = cached_response_json
        else:
            try:
                # request the elevations from the API
                log('Requesting node elevations: {}'.format(url))
                time.sleep(pause_duration)
                response = requests.get(url)
                response_json = response.json()
                save_to_cache(url, response_json)
            except Exception as e:
                log(e)
                log('Server responded with {}: {}'.format(response.status_code, response.reason))

        # append these elevation results to the list of all results
        results.extend(response_json['results'])

    # sanity check that all our vectors have the same number of elements
    if not (len(results) == len(G.nodes()) == len(node_points)):
        raise Exception('Graph has {} nodes but we received {} results from the elevation API.'.format(len(G.nodes()), len(results)))
    else:
        log('Graph has {} nodes and we received {} results from the elevation API.'.format(len(G.nodes()), len(results)))

    # add elevation as an attribute to the nodes
    df = pd.DataFrame(node_points, columns=['node_points'])
    df['elevation'] = [result['elevation'] for result in results]
    df['elevation'] = df['elevation'].round(3) # round to millimeter
    nx.set_node_attributes(G, name='elevation', values=df['elevation'].to_dict())
    log('Added elevation data to all nodes.')

    return G



def add_edge_grades(G, add_absolute=True): # pragma: no cover
    """
    Get the directed grade (ie, rise over run) for each edge in the network and
    add it to the edge as an attribute. Nodes must have elevation attributes to
    use this function.

    Parameters
    ----------
    G : networkx multidigraph
    add_absolute : bool
        if True, also add the absolute value of the grade as an edge attribute

    Returns
    -------
    G : networkx multidigraph
    """

    # for each edge, calculate the difference in elevation from origin to
    # destination, then divide by edge length
    for u, v, data in G.edges(keys=False, data=True):
        elevation_change = G.nodes[v]['elevation'] - G.nodes[u]['elevation']

        # round to ten-thousandths decimal place
        try:
            grade = round(elevation_change / data['length'], 4)
        except ZeroDivisionError:
            grade = None

        # add grade and (optionally) grade absolute value to the edge data
        data['grade'] = grade
        if add_absolute:
            data['grade_abs'] = abs(grade)

    log('Added grade data to all edges.')
    return G
