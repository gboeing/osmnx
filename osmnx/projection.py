"""Project spatial geometries and street networks."""

import math

import geopandas as gpd
import networkx as nx
from pyproj import CRS
from shapely.geometry import Point

from . import settings
from . import utils


def _is_crs_utm(crs):
    """
    Determine if a CRS is a UTM CRS.

    Parameters
    ----------
    crs : dict or string or pyproj.CRS
        a coordinate reference system

    Returns
    -------
    bool
        True if crs is UTM, False otherwise
    """
    if not crs:
        return False
    crs_obj = CRS.from_user_input(crs)
    if crs_obj.coordinate_operation and crs_obj.coordinate_operation.name.upper().startswith("UTM"):
        return True
    return False


def project_geometry(geometry, crs=None, to_crs=None, to_latlong=False):
    """
    Project a shapely (Multi)Polygon from lat-lng to UTM, or vice-versa.

    Parameters
    ----------
    geometry : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        the geometry to project
    crs : dict or string or pyproj.CRS
        the starting coordinate reference system of the passed-in geometry,
        default value (None) will set settings.default_crs as the CRS
    to_crs : dict or string or pyproj.CRS
        if not None, just project to this CRS instead of to UTM
    to_latlong : bool
        if True, project from crs to lat-lng, if False, project from crs to
        local UTM zone

    Returns
    -------
    geometry_proj, crs : tuple
        the projected shapely geometry and the crs of the projected geometry
    """
    if crs is None:
        crs = settings.default_crs

    gdf = gpd.GeoDataFrame()
    gdf.crs = crs
    gdf["geometry"] = None
    gdf.loc[0, "geometry"] = geometry
    gdf_proj = project_gdf(gdf, to_crs=to_crs, to_latlong=to_latlong)
    geometry_proj = gdf_proj["geometry"].iloc[0]
    return geometry_proj, gdf_proj.crs


def project_gdf(gdf, to_crs=None, to_latlong=False):
    """
    Project a GeoDataFrame to UTM.

    Automatically chooses the UTM zone appropriate for its geometries'
    centroid. The simple calculation in this function works well for most
    latitudes, but won't work for some far northern locations like Svalbard
    and parts of far northern Norway.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        the gdf to be projected
    to_crs : dict or string or pyproj.CRS
        if not None, project to this CRS instead of to UTM
    to_latlong : bool
        if True, projects to settings.default_crs instead of to UTM

    Returns
    -------
    gdf_proj : geopandas.GeoDataFrame
        the projected GeoDataFrame
    """
    if len(gdf) < 1:
        raise ValueError("Cannot project an empty GeoDataFrame")

    # if to_crs was passed-in, use this value to project the gdf
    if to_crs is not None:
        gdf_proj = gdf.to_crs(to_crs)

    # if to_crs was not passed-in, calculate the centroid of the geometry to
    # determine UTM zone
    else:
        if to_latlong:
            # if to_latlong is True, project the gdf to latlong
            latlong_crs = settings.default_crs
            gdf_proj = gdf.to_crs(latlong_crs)
            utils.log("Projected GeoDataFrame to settings.default_crs")
        else:
            # else, project the gdf to UTM
            # if GeoDataFrame is already in UTM, just return it
            if _is_crs_utm(gdf.crs):
                return gdf

            # calculate the centroid of the union of all the geometries in the
            # GeoDataFrame
            avg_longitude = gdf["geometry"].unary_union.centroid.x

            # calculate the UTM zone from this avg longitude and define the UTM
            # CRS to project
            utm_zone = int(math.floor((avg_longitude + 180) / 6.0) + 1)
            utm_crs = f"+proj=utm +zone={utm_zone} +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

            # project the GeoDataFrame to the UTM CRS
            gdf_proj = gdf.to_crs(utm_crs)
            utils.log(f"Projected GeoDataFrame to UTM-{utm_zone}")

    return gdf_proj


def project_graph(G, to_crs=None):
    """
    Project graph from lat-lng to UTM zone appropriate for its centroid.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        the networkx graph to be projected
    to_crs : dict or string or pyproj.CRS
        if not None, just project to this CRS instead of to UTM

    Returns
    -------
    G_proj : networkx.MultiDiGraph
        the projected graph
    """
    G_proj = G.copy()

    # create a GeoDataFrame of the nodes, name it, convert osmid to str
    nodes, data = zip(*G_proj.nodes(data=True))
    gdf_nodes = gpd.GeoDataFrame(list(data), index=nodes)
    gdf_nodes.crs = G_proj.graph["crs"]

    # create new lat-lng columns just to save that data for later reference
    # if they do not already exist (i.e., don't overwrite in subsequent re-projections)
    if "lon" not in gdf_nodes.columns or "lat" not in gdf_nodes.columns:
        gdf_nodes["lon"] = gdf_nodes["x"]
        gdf_nodes["lat"] = gdf_nodes["y"]

    # create a geometry column from x/y columns
    gdf_nodes["geometry"] = gdf_nodes.apply(lambda row: Point(row["x"], row["y"]), axis=1)
    gdf_nodes.set_geometry("geometry", inplace=True)
    utils.log("Created a GeoDataFrame from graph")

    # project the nodes GeoDataFrame to UTM
    gdf_nodes_utm = project_gdf(gdf_nodes, to_crs=to_crs)

    # extract data for all edges that have geometry attribute
    edges_with_geom = []
    for u, v, key, data in G_proj.edges(keys=True, data=True):
        if "geometry" in data:
            edges_with_geom.append({"u": u, "v": v, "key": key, "geometry": data["geometry"]})

    # create an edges GeoDataFrame and project to UTM, if there were any edges
    # with a geometry attribute. geom attr only exists if graph has been
    # simplified, otherwise you don't have to project anything for the edges
    # because the nodes still contain all spatial data
    if len(edges_with_geom) > 0:
        gdf_edges = gpd.GeoDataFrame(edges_with_geom)
        gdf_edges.crs = G_proj.graph["crs"]
        gdf_edges_utm = project_gdf(gdf_edges, to_crs=to_crs)

    # extract projected x and y values from the nodes' geometry column
    gdf_nodes_utm["x"] = gdf_nodes_utm["geometry"].map(lambda point: point.x)
    gdf_nodes_utm["y"] = gdf_nodes_utm["geometry"].map(lambda point: point.y)
    gdf_nodes_utm = gdf_nodes_utm.drop("geometry", axis=1)
    utils.log("Extracted projected node geometries from GeoDataFrame")

    # clear the graph to make it a blank slate for the projected data
    edges = list(G_proj.edges(keys=True, data=True))
    G_proj.clear()

    # add the projected nodes and all their attributes to the graph
    G_proj.add_nodes_from(gdf_nodes_utm.index)
    attributes = gdf_nodes_utm.to_dict()
    for label in gdf_nodes_utm.columns:
        nx.set_node_attributes(G_proj, name=label, values=attributes[label])

    # add the edges and all their attributes (including reconstructed geometry,
    # when it exists) to the graph
    for u, v, key, attributes in edges:
        if "geometry" in attributes:
            mask = (
                (gdf_edges_utm["u"] == u)
                & (gdf_edges_utm["v"] == v)
                & (gdf_edges_utm["key"] == key)
            )
            row = gdf_edges_utm[mask]
            attributes["geometry"] = row["geometry"].iloc[0]

        # attributes dict contains key, so we don't need to explicitly pass it here
        G_proj.add_edge(u, v, **attributes)

    # set the graph's CRS attribute to the new, projected CRS and return the
    # projected graph
    G_proj.graph["crs"] = gdf_nodes_utm.crs
    if "streets_per_node" in G.graph:
        G_proj.graph["streets_per_node"] = G.graph["streets_per_node"]
    utils.log("Rebuilt projected graph")
    return G_proj
