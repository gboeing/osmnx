"""Project spatial geometries and street networks."""

import math

import geopandas as gpd
from pyproj import CRS

from . import settings
from . import simplification
from . import utils
from . import utils_graph


def project_geometry(geometry, crs=None, to_crs=None, to_latlong=False):
    """
    Project a shapely geometry from its current CRS to another.

    If to_crs is None, project to the UTM CRS for the UTM zone in which the
    geometry's centroid lies. Otherwise project to the CRS defined by to_crs.

    Parameters
    ----------
    geometry : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        the geometry to project
    crs : dict or string or pyproj.CRS
        the starting CRS of the passed-in geometry. if None, it will be set to
        settings.default_crs
    to_crs : dict or string or pyproj.CRS
        if None, project to UTM zone in which geometry's centroid lies,
        otherwise project to this CRS
    to_latlong : bool
        if True, project to settings.default_crs and ignore to_crs

    Returns
    -------
    geometry_proj, crs : tuple
        the projected geometry and its new CRS
    """
    if crs is None:
        crs = settings.default_crs

    gdf = gpd.GeoDataFrame(geometry=[geometry], crs=crs)
    gdf_proj = project_gdf(gdf, to_crs=to_crs, to_latlong=to_latlong)
    geometry_proj = gdf_proj["geometry"].iloc[0]
    return geometry_proj, gdf_proj.crs


def project_gdf(gdf, to_crs=None, to_latlong=False):
    """
    Project a GeoDataFrame from its current CRS to another.

    If to_crs is None, project to the UTM CRS for the UTM zone in which the
    GeoDataFrame's centroid lies. Otherwise project to the CRS defined by
    to_crs. The simple UTM zone calculation in this function works well for
    most latitudes, but may not work for some extreme northern locations like
    Svalbard or far northern Norway.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        the GeoDataFrame to be projected
    to_crs : dict or string or pyproj.CRS
        if None, project to UTM zone in which gdf's centroid lies, otherwise
        project to this CRS
    to_latlong : bool
        if True, project to settings.default_crs and ignore to_crs

    Returns
    -------
    gdf_proj : geopandas.GeoDataFrame
        the projected GeoDataFrame
    """
    if gdf.crs is None or len(gdf) < 1:
        raise ValueError("GeoDataFrame must have a valid CRS and cannot be empty")

    # if to_latlong is True, project the gdf to latlong
    if to_latlong:
        gdf_proj = gdf.to_crs(settings.default_crs)
        utils.log(f"Projected GeoDataFrame to {settings.default_crs}")

    # else if to_crs was passed-in, project gdf to this CRS
    elif to_crs is not None:
        gdf_proj = gdf.to_crs(to_crs)
        utils.log(f"Projected GeoDataFrame to {to_crs}")

    # otherwise, automatically project the gdf to UTM
    else:
        if CRS.from_user_input(gdf.crs).is_projected:
            raise ValueError("Geometry must be unprojected to calculate UTM zone")

        # calculate longitude of centroid of union of all geometries in gdf
        avg_lng = gdf["geometry"].unary_union.centroid.x

        # calculate UTM zone from avg longitude to define CRS to project to
        utm_zone = int(math.floor((avg_lng + 180) / 6.0) + 1)
        utm_crs = f"+proj=utm +zone={utm_zone} +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

        # project the GeoDataFrame to the UTM CRS
        gdf_proj = gdf.to_crs(utm_crs)
        utils.log(f"Projected GeoDataFrame to {gdf_proj.crs}")

    return gdf_proj


def project_graph(G, to_crs=None):
    """
    Project graph from its current CRS to another.

    If to_crs is None, project the graph to the UTM CRS for the UTM zone in
    which the graph's centroid lies. Otherwise, project the graph to the CRS
    defined by to_crs.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        the graph to be projected
    to_crs : dict or string or pyproj.CRS
        if None, project graph to UTM zone in which graph centroid lies,
        otherwise project graph to this CRS

    Returns
    -------
    G_proj : networkx.MultiDiGraph
        the projected graph
    """
    # STEP 1: PROJECT THE NODES
    gdf_nodes = utils_graph.graph_to_gdfs(G, edges=False)

    # create new lat/lng columns to preserve lat/lng for later reference if
    # cols do not already exist (ie, don't overwrite in later re-projections)
    if "lon" not in gdf_nodes.columns or "lat" not in gdf_nodes.columns:
        gdf_nodes["lon"] = gdf_nodes["x"]
        gdf_nodes["lat"] = gdf_nodes["y"]

    # project the nodes GeoDataFrame and extract the projected x/y values
    gdf_nodes_proj = project_gdf(gdf_nodes, to_crs=to_crs)
    gdf_nodes_proj["x"] = gdf_nodes_proj["geometry"].x
    gdf_nodes_proj["y"] = gdf_nodes_proj["geometry"].y
    gdf_nodes_proj = gdf_nodes_proj.drop(columns=["geometry"])

    # STEP 2: PROJECT THE EDGES
    if simplification._is_simplified(G):
        # if graph has previously been simplified, project the edge geometries
        gdf_edges = utils_graph.graph_to_gdfs(G, nodes=False, fill_edge_geometry=True)
        gdf_edges_proj = project_gdf(gdf_edges, to_crs=gdf_nodes_proj.crs)
    else:
        # if not, you don't have to project these edges because the nodes
        # contain all the spatial data in the graph (unsimplified edges have
        # no geometry attributes)
        gdf_edges_proj = utils_graph.graph_to_gdfs(G, nodes=False, fill_edge_geometry=False).drop(
            columns=["geometry"]
        )

    # STEP 3: REBUILD GRAPH
    # turn projected node/edge gdfs into a graph and update its CRS attribute
    G_proj = utils_graph.graph_from_gdfs(gdf_nodes_proj, gdf_edges_proj, G.graph)
    G_proj.graph["crs"] = gdf_nodes_proj.crs

    utils.log(f"Projected graph with {len(G)} nodes and {len(G.edges)} edges")
    return G_proj
