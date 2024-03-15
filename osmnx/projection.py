"""Project a graph, GeoDataFrame, or geometry to a different CRS."""

import geopandas as gpd

from . import convert
from . import settings
from . import utils


def is_projected(crs):
    """
    Determine if a coordinate reference system is projected or not.

    Parameters
    ----------
    crs : string or pyproj.CRS
        the identifier of the coordinate reference system, which can be
        anything accepted by `pyproj.CRS.from_user_input()` such as an
        authority string or a WKT string

    Returns
    -------
    projected : bool
        True if crs is projected, otherwise False
    """
    return gpd.GeoSeries(crs=crs).crs.is_projected


def project_geometry(geometry, crs=None, to_crs=None, to_latlong=False):
    """
    Project a Shapely geometry from its current CRS to another.

    If `to_latlong` is `True`, this projects the GeoDataFrame to the CRS
    defined by `settings.default_crs`, otherwise it projects it to the CRS
    defined by `to_crs`. If `to_crs` is `None`, it projects it to the CRS of
    an appropriate UTM zone given `geometry`'s bounds.

    Parameters
    ----------
    geometry : shapely geometry
        the geometry to be projected
    crs : string or pyproj.CRS
        the initial CRS of `geometry`. if None, it will be set to
        `settings.default_crs`
    to_crs : string or pyproj.CRS
        if None, project to an appropriate UTM zone, otherwise project to
        this CRS
    to_latlong : bool
        if True, project to `settings.default_crs` and ignore `to_crs`

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

    If `to_latlong` is `True`, this projects the GeoDataFrame to the CRS
    defined by `settings.default_crs`, otherwise it projects it to the CRS
    defined by `to_crs`. If `to_crs` is `None`, it projects it to the CRS of
    an appropriate UTM zone given `gdf`'s bounds.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        the GeoDataFrame to be projected
    to_crs : string or pyproj.CRS
        if None, project to an appropriate UTM zone, otherwise project to
        this CRS
    to_latlong : bool
        if True, project to `settings.default_crs` and ignore `to_crs`

    Returns
    -------
    gdf_proj : geopandas.GeoDataFrame
        the projected GeoDataFrame
    """
    if gdf.crs is None or len(gdf) < 1:  # pragma: no cover
        msg = "GeoDataFrame must have a valid CRS and cannot be empty"
        raise ValueError(msg)

    # if to_latlong is True, project the gdf to the default_crs
    if to_latlong:
        to_crs = settings.default_crs

    # else if to_crs is None, project gdf to an appropriate UTM zone
    elif to_crs is None:
        to_crs = gdf.estimate_utm_crs()

    # project the gdf
    gdf_proj = gdf.to_crs(to_crs)
    crs_desc = f"{gdf_proj.crs.to_string()} / {gdf_proj.crs.name}"
    utils.log(f"Projected GeoDataFrame to {crs_desc!r}")
    return gdf_proj


def project_graph(G, to_crs=None, to_latlong=False):
    """
    Project a graph from its current CRS to another.

    If `to_latlong` is `True`, this projects the GeoDataFrame to the CRS
    defined by `settings.default_crs`, otherwise it projects it to the CRS
    defined by `to_crs`. If `to_crs` is `None`, it projects it to the CRS of
    an appropriate UTM zone given `G`'s bounds.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        the graph to be projected
    to_crs : string or pyproj.CRS
        if None, project to an appropriate UTM zone, otherwise project to
        this CRS
    to_latlong : bool
        if True, project to `settings.default_crs` and ignore `to_crs`

    Returns
    -------
    G_proj : networkx.MultiDiGraph
        the projected graph
    """
    if to_latlong:
        to_crs = settings.default_crs

    # STEP 1: PROJECT THE NODES
    gdf_nodes = convert.graph_to_gdfs(G, edges=False)

    # create new lat/lon columns to preserve lat/lon for later reference if
    # cols do not already exist (ie, don't overwrite in later re-projections)
    if "lon" not in gdf_nodes.columns or "lat" not in gdf_nodes.columns:
        gdf_nodes["lon"] = gdf_nodes["x"]
        gdf_nodes["lat"] = gdf_nodes["y"]

    # project the nodes GeoDataFrame and extract the projected x/y values
    gdf_nodes_proj = project_gdf(gdf_nodes, to_crs=to_crs)
    gdf_nodes_proj["x"] = gdf_nodes_proj["geometry"].x
    gdf_nodes_proj["y"] = gdf_nodes_proj["geometry"].y
    to_crs = gdf_nodes_proj.crs
    gdf_nodes_proj = gdf_nodes_proj.drop(columns=["geometry"])

    # STEP 2: PROJECT THE EDGES
    if "simplified" in G.graph and G.graph["simplified"]:
        # if graph has previously been simplified, project the edge geometries
        gdf_edges = convert.graph_to_gdfs(G, nodes=False, fill_edge_geometry=False)
        gdf_edges_proj = project_gdf(gdf_edges, to_crs=to_crs)
    else:
        # if not, you don't have to project these edges because the nodes
        # contain all the spatial data in the graph (unsimplified edges have
        # no geometry attributes)
        gdf_edges_proj = convert.graph_to_gdfs(G, nodes=False, fill_edge_geometry=False).drop(
            columns=["geometry"]
        )

    # STEP 3: REBUILD GRAPH
    # turn projected node/edge gdfs into a graph and update its CRS attribute
    G_proj = convert.graph_from_gdfs(gdf_nodes_proj, gdf_edges_proj, G.graph)
    G_proj.graph["crs"] = to_crs

    utils.log(f"Projected graph with {len(G)} nodes and {len(G.edges)} edges")
    return G_proj
