"""Project a graph, GeoDataFrame, or geometry to a different CRS."""

from __future__ import annotations

import logging as lg
from typing import TYPE_CHECKING
from typing import Any

import geopandas as gpd

from . import convert
from . import settings
from . import utils

if TYPE_CHECKING:
    import networkx as nx
    from shapely import Geometry


def is_projected(crs: Any) -> bool:  # noqa: ANN401
    """
    Determine if a coordinate reference system is projected or not.

    Parameters
    ----------
    crs
        The identifier of the coordinate reference system. This can be
        anything accepted by `pyproj.CRS.from_user_input()`, such as an
        authority string or a WKT string.

    Returns
    -------
    projected
        True if `crs` is projected, otherwise False
    """
    return bool(gpd.GeoSeries(crs=crs).crs.is_projected)


def project_geometry(
    geom: Geometry,
    *,
    crs: Any | None = None,  # noqa: ANN401
    to_crs: Any | None = None,  # noqa: ANN401
    to_latlong: bool = False,
) -> tuple[Geometry, Any]:
    """
    Project a Shapely geometry from its current CRS to another.

    If `to_latlong` is True, this projects the geometry to the coordinate
    reference system defined by `settings.default_crs`. Otherwise it projects
    it to the CRS defined by `to_crs`. If `to_crs` is `None`, it projects it
    to the CRS of an appropriate UTM zone given `geometry`'s bounds.

    Parameters
    ----------
    geom
        The geometry to be projected.
    crs
        The initial CRS of `geometry`. If None, it will be set to
        `settings.default_crs`.
    to_crs
        If None, project to an appropriate UTM zone. Otherwise project to this
        CRS.
    to_latlong
        If True, project to `settings.default_crs` and ignore `to_crs`.

    Returns
    -------
    geom_proj, crs
        The projected geometry and its new CRS.
    """
    if crs is None:
        crs = settings.default_crs

    gdf = gpd.GeoDataFrame(geometry=[geom], crs=crs)
    gdf_proj = project_gdf(gdf, to_crs=to_crs, to_latlong=to_latlong)
    geom_proj = gdf_proj["geometry"].iloc[0]
    return geom_proj, gdf_proj.crs


def project_gdf(
    gdf: gpd.GeoDataFrame,
    *,
    to_crs: Any | None = None,  # noqa: ANN401
    to_latlong: bool = False,
) -> gpd.GeoDataFrame:
    """
    Project a GeoDataFrame from its current CRS to another.

    If `to_latlong` is True, this projects the GeoDataFrame to the coordinate
    reference system defined by `settings.default_crs`. Otherwise it projects
    it to the CRS defined by `to_crs`. If `to_crs` is `None`, it projects it
    to the CRS of an appropriate UTM zone given `gdf`'s bounds.

    Parameters
    ----------
    gdf
        The GeoDataFrame to be projected.
    to_crs
        If None, project to an appropriate UTM zone. Otherwise project to
        this CRS.
    to_latlong
        If True, project to `settings.default_crs` and ignore `to_crs`.

    Returns
    -------
    gdf_proj
        The projected GeoDataFrame.
    """
    if gdf.crs is None or len(gdf) == 0:  # pragma: no cover
        msg = "`gdf` must have a valid CRS and cannot be empty."
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

    msg = f"Projected GeoDataFrame to {crs_desc!r}"
    utils.log(msg, level=lg.INFO)
    return gdf_proj


def project_graph(
    G: nx.MultiDiGraph,
    *,
    to_crs: Any | None = None,  # noqa: ANN401
    to_latlong: bool = False,
) -> nx.MultiDiGraph:
    """
    Project a graph from its current CRS to another.

    If `to_latlong` is True, this projects the graph to the coordinate
    reference system defined by `settings.default_crs`. Otherwise it projects
    it to the CRS defined by `to_crs`. If `to_crs` is `None`, it projects it
    to the CRS of an appropriate UTM zone given `geometry`'s bounds.

    Parameters
    ----------
    G
        The graph to be projected.
    to_crs
        If None, project to an appropriate UTM zone. Otherwise project to
        this CRS.
    to_latlong
        If True, project to `settings.default_crs` and ignore `to_crs`.

    Returns
    -------
    G_proj
        The projected graph.
    """
    if to_latlong:
        to_crs = settings.default_crs

    # STEP 1: PROJECT THE NODES
    gdf_nodes = convert.graph_to_gdfs(G, edges=False)

    # project the nodes GeoDataFrame and extract the projected x/y values
    gdf_nodes_proj = project_gdf(gdf_nodes, to_crs=to_crs)
    gdf_nodes_proj["x"] = gdf_nodes_proj["geometry"].x
    gdf_nodes_proj["y"] = gdf_nodes_proj["geometry"].y
    to_crs = gdf_nodes_proj.crs

    # STEP 2: PROJECT THE EDGES
    if G.graph.get("simplified"):
        # if graph has previously been simplified, project the edge geometries
        gdf_edges = convert.graph_to_gdfs(G, nodes=False, fill_edge_geometry=False)
        gdf_edges_proj = project_gdf(gdf_edges, to_crs=to_crs)
    else:
        # if not, you don't have to project these edges because the nodes
        # contain all the spatial data in the graph (unsimplified edges have
        # no geometry attributes)
        gdf_edges_proj = convert.graph_to_gdfs(G, nodes=False, fill_edge_geometry=False)

    # STEP 3: REBUILD GRAPH
    # turn projected node/edge gdfs into a graph and update its CRS attribute
    G_proj = convert.graph_from_gdfs(gdf_nodes_proj, gdf_edges_proj, graph_attrs=G.graph)
    G_proj.graph["crs"] = to_crs

    msg = f"Projected graph with {len(G)} nodes and {len(G.edges)} edges"
    utils.log(msg, level=lg.INFO)
    return G_proj
