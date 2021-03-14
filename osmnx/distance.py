"""Calculate distances and shortest paths and find nearest node/edge(s) to point(s)."""

import itertools

import networkx as nx
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from shapely.geometry import Point

from . import projection
from . import utils
from . import utils_geo
from . import utils_graph

# scikit-learn is optional dependency for unprojected nearest-neighbor search
try:
    from sklearn.neighbors import BallTree
except ImportError:
    BallTree = None


def great_circle_vec(lat1, lng1, lat2, lng2, earth_radius=6_371_009):
    """
    Calculate great-circle distances between points.

    Vectorized function to calculate the great-circle distance between two
    points' coordinates or between arrays of points' coordinates using the
    haversine formula. Expects coordinates in decimal degrees.

    Parameters
    ----------
    lat1 : float or np.array of float
        first point's latitude coordinate
    lng1 : float or np.array of float
        first point's longitude coordinate
    lat2 : float or np.array of float
        second point's latitude coordinate
    lng2 : float or np.array of float
        second point's longitude coordinate
    earth_radius : int or float
        radius of earth in units in which distance will be returned
        (default is meters)

    Returns
    -------
    dist : float or np.array
        distance or array of distances from (lat1, lng1) to (lat2, lng2) in
        units of earth_radius
    """
    phi1 = np.deg2rad(lat1)
    phi2 = np.deg2rad(lat2)
    d_phi = phi2 - phi1

    theta1 = np.deg2rad(lng1)
    theta2 = np.deg2rad(lng2)
    d_theta = theta2 - theta1

    h = np.sin(d_phi / 2) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(d_theta / 2) ** 2
    h = np.minimum(1.0, h)  # protect against floating point errors

    arc = 2 * np.arcsin(np.sqrt(h))

    # return distance in units of earth_radius
    return arc * earth_radius


def euclidean_dist_vec(y1, x1, y2, x2):
    """
    Calculate Euclidean distances between points.

    Vectorized function to calculate the Euclidean distance between two
    points' coordinates or between arrays of points' coordinates. For most
    accurate results, use projected coordinates rather than decimal degrees.

    Parameters
    ----------
    y1 : float or np.array of float
        first point's y coordinate
    x1 : float or np.array of float
        first point's x coordinate
    y2 : float or np.array of float
        second point's y coordinate
    x2 : float or np.array of float
        second point's x coordinate

    Returns
    -------
    dist : float or np.array of float
        distance or array of distances from (x1, y1) to (x2, y2) in
        coordinates' units
    """
    # pythagorean theorem
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5


def nearest_nodes(G, X, Y, return_dist=False):
    """
    Find the nearest node(s) to some point(s).

    If the graph is projected, this will use a k-d tree for fast euclidean
    nearest-neighbor search. If it is unprojected, this will use a ball tree
    for fast haversine nearest-neighbor search, which requires that
    scikit-learn is installed as an optional dependency.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        graph in which to find nearest nodes
    X : list
        points' x or longitude coordinates, in same CRS and units as graph
    Y : list
        points' y or latitude coordinates, in same CRS and units as graph
    return_dist : bool
        optionally also return distance between points and nearest nodes

    Returns
    -------
    nn or (nn, dist): list or tuple of list
        nearest node IDs or optionally tuple of lists where `dist` contains
        distances between the points and their nearest nodes
    """
    nodes = utils_graph.graph_to_gdfs(G, edges=False, node_geometry=False)[["x", "y"]]

    if projection.is_projected(G.graph["crs"]):
        # if projected, use k-d tree for euclidean nearest-neighbor search
        dist, pos = cKDTree(nodes).query(np.array([X, Y]).T, k=1)
        nn = nodes.index[pos].tolist()

    else:
        # if unprojected, use ball tree for haversine nearest-neighbor search
        if BallTree is None:
            raise ImportError("scikit-learn must be installed to search an unprojected graph")
        # haversine requires lat, lng coords in radians
        nodes_rad = np.deg2rad(nodes[["y", "x"]])
        points_rad = np.deg2rad(np.array([Y, X]).T)
        dist, pos = BallTree(nodes_rad, metric="haversine").query(points_rad, k=1)
        nn = nodes.index[pos[:, 0]].tolist()

    if return_dist:
        return nn, dist
    else:
        return nn


def nearest_edges(G, X, Y, spacing=None, return_dist=False):
    """
    Find the nearest edge(s) to some point(s).

    If spacing is None, search for the nearest edge to each point, one at a
    time, using an r-tree and minimizing the euclidean distances from the
    point to the possible matches. Can be slow if searching many points. For
    best results, use a projected graph and points.

    Otherwise, if spacing is not None: interpolate points along edges to use a
    k-d tree or ball tree for search. if the graph is projected, this will
    use a k-d tree for fast euclidean nearest-neighbor search. If graph is
    unprojected, this will use a ball tree for fast haversine nearest-neighbor
    search, which requires that scikit-learn is installed as an optional
    dependency.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        graph in which to find nearest edges
    X : list
        points' x or longitude coordinates, in same CRS and units as graph
    Y : list
        points' y or latitude coordinates, in same CRS and units as graph
    spacing : float
        spacing between interpolated points, in same units as graph. smaller
        values generate more points.
    return_dist : bool
        optionally also return distance between points and nearest edges

    Returns
    -------
    ne or (ne, dist): list or tuple of list
        nearest edges as (u, v, key) or optionally tuple of lists where `dist`
        contains distances between the points and their nearest edges
    """
    geoms = ox.utils_graph.graph_to_gdfs(G, nodes=False)["geometry"]

    # if no interpolation spacing is provided, use an r-tree to find possible
    # matches, then minimize euclidean distance from point to possible matches
    if not spacing:
        ne = list()
        dist = list()
        for xy in zip(X, Y):
            distances = geoms.iloc[list(sindex.nearest(xy))].distance(Point(xy))
            ne.append(distances.idxmin())
            dist.append(distances.min())

    # otherwise, interpolate points along edges for k-d tree or ball tree
    else:
        uvk_xy = list()
        for uvk, geom in zip(geoms.index, geoms.values):
            uvk_xy.extend((uvk, xy) for xy in ox.utils_geo.interpolate_points(geom, spacing))
        labels, xy = zip(*uvk_xy)
        vertices = pd.DataFrame(xy, index=labels, columns=["x", "y"])

        if ox.projection.is_projected(G.graph["crs"]):
            # if projected, use k-d tree for euclidean nearest-neighbor search
            dist, pos = cKDTree(vertices).query(np.array([X, Y]).T, k=1)
            ne = vertices.index[pos].tolist()

        else:
            # if unprojected, use ball tree for haversine nearest-neighbor search
            if BallTree is None:
                raise ImportError("scikit-learn must be installed to search an unprojected graph")
            # haversine requires lat, lng coords in radians
            vertices_rad = np.deg2rad(vertices[["y", "x"]])
            points_rad = np.deg2rad(np.array([Y, X]).T)
            dist, pos = BallTree(vertices_rad, metric="haversine").query(points_rad, k=1)
            ne = vertices.index[pos[:, 0]].tolist()

    if return_dist:
        return ne, dist
    else:
        return ne


def get_nearest_node(G, point, method="haversine", return_dist=False):
    """
    Find the nearest node to a point.

    Return the graph node nearest to some (lat, lng) or (y, x) point and
    optionally the distance between the node and the point. This function can
    use either the haversine formula or Euclidean distance.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    point : tuple
        The (lat, lng) or (y, x) point for which we will find the nearest node
        in the graph
    method : string {'haversine', 'euclidean'}
        Which method to use for calculating distances to find nearest node.
        If 'haversine', graph nodes' coordinates must be in units of decimal
        degrees. If 'euclidean', graph nodes' coordinates must be projected.
    return_dist : bool
        Optionally also return the distance (in meters if haversine, or graph
        node coordinate units if euclidean) between the point and the nearest
        node

    Returns
    -------
    int or tuple of (int, float)
        Nearest node ID or optionally a tuple of (node ID, dist), where dist
        is the distance (in meters if haversine, or graph node coordinate
        units if euclidean) between the point and nearest node
    """
    if not G:
        raise ValueError("G must contain at least one node")

    # dump graph node coordinates into a pandas dataframe indexed by node id
    # with x and y columns
    coords = ((n, d["x"], d["y"]) for n, d in G.nodes(data=True))
    df = pd.DataFrame(coords, columns=["node", "x", "y"]).set_index("node")

    # add columns to df for the (constant) coordinates of reference point
    df["ref_y"] = point[0]
    df["ref_x"] = point[1]

    # calculate the distance between each node and the reference point
    if method == "haversine":
        # calculate distances using haversine for spherical lat-lng geometries
        dists = great_circle_vec(lat1=df["ref_y"], lng1=df["ref_x"], lat2=df["y"], lng2=df["x"])

    elif method == "euclidean":
        # calculate distances using euclid's formula for projected geometries
        dists = euclidean_dist_vec(y1=df["ref_y"], x1=df["ref_x"], y2=df["y"], x2=df["x"])

    else:
        raise ValueError('method argument must be either "haversine" or "euclidean"')

    # nearest node's ID is the index label of the minimum distance
    nearest_node = dists.idxmin()
    utils.log(f"Found nearest node ({nearest_node}) to point {point}")

    # if caller requested return_dist, return distance between the point and the
    # nearest node as well
    if return_dist:
        return nearest_node, dists.loc[nearest_node]
    else:
        return nearest_node


def get_nearest_edge(G, point, return_geom=False, return_dist=False):
    """
    Find the nearest edge to a point by minimum Euclidean distance.

    For best results, both G and point should be projected.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    point : tuple
        the (lat, lng) or (y, x) point for which we will find the nearest edge
        in the graph
    return_geom : bool
        Optionally return the geometry of the nearest edge
    return_dist : bool
        Optionally return the distance in graph's coordinates' units between
        the point and the nearest edge

    Returns
    -------
    tuple
        Graph edge unique identifier as a tuple of (u, v, key).
        Or a tuple of (u, v, key, geom) if return_geom is True.
        Or a tuple of (u, v, key, dist) if return_dist is True.
        Or a tuple of (u, v, key, geom, dist) if return_geom and return_dist are True.
    """
    # convert lat,lng (y,x) point to x,y for shapely distance operation
    xy_point = Point(reversed(point))

    # calculate euclidean distance from each edge's geometry to this point
    gs_edges = utils_graph.graph_to_gdfs(G, nodes=False)["geometry"]
    uvk_geoms = zip(gs_edges.index, gs_edges.values)
    distances = ((uvk, geom, xy_point.distance(geom)) for uvk, geom in uvk_geoms)

    # the nearest edge minimizes the distance to the point
    (u, v, key), geom, dist = min(distances, key=lambda x: x[2])
    utils.log(f"Found nearest edge ({u, v, key}) to point {point}")

    # return results requested by caller
    if return_dist and return_geom:
        return u, v, key, geom, dist
    elif return_dist:
        return u, v, key, dist
    elif return_geom:
        return u, v, key, geom
    else:
        return u, v, key


def get_nearest_nodes(G, X, Y, method=None, return_dist=False):
    """
    Find the nearest node to each point in a list of points.

    Pass in points as separate lists of X and Y coordinates. The 'kdtree'
    method is by far the fastest with large data sets, but only finds
    approximate nearest nodes if working in unprojected coordinates like
    lat-lng (it precisely finds the nearest node if working in projected
    coordinates). The 'balltree' method is second fastest with large data
    sets but it is precise if working in unprojected coordinates like lat-lng.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    X : list-like
        the longitudes or x coordinates for which we will find the nearest
        node in the graph
    Y : list-like
        the latitudes or y coordinates for which we will find the nearest node
        in the graph
    method : string {None, 'kdtree', 'balltree'}
        Which method to use for finding the nearest node to each point. If
        None, we manually find each node one at a time using
        utils.get_nearest_node and haversine. If 'kdtree' we use
        scipy.spatial.cKDTree for very fast euclidean search. If
        'balltree', we use sklearn.neighbors.BallTree for fast
        haversine search.
    return_dist : bool
        Optionally also return the distance (in meters if haversine, or graph
        node coordinate units if euclidean) between the point and the nearest
        node

    Returns
    -------
    nn or (dist, nn): np.array or tuple of (np.array<float>, np.array<int>)
        array of node IDs representing the node nearest to each point in the
        passed-in list of points, or optionally tuple of arrays (dist,nn) where dist
        is an array of the distances (in meters if haversine, or graph node coordinate
        units if euclidean) between the point and nearest node
    """
    if method is None:
        # calculate nearest node (and optionally dist) one at a time for each point
        result = [
            get_nearest_node(G, (y, x), method="haversine", return_dist=return_dist)
            for x, y in zip(X, Y)
        ]
        if return_dist:
            nn, dist = map(list, zip(*result))
        else:
            nn = result
    elif method == "kdtree":

        # build a k-d tree for euclidean nearest node search
        nodes = pd.DataFrame(
            {"x": nx.get_node_attributes(G, "x"), "y": nx.get_node_attributes(G, "y")}
        )
        tree = cKDTree(data=nodes[["x", "y"]], compact_nodes=True, balanced_tree=True)

        # query the tree for nearest node to each point
        points = np.array([X, Y]).T
        dist, idx = tree.query(points, k=1)
        nn = nodes.iloc[idx].index

    elif method == "balltree":

        # check if we were able to import sklearn.neighbors.BallTree
        if BallTree is None:
            raise ImportError("scikit-learn must be installed to use this optional feature")

        # haversine requires data in form of [lat, lng] and inputs/outputs in
        # units of radians
        nodes = pd.DataFrame(
            {"x": nx.get_node_attributes(G, "x"), "y": nx.get_node_attributes(G, "y")}
        )
        nodes_rad = np.deg2rad(nodes[["y", "x"]].astype(np.float))
        points = np.array([Y.astype(np.float), X.astype(np.float)]).T
        points_rad = np.deg2rad(points)

        # build a ball tree for haversine nearest node search
        tree = BallTree(nodes_rad, metric="haversine")

        # query the tree for nearest node to each point
        if return_dist:
            dist, idx = tree.query(points_rad, k=1, return_distance=True)
            nn = nodes.iloc[idx[:, 0]].index
        else:
            idx = tree.query(points_rad, k=1, return_distance=False)
            nn = nodes.iloc[idx[:, 0]].index

    else:
        raise ValueError("You must pass a valid method name, or None.")

    utils.log(f"Found nearest nodes to {len(X)} points")

    if return_dist:
        return np.array(nn), np.array(dist)
    else:
        return np.array(nn)


def get_nearest_edges(G, X, Y, method=None, dist=0.0001):
    """
    Find the nearest edge to each point in a list of points.

    Pass in points as separate lists of X and Y coordinates. The 'kdtree'
    method is by far the fastest with large data sets, but only finds
    approximate nearest edges if working in unprojected coordinates like
    lat-lng (it precisely finds the nearest edge if working in projected
    coordinates). The 'balltree' method is second fastest with large data
    sets, but it is precise if working in unprojected coordinates like
    lat-lng. As a rule of thumb, if you have a small graph just use
    method=None. If you have a large graph with lat-lng coordinates, use
    method='balltree'. If you have a large graph with projected coordinates,
    use method='kdtree'. Note that if you are working in units of lat-lng,
    the X vector corresponds to longitude and the Y vector corresponds
    to latitude. The method creates equally distanced points along the edges
    of the network. Then, these points are used in a kdTree or BallTree search
    to identify which is nearest. Note that this method will not give exact
    perpendicular point along the edge, but the smaller the *dist* parameter,
    the closer (but slower) the solution will be.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    X : list-like
        the longitudes or x coordinates for which we will find the nearest
        edge in the graph. For projected graphs use the projected coordinates,
        usually in meters.
    Y : list-like
        the latitudes or y coordinates for which we will find the nearest edge
        in the graph. For projected graphs use the projected coordinates,
        usually in meters.
    method : string {None, 'kdtree', 'balltree'}
        Which method to use for finding nearest edge to each point. If None,
        we manually find each edge one at a time using get_nearest_edge. If
        'kdtree' we use scipy.spatial.cKDTree for very fast euclidean search.
        Recommended for projected graphs. If 'balltree', we use
        sklearn.neighbors.BallTree for fast haversine search. Recommended for
        unprojected graphs.
    dist : float
        spacing length along edges. Units are the same as the graph's
        geometries. The smaller the value, the more points are created.

    Returns
    -------
    ne : np.array
        array of edge IDs representing the edge nearest to each point in the
        passed-in list of points. Edge IDs are represented by u, v, key where
        u and v the node IDs of the nodes the edge links.
    """
    if method is None:
        # calculate nearest edge one at a time for each (y, x) point
        ne = [get_nearest_edge(G, (y, x)) for x, y in zip(X, Y)]

    elif method == "kdtree":

        # transform graph into DataFrame
        edges = utils_graph.graph_to_gdfs(G, nodes=False).reset_index()

        # transform edges into evenly spaced points
        edges["points"] = edges.apply(
            lambda x: utils_geo.redistribute_vertices(x.geometry, dist), axis=1
        )

        # develop edges data for each created points
        extended = (
            edges["points"]
            .apply([pd.Series])
            .stack()
            .reset_index(level=1, drop=True)
            .join(edges)
            .reset_index()
        )

        # Prepare btree arrays
        nbdata = np.array(
            list(
                zip(
                    extended["Series"].apply(lambda x: x.x), extended["Series"].apply(lambda x: x.y)
                )
            )
        )

        # build a k-d tree for euclidean nearest node search
        btree = cKDTree(data=nbdata, compact_nodes=True, balanced_tree=True)

        # query the tree for nearest node to each point
        points = np.array([X, Y]).T
        dist, idx = btree.query(points, k=1)  # Returns ids of closest point
        eidx = extended.loc[idx, "index"]
        ne = edges.loc[eidx, ["u", "v", "key"]]

    elif method == "balltree":

        # check if we were able to import sklearn.neighbors.BallTree successfully
        if BallTree is None:
            raise ImportError("scikit-learn must be installed to use this optional feature")

        # transform graph into DataFrame
        edges = utils_graph.graph_to_gdfs(G, nodes=False).reset_index()

        # transform edges into evenly spaced points
        edges["points"] = edges.apply(
            lambda x: utils_geo.redistribute_vertices(x.geometry, dist), axis=1
        )

        # develop edges data for each created points
        extended = (
            edges["points"]
            .apply([pd.Series])
            .stack()
            .reset_index(level=1, drop=True)
            .join(edges)
            .reset_index()
        )

        # haversine requires data in form of [lat, lng] and inputs/outputs in units of radians
        nodes = pd.DataFrame(
            {
                "x": extended["Series"].apply(lambda x: x.x),
                "y": extended["Series"].apply(lambda x: x.y),
            }
        )
        nodes_rad = np.deg2rad(nodes[["y", "x"]].values.astype(np.float))
        points = np.array([Y, X]).T
        points_rad = np.deg2rad(points)

        # build a ball tree for haversine nearest node search
        tree = BallTree(nodes_rad, metric="haversine")

        # query the tree for nearest node to each point
        idx = tree.query(points_rad, k=1, return_distance=False)
        eidx = extended.loc[idx[:, 0], "index"]
        ne = edges.loc[eidx, ["u", "v", "key"]]

    else:
        raise ValueError("You must pass a valid method name, or None.")

    utils.log(f"Found nearest edges to {len(X)} points")

    return np.array(ne)


def shortest_path(G, orig, dest, weight="length"):
    """
    Get shortest path from origin node to destination node.

    See also `k_shortest_paths` to get multiple shortest paths.

    This function is a convenience wrapper around networkx.shortest_path. For
    more functionality or different algorithms, use networkx directly.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    orig : int
        origin node ID
    dest : int
        destination node ID
    weight : string
        edge attribute to minimize when solving shortest path. default is edge
        length in meters.

    Returns
    -------
    path : list
        list of node IDs consituting the shortest path
    """
    path = nx.shortest_path(G, orig, dest, weight=weight)
    return path


def k_shortest_paths(G, orig, dest, k, weight="length"):
    """
    Get k shortest paths from origin node to destination node.

    See also `shortest_path` to get just the one shortest path.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    orig : int
        origin node ID
    dest : int
        destination node ID
    k : int
        number of shortest paths to get
    weight : string
        edge attribute to minimize when solving shortest paths. default is
        edge length in meters.

    Returns
    -------
    generator
        a generator of k shortest paths ordered by total weight. each path is
        a list of node IDs.
    """
    paths_gen = nx.shortest_simple_paths(utils_graph.get_digraph(G, weight), orig, dest, weight)
    for path in itertools.islice(paths_gen, 0, k):
        yield path
