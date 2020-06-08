"""Functions to calculate distances and find nearest node/edge(s) to point(s)."""

import networkx as nx
import numpy as np
import pandas as pd
from shapely.geometry import Point

from . import utils
from . import utils_geo
from . import utils_graph

# scipy and sklearn are optional dependencies for faster nearest node search
try:
    from scipy.spatial import cKDTree
except ImportError:
    cKDTree = None
try:
    from sklearn.neighbors import BallTree
except ImportError:
    BallTree = None


def great_circle_vec(lat1, lng1, lat2, lng2, earth_radius=6371009):
    """
    Calculate great-circle distances.

    Vectorized function to calculate the great-circle distance between two
    points or between vectors of points, using haversine.

    Parameters
    ----------
    lat1 : float or array of float
        first lat coord
    lng1 : float or array of float
        first lng coord
    lat2 : float or array of float
        second lat coord
    lng2 : float or array of float
        second lng coord
    earth_radius : numeric
        radius of earth in units in which distance will be returned
        (default is meters)

    Returns
    -------
    dist : float or np.array of floats
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
    dist = arc * earth_radius
    return dist


def euclidean_dist_vec(y1, x1, y2, x2):
    """
    Calculate euclidean distances.

    Vectorized function to calculate the euclidean distance between two points
    or between arrays of points.

    Parameters
    ----------
    y1 : float or np.array of float
        first y coord
    x1 : float or np.array of float
        first x coord
    y2 : float or np.array of float
        second y coord
    x2 : float or np.array of float
        second x coord

    Returns
    -------
    dist : float or np.array of float
        distance or vector of distances from (x1, y1) to (x2, y2) in graph
        units
    """
    # euclid's formula
    dist = ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5
    return dist


def get_nearest_node(G, point, method="haversine", return_dist=False):
    """
    Find node nearest to a point.

    Return the graph node nearest to some specified (lat, lng) or (y, x) point
    and optionally the distance between the node and the point. This function
    can use either a haversine or euclidean distance calculator.

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
    if len(G) < 1:
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
    Return the nearest edge to a point, by minimum euclidean distance.

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
    # get u, v, key, geom from all the graph edges
    gdf_edges = utils_graph.graph_to_gdfs(G, nodes=False, fill_edge_geometry=True)
    edges = gdf_edges[["u", "v", "key", "geometry"]].values

    # convert lat/lng point to x/y for shapely distance operation
    xy_point = Point(reversed(point))

    # calculate euclidean distance from each edge's geometry to this point
    edge_distances = [(edge, xy_point.distance(edge[3])) for edge in edges]

    # the nearest edge minimizes the distance to the point
    (u, v, key, geom), dist = min(edge_distances, key=lambda x: x[1])
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


def get_nearest_nodes(G, X, Y, method=None):
    """
    Return the graph nodes nearest to a list of points.

    Pass in points as separate vectors of X and Y coordinates. The 'kdtree'
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
        The vector of longitudes or x's for which we will find the nearest
        node in the graph
    Y : list-like
        The vector of latitudes or y's for which we will find the nearest
        node in the graph
    method : string {None, 'kdtree', 'balltree'}
        Which method to use for finding nearest node to each point.
        If None, we manually find each node one at a time using
        utils.get_nearest_node and haversine. If 'kdtree' we use
        scipy.spatial.cKDTree for very fast euclidean search. If
        'balltree', we use sklearn.neighbors.BallTree for fast
        haversine search.

    Returns
    -------
    nn : np.array
        list of nearest node IDs
    """
    if method is None:

        # calculate nearest node one at a time for each point
        nn = [get_nearest_node(G, (y, x), method="haversine") for x, y in zip(X, Y)]

    elif method == "kdtree":

        # check if we were able to import scipy.spatial.cKDTree
        if not cKDTree:
            raise ImportError("The scipy package must be installed to use this optional feature.")

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
        if not BallTree:
            raise ImportError(
                "The scikit-learn package must be installed to use this optional feature."
            )

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
        idx = tree.query(points_rad, k=1, return_distance=False)
        nn = nodes.iloc[idx[:, 0]].index

    else:
        raise ValueError("You must pass a valid method name, or None.")

    utils.log(f"Found nearest nodes to {len(X)} points")
    return np.array(nn)


def get_nearest_edges(G, X, Y, method=None, dist=0.0001):
    """
    Return the graph edges nearest to a list of points.

    Pass in points as separate vectors of X and Y coordinates. The 'kdtree'
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
    to identify which is nearest.Note that this method will not give the exact
    perpendicular point along the edge, but the smaller the *dist* parameter,
    the closer the solution will be.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    X : list-like
        The vector of longitudes or x's for which we will find the nearest
        edge in the graph. For projected graphs, use the projected coordinates,
        usually in meters.
    Y : list-like
        The vector of latitudes or y's for which we will find the nearest
        edge in the graph. For projected graphs, use the projected coordinates,
        usually in meters.
    method : string {None, 'kdtree', 'balltree'}
        Which method to use for finding nearest edge to each point.
        If None, we manually find each edge one at a time using
        get_nearest_edge. If 'kdtree' we use
        scipy.spatial.cKDTree for very fast euclidean search. Recommended for
        projected graphs. If 'balltree', we use sklearn.neighbors.BallTree for
        fast haversine search. Recommended for unprojected graphs.

    dist : float
        spacing length along edges. Units are the same as the geom; Degrees for
        unprojected geometries and meters for projected geometries. The smaller
        the value, the more points are created.

    Returns
    -------
    ne : np.array
        array of nearest edges represented by u and v (the IDs of the nodes
        they link) and key
    """
    if method is None:
        # calculate nearest edge one at a time for each (y, x) point
        ne = [get_nearest_edge(G, (y, x)) for x, y in zip(X, Y)]

    elif method == "kdtree":

        # check if we were able to import scipy.spatial.cKDTree successfully
        if not cKDTree:
            raise ImportError("The scipy package must be installed to use this optional feature.")

        # transform graph into DataFrame
        edges = utils_graph.graph_to_gdfs(G, nodes=False, fill_edge_geometry=True)

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
        if not BallTree:
            raise ImportError(
                "The scikit-learn package must be installed to use this optional feature."
            )

        # transform graph into DataFrame
        edges = utils_graph.graph_to_gdfs(G, nodes=False, fill_edge_geometry=True)

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
