"""Shared constants and builders for OSMnx tests."""

from __future__ import annotations

import json
from typing import TypeAlias

import matplotlib as mpl

mpl.use("Agg")

import networkx as nx
import requests
from shapely import LineString

import osmnx as ox

LOCATION_POINT = (37.791427, -122.410018)
ADDRESS = "Transamerica Pyramid, 600 Montgomery Street, San Francisco, California, USA"
PLACE = {"city": "Piedmont", "state": "California", "country": "USA"}
TAGS: dict[str, bool | str | list[str]] = {
    "landuse": True,
    "building": True,
    "highway": True,
    "amenity": True,
}

ResponseJson: TypeAlias = dict[str, object] | list[dict[str, object]]
HTTP_OK = 200
HTTP_ERROR = 500


def drive_graph() -> nx.MultiDiGraph:
    """
    Return the committed-cache drive graph used across offline tests.

    Returns
    -------
    nx.MultiDiGraph
        Cached fixture graph.
    """
    return ox.graph_from_point(
        LOCATION_POINT,
        dist=500,
        network_type="drive",
        simplify=False,
        retain_all=True,
    )


def toy_graph(*, crs: str = "epsg:4326") -> nx.MultiDiGraph:
    """
    Return a tiny graph with enough attrs for stats, routing, and plotting tests.

    Parameters
    ----------
    crs
        Graph CRS.

    Returns
    -------
    nx.MultiDiGraph
        Synthetic graph fixture.
    """
    G = nx.MultiDiGraph(crs=crs)
    G.add_node(1, x=0.0, y=0.0, street_count=1, elevation=0.0)
    G.add_node(2, x=1.0, y=0.0, street_count=2, elevation=10.0)
    G.add_node(3, x=2.0, y=0.0, street_count=1, elevation=20.0)
    G.add_edge(
        1,
        2,
        osmid=10,
        length=1.0,
        highway="residential",
        maxspeed="25 mph",
        geometry=LineString([(0, 0), (1, 0)]),
    )
    G.add_edge(
        2,
        3,
        osmid=11,
        length=1.0,
        highway=["primary", "secondary"],
        maxspeed=["30 mph", "50"],
        geometry=LineString([(1, 0), (2, 0)]),
    )
    return G


class Response(requests.Response):
    """
    Minimal response object for HTTP parser tests.

    Parameters
    ----------
    payload
        JSON payload to return.
    ok
        Whether the response should be treated as successful.
    status_code
        HTTP status code to expose.
    """

    def __init__(self, payload: ResponseJson, *, ok: bool = True, status_code: int = 200) -> None:
        """
        Instantiate a minimal response object.

        Parameters
        ----------
        payload
            JSON payload to return.
        ok
            Whether the response should be treated as successful.
        status_code
            HTTP status code to expose.
        """
        super().__init__()
        self._payload = payload
        self.status_code = status_code if ok or status_code != HTTP_OK else HTTP_ERROR
        self.reason = "OK" if ok else "Error"
        self._content = json.dumps(payload).encode()
        self.url = "https://example.com/api"

    def json(self, **kwargs: object) -> ResponseJson:
        """
        Return the configured JSON payload.

        Parameters
        ----------
        **kwargs
            Ignored JSON decoder keyword arguments.

        Returns
        -------
        ResponseJson
            The configured JSON payload.
        """
        del kwargs
        return self._payload
