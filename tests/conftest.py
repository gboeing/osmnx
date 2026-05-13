#!/usr/bin/env python
"""Shared pytest fixtures for deterministic OSMnx tests."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import TYPE_CHECKING
from typing import Any
from typing import TypeAlias

import matplotlib as mpl
import networkx as nx
import pytest
import requests
from shapely import LineString

mpl.use("Agg")

if TYPE_CHECKING:
    from collections.abc import Iterator

# Shared test data and builders.
LOCATION_POINT = (37.791427, -122.410018)
ADDRESS = "Transamerica Pyramid, 600 Montgomery Street, San Francisco, California, USA"
PLACE = {"city": "Piedmont", "state": "California", "country": "USA"}
TAGS: dict[str, bool | str | list[str]] = {
    "landuse": True,
    "building": True,
    "highway": True,
    "amenity": True,
}

_ResponseJson: TypeAlias = dict[str, object] | list[dict[str, object]]
HTTP_OK = 200
HTTP_ERROR = 500


def _drive_graph() -> nx.MultiDiGraph:
    import osmnx as ox  # noqa: PLC0415

    return ox.graph_from_point(
        LOCATION_POINT,
        dist=500,
        network_type="drive",
        simplify=False,
        retain_all=True,
    )


def _toy_graph(*, crs: str = "epsg:4326") -> nx.MultiDiGraph:
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


class _Response(requests.Response):
    def __init__(self, payload: _ResponseJson, *, ok: bool = True, status_code: int = 200) -> None:
        super().__init__()
        self._payload = payload
        self.status_code = status_code if ok or status_code != HTTP_OK else HTTP_ERROR
        self.reason = "OK" if ok else "Error"
        self._content = json.dumps(payload).encode()
        self.url = "https://example.com/api"

    def json(self, **kwargs: object) -> _ResponseJson:
        del kwargs
        return self._payload


HTTP_CACHE_DIR = Path(__file__).parent / "input_data" / "http_cache"


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """
    Skip live online tests unless explicitly requested.

    Parameters
    ----------
    config : pytest.Config
        Pytest configuration object (unused).
    items : list of pytest.Item
        Collected test items.
    """
    del config

    if os.environ.get("OSMNX_RUN_ONLINE_TESTS"):
        return

    skip_online = pytest.mark.skip(reason="set OSMNX_RUN_ONLINE_TESTS=1 to run online tests")
    for item in items:
        if item.get_closest_marker("online"):
            item.add_marker(skip_online)


@pytest.fixture(autouse=True)
def _isolate_settings(tmp_path: Path) -> Iterator[None]:
    """
    Restore global settings after each test and isolate generated files.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory unique to the test invocation.

    Yields
    ------
    None
        Control is yielded to the test, after which original settings are restored.
    """
    import osmnx as ox  # noqa: PLC0415

    original_settings = {
        name: getattr(ox.settings, name)
        for name in dir(ox.settings)
        if not name.startswith("_") and name.islower()
    }

    ox.settings.data_folder = tmp_path / "data"
    ox.settings.logs_folder = tmp_path / "logs"
    ox.settings.imgs_folder = tmp_path / "imgs"
    ox.settings.cache_folder = tmp_path / "cache"
    ox.settings.log_console = False
    ox.settings.log_file = False
    ox.settings.use_cache = True

    yield

    for name, value in original_settings.items():
        setattr(ox.settings, name, value)


@pytest.fixture(autouse=True)
def _block_network(
    monkeypatch: pytest.MonkeyPatch,
    request: pytest.FixtureRequest,
) -> None:
    """
    Prevent accidental live HTTP calls in the default offline suite.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture for safely modifying objects during tests.
    request : pytest.FixtureRequest
        Provides access to the requesting test context.

    Returns
    -------
    None
        This fixture modifies global state and does not return a value.
    """
    if request.node.get_closest_marker("online") and os.environ.get("OSMNX_RUN_ONLINE_TESTS"):
        return

    def _blocked_request(*args: Any, **kwargs: Any) -> None:  # noqa: ANN401, ARG001
        msg = (
            "Network access is blocked in offline tests. Mark the test with "
            "`@pytest.mark.online` and set OSMNX_RUN_ONLINE_TESTS=1 to allow it."
        )
        raise AssertionError(msg)

    monkeypatch.setattr(requests, "get", _blocked_request)
    monkeypatch.setattr(requests, "post", _blocked_request)


@pytest.fixture
def http_cache() -> Path:
    """
    Point OSMnx cache lookups at committed raw HTTP response fixtures.

    Returns
    -------
    pathlib.Path
        Path to the HTTP cache directory used for tests.
    """
    import osmnx as ox  # noqa: PLC0415

    ox.settings.cache_folder = HTTP_CACHE_DIR
    ox.settings.use_cache = True
    ox.settings.overpass_rate_limit = False
    return HTTP_CACHE_DIR
