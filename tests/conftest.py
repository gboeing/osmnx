# ruff: noqa: PLC0415, TC003, NPY002
"""Shared pytest fixtures for deterministic OSMnx tests."""

from __future__ import annotations

import os
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import requests

HTTP_CACHE_DIR = Path(__file__).parent / "input_data" / "http_cache"


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """
    Skip live network tests unless explicitly requested.

    Parameters
    ----------
    config : pytest.Config
        Pytest configuration object (unused).
    items : list of pytest.Item
        Collected test items.
    """
    del config

    if os.environ.get("OSMNX_RUN_NETWORK_TESTS"):
        return

    skip_network = pytest.mark.skip(reason="set OSMNX_RUN_NETWORK_TESTS=1 to run network tests")
    for item in items:
        if item.get_closest_marker("network"):
            item.add_marker(skip_network)


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
    import osmnx as ox

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

    np.random.seed(0)

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
    if request.node.get_closest_marker("network") and os.environ.get("OSMNX_RUN_NETWORK_TESTS"):
        return

    def _blocked_request(*args: Any, **kwargs: Any) -> None:  # noqa: ANN401, ARG001
        msg = (
            "Network access is blocked in offline tests. Mark the test with "
            "`@pytest.mark.network` and set OSMNX_RUN_NETWORK_TESTS=1 to allow it."
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
    import osmnx as ox

    ox.settings.cache_folder = HTTP_CACHE_DIR
    ox.settings.use_cache = True
    ox.settings.overpass_rate_limit = False
    return HTTP_CACHE_DIR
