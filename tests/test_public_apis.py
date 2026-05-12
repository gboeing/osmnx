# ruff: noqa: PLR2004, S101
# numpydoc ignore=PR01,RT01
"""Live public web API compatibility tests."""

from __future__ import annotations

import networkx as nx
import pytest

import osmnx as ox

pytestmark = pytest.mark.network

LOCATION_POINT = (37.791427, -122.410018)
ADDRESS = "Transamerica Pyramid, 600 Montgomery Street, San Francisco, California, USA"
PLACE = {"city": "Piedmont", "state": "California", "country": "USA"}
TAGS: dict[str, bool | str | list[str]] = {"building": True, "amenity": True, "highway": "bus_stop"}


@pytest.fixture(autouse=True)
def _configure_live_api_tests() -> None:
    """Disable cache so these tests exercise current public API behavior."""
    ox.settings.use_cache = False
    ox.settings.overpass_rate_limit = True
    ox.settings.requests_timeout = 180


def _assert_valid_graph(G: nx.MultiDiGraph) -> None:
    ox.convert.validate_graph(G)
    assert len(G) > 0
    assert len(G.edges) > 0
    assert G.graph["crs"] == ox.settings.default_crs


def test_live_nominatim_downloaders() -> None:
    """Smoke-test public Nominatim geocoding endpoints."""
    point = ox.geocode(ADDRESS)
    gdf_place = ox.geocode_to_gdf(PLACE, which_result=1)
    gdf_osmid = ox.geocode_to_gdf("R2999176", by_osmid=True)

    assert len(point) == 2
    assert len(gdf_place) == 1
    assert len(gdf_osmid) == 1
    assert gdf_place.crs == ox.settings.default_crs
    assert not gdf_place.geometry.is_empty.any()


def test_live_overpass_graph_downloaders() -> None:
    """Smoke-test public Overpass graph downloader entry points."""
    bbox = ox.utils_geo.bbox_from_point(LOCATION_POINT, dist=250)
    polygon = ox.geocode_to_gdf(PLACE, which_result=1).geometry.iloc[0]

    graphs = [
        ox.graph_from_bbox(bbox, network_type="drive"),
        ox.graph_from_point(LOCATION_POINT, dist=250, network_type="drive"),
        ox.graph_from_address(ADDRESS, dist=250, network_type="drive"),
        ox.graph_from_place(PLACE, network_type="drive", which_result=1),
        ox.graph_from_polygon(polygon, network_type="drive"),
    ]

    for G in graphs:
        _assert_valid_graph(G)


def test_live_overpass_feature_downloaders() -> None:
    """Smoke-test public Overpass feature downloader entry points."""
    bbox = ox.utils_geo.bbox_from_point(LOCATION_POINT, dist=250)
    polygon = ox.geocode_to_gdf(PLACE, which_result=1).geometry.iloc[0]

    gdfs = [
        ox.features_from_bbox(bbox, tags=TAGS),
        ox.features_from_point(LOCATION_POINT, tags=TAGS, dist=250),
        ox.features_from_address(ADDRESS, tags=TAGS, dist=250),
        ox.features_from_place(PLACE, tags=TAGS, which_result=1),
        ox.features_from_polygon(polygon, tags=TAGS),
    ]

    for gdf in gdfs:
        ox.convert.validate_features_gdf(gdf)
        assert len(gdf) > 0
        assert gdf.crs == ox.settings.default_crs


def test_live_elevation_downloader() -> None:
    """Smoke-test a public Google-compatible elevation endpoint."""
    G = nx.MultiDiGraph(crs=ox.settings.default_crs)
    G.add_node(1, x=-122.410018, y=37.791427, street_count=1)
    G.add_node(2, x=-122.409018, y=37.792427, street_count=1)
    G.add_edge(1, 2, osmid=1, length=100.0)

    ox.settings.elevation_url_template = (
        "https://api.opentopodata.org/v1/aster30m?locations={locations}&key={key}"
    )
    G = ox.elevation.add_node_elevations_google(G, batch_size=100, pause=1)

    elevations = dict(G.nodes(data="elevation"))
    assert set(elevations) == {1, 2}
    assert all(isinstance(value, float) for value in elevations.values())
