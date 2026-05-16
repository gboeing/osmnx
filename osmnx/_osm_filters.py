"""Reusable OSM network tag filters for Overpass queries and local data."""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from typing import Final
from typing import Literal
from typing import TypeAlias

from . import settings

_FilterOperator: TypeAlias = Literal["exists", "=", "!=", "~", "!~"]
_Tags: TypeAlias = Mapping[str, str]

_FILTER_OPERATORS: Final[dict[str, _FilterOperator]] = {
    "=": "=",
    "!=": "!=",
    "~": "~",
    "!~": "!~",
}
_OVERPASS_FILTER_PATTERN: Final = re.compile(
    r"""
    \[
    (?:"(?P<quoted_key>[^"]+)"|(?P<unquoted_key>[A-Za-z0-9_:-]+))
    (?:
        (?P<operator>!=|!~|=|~)
        (?:"(?P<quoted_value>[^"]*)"|(?P<unquoted_value>[^\]\s]+))
    )?
    \]
    """,
    re.VERBOSE,
)


@dataclass(frozen=True)
class _TagFilter:
    """
    Represent one OSM tag predicate.

    This models the subset of Overpass QL tag filters used by OSMnx's built-in
    network filters: key existence, exact equality/inequality, and regular
    expression matching/non-matching. Keeping this representation separate from
    the Overpass string lets local sources, such as PBF files, evaluate the
    same predicates against a way's tags without parsing query text.

    Parameters
    ----------
    key
        OSM tag key to inspect.
    operator
        Predicate operator. The default, "exists", means the key must be
        present regardless of its value.
    value
        Value or regular expression pattern for all non-existence operators.
    """

    key: str
    operator: _FilterOperator = "exists"
    value: str | None = None

    def matches(self, tags: _Tags) -> bool:
        """
        Return True if this predicate matches a way's tags.

        Overpass negative filters such as ["access"!~"private"] match elements
        where the key is absent as well as elements where the key's value does
        not match the expression. The local evaluator mirrors that behavior so
        PBF filtering follows the same semantics as Overpass queries.
        """
        if self.operator == "exists":
            return self.key in tags

        if self.value is None:  # pragma: no cover
            msg = "A tag filter value is required unless operator='exists'."
            raise ValueError(msg)

        tag_value = tags.get(self.key)

        if self.operator == "=":
            return tag_value == self.value
        if self.operator == "!=":
            return tag_value != self.value
        if tag_value is None:
            return self.operator == "!~"

        regex_matches = re.search(self.value, tag_value) is not None
        if self.operator == "~":
            return regex_matches
        return not regex_matches

    def to_overpass(self) -> str:
        """
        Convert this predicate to its Overpass QL tag filter form.

        The returned string is intended to be concatenated after an Overpass
        element selector, for example `way["highway"]["area"!~"yes"]`.
        """
        if self.operator == "exists":
            return f'["{self.key}"]'
        return f'["{self.key}"{self.operator}"{self.value}"]'


@dataclass(frozen=True)
class _DefaultAccessFilter:
    """
    Represent the configurable default access predicate.

    `settings.default_access` is a user-configurable Overpass filter string, so
    it cannot be frozen into `_NETWORK_FILTERS` at import time. This lightweight
    sentinel marks where the current setting should be compiled into a concrete
    `_TagFilter` when building a `_WayFilter`.
    """


_FilterClause: TypeAlias = _TagFilter | _DefaultAccessFilter


@dataclass(frozen=True)
class _WayFilter:
    """
    Represent a complete way filter as ORs of ANDed tag predicates.

    Each inner tuple is an intersection of tag predicates. The outer tuple is
    their union. This matches OSMnx's existing `custom_filter` semantics: a
    string is one intersection, while a list of strings is a union.

    Parameters
    ----------
    clauses
        Groups of tag filters to match or render.
    overpass_filters
        Raw Overpass QL filters to preserve for rendering only.
    """

    clauses: tuple[tuple[_TagFilter, ...], ...]
    overpass_filters: tuple[str, ...] | None = None

    def matches(self, tags: _Tags) -> bool:
        """Return True if any filter clause group matches tags."""
        if self.overpass_filters is not None:
            msg = "Cannot locally match unparsed custom Overpass filters."
            raise ValueError(msg)
        return any(all(filter_.matches(tags) for filter_ in group) for group in self.clauses)

    def to_overpass_filters(self) -> list[str]:
        """Render each clause group as an Overpass QL tag filter string."""
        if self.overpass_filters is not None:
            return list(self.overpass_filters)
        return ["".join(filter_.to_overpass() for filter_ in group) for group in self.clauses]

_DEFAULT_ACCESS_FILTER = _DefaultAccessFilter()
_HIGHWAY_EXISTS = _TagFilter("highway")
_AREA_NOT_YES = _TagFilter("area", "!~", "yes")

# Each network type is an intersection of tag predicates: a way must satisfy
# every clause in the tuple to be included. The tuple order is preserved when
# rendering Overpass QL so existing query strings remain stable.
_NETWORK_FILTERS: Final[dict[str, tuple[_FilterClause, ...]]] = {
    "drive": (
        _HIGHWAY_EXISTS,
        _AREA_NOT_YES,
        _DEFAULT_ACCESS_FILTER,
        _TagFilter(
            "highway",
            "!~",
            "abandoned|bridleway|bus_guideway|construction|corridor|cycleway|elevator|"
            "escalator|footway|no|path|pedestrian|planned|platform|proposed|raceway|"
            "razed|rest_area|service|services|steps|track",
        ),
        _TagFilter("motor_vehicle", "!~", "no"),
        _TagFilter("motorcar", "!~", "no"),
        _TagFilter(
            "service",
            "!~",
            "alley|driveway|emergency_access|parking|parking_aisle|private",
        ),
    ),
    "drive_service": (
        _HIGHWAY_EXISTS,
        _AREA_NOT_YES,
        _DEFAULT_ACCESS_FILTER,
        _TagFilter(
            "highway",
            "!~",
            "abandoned|bridleway|bus_guideway|construction|corridor|cycleway|elevator|"
            "escalator|footway|no|path|pedestrian|planned|platform|proposed|raceway|"
            "razed|rest_area|services|steps|track",
        ),
        _TagFilter("motor_vehicle", "!~", "no"),
        _TagFilter("motorcar", "!~", "no"),
        _TagFilter("service", "!~", "emergency_access|parking|parking_aisle|private"),
    ),
    "walk": (
        _HIGHWAY_EXISTS,
        _AREA_NOT_YES,
        _DEFAULT_ACCESS_FILTER,
        _TagFilter(
            "highway",
            "!~",
            "abandoned|bus_guideway|construction|cycleway|motor|no|planned|platform|"
            "proposed|raceway|razed|rest_area|services",
        ),
        _TagFilter("foot", "!~", "no"),
        _TagFilter("service", "!~", "private"),
        _TagFilter("sidewalk", "!~", "separate"),
        _TagFilter("sidewalk:both", "!~", "separate"),
        _TagFilter("sidewalk:left", "!~", "separate"),
        _TagFilter("sidewalk:right", "!~", "separate"),
    ),
    "bike": (
        _HIGHWAY_EXISTS,
        _AREA_NOT_YES,
        _DEFAULT_ACCESS_FILTER,
        _TagFilter(
            "highway",
            "!~",
            "abandoned|bus_guideway|construction|corridor|elevator|escalator|footway|"
            "motor|no|planned|platform|proposed|raceway|razed|rest_area|services|steps",
        ),
        _TagFilter("bicycle", "!~", "no"),
        _TagFilter("service", "!~", "private"),
    ),
    "all_public": (
        _HIGHWAY_EXISTS,
        _AREA_NOT_YES,
        _DEFAULT_ACCESS_FILTER,
        _TagFilter(
            "highway",
            "!~",
            "abandoned|construction|no|planned|platform|proposed|raceway|razed|rest_area|services",
        ),
        _TagFilter("service", "!~", "private"),
    ),
    "all": (
        _HIGHWAY_EXISTS,
        _AREA_NOT_YES,
        _TagFilter(
            "highway",
            "!~",
            "abandoned|construction|no|planned|platform|proposed|raceway|razed|rest_area|services",
        ),
    ),
}


def _get_way_filter(
    network_type: str,
    custom_filter: str | list[str] | None,
    *,
    parse_custom_filter: bool = True,
) -> _WayFilter:
    """
    Create a reusable way filter from a network type or custom filter.

    Parameters
    ----------
    network_type
        {"all", "all_public", "bike", "drive", "drive_service", "walk"}
        What type of street network to retrieve if `custom_filter` is None.
    custom_filter
        A custom ways filter to use instead of the `network_type` presets. If
        a string, the intersection of its tag filters will be used. If a list,
        the union of the list items will be used.
    parse_custom_filter
        If True, parse `custom_filter` into locally matchable predicates. If
        False, preserve it as raw Overpass QL for rendering only.

    Returns
    -------
    way_filter
        Reusable way filter for rendering Overpass QL or matching local tags.
    """
    if isinstance(custom_filter, list):
        if parse_custom_filter:
            return _WayFilter(tuple(_parse_overpass_filter(filter_) for filter_ in custom_filter))
        return _WayFilter((), overpass_filters=tuple(custom_filter))
    if isinstance(custom_filter, str):
        if parse_custom_filter:
            return _WayFilter((_parse_overpass_filter(custom_filter),))
        return _WayFilter((), overpass_filters=(custom_filter,))

    return _WayFilter((_get_network_filter_clauses(network_type),))


def _get_network_filter(network_type: str) -> str:
    """
    Create an Overpass QL filter for the specified network type.

    Parameters
    ----------
    network_type
        {"all", "all_public", "bike", "drive", "drive_service", "walk"}
        What type of street network to retrieve.

    Returns
    -------
    way_filter
        The Overpass query filter.
    """
    return _get_way_filter(network_type, custom_filter=None).to_overpass_filters()[0]


def _network_filter_matches(network_type: str, tags: _Tags) -> bool:
    """
    Return True if tags match the specified network type's OSM tag filters.

    Parameters
    ----------
    network_type
        {"all", "all_public", "bike", "drive", "drive_service", "walk"}
        What type of street network to match.
    tags
        The way's OSM tags.

    Returns
    -------
    matches
        True if the tags match this network type's filter.
    """
    return _get_way_filter(network_type, custom_filter=None).matches(tags)


def _get_network_filter_clauses(network_type: str) -> tuple[_TagFilter, ...]:
    """
    Compile a network type into concrete tag filters.

    This resolves `settings.default_access` once per `_WayFilter` construction
    so local PBF scans do not repeatedly parse it for every way.

    Parameters
    ----------
    network_type
        {"all", "all_public", "bike", "drive", "drive_service", "walk"}
        What type of street network to retrieve.

    Returns
    -------
    tag_filters
        Concrete tag filters.
    """
    if network_type not in _NETWORK_FILTERS:  # pragma: no cover
        msg = f"Unrecognized network_type {network_type!r}."
        raise ValueError(msg)

    tag_filters: list[_TagFilter] = []
    for filter_ in _NETWORK_FILTERS[network_type]:
        if isinstance(filter_, _DefaultAccessFilter):
            tag_filters.extend(_parse_overpass_filter(settings.default_access))
        else:
            tag_filters.append(filter_)

    return tuple(tag_filters)


@cache
def _parse_overpass_filter(overpass_filter: str) -> tuple[_TagFilter, ...]:
    """
    Parse a simple Overpass QL tag filter string into reusable tag filters.

    This parser intentionally accepts only the compact tag-filter subset used
    by `settings.default_access` and the built-in network filters. Unsupported
    syntax raises ValueError instead of guessing at broader Overpass QL
    semantics.

    Parameters
    ----------
    overpass_filter
        Overpass QL tag filter string.

    Returns
    -------
    tag_filters
        Parsed tag filters.
    """
    position = 0
    filters: list[_TagFilter] = []
    for match in _OVERPASS_FILTER_PATTERN.finditer(overpass_filter):
        if overpass_filter[position : match.start()].strip():
            msg = f"Unsupported Overpass filter syntax: {overpass_filter!r}."
            raise ValueError(msg)

        key = match.group("quoted_key") or match.group("unquoted_key")
        if key is None:  # pragma: no cover
            msg = f"Overpass filter is missing a key: {overpass_filter!r}."
            raise ValueError(msg)

        quoted_value = match.group("quoted_value")
        value = quoted_value if quoted_value is not None else match.group("unquoted_value")
        filters.append(
            _make_tag_filter(key=key, operator=match.group("operator"), value=value),
        )
        position = match.end()

    if overpass_filter[position:].strip():
        msg = f"Unsupported Overpass filter syntax: {overpass_filter!r}."
        raise ValueError(msg)

    return tuple(filters)


def _make_tag_filter(key: str, operator: str | None, value: str | None) -> _TagFilter:
    """
    Create a `_TagFilter` from parsed Overpass QL components.

    Parameters
    ----------
    key
        OSM tag key.
    operator
        Parsed Overpass QL operator, or None for a key-existence filter.
    value
        Parsed value or regular expression pattern.

    Returns
    -------
    tag_filter
        A reusable tag filter.
    """
    if operator is None:
        return _TagFilter(key)

    if value is None:  # pragma: no cover
        msg = f"Overpass filter for {key!r} is missing a value."
        raise ValueError(msg)

    filter_operator = _validate_filter_operator(operator)
    return _TagFilter(key=key, operator=filter_operator, value=value)


def _validate_filter_operator(operator: str) -> _FilterOperator:
    """
    Validate and narrow a parsed Overpass operator.

    Parameters
    ----------
    operator
        Parsed operator string.

    Returns
    -------
    filter_operator
        Operator narrowed to the internal `_FilterOperator` type.
    """
    if operator in _FILTER_OPERATORS:
        return _FILTER_OPERATORS[operator]

    msg = f"Unsupported Overpass filter operator {operator!r}."
    raise ValueError(msg)
