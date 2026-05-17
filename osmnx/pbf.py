"""
Create graphs from local OpenStreetMap PBF files.

This module requires the optional ``osmium`` dependency. Install it with
``osmnx[pbf]``.
"""

from __future__ import annotations

from .graph import graph_from_pbf

__all__ = ["graph_from_pbf"]
