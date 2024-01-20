"""
Configuration file for the Sphinx documentation builder.

For the full list of built-in configuration values, see the documentation:
https://www.sphinx-doc.org/en/master/usage/configuration.html
"""

import sys
from pathlib import Path

# project info
author = "Geoff Boeing"
copyright = "2016-2024, Geoff Boeing"  # noqa: A001
project = "OSMnx"

# go up two levels from current working dir (/docs/source) to package root
pkg_root_path = str(Path.cwd().parent.parent)
sys.path.insert(0, pkg_root_path)

# dynamically load version from /osmnx/_version.py
with Path.open(Path("../../osmnx/_version.py")) as f:
    version = release = f.read().split(" = ")[1].replace('"', "")

# mock import all required + optional dependency packages because readthedocs
# does not have them installed
autodoc_mock_imports = [
    "geopandas",
    "matplotlib",
    "networkx",
    "numpy",
    "osgeo",
    "pandas",
    "rasterio",
    "requests",
    "scipy",
    "shapely",
    "sklearn",
]

autodoc_typehints = "none"

# general configuration and options for HTML output
# see https://www.sphinx-doc.org/en/master/usage/configuration.html
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon"]
html_static_path: list[str] = []
html_theme = "furo"
language = "en"
needs_sphinx = "7"  # same value as pinned in /docs/requirements.txt
root_doc = "index"
source_suffix = ".rst"
templates_path: list[str] = []
