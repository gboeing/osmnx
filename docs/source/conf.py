# ruff: noqa: INP001
"""
Configuration file for the Sphinx documentation builder.

For the full list of built-in configuration values, see the documentation:
https://www.sphinx-doc.org/en/master/usage/configuration.html
"""

import re
import sys
from pathlib import Path

# project info
author = "Geoff Boeing"
copyright = "2016-2025, Geoff Boeing"  # noqa: A001
project = "OSMnx"

# go up two levels from current working dir (/docs/source) to package root
pkg_root_path = str(Path.cwd().parent.parent)
sys.path.insert(0, pkg_root_path)

# dynamically load version from __version__ variable in /osmnx/_version.py
file_text = Path("../../osmnx/_version.py").read_text(encoding="utf-8")
regex = re.compile("^__version__ = ['\"]([^'\"]*)['\"]", re.MULTILINE)
match = re.search(regex, file_text)
if match:
    version = release = match.group(1)
else:
    msg = "Unable to find version string in file."
    raise ValueError(msg)

# mock import all required + optional dependency packages because readthedocs
# does not have them installed
autodoc_mock_imports = [
    "geopandas",
    "matplotlib",
    "networkx",
    "numpy",
    "pandas",
    "rasterio",
    "requests",
    "rio-vrt",
    "scipy",
    "shapely",
    "sklearn",
]

# linkcheck for stackoverflow gets HTTP 403 in CI environment
linkcheck_ignore = [r"https://stackoverflow\.com/.*"]

# type annotations configuration
autodoc_typehints = "description"
napoleon_use_param = True
napoleon_use_rtype = False
typehints_document_rtype = True
typehints_use_rtype = False
typehints_fully_qualified = False

# general configuration and options for HTML output
# see https://www.sphinx-doc.org/en/master/usage/configuration.html
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon", "sphinx_autodoc_typehints"]
html_static_path: list[str] = []
html_theme = "furo"
language = "en"
needs_sphinx = "7"  # match version from pyproject.toml optional-dependencies
root_doc = "index"
source_suffix = ".rst"
templates_path: list[str] = []
