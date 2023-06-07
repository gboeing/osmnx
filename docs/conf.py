"""
Configuration file for the Sphinx documentation builder.

For the full list of built-in configuration values, see the documentation:
https://www.sphinx-doc.org/en/master/usage/configuration.html
"""

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import sys
from pathlib import Path

# go up one level from /docs/ to the package root
sys.path.insert(0, str(Path().resolve().parent))

author = "Geoff Boeing"
copyright = "2016–2023, Geoff Boeing"
project = "OSMnx"
version = release = "1.4.0-alpha"

# mock import these packages because readthedocs doesn't have them installed
autodoc_mock_imports = [
    "geopandas",
    "matplotlib",
    "networkx",
    "numpy",
    "osgeo",
    "pandas",
    "pyproj",
    "rasterio",
    "requests",
    "scipy",
    "shapely",
    "sklearn",
]

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon"]
language = "en"
needs_sphinx = "6"
root_doc = "index"
source_suffix = ".rst"
templates_path = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_static_path = []
html_theme = "sphinx_rtd_theme"
