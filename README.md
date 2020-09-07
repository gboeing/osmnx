[![PyPI version](https://badge.fury.io/py/osmnx.svg)](https://badge.fury.io/py/osmnx)
[![PyPI Downloads](https://img.shields.io/pypi/dm/osmnx.svg)](https://badge.fury.io/py/osmnx)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/osmnx/badges/downloads.svg)](https://anaconda.org/conda-forge/osmnx)
[![Documentation Status](https://readthedocs.org/projects/osmnx/badge/?version=latest)](https://osmnx.readthedocs.io/)
[![Build Status](https://travis-ci.org/gboeing/osmnx.svg?branch=master)](https://travis-ci.org/gboeing/osmnx)
[![Coverage Status](https://coveralls.io/repos/github/gboeing/osmnx/badge.svg?branch=master)](https://coveralls.io/github/gboeing/osmnx?branch=master)


# OSMnx

**Python for street networks**

Retrieve, model, analyze, and visualize OpenStreetMap street networks and other spatial data.

**Citation info**: Boeing, G. 2017. "[OSMnx: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks](https://geoffboeing.com/publications/osmnx-complex-street-networks/)." *Computers, Environment and Urban Systems* 65, 126-139. doi:10.1016/j.compenvurbsys.2017.05.004



## Features

**OSMnx** is a Python package that lets you download spatial geometries and model, project, visualize, and analyze real-world street networks from OpenStreetMap's APIs. Users can download and model walkable, drivable, or bikeable urban networks with a single line of Python code, and then easily analyze and visualize them. You can just as easily download and work with amenities/points of interest, building footprints, elevation data, street bearings/orientations, speed/travel time, and network routing.

OSMnx is built on top of geopandas, networkx, and matplotlib and interacts with OpenStreetMap's APIs to:

  * Download and model street networks or other networked infrastructure anywhere in the world with a single line of code
  * Download any other spatial geometries, place boundaries, building footprints, or points of interest as a GeoDataFrame
  * Download by city name, polygon, bounding box, or point/address + network distance
  * Download drivable, walkable, bikeable, or all street networks
  * Download node elevations and calculate edge grades (inclines)
  * Impute missing speeds and calculate graph edge travel times
  * Simplify and correct the network's topology to clean-up nodes and consolidate intersections
  * Fast map-matching of points, routes, or trajectories to nearest graph edges or nodes
  * Save networks to disk as shapefiles, geopackages, and GraphML
  * Save/load street network to/from a local .osm xml file
  * Conduct topological and spatial analyses to automatically calculate dozens of indicators
  * Calculate and visualize street bearings and orientations
  * Calculate and visualize shortest-path routes that minimize distance, travel time, elevation, etc
  * Visualize street networks as a static map or interactive leaflet web map
  * Visualize travel distance and travel time with isoline and isochrone maps
  * Plot figure-ground diagrams of street networks and building footprints

Examples and demonstrations of these features are in the [examples repo](https://github.com/gboeing/osmnx-examples). More feature development details are in the change log.



## Installation

If you have any trouble with the installation, read the [docs](https://osmnx.readthedocs.io/) for more info.

Install OSMnx in a clean [conda](https://anaconda.org/conda-forge/osmnx) environment:

```
conda config --prepend channels conda-forge
conda create -n ox --strict-channel-priority osmnx
```

Alternatively, you can run OSMnx + Jupyter directly from its official [docker container](https://hub.docker.com/r/gboeing/osmnx).



## Documentation and Usage

Documentation available on [readthedocs](https://osmnx.readthedocs.io/).

"How do I use OSMnx?" Usage examples and tutorials available in the [examples repo](https://github.com/gboeing/osmnx-examples).

Examples of projects and blog posts [using OSMnx](https://geoffboeing.com/2018/03/osmnx-features-roundup/).

If you use OSMnx in your work, please cite the [journal article](https://geoffboeing.com/publications/osmnx-complex-street-networks/).
