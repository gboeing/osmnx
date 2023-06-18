[![PyPI Version](https://badge.fury.io/py/osmnx.svg)](https://badge.fury.io/py/osmnx)
[![PyPI Downloads](https://static.pepy.tech/personalized-badge/osmnx?period=total&units=international_system&left_color=grey&right_color=brightgreen&left_text=downloads)](https://pepy.tech/project/osmnx)
[![Anaconda Downloads](https://anaconda.org/conda-forge/osmnx/badges/downloads.svg)](https://anaconda.org/conda-forge/osmnx)
[![Documentation Status](https://readthedocs.org/projects/osmnx/badge/?version=latest)](https://osmnx.readthedocs.io/)
[![Build Status](https://github.com/gboeing/osmnx/workflows/tests/badge.svg?branch=main)](https://github.com/gboeing/osmnx/actions?query=workflow%3A%22tests%22)
[![Coverage Status](https://codecov.io/gh/gboeing/osmnx/branch/main/graph/badge.svg)](https://codecov.io/gh/gboeing/osmnx)

# OSMnx

**OSMnx** is a Python package that lets you download geospatial data from OpenStreetMap and model, project, visualize, and analyze real-world street networks and any other geospatial geometries. You can download and model walking, driving, or biking networks with a single line of code then easily analyze and visualize them. You can just as easily download and work with other infrastructure types, amenities/points of interest, building footprints, elevation data, street bearings/orientations, speed/travel time, and routing.

If you use OSMnx in your work, please cite the journal article.

**Citation info**: Boeing, G. 2017. "[OSMnx: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks](https://geoffboeing.com/publications/osmnx-complex-street-networks/)." *Computers, Environment and Urban Systems* 65, 126–139.


## Getting Started

*How do I install OSMnx?* Follow the [installation](https://osmnx.readthedocs.io/en/stable/#installation) instructions.

*How do I use OSMnx?* Check out the usage examples/tutorials in the [examples](https://github.com/gboeing/osmnx-examples) repo.

*How does this or that function work?* Read its [documentation](https://osmnx.readthedocs.io/).

*What can I do with OSMnx?* Check out recent [projects](https://geoffboeing.com/2018/03/osmnx-features-roundup/) and blog posts that use OSMnx.

*I have a usage question.* Please ask it on [StackOverflow](https://stackoverflow.com/search?q=osmnx).


## Features

OSMnx is pronounced as the initialism: "oh-ess-em-en-ex". It is built on top of [NetworkX](https://networkx.org/) and [GeoPandas](https://geopandas.org/), and interacts with [OpenStreetMap](https://www.openstreetmap.org/) APIs to:

  * Download and model street networks or other networked infrastructure anywhere in the world with a single line of code
  * Download any other spatial geometries (e.g., political boundaries, building footprints, grocery stores, transit stops) as a GeoDataFrame
  * Download by city name, polygon, bounding box, or point/address + distance
  * Model driving, walking, biking, and other travel modes
  * Download node elevations and calculate edge grades (inclines)
  * Impute missing speeds and calculate graph edge travel times
  * Simplify and correct the network's topology to clean-up nodes and consolidate complex intersections
  * Fast map-matching of points, routes, or trajectories to nearest graph edges or nodes
  * Save networks to disk as GeoPackages or GraphML files
  * Save/load a street network to/from a .osm XML file
  * Conduct topological and spatial analyses to automatically calculate dozens of indicators
  * Calculate and visualize street bearings and orientations
  * Calculate and visualize shortest-path routes that minimize distance, travel time, elevation, etc
  * Explore street networks as a static map or interactive web map
  * Visualize travel distance and travel time with isoline and isochrone maps
  * Plot figure-ground diagrams of street networks and building footprints

All of these features are demonstrated in the [examples](https://github.com/gboeing/osmnx-examples) repo and documented in the [documentation](https://osmnx.readthedocs.io/). Feature development details are in the [change log](CHANGELOG.md). Read the [journal article](https://geoffboeing.com/publications/osmnx-complex-street-networks/) for further technical details.

## License

OSMnx is open source and licensed under the MIT license. OpenStreetMap's open data [license](https://www.openstreetmap.org/copyright/) requires that derivative works provide proper attribution.
