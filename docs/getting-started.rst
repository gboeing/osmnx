Getting Started
===============

Preliminaries
-------------

You can install OSMnx by following the :doc:`installation` guide.

After it's installed, get started by reading the :doc:`osmnx` and working through the `OSMnx Examples`_ repository for introductory usage demonstrations and sample code.

Make sure you have read the `NetworkX`_ and `GeoPandas`_ user guides if you're not already familiar with these packages, as OSMnx uses their data structures and functionality.

What Is OSMnx?
--------------

OSMnx is pronounced as the initialism: "oh-ess-em-en-ex". It is built on top of `NetworkX`_ and `GeoPandas`_, and interacts with `OpenStreetMap`_ APIs to:

* Download and model street networks or other networked infrastructure anywhere in the world with a single line of code
* Download any other spatial geometries (e.g., political boundaries, building footprints, amenities, transit lines or stops) as a GeoDataFrame
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

Usage examples and demonstrations of all these features are in the `OSMnx Examples`_ repository.

How Does It Work?
-----------------

OSMnx geocodes place names and addresses with the OpenStreetMap Nominatim API. Using OSMnx's :code:`geometries` module, you can retrieve any geospatial objects (such as building footprints, grocery stores, schools, public parks, transit stops, etc) from the OpenStreetMap Overpass API as a GeoPandas GeoDataFrame.

Using OSMnx's :code:`graph` module, you can retrieve any spatial network data (such as streets, paths, canals, etc) from the Overpass API and model them as NetworkX MultiDiGraphs. MultiDiGraphs are nonplanar directed graphs with possible self-loops and parallel edges. Thus, a one-way street will be represented with a single directed edge from node u to node v, but a bidirectional street will be represented with two reciprocal directed edges (with identical geometries): one from node u to node v and another from v to u, to represent both possible directions of flow.

OSMnx automatically processes network topology from the original raw OpenStreetMap data such that nodes represent intersections/dead-ends and edges represent the street segments that link them. Because they are nonplanar, they correctly model the topology of interchanges, bridges, and tunnels. That is, edge crossings in a two-dimensional plane are not intersections in an OSMnx model unless they are truly intersections in the real world.

OSMnx can convert a MultiDiGraph to a MultiGraph if you prefer an undirected representation of the network. It can also convert a MultiDiGraph to/from GeoPandas node and edge GeoDataFrames. This allows you to load any node/edge ShapeFiles or GeoPackage layers as GeoDataFrames then convert them to a MultiDiGraph for graph analysis.

More feature development details are in the `change log`_. Read the `journal article`_ for further technical details. Package usage is detailed in the :doc:`osmnx`.

.. _OSMnx Examples: https://github.com/gboeing/osmnx-examples
.. _GeoPandas: https://geopandas.org/
.. _NetworkX: https://networkx.org/
.. _OpenStreetMap: https://www.openstreetmap.org/
.. _journal article: https://geoffboeing.com/publications/osmnx-complex-street-networks/
.. _change log: https://github.com/gboeing/osmnx/blob/main/CHANGELOG.md
