Getting Started
===============

Preliminaries
-------------

You can install OSMnx by following the :doc:`installation` guide.

After it's installed, get started by reading the :doc:`osmnx` and working through the step-by-step `OSMnx Examples`_ repository for introductory usage demonstrations and sample code.

Make sure you have read the `NetworkX`_ and `GeoPandas`_ user guides if you're not already familiar with these packages, as OSMnx uses their data structures and functionality.

Features
--------

OSMnx is pronounced as the initialism: "oh-ess-em-en-ex". It is built on top of NetworkX and GeoPandas, and interacts with `OpenStreetMap`_ APIs to:

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

Examples and demonstrations of all these features are in the `OSMnx Examples`_ repository and package usage is detailed in the :doc:`osmnx`.

Using OSMnx
-----------

Querying
^^^^^^^^

OSMnx geocodes place names and addresses with the OpenStreetMap Nominatim API.

Urban Amenities
^^^^^^^^^^^^^^^

Using OSMnx's :code:`geometries` module, you can retrieve any geospatial objects (such as building footprints, grocery stores, schools, public parks, transit stops, etc) from the OpenStreetMap Overpass API as a GeoPandas GeoDataFrame.

Modeling a Network
^^^^^^^^^^^^^^^^^^

Using OSMnx's :code:`graph` module, you can retrieve any spatial network data (such as streets, paths, canals, etc) from the Overpass API and model them as NetworkX MultiDiGraphs. MultiDiGraphs are nonplanar directed graphs with possible self-loops and parallel edges.

Thus, a one-way street will be represented with a single directed edge from node *u* to node *v*, but a bidirectional street will be represented with two reciprocal directed edges (with identical geometries): one from node *u* to node *v* and another from *v* to *u*, to represent both possible directions of flow. Because these graphs are nonplanar, they correctly model the topology of interchanges, bridges, and tunnels. That is, edge crossings in a two-dimensional plane are not intersections in an OSMnx model unless they represent true junctions in the three-dimensional real world.

Topology Clean-Up
^^^^^^^^^^^^^^^^^

OSMnx automatically processes network topology from the original raw OpenStreetMap data such that nodes represent intersections/dead-ends and edges represent the street segments that link them.

Converting Graphs
^^^^^^^^^^^^^^^^^

OSMnx can convert a MultiDiGraph to a MultiGraph if you prefer an undirected representation of the network, or to a DiGraph if you prefer a directed representation without any parallel edges.

It can also convert a MultiDiGraph to/from GeoPandas node and edge GeoDataFrames. This allows you to load arbitrary node/edge ShapeFiles or GeoPackage layers as GeoDataFrames then model them as a MultiDiGraph for graph analysis.

Projection.

You can save your OSMnx graph to disk as a GraphML file, GeoPackage, or .osm XML file.

Working with Elevation
^^^^^^^^^^^^^^^^^^^^^^

Local raster file or web service.

Network Statistics
^^^^^^^^^^^^^^^^^^

Measures and bearing stats.

Routing
^^^^^^^

Speed and travel time. Nearest nodes/edges. Shortest paths.

Visualization
^^^^^^^^^^^^^

Plot graphs, routes, figure-ground diagrams, building footprints, orientation rose diagrams, interactive web maps.

More Info
^^^^^^^^^

All of this functionality is demonstrated step-by-step in the `OSMnx Examples`_ repository, and usage is detailed in the :doc:`osmnx`.

More feature development details are in the `Change Log`_. Consult the :doc:`further-reading` resources for additional technical details and research.

.. _OSMnx Examples: https://github.com/gboeing/osmnx-examples
.. _GeoPandas: https://geopandas.org/
.. _NetworkX: https://networkx.org/
.. _OpenStreetMap: https://www.openstreetmap.org/
.. _Change Log: https://github.com/gboeing/osmnx/blob/main/CHANGELOG.md
