Getting Started
===============

Get Started in 4 Steps
----------------------

1. Install OSMnx by following the :doc:`installation` guide.

2. Read ":ref:`Introducing OSMnx`" below on this page.

3. Work through the `OSMnx Examples`_ gallery for step-by-step tutorials and sample code.

4. Consult the :doc:`user-reference` for complete details on using the package.

Finally, if you're not already familiar with `NetworkX`_ and `GeoPandas`_, make sure you read their user guides as OSMnx uses their data structures and functionality.

.. _Introducing OSMnx:

Introducing OSMnx
-----------------

This quick introduction explains key concepts and the basic functionality of OSMnx.

Overview
^^^^^^^^

OSMnx is pronounced as the initialism: "oh-ess-em-en-ex". It is built on top of NetworkX and GeoPandas, and interacts with `OpenStreetMap`_ APIs to:

* Download and model street networks or other infrastructure anywhere in the world with a single line of code
* Download geospatial features (e.g., political boundaries, building footprints, grocery stores, transit stops) as a GeoDataFrame
* Query by city name, polygon, bounding box, or point/address + distance
* Model driving, walking, biking, and other travel modes
* Attach node elevations from a local raster file or web service and calculate edge grades
* Impute missing speeds and calculate graph edge travel times
* Simplify and correct the network's topology to clean-up nodes and consolidate complex intersections
* Fast map-matching of points, routes, or trajectories to nearest graph edges or nodes
* Save/load network to/from disk as GraphML, GeoPackage, or .osm XML file
* Conduct topological and spatial analyses to automatically calculate dozens of indicators
* Calculate and visualize street bearings and orientations
* Calculate and visualize shortest-path routes that minimize distance, travel time, elevation, etc
* Explore street networks and geospatial features as a static map or interactive web map
* Visualize travel distance and travel time with isoline and isochrone maps
* Plot figure-ground diagrams of street networks and building footprints

The `OSMnx Examples`_ gallery contains tutorials and demonstrations of all these features, and package usage is detailed in the :doc:`user-reference`.

Configuration
^^^^^^^^^^^^^

You can configure OSMnx using the :code:`settings` module. Here you can adjust logging behavior, caching, server endpoints, and more. You can also configure OSMnx to retrieve historical snapshots of OpenStreetMap data as of a certain date.

Geocoding and Querying
^^^^^^^^^^^^^^^^^^^^^^

OSMnx geocodes place names and addresses with the OpenStreetMap Nominatim API. You can use the :code:`geocoder` module to geocode place names or addresses to lat-lng coordinates. Or, you can retrieve places or any other geospatial features by name or OpenStreetMap ID.

Using the :code:`features` and :code:`graph` modules, as described below, you can download data by lat-lng point, address, bounding box, bounding polygon, or place name (e.g., neighborhood, city, county, etc).

Urban Amenities
^^^^^^^^^^^^^^^

Using OSMnx's :code:`features` module, you can search for and download any geospatial `features`_ (such as building footprints, grocery stores, schools, public parks, transit stops, etc) from the OpenStreetMap Overpass API as a GeoPandas GeoDataFrame. This uses OpenStreetMap `tags`_ to search for matching features.

Modeling a Network
^^^^^^^^^^^^^^^^^^

Using OSMnx's :code:`graph` module, you can retrieve any spatial network data (such as streets, paths, canals, etc) from the Overpass API and model them as NetworkX MultiDiGraphs. MultiDiGraphs are nonplanar directed graphs with possible self-loops and parallel edges.

Thus, a one-way street will be represented with a single directed edge from node *u* to node *v*, but a bidirectional street will be represented with two reciprocal directed edges (with identical geometries): one from node *u* to node *v* and another from *v* to *u*, to represent both possible directions of flow. Because these graphs are nonplanar, they correctly model the topology of interchanges, bridges, and tunnels. That is, edge crossings in a two-dimensional plane are not intersections in an OSMnx model unless they represent true junctions in the three-dimensional real world.

Topology Clean-Up
^^^^^^^^^^^^^^^^^

OSMnx's :code:`simplification` module automatically processes network topology from the original raw OpenStreetMap data such that nodes represent intersections/dead-ends and edges represent the street segments that link them. This takes two primary forms: graph simplification and intersection consolidation.

Graph simplification cleans up the graph's topology so that nodes represent intersections or dead-ends and edges represent street segments. This is important because in OpenStreetMap raw data, ways comprise sets of straight-line segments between nodes: that is, nodes are vertices for streets' curving line geometries, not just intersections and dead-ends. By default, OSMnx simplifies this topology by discarding non-intersection/dead-end nodes, while retaining the complete true edge geometry as an edge attribute.

Intersection consolidation is important because many real-world street networks feature complex intersections and traffic circles, resulting in a cluster of graph nodes where there is really just one true intersection as we would think of it in transportation or urban design. Similarly, divided roads are often represented by separate centerline edges: the intersection of two divided roads thus creates 4 nodes, representing where each edge intersects a perpendicular edge, but these 4 nodes represent a single intersection in the real world. OSMnx can consolidate such complex intersections and optionally rebuild the graph's edge topology accordingly.

Converting Graphs
^^^^^^^^^^^^^^^^^

OSMnx can convert a MultiDiGraph to a MultiGraph if you prefer an undirected representation of the network, or to a DiGraph if you prefer a directed representation without any parallel edges.

It can also convert a MultiDiGraph to/from GeoPandas node and edge GeoDataFrames. This allows you to load arbitrary node/edge ShapeFiles or GeoPackage layers as GeoDataFrames then model them as a MultiDiGraph for graph analysis.

You can easily project your graphs to different coordinate reference systems using the :code:`projection` module. If you're unsure which CRS you want to project to, OSMnx can automatically determine an appropriate UTM CRS for you.

You can save your OSMnx graph to disk as a GraphML file or GeoPackage using the :code:`io` module. Use the GraphML format whenever saving a graph for later work with OSMnx.

Working with Elevation
^^^^^^^^^^^^^^^^^^^^^^

Using the :code:`elevation` module, you can automatically attach elevations to the graph's nodes from a local raster file or web service such as the Google Maps Elevation API. You can also calculate edge grades (i.e., rise-over-run).

Network Statistics
^^^^^^^^^^^^^^^^^^

You can use the :code:`stats` module to calculate a variety of geometric and topological measures as well as street network bearing/orientation statistics. These measures define streets as the edges in an undirected representation of the graph to prevent double-counting bidirectional edges of a two-way street.

You can use NetworkX directly to calculate additional topological network measures.

Routing
^^^^^^^

The :code:`speed` module can impute missing speeds to the graph edges. This imputation can obviously be imprecise, and the user can override it by passing in arguments that define local speed limits. It can also calculate free-flow travel times for each edge.

The :code:`distance` module can find the nearest node(s) or edge(s) to coordinates using a fast spatial index. It can also solve shortest paths for network routing, parallelized with multiprocessing, using different weights (e.g., distance, travel time, elevation change, etc).

Visualization
^^^^^^^^^^^^^

You can plot graphs, routes, network figure-ground diagrams, building footprints, and street network orientation rose diagrams (polar histograms) with the :code:`plot` module. You can also explore street networks, routes, or geospatial features as interactive folium/leaflet web maps.

More Info
---------

All of this functionality is demonstrated step-by-step in the `OSMnx Examples`_ gallery, and usage is detailed in the :doc:`user-reference`. More feature development details are in the `Change Log`_. Consult the :doc:`further-reading` resources for additional technical details and research.

Frequently Asked Questions
--------------------------

*How do I install OSMnx?* Follow the :doc:`installation` guide.

*How do I use OSMnx?* Check out the step-by-step tutorials in the `OSMnx Examples`_ gallery.

*How does this or that function work?* Consult the :doc:`user-reference`.

*What can I do with OSMnx?* Check out recent `projects`_ that use OSMnx.

*I have a usage question.* Please ask it on `StackOverflow`_.

.. _OSMnx Examples: https://github.com/gboeing/osmnx-examples
.. _GeoPandas: https://geopandas.org/
.. _NetworkX: https://networkx.org/
.. _OpenStreetMap: https://www.openstreetmap.org/
.. _features: https://wiki.openstreetmap.org/wiki/Map_features
.. _tags: https://wiki.openstreetmap.org/wiki/Tags
.. _Change Log: https://github.com/gboeing/osmnx/blob/main/CHANGELOG.md
.. _projects: https://geoffboeing.com/2018/03/osmnx-features-roundup/
.. _StackOverflow: https://stackoverflow.com/search?q=osmnx
