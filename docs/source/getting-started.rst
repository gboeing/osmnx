Getting Started
===============

Get Started in 4 Steps
----------------------

1. Install OSMnx by following the :doc:`installation` guide.

2. Read the :ref:`introducing-osmnx` section on this page.

3. Work through the OSMnx `Examples Gallery`_ for step-by-step tutorials and sample code.

4. Consult the :doc:`user-reference` for complete details on using the package.

Finally, if you're not already familiar with `NetworkX`_ and `GeoPandas`_, make sure you read their user guides as OSMnx uses their data structures.

.. _introducing-osmnx:

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
* Save/load network to/from disk as GraphML, GeoPackage, or OSM XML file
* Conduct topological and spatial analyses to automatically calculate dozens of indicators
* Calculate and visualize street bearings and orientations
* Calculate and visualize shortest-path routes that minimize distance, travel time, elevation, etc
* Explore street networks and geospatial features as a static map or interactive web map
* Visualize travel distance and travel time with isoline and isochrone maps
* Plot figure-ground diagrams of street networks and building footprints

The OSMnx `Examples Gallery`_ contains tutorials and demonstrations of all these features, and package usage is detailed in the :doc:`user-reference`.

Configuration
^^^^^^^^^^^^^

You can configure OSMnx using the :code:`settings` module. Here you can adjust logging behavior, caching, server endpoints, and more. You can also configure OSMnx to retrieve historical snapshots of OpenStreetMap data as of a certain date.

Read more about the :ref:`settings <osmnx-settings-module>` module in the User Reference.

Geocoding and Querying
^^^^^^^^^^^^^^^^^^^^^^

OSMnx geocodes place names and addresses with the OpenStreetMap `Nominatim`_ API. You can use the :code:`geocoder` module to geocode place names or addresses to lat-lon coordinates. Or, you can retrieve place boundaries or any other OpenStreetMap elements by name or ID. Read more about the :ref:`geocoder <osmnx-geocoder-module>` module in the User Reference.

Using the :code:`features` and :code:`graph` modules, as described below, you can download data by lat-lon point, address, bounding box, bounding polygon, or place name (e.g., neighborhood, city, county, etc).

Urban Amenities
^^^^^^^^^^^^^^^

Using OSMnx's :code:`features` module, you can search for and download any geospatial `features`_ (such as building footprints, grocery stores, schools, public parks, transit stops, etc) from the OpenStreetMap `Overpass`_ API as a GeoPandas GeoDataFrame. This uses OpenStreetMap `tags`_ to search for matching `elements`_.

Read more about the :ref:`features <osmnx-features-module>` module in the User Reference.

Modeling a Network
^^^^^^^^^^^^^^^^^^

Using OSMnx's :code:`graph` module, you can retrieve any spatial network data (such as streets, paths, rail, canals, etc) from the Overpass API and model them as NetworkX `MultiDiGraphs`_.

In short, MultiDiGraphs are nonplanar directed graphs with possible self-loops and parallel edges. Thus, a one-way street will be represented with a single directed edge from node *u* to node *v*, but a bidirectional street will be represented with two reciprocal directed edges (with identical geometries): one from node *u* to node *v* and another from *v* to *u*, to represent both possible directions of flow. Because these graphs are nonplanar, they correctly model the topology of interchanges, bridges, and tunnels. That is, edge crossings in a two-dimensional plane are not intersections in an OSMnx model unless they represent true junctions in the three-dimensional real world.

The :code:`graph` module uses filters to query the Overpass API: you can either specify a built-in network type or provide your own custom filter with `Overpass QL`_. Under the hood, OSMnx does several things to generate the best possible model. It initially creates a 500m-buffered graph before truncating it to your desired query area, to ensure accurate streets-per-node stats and to attenuate graph perimeter effects. By default, it returns the largest weakly connected component. It also simplifies the graph topology as discussed below.

Read more about the :ref:`graph <osmnx-graph-module>` module in the User Reference and refer to the official reference paper at the :doc:`further-reading` page for complete modeling details.

Topology Clean-Up
^^^^^^^^^^^^^^^^^

The :code:`simplification` module automatically processes the network's topology from the original raw OpenStreetMap data, such that nodes represent intersections/dead-ends and edges represent the street segments that link them. This takes two primary forms: graph simplification and intersection consolidation.

**Graph simplification** cleans up the graph's topology so that nodes represent intersections or dead-ends and edges represent street segments. This is important because in OpenStreetMap raw data, ways comprise sets of straight-line segments between nodes: that is, nodes are vertices for streets' curving line geometries, not just intersections and dead-ends. By default, OSMnx simplifies this topology by discarding non-intersection/dead-end nodes while retaining the complete true edge geometry as an edge attribute. When multiple OpenStreetMap ways are merged into a single graph edge, the ways' attribute values can be aggregated into a single value.

**Intersection consolidation** is important because many real-world street networks feature complex intersections and traffic circles, resulting in a cluster of graph nodes where there is really just one true intersection as we would think of it in transportation or urban design. Similarly, divided roads are often represented by separate centerline edges: the intersection of two divided roads thus creates 4 nodes, representing where each edge intersects a perpendicular edge, but these 4 nodes represent a single intersection in the real world. OSMnx can consolidate such complex intersections into a single node and optionally rebuild the graph's edge topology accordingly. When multiple OpenStreetMap nodes are merged into a single graph node, the nodes' attribute values can be aggregated into a single value.

Read more about the :ref:`simplification <osmnx-simplification-module>` module in the User Reference.

Model Attributes
^^^^^^^^^^^^^^^^

An OSMnx model has some standard required attributes, plus some optional attributes. The latter are sometimes present based on the source OSM data's tagging, the :code:`settings` module configuration, and any processing you may have done to add additional attributes (as noted in various functions' documentation).

As a NetworkX `MultiDiGraph`_ object, it has top-level :code:`graph`, :code:`nodes`, and :code:`edges` attributes. The :code:`graph` attribute dictionary must contain a "crs" key defining its coordinate reference system. The :code:`nodes` are identified by OSM ID and each must contain a :code:`data` attribute dictionary that must have "x" and "y" keys defining its coordinates and a "street_count" key defining how many physical streets are incident to it. The :code:`edges` are identified by a 3-tuple of "u" (source node ID), "v" (target node ID), and "key" (to differentiate parallel edges), and each must contain a :code:`data` attribute dictionary that must have an "osmid" key defining its OSM ID and a "length" key defining its length in meters.

The OSMnx :code:`graph` module automatically creates MultiDiGraphs with these required attributes, plus additional optional attributes based on the :code:`settings` module configuration. If you instead manually create your own graph model, make sure it has these required attributes at a minimum.

Convert, Project, Save
^^^^^^^^^^^^^^^^^^^^^^

OSMnx's :code:`convert` module can convert a MultiDiGraph to a `DiGraph`_ if you prefer a directed representation of the network without any parallel edges, or to a `MultiGraph`_ if you need an undirected representation for use with functions or algorithms that only accept a MultiGraph object. If you just want a fully bidirectional graph (such as for a walking network), just configure the `settings` module's `bidirectional_network_types` before creating your graph.

The :code:`convert` module can also convert a MultiDiGraph to/from GeoPandas node and edge `GeoDataFrames`_. The nodes GeoDataFrame is indexed by OSM ID and the edges GeoDataFrame is multi-indexed by :code:`u, v, key` just like a NetworkX edge. This allows you to load arbitrary node/edge ShapeFiles or GeoPackage layers as GeoDataFrames then model them as a MultiDiGraph for graph analysis. Read more about the :ref:`convert <osmnx-convert-module>` module in the User Reference.

You can easily project your graph to different coordinate reference systems using the :code:`projection` module. If you're unsure which `CRS`_ you want to project to, OSMnx can automatically determine an appropriate UTM CRS for you. Read more about the :ref:`projection <osmnx-projection-module>` module in the User Reference.

Using the :code:`io` module, you can save your graph to disk as a GraphML file (to load into other network analysis software), a GeoPackage (to load into other GIS software), or an OSM XML file. Use the GraphML format whenever saving a graph for later work with OSMnx. Read more about the :ref:`io <osmnx-io-module>` module in the User Reference.

Network Measures
^^^^^^^^^^^^^^^^

You can use the :code:`stats` module to calculate a variety of geometric and topological measures as well as street network bearing and orientation statistics. These measures define streets as the edges in an undirected representation of the graph to prevent double-counting bidirectional edges of a two-way street. You can easily generate common stats in transportation studies, urban design, and network science, including intersection density, circuity, average node degree (connectedness), betweenness centrality, and much more. Read more about the :ref:`stats <osmnx-stats-module>` module in the User Reference.

You can also use NetworkX directly to calculate additional topological network measures.

Working with Elevation
^^^^^^^^^^^^^^^^^^^^^^

The :code:`elevation` module lets you automatically attach elevations to the graph's nodes from a local raster file or the Google Maps `Elevation API`_ (or equivalent web API with a compatible interface). You can also calculate edge grades (i.e., rise-over-run) and analyze the steepness of certain streets or routes.

Read more about the :ref:`elevation <osmnx-elevation-module>` module in the User Reference.

Routing
^^^^^^^

The :code:`distance` module can find the nearest node(s) or edge(s) to coordinates using a fast spatial index. The :code:`routing` module can solve shortest paths for network routing, parallelized with multiprocessing, using different weights (e.g., distance, travel time, elevation change, etc). It can also impute missing speeds to the graph edges. This imputation can obviously be imprecise, so the user can override it by passing in arguments that define local speed limits. It can also calculate free-flow travel times for each edge.

Read more about the :ref:`distance <osmnx-distance-module>` and :ref:`routing <osmnx-routing-module>` modules in the User Reference.

Visualization
^^^^^^^^^^^^^

You can plot graphs, routes, network figure-ground diagrams, building footprints, and street network orientation rose diagrams (aka, polar histograms) with the :code:`plot` module. You can also explore street networks, routes, or geospatial features as interactive `Folium`_ web maps.

Read more about the :ref:`plot <osmnx-plot-module>` module in the User Reference.

Usage Limits
^^^^^^^^^^^^

Refer to the `Nominatim Usage Policy`_ and `Overpass Commons`_ documentation for API usage limits and restrictions to which you must adhere. If you configure OSMnx to use an alternative API instance, ensure you understand and follow their policies. If you feel you need to exceed these limits, consider installing your own hosted instance and setting OSMnx to use it.

More Info
---------

All of this functionality is demonstrated step-by-step in the OSMnx `Examples Gallery`_, and usage is detailed in the :doc:`user-reference`. Feature development details are in the `Changelog`_. Consult the :doc:`further-reading` resources for additional technical details and research.

Frequently Asked Questions
--------------------------

*How do I install OSMnx?* Follow the :doc:`installation` guide.

*How do I use OSMnx?* Check out the step-by-step tutorials in the OSMnx `Examples Gallery`_.

*How does this or that function work?* Consult the :doc:`user-reference`.

*What can I do with OSMnx?* Check out recent `projects`_ that use OSMnx.

*I have a usage question.* Please ask it on `StackOverflow`_.


.. _Changelog: https://github.com/gboeing/osmnx/blob/main/CHANGELOG.md
.. _CRS: https://en.wikipedia.org/wiki/Coordinate_reference_system
.. _DiGraph: https://networkx.org/documentation/stable/reference/classes/digraph.html
.. _elements: https://wiki.openstreetmap.org/wiki/Elements
.. _Elevation API: https://developers.google.com/maps/documentation/elevation
.. _Examples Gallery: https://github.com/gboeing/osmnx-examples
.. _features: https://wiki.openstreetmap.org/wiki/Map_features
.. _Folium: https://python-visualization.github.io/folium/
.. _GeoDataFrames: https://geopandas.org/en/stable/docs/reference/geodataframe.html
.. _GeoPandas: https://geopandas.org
.. _MultiDiGraph: https://networkx.org/documentation/stable/reference/classes/multidigraph.html
.. _MultiDiGraphs: https://networkx.org/documentation/stable/reference/classes/multidigraph.html
.. _MultiGraph: https://networkx.org/documentation/stable/reference/classes/multigraph.html
.. _NetworkX: https://networkx.org
.. _Nominatim: https://nominatim.org
.. _Nominatim Usage Policy: https://operations.osmfoundation.org/policies/nominatim/
.. _OpenStreetMap: https://www.openstreetmap.org
.. _Overpass: https://wiki.openstreetmap.org/wiki/Overpass_API
.. _Overpass Commons: https://dev.overpass-api.de/overpass-doc/en/preface/commons.html
.. _Overpass QL: https://wiki.openstreetmap.org/wiki/Overpass_API/Overpass_QL
.. _projects: https://geoffboeing.com/2018/03/osmnx-features-roundup
.. _StackOverflow: https://stackoverflow.com/search?q=osmnx
.. _tags: https://wiki.openstreetmap.org/wiki/Tags
