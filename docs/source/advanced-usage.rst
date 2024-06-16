Advanced Usage
==============

This guide offers further details on advanced use cases. Refer to the :doc:`getting-started` guide for an introduction to the package.

Model Attributes
----------------

An OSMnx model has some standard required attributes, plus some optional attributes. The latter are sometimes present based on the source OSM data's tagging, the :code:`settings` module configuration, and any processing you may have done to add additional attributes (as noted in various functions' documentation).

As a NetworkX `MultiDiGraph`_ object, it has top-level :code:`graph`, :code:`nodes`, and :code:`edges` attributes. The :code:`graph` attribute dictionary must contain a "crs" key defining its coordinate reference system. The :code:`nodes` are identified by OSM ID and each must contain a :code:`data` attribute dictionary that must have "x" and "y" keys defining its coordinates and a "street_count" key defining how many physical streets are incident to it. The :code:`edges` are identified by a 3-tuple of "u" (source node ID), "v" (target node ID), and "key" (to differentiate parallel edges), and each must contain a :code:`data` attribute dictionary that must have an "osmid" key defining its OSM ID and a "length" key defining its length in meters.

The OSMnx :code:`graph` module automatically creates MultiDiGraphs with these required attributes, plus additional optional attributes based on the :code:`settings` module configuration. If you instead manually create your own graph model, make sure it has these required attributes at a minimum.

Topology Clean-Up
-----------------

The :code:`simplification` module automatically processes the network's topology from the original raw OpenStreetMap data, such that nodes represent intersections/dead-ends and edges represent the street segments that link them. This takes two primary forms: graph simplification and intersection consolidation.

**Graph simplification** cleans up the graph's topology so that nodes represent intersections or dead-ends and edges represent street segments. This is important because in OpenStreetMap raw data, ways comprise sets of straight-line segments between nodes: that is, nodes are vertices for streets' curving line geometries, not just intersections and dead-ends. By default, OSMnx simplifies this topology by discarding non-intersection/dead-end nodes while retaining the complete true edge geometry as an edge attribute. When multiple OpenStreetMap ways are merged into a single graph edge, the ways' attribute values can be aggregated into a single value.

**Intersection consolidation** is important because many real-world street networks feature complex intersections and traffic circles, resulting in a cluster of graph nodes where there is really just one true intersection as we would think of it in transportation or urban design. Similarly, divided roads are often represented by separate centerline edges: the intersection of two divided roads thus creates 4 nodes, representing where each edge intersects a perpendicular edge, but these 4 nodes represent a single intersection in the real world. OSMnx can consolidate such complex intersections into a single node and optionally rebuild the graph's edge topology accordingly. When multiple OpenStreetMap nodes are merged into a single graph node, the nodes' attribute values can be aggregated into a single value.

Convert, Project, Save
----------------------

OSMnx's :code:`convert` module can convert a MultiDiGraph to a `MultiGraph`_ if you prefer an undirected representation of the network, or to a `DiGraph`_ if you prefer a directed representation without any parallel edges. It can also convert a MultiDiGraph to/from GeoPandas node and edge `GeoDataFrames`_. The nodes GeoDataFrame is indexed by OSM ID and the edges GeoDataFrame is multi-indexed by :code:`u, v, key` just like a NetworkX edge. This allows you to load arbitrary node/edge ShapeFiles or GeoPackage layers as GeoDataFrames then model them as a MultiDiGraph for graph analysis.

You can easily project your graph to different coordinate reference systems using the :code:`projection` module. If you're unsure which `CRS`_ you want to project to, OSMnx can automatically determine an appropriate UTM CRS for you.

Using the :code:`io` module, you can save your graph to disk as a GraphML file (to load into other network analysis software), a GeoPackage (to load into other GIS software), or an OSM XML file. Use the GraphML format whenever saving a graph for later work with OSMnx.

Working with Elevation
----------------------

The :code:`elevation` module lets you automatically attach elevations to the graph's nodes from a local raster file or a web service like the Google Maps `Elevation API`_. You can also calculate edge grades (i.e., rise-over-run) and analyze the steepness of certain streets or routes.


Usage Limits
------------

Refer to the `Nominatim Usage Policy`_ and `Overpass Commons`_ documentation for API usage limits and restrictions to which you must adhere. If you configure OSMnx to use an alternative API instance, ensure you understand and follow their policies. If you feel you need to exceed these limits, consider installing your own hosted instance and setting OSMnx to use it.

More Info
---------

All of this functionality is demonstrated step-by-step in the OSMnx `Examples Gallery`_, and usage is detailed in the :doc:`user-reference`. More feature development details are in the `Changelog`_. Consult the :doc:`further-reading` resources for additional technical details and research.


.. _Examples Gallery: https://github.com/gboeing/osmnx-examples
.. _Changelog: https://github.com/gboeing/osmnx/blob/main/CHANGELOG.md
.. _MultiDiGraph: https://networkx.org/documentation/stable/reference/classes/multidigraph.html
.. _MultiGraph: https://networkx.org/documentation/stable/reference/classes/multigraph.html
.. _DiGraph: https://networkx.org/documentation/stable/reference/classes/digraph.html
.. _GeoDataFrames: https://geopandas.org/en/stable/docs/reference/geodataframe.html
.. _CRS: https://en.wikipedia.org/wiki/Coordinate_reference_system
.. _Elevation API: https://developers.google.com/maps/documentation/elevation
.. _Nominatim Usage Policy: https://operations.osmfoundation.org/policies/nominatim/
.. _Overpass Commons: https://dev.overpass-api.de/overpass-doc/en/preface/commons.html
