User Reference
==============

User reference for the OSMnx package.

This guide covers usage of all public modules and functions. Every function can be accessed via `ox.module_name.function_name()` and most can also be accessed directly via `ox.function_name()` as a shortcut. Less-common functions are accessible only via `ox.module_name.function_name()`.

OSMnx geocodes place names and addresses with the OpenStreetMap Nominatim API. Using OSMnx's :code:`geometries` module, you can retrieve any geospatial objects (such as building footprints, grocery stores, schools, public parks, transit stops, etc) from the OpenStreetMap Overpass API as a GeoPandas GeoDataFrame. Using OSMnx's :code:`graph` module, you can retrieve any spatial network data (such as streets, paths, canals, etc) from the Overpass API and model them as NetworkX MultiDiGraphs.

OSMnx automatically processes network topology from the original raw OpenStreetMap data such that nodes represent intersections/dead-ends and edges represent the street segments that link them. MultiDiGraphs are nonplanar directed graphs with possible self-loops and parallel edges. Thus, a one-way street will be represented with a single directed edge from node u to node v, but a bidirectional street will be represented with two reciprocal directed edges (with identical geometries): one from node u to node v and another from v to u, to represent both possible directions of flow. OSMnx can convert a MultiDiGraph to a MultiGraph if you prefer an undirected representation of the network. It can also convert a MultiDiGraph to/from GeoPandas node and edge GeoDataFrames.

osmnx.bearing module
--------------------

.. automodule:: osmnx.bearing
    :members:

osmnx.distance module
---------------------

.. automodule:: osmnx.distance
    :members:

osmnx.downloader module
-----------------------

.. automodule:: osmnx.downloader
    :members:

osmnx.elevation module
----------------------

.. automodule:: osmnx.elevation
    :members:

osmnx.folium module
-------------------

.. automodule:: osmnx.folium
    :members:

osmnx.geocoder module
---------------------

.. automodule:: osmnx.geocoder
    :members:

osmnx.geometries module
-----------------------

.. automodule:: osmnx.geometries
    :members:

osmnx.graph module
------------------

.. automodule:: osmnx.graph
    :members:

osmnx.io module
---------------

.. automodule:: osmnx.io
    :members:

osmnx.osm_xml module
--------------------

.. automodule:: osmnx.osm_xml
    :members:

osmnx.plot module
-----------------

.. automodule:: osmnx.plot
    :members:

osmnx.projection module
-----------------------

.. automodule:: osmnx.projection
    :members:

osmnx.settings module
---------------------

.. automodule:: osmnx.settings
    :members:

osmnx.simplification module
---------------------------

.. automodule:: osmnx.simplification
    :members:

osmnx.speed module
------------------

.. automodule:: osmnx.speed
    :members:

osmnx.stats module
------------------

.. automodule:: osmnx.stats
    :members:

osmnx.truncate module
---------------------

.. automodule:: osmnx.truncate
    :members:

osmnx.utils module
------------------

.. automodule:: osmnx.utils
    :members:

osmnx.utils_geo module
----------------------

.. automodule:: osmnx.utils_geo
    :members:

osmnx.utils_graph module
------------------------

.. automodule:: osmnx.utils_graph
    :members:
