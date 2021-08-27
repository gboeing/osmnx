OSMnx |version|
===============

OSMnx is a Python package that lets you download geospatial data from OpenStreetMap and model, project, visualize, and analyze real-world street networks and any other geospatial geometries. You can download and model walkable, drivable, or bikeable urban networks with a single line of Python code then easily analyze and visualize them. You can just as easily download and work with other infrastructure types, amenities/points of interest, building footprints, elevation data, street bearings/orientations, and speed/travel time.

If you use OSMnx in your work, please cite the journal article:

Boeing, G. 2017. `OSMnx: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks`_. *Computers, Environment and Urban Systems* 65, 126-139. doi:10.1016/j.compenvurbsys.2017.05.004

.. _OSMnx\: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks: https://geoffboeing.com/publications/osmnx-complex-street-networks/



Installation
------------

You can install OSMnx with `conda`_:

.. code-block:: shell

    conda config --prepend channels conda-forge
    conda create -n ox --strict-channel-priority osmnx

If you want other packages, such as :code:`jupyterlab`, installed in this environment as well, just add their names after :code:`osmnx` above. See the conda documentation for further details. To upgrade OSMnx to a newer release, remove the conda environment you created and then create a new one again, as above. Don't just run "conda update" or you could get package conflicts.

You can also run OSMnx + Jupyter directly from its official `Docker container`_, or you can install OSMnx via `pip`_ if you already have all of its dependencies installed and fully tested on your system. Note: installing the dependencies with pip is nontrivial. If you don't know exactly what you're doing, just use conda as described above.

.. _conda: https://conda.io/
.. _Docker container: https://hub.docker.com/r/gboeing/osmnx
.. _pip: https://pypi.org/project/osmnx/



Usage
-----

To get started with OSMnx, read its `user reference`_ and work through its `examples`_ repo for introductory usage demonstrations and sample code. Make sure you have read the `GeoPandas`_ and `NetworkX`_ user guides if you're not already familiar with these packages, as OSMnx uses their data structures and functionality.

OSMnx is built on top of GeoPandas, NetworkX, and matplotlib and interacts with OpenStreetMap's APIs to:

  * Download and model street networks or other networked infrastructure anywhere in the world with a single line of code
  * Download any other spatial geometries, place boundaries, building footprints, or points of interest as a GeoDataFrame
  * Download by city name, polygon, bounding box, or point/address + network distance
  * Download drivable, walkable, bikeable, or all street networks
  * Download node elevations and calculate edge grades (inclines)
  * Impute missing speeds and calculate graph edge travel times
  * Simplify and correct the network's topology to clean-up nodes and consolidate intersections
  * Fast map-matching of points, routes, or trajectories to nearest graph edges or nodes
  * Save networks to disk as shapefiles, GeoPackages, and GraphML
  * Save/load street network to/from a local .osm XML file
  * Conduct topological and spatial analyses to automatically calculate dozens of indicators
  * Calculate and visualize street bearings and orientations
  * Calculate and visualize shortest-path routes that minimize distance, travel time, elevation, etc
  * Visualize street networks as a static map or interactive Leaflet web map
  * Visualize travel distance and travel time with isoline and isochrone maps
  * Plot figure-ground diagrams of street networks and building footprints

OSMnx geocodes place names and addresses with the OpenStreetMap Nominatim API. Using OSMnx's :code:`geometries` module, you can retrieve any geospatial objects (such as building footprints, grocery stores, schools, public parks, transit stops, etc) from the OpenStreetMap Overpass API as a GeoPandas GeoDataFrame. Using OSMnx's :code:`graph` module, you can retrieve any spatial network data (such as streets, paths, canals, etc) from the Overpass API and model them as NetworkX MultiDiGraphs.

OSMnx automatically processes network topology from the original raw OpenStreetMap data such that nodes represent intersections/dead-ends and edges represent the street segments that link them. MultiDiGraphs are nonplanar directed graphs with possible self-loops and parallel edges. Thus, a one-way street will be represented with a single directed edge from node u to node v, but a bidirectional street will be represented with two reciprocal directed edges (with identical geometries): one from node u to node v and another from v to u, to represent both possible directions of flow. OSMnx can convert a MultiDiGraph to a MultiGraph if you prefer an undirected representation of the network. It can also convert a MultiDiGraph to/from GeoPandas node and edge GeoDataFrames.

Usage examples and demonstrations of these features are in the `examples`_ GitHub repo. More feature development details are in the `change log`_. Read the `journal article`_ for further technical details. Package usage is detailed in the `user reference`_.

.. _user reference: osmnx.html
.. _examples: https://github.com/gboeing/osmnx-examples
.. _GeoPandas: https://geopandas.org/
.. _NetworkX: https://networkx.org/
.. _journal article: https://geoffboeing.com/publications/osmnx-complex-street-networks/
.. _change log: https://github.com/gboeing/osmnx/blob/master/CHANGELOG.md



User reference
--------------

.. toctree::
   :maxdepth: 2

   osmnx

.. toctree::
   :maxdepth: 1

   internals



Support
-------

If you have a usage question, please ask it on `StackOverflow`_. If you've discovered a bug in OSMnx, please open an issue at the OSMnx `GitHub`_ repo.

.. _GitHub: https://github.com/gboeing/osmnx
.. _StackOverflow: https://stackoverflow.com/search?q=osmnx



License
-------

OSMnx is licensed under the MIT license. OpenStreetMap's open data `license`_ requires that derivative works provide proper attribution.

.. _license: https://www.openstreetmap.org/copyright



Indices
-------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
