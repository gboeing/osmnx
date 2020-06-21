
OSMnx documentation
===================

.. toctree::
   :maxdepth: 2
   :caption: Reference

   osmnx

**OSMnx**: retrieve, model, analyze, and visualize street networks from OpenStreetMap. OSMnx is a Python package that lets you download spatial geometries and model, project, visualize, and analyze real-world street networks from OpenStreetMap's APIs. Users can download and model walkable, drivable, or bikeable urban networks with a single line of Python code, and then easily analyze and visualize them. You can just as easily download and work with amenities/points of interest, building footprints, elevation data, street bearings/orientations, speed/travel time, and network routing.


Citation info
-------------

If you use OSMnx in your work, please cite the journal article:

Boeing, G. 2017. "`OSMnx: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks`_." *Computers, Environment and Urban Systems* 65, 126-139. doi:10.1016/j.compenvurbsys.2017.05.004


Features
--------

OSMnx is built on top of geopandas, networkx, and matplotlib and interacts with OpenStreetMap's APIs to:

  * Download street networks anywhere in the world with a single line of code
  * Download other infrastructure types, place boundaries, building footprints, and points of interest
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
  * Visualize street network as a static map or interactive leaflet web map
  * Visualize travel distance and travel time with isoline and isochrone maps
  * Plot figure-ground diagrams of street networks and/or building footprints

Examples and demonstrations of these features are in the GitHub repo (see below). More feature development details are in the `change log`_.


Installation
------------

You can install OSMnx with conda:

.. code-block:: shell

    conda config --prepend channels conda-forge
    conda create -n ox --strict-channel-priority osmnx

Alternatively, you can run OSMnx + Jupyter directly from its official `docker container`_, or you can install OSMnx via `pip`_ if you already have all of its dependencies installed on your system.


Examples
--------

For code and usage examples/demos, see the `examples`_ GitHub repo.


Support
-------

If you've discovered a bug in OSMnx, please open an `issue`_ at the `OSMnx GitHub repo`_  documenting what is broken in the package. Alternatively, if you have a usage question, please ask it on `StackOverflow`_.


License
-------

The project is licensed under the MIT license.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. toctree::
   :maxdepth: 1
   :caption: Complete reference, including private internal functions:

   osmnx-complete

.. _change log: https://github.com/gboeing/osmnx/blob/master/CHANGELOG.md
.. _docker container: https://hub.docker.com/r/gboeing/osmnx
.. _conda environment: https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
.. _examples: https://github.com/gboeing/osmnx-examples
.. _issue: https://github.com/gboeing/osmnx/issues
.. _OSMnx GitHub repo: https://github.com/gboeing/osmnx
.. _pip: https://pypi.python.org/pypi/OSMnx
.. _StackOverflow: https://stackoverflow.com/search?q=osmnx
.. _OSMnx\: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks: http://geoffboeing.com/publications/osmnx-complex-street-networks/
