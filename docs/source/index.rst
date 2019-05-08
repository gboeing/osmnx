
OSMnx documentation
===================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   osmnx

**OSMnx**: retrieve, model, analyze, and visualize street networks from OpenStreetMap. OSMnx is a Python package that lets you download spatial geometries and model, project, visualize, and analyze street networks and  other spatial data from OpenStreetMap's APIs. Users can download and model walkable, drivable, or bikable urban networks with a single line of Python code, and then easily analyze and visualize them.


Citation info
-------------

Boeing, G. 2017. "`OSMnx: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks`_." *Computers, Environment and Urban Systems* 65, 126-139. doi:10.1016/j.compenvurbsys.2017.05.004


Features
--------

OSMnx is built on top of geopandas, networkx, and matplotlib and works with OpenStreetMap's APIs to:

  * Download street networks anywhere in the world with a single line of code
  * Download other infrastructure network types, place polygons, building footprints, and points of interest
  * Download by city name, polygon, bounding box, or point/address + network distance
  * Download drivable, walkable, bikeable, or all street networks
  * Load street network from a local .osm file
  * Visualize street network as a static image or interactive leaflet web map
  * Simplify and correct the network's topology to clean and consolidate intersections
  * Save networks to disk as shapefiles or GraphML
  * Conduct topological and spatial analyses to automatically calculate dozens of indicators
  * Calculate and plot shortest-path routes as a static image or leaflet web map
  * Fast map-matching of points, routes, or trajectories to nearest graph edges or nodes
  * Plot figure-ground diagrams of street networks and/or building footprints
  * Download node elevations and calculate edge grades
  * Visualize travel distance and travel time with isoline and isochrone maps
  * Calculate and visualize street bearings and orientations

Examples and demonstrations of these features are in the GitHub repo (see below). More feature development details are in the `change log`_.


Installation
------------

You can install OSMnx with conda:

.. code-block:: shell

    conda create -n ox -c conda-forge python=3 osmnx

Alternatively, you can run OSMnx + Jupyter directly from this `docker container`_, or you can install OSMnx via pip (if you already have all of its dependencies installed on your system):

.. code-block:: shell

    pip install osmnx

If you have any trouble with the installation, try installing OSMnx in a new, clean `conda environment`_ using conda-forge and strict channel priority:

.. code-block:: shell

    conda config --prepend channels conda-forge
    conda create -n ox --strict-channel-priority python=3 osmnx


Examples
--------

For examples and demos, see the `examples`_ GitHub repo.


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


.. _change log: https://github.com/gboeing/osmnx/blob/master/CHANGELOG.md
.. _docker container: https://hub.docker.com/r/gboeing/osmnx
.. _conda environment: https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
.. _examples: https://github.com/gboeing/osmnx-examples
.. _issue: https://github.com/gboeing/osmnx/issues
.. _OSMnx GitHub repo: https://github.com/gboeing/osmnx
.. _StackOverflow: https://stackoverflow.com/search?q=osmnx
.. _OSMnx\: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks: http://geoffboeing.com/publications/osmnx-complex-street-networks/
