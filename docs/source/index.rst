
OSMnx documentation
===================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   osmnx

**OSMnx**: retrieve, construct, analyze, and visualize street networks from
OpenStreetMap. OSMnx is a Python package that lets you download spatial
geometries and construct, project, visualize, and analyze street networks from
OpenStreetMap's APIs. Users can download and construct walkable, drivable,
or bikable urban networks with a single line of Python code, and then easily
analyze and visualize them.


Citation info
-------------

Boeing, G. 2017. "`OSMnx: New Methods for Acquiring, Constructing, Analyzing,
and Visualizing Complex Street Networks`_." *Computers, Environment and Urban
Systems* 65, 126-139. doi:10.1016/j.compenvurbsys.2017.05.004


Features
--------

OSMnx is built on top of geopandas, networkx, and matplotlib and works with
OpenStreetMap's APIs to:

  * Download street networks anywhere in the world with a single line of code
  * Download other infrastructure network types, place polygons, or building footprints
  * Download by city name, polygon, bounding box, or point/address + network distance
  * Download drivable, walkable, bikeable, or all street networks
  * Load street network from a local .osm file
  * Visualize street network as a static image or leaflet web map
  * Simplify and correct the network's topology to clean and consolidate intersections
  * Save networks to disk as shapefiles or GraphML
  * Conduct topological and spatial analyses to automatically calculate dozens of indicators
  * Calculate and plot shortest-path routes as a static image or leaflet web map
  * Plot figure-ground diagrams of street networks and/or building footprints
  * Download node elevations and calculate edge grades
  * Visualize travel distance and travel time with isoline and isochrone maps
  * Calculate and visualize street bearings and orientations

Examples and demonstrations of these features are in the GitHub repo (see below).
More feature development details are in the `change log`_.


Installation
------------

Install OSMnx with conda (easiest):

.. code-block:: shell

    conda install -c conda-forge osmnx

or with pip:

.. code-block:: shell

    pip install osmnx

If you have any trouble with the installation, try installing OSMnx in a new,
clean `virtual environment`_ using conda and conda-forge:

.. code-block:: shell

    conda create --override-channels -c conda-forge -n OSMNX python=3 osmnx
    source activate OSMNX


Examples
--------

For examples and demos, see the `examples`_ GitHub repo.


Support
-------

The `issue tracker`_ is at the `OSMnx GitHub repo`_.


License
-------

The project is licensed under the MIT license.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _change log: https://github.com/gboeing/osmnx/blob/master/CHANGELOG.md
.. _virtual environment: https://conda.io/docs/using/envs.html
.. _examples: https://github.com/gboeing/osmnx-examples
.. _issue tracker: https://github.com/gboeing/osmnx/issues
.. _OSMnx GitHub repo: https://github.com/gboeing/osmnx
.. _OSMnx\: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks: http://geoffboeing.com/publications/osmnx-complex-street-networks/
