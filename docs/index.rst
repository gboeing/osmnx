OSMnx |version|
===============

OSMnx is a Python package that lets you download geospatial data from OpenStreetMap and model, project, visualize, and analyze real-world street networks and any other geospatial geometries. You can download and model walking, driving, or biking networks with a single line of code then easily analyze and visualize them. You can just as easily download and work with other infrastructure types, amenities/points of interest, building footprints, elevation data, street bearings/orientations, speed/travel time, and routing.

Citation
--------

If you use OSMnx in your work, please cite the journal article:

Boeing, G. 2017. `OSMnx: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks`_. *Computers, Environment and Urban Systems* 65, 126-139.

.. _OSMnx\: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks: https://geoffboeing.com/publications/osmnx-complex-street-networks/


Installation
------------

You can install OSMnx with conda:

.. code-block:: shell

    conda config --prepend channels conda-forge
    conda create -n ox --strict-channel-priority osmnx

For more options and details, read the :doc:`installation` guide.


Getting Started
---------------

Read the :doc:`getting-started` guide for an introduction to the OSMnx package, what you can do with it, and how to use it.


Support
-------

If you have a "how-to" or usage question, please ask it on `StackOverflow`_. If you've discovered a bug in OSMnx, please open an issue at the OSMnx `GitHub`_ repo.

.. _GitHub: https://github.com/gboeing/osmnx
.. _StackOverflow: https://stackoverflow.com/search?q=osmnx


License
-------

OSMnx is licensed under the MIT license. OpenStreetMap's open data `license`_ requires that derivative works provide proper attribution.

.. _license: https://www.openstreetmap.org/copyright


Documentation
-------------

.. toctree::
   :maxdepth: 1

   installation

.. toctree::
   :maxdepth: 1

   getting-started

.. toctree::
   :maxdepth: 1

   osmnx

.. toctree::
   :maxdepth: 1

   internals

.. toctree::
   :maxdepth: 1

   further-reading



Indices
-------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
