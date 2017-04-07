

OSMnx documentation
===================

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   osmnx

**OSMnx**: retrieve, construct, analyze, and visualize street networks from OpenStreetMap. 
OSMnx is a Python package that lets you download spatial geometries and construct, project, visualize, 
and analyze street networks from OpenStreetMap's APIs. Users can download and construct walkable, drivable, 
or bikable urban networks with a single line of Python code, and then easily analyze and visualize them.


Installation
------------

Install OSMnx with pip by running:

.. code-block:: shell
    
    pip install osmnx
    
or with conda by running:
    
.. code-block:: shell
    
    conda install -c conda-forge osmnx
    
If you have any trouble with the installation, try installing OSMnx in a new, clean `virtual environment`_:

.. code-block:: shell
    
    conda create --yes -c conda-forge -n OSMNX python=3 osmnx
    source activate OSMNX
    
Examples
--------

For examples and demos, see the GitHub repo: `https://github.com/gboeing/osmnx`_
    
Support
-------

Issue tracker: `https://github.com/gboeing/osmnx/issues`_

License
-------

The project is licensed under the MIT license.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _virtual environment: https://conda.io/docs/using/envs.html
.. _https://github.com/gboeing/osmnx: https://github.com/gboeing/osmnx
.. _https://github.com/gboeing/osmnx/issues: https://github.com/gboeing/osmnx/issues
