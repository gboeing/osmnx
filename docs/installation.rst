Installation
============

Conda
-----

The official supported way to install OSMnx is with conda:

.. code-block:: shell

    conda create -n ox -c conda-forge --strict-channel-priority osmnx

This creates a new conda environment and installs OSMnx into it, via the conda-forge channel. If you want other packages, such as :code:`jupyterlab`, installed in this environment as well, just add their names after :code:`osmnx` above.

To upgrade OSMnx to a newer release, remove the conda environment you created and then create a new one again, as above. Don't just run "conda update" or you could get package conflicts. See the `conda`_ and `conda-forge`_ documentation for more details.

Docker
------

You can run OSMnx + JupyterLab directly from the official OSMnx `Docker`_ image.

Pip
---

If you already have all of its dependencies installed and tested on your system, you can install OSMnx with `pip`_. OSMnx is written in pure Python and its installation alone is thus trivially simple. However, it depends on other packages written in C. Installing those dependencies with pip can be challenging depending on your system's configuration. Therefore, if you don't know exactly what you're doing, just follow the conda instructions above to avoid installation problems.

.. _conda: https://conda.io/
.. _conda-forge: https://conda-forge.org/
.. _Docker: https://hub.docker.com/r/gboeing/osmnx
.. _pip: https://pypi.org/project/osmnx/
