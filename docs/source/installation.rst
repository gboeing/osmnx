Installation
============

Conda
-----

The official supported way to install OSMnx is with `conda`_ or `mamba`_:

.. code-block:: shell

    conda create -n ox -c conda-forge --strict-channel-priority osmnx

This creates a new conda environment and installs OSMnx into it, via the conda-forge channel. If you want other packages, such as :code:`jupyterlab`, installed in this environment as well, just add their names after :code:`osmnx` above.

To upgrade OSMnx to a newer release, remove the conda environment you created and then create a new one again, as above. Don't just run "conda update" or you could get package conflicts. See the `conda`_ and `conda-forge`_ documentation for more details.

Docker
------

You can run OSMnx + JupyterLab directly from the official OSMnx `Docker`_ image.

Pip
---

You should be able to install OSMnx with `uv`_ or `pip`_ (into a virtual environment) but this is not officially supported.

OSMnx is written in pure Python and its installation alone is thus trivially simple if you already have its dependencies installed and tested on your system. However, OSMnx depends on other packages that in turn depend on compiled C/C++ libraries, and installing those dependencies with pip is sometimes challenging depending on your specific system's configuration. If precompiled binaries are not available for your system, you may need to compile and configure those dependencies manually by following their installation instructions.

So, if you're not sure what you're doing, just follow the conda instructions above to avoid installation problems.

.. _conda: https://conda.io/
.. _conda-forge: https://conda-forge.org/
.. _Docker: https://hub.docker.com/r/gboeing/osmnx
.. _mamba: https://mamba.readthedocs.io/
.. _pip: https://pip.pypa.io/
.. _uv: https://docs.astral.sh/uv/
