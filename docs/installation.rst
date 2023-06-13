Installation
============

Conda
-----

The official supported way to install OSMnx is with conda:

.. code-block:: shell

    conda config --prepend channels conda-forge
    conda create -n ox --strict-channel-priority osmnx

If you want other packages, such as :code:`jupyterlab`, installed in this environment as well, just add their names after :code:`osmnx` above. See the `conda`_ documentation for further details.

To upgrade OSMnx to a newer release, remove the conda environment you created and then create a new one again, as above. Don't just run "conda update" or you could get package conflicts.

You can equivalently use mamba as a drop-in replacement for conda.

Docker
------

You can run OSMnx + Jupyter directly with the official OSMnx `Docker`_ image.

Pip
---

If you already have all of its dependencies installed and tested on your system, you can install OSMnx with `pip`_.

The OSMnx package is pure Python and its installation alone is thus trivially simple. However, OSMnx depends on other packages which themselves may rely on other tricky dependencies or C extensions. Installing these dependencies with pip can be complicated depending on your system's configuration. Conda makes the complete installation seamless: if you don't know exactly what you're doing, just use conda as described above.

.. _conda: https://conda.io/
.. _Docker: https://hub.docker.com/r/gboeing/osmnx
.. _pip: https://pypi.org/project/osmnx/
