"""
OSMnx setup script.

See license in LICENSE.txt.
"""

import os

from setuptools import setup

# provide a long description using reStructuredText
LONG_DESCRIPTION = r"""
OSMnx is a Python package that lets you download geospatial data from
OpenStreetMap and model, project, visualize, and analyze real-world street
networks and any other geospatial geometries. You can download and model
walking, driving, or biking networks with a single line of code then easily
analyze and visualize them. You can just as easily download and work with
other infrastructure types, amenities/points of interest, building footprints,
elevation data, street bearings/orientations, speed/travel time, and routing.

Citation info: Boeing, G. 2017. `OSMnx: New Methods for Acquiring,
Constructing, Analyzing, and Visualizing Complex Street Networks`_.
*Computers, Environment and Urban Systems* 65, 126-139.
doi:10.1016/j.compenvurbsys.2017.05.004

To get started, read the `OSMnx Documentation`_ and work through the
step-by-step `OSMnx Examples`_ repository for introductory usage
demonstrations and sample code.

.. _OSMnx\: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks: http://geoffboeing.com/publications/osmnx-complex-street-networks/
.. _OSMnx Examples: https://github.com/gboeing/osmnx-examples
.. _OSMnx Documentation: https://osmnx.readthedocs.io/
"""

# list of classifiers from the PyPI classifiers trove
CLASSIFIERS = [
    "Development Status :: 5 - Production/Stable",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3 :: Only",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Topic :: Scientific/Engineering :: GIS",
    "Topic :: Scientific/Engineering :: Information Analysis",
    "Topic :: Scientific/Engineering :: Mathematics",
    "Topic :: Scientific/Engineering :: Physics",
    "Topic :: Scientific/Engineering :: Visualization",
]

DESC = (
    "Retrieve, model, analyze, and visualize OpenStreetMap street networks and other spatial data"
)

# only specify install_requires if not in RTD environment
if os.getenv("READTHEDOCS") == "True":
    INSTALL_REQUIRES = []
else:
    with open("requirements.txt") as f:
        INSTALL_REQUIRES = [line.strip() for line in f.readlines()]

# now call setup
setup(
    name="osmnx",
    version="1.4.1dev",
    description=DESC,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/x-rst",
    classifiers=CLASSIFIERS,
    url="https://github.com/gboeing/osmnx",
    author="Geoff Boeing",
    author_email="boeing@usc.edu",
    license="MIT",
    platforms="any",
    packages=["osmnx"],
    python_requires=">=3.8",
    install_requires=INSTALL_REQUIRES,
    extras_require={
        "all": ["gdal", "matplotlib", "rasterio", "scikit-learn", "scipy"],
        "entropy": ["scipy"],
        "nearest_neighbor": ["scikit-learn", "scipy"],
        "visualization": ["matplotlib"],
        "raster": ["gdal", "rasterio"],
    },
)
