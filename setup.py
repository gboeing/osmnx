"""
OSMnx setup script.

See license in LICENSE.txt.
"""

import os

from setuptools import setup

# provide a long description using reStructuredText
LONG_DESCRIPTION = r"""
**OSMnx** is a Python package that lets you download spatial geometries and
model, project, visualize, and analyze real-world street networks from
OpenStreetMap's APIs. Users can download and model walkable, drivable, or
bikeable urban networks with a single line of Python code, and then easily
analyze and visualize them. You can just as easily download and work with
amenities/points of interest, building footprints, elevation data, street
bearings/orientations, speed/travel time, and network routing.

Citation info: Boeing, G. 2017. "`OSMnx: New Methods for Acquiring,
Constructing, Analyzing, and Visualizing Complex Street Networks`_."
*Computers, Environment and Urban Systems* 65, 126-139.
doi:10.1016/j.compenvurbsys.2017.05.004

Read the `docs`_ or see usage examples and demos on `GitHub`_.

.. _GitHub: https://github.com/gboeing/osmnx-examples
.. _docs: https://osmnx.readthedocs.io
.. _OSMnx\: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks: http://geoffboeing.com/publications/osmnx-complex-street-networks/
"""

# list of classifiers from the PyPI classifiers trove
CLASSIFIERS = [
    "Development Status :: 5 - Production/Stable",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: GIS",
    "Topic :: Scientific/Engineering :: Visualization",
    "Topic :: Scientific/Engineering :: Physics",
    "Topic :: Scientific/Engineering :: Mathematics",
    "Topic :: Scientific/Engineering :: Information Analysis",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
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
    version="0.16.2",
    description=DESC,
    long_description=LONG_DESCRIPTION,
    classifiers=CLASSIFIERS,
    url="https://github.com/gboeing/osmnx",
    author="Geoff Boeing",
    author_email="boeing@usc.edu",
    license="MIT",
    platforms="any",
    packages=["osmnx"],
    python_requires=">=3.6",
    install_requires=INSTALL_REQUIRES,
    extras_require={
        "folium": ["folium>=0.11"],
        "kdtree": ["scipy>=1.5"],
        "balltree": ["scikit-learn>=0.23"],
    },
)
