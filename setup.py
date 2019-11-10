# OSMnx
# See full license in LICENSE.txt

from setuptools import setup

# provide a long description using reStructuredText
long_description = """
**OSMnx** is a package to easily download, model, project, visualize, and
analyze complex street networks from OpenStreetMap in Python with NetworkX.

You can get a city's or neighborhood's walking, driving, or biking network with
a single line of Python code. Then you can easily visualize cul-de-sacs or
one-way streets, plot shortest-path routes, or calculate stats like intersection
density, average node connectivity, betweenness centrality, etc.

Citation info: Boeing, G. 2017. "`OSMnx: New Methods for Acquiring, Constructing, Analyzing,
and Visualizing Complex Street Networks`_." *Computers, Environment and Urban
Systems* 65, 126-139. doi:10.1016/j.compenvurbsys.2017.05.004

See the examples and demos on `GitHub`_ or read more about `OSMnx`_.

.. _GitHub: https://github.com/gboeing/osmnx
.. _OSMnx: http://geoffboeing.com/2016/11/osmnx-python-street-networks/
.. _OSMnx\: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks: http://geoffboeing.com/publications/osmnx-complex-street-networks/
"""

# list of classifiers from the PyPI classifiers trove
classifiers = ['Development Status :: 5 - Production/Stable',
               'License :: OSI Approved :: MIT License',
               'Operating System :: OS Independent',
               'Intended Audience :: Science/Research',
               'Topic :: Scientific/Engineering :: GIS',
               'Topic :: Scientific/Engineering :: Visualization',
               'Topic :: Scientific/Engineering :: Physics',
               'Topic :: Scientific/Engineering :: Mathematics',
               'Topic :: Scientific/Engineering :: Information Analysis',
               'Programming Language :: Python',
               'Programming Language :: Python :: 2',
               'Programming Language :: Python :: 2.7',
               'Programming Language :: Python :: 3',
               'Programming Language :: Python :: 3.4',
               'Programming Language :: Python :: 3.5',
               'Programming Language :: Python :: 3.6',
               'Programming Language :: Python :: 3.7']

with open('requirements.txt') as f:
    requirements_lines = f.readlines()
install_requires = [r.strip() for r in requirements_lines]

# now call setup
setup(name='osmnx',
      version='0.11dev',
      description='Retrieve, model, analyze, and visualize OpenStreetMap street networks and other spatial data',
      long_description=long_description,
      classifiers=classifiers,
      url='https://github.com/gboeing/osmnx',
      author='Geoff Boeing',
      author_email='boeing@usc.edu',
      license='MIT',
      platforms='any',
      packages=['osmnx'],
      install_requires=install_requires,
      extras_require={'folium':['folium>=0.6'],
                      'kdtree':['scipy>=1.1'],
                      'balltree':['scikit-learn>=0.19']})
