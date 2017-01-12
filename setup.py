# OSMnx
# See full license in LICENSE.txt

from setuptools import setup

# provide a long description using reStructuredText
long_description = """
**OSMnx** is a package to easily download, construct, project, visualize, and analyze complex street 
networks from OpenStreetMap in Python with NetworkX. 

You can get a city's or neighborhood's walking, driving, or biking network with a single line of Python 
code. Then you can easily visualize cul-de-sacs or one-way streets, plot shortest-path routes, or 
calculate stats like intersection density, average node connectivity, betweenness centrality, etc.

See the examples and demos on `GitHub`_.

|ModenaImageLink|_

.. _GitHub: https://github.com/gboeing/osmnx
.. _ModenaImageLink: http://geoffboeing.com/2016/11/osmnx-python-street-networks/
.. |ModenaImageLink| image:: https://i0.wp.com/geoffboeing.com/wp-content/uploads/2016/10/modena-italy-street-network.png?w=400
                     :alt: OSMnx Modena Italy urban city street network with Python matplotlib networkx geopandas
"""

# list of classifiers from the PyPI classifiers trove
classifiers=['Development Status :: 5 - Production/Stable',
             'License :: OSI Approved :: MIT License',
             'Operating System :: OS Independent',
             'Intended Audience :: Science/Research',
             'Topic :: Scientific/Engineering :: GIS',
             'Topic :: Scientific/Engineering :: Visualization',
             'Topic :: Scientific/Engineering :: Physics',
             'Topic :: Scientific/Engineering :: Information Analysis',
             'Programming Language :: Python',
             'Programming Language :: Python :: 2',
             'Programming Language :: Python :: 3',
             'Programming Language :: Python :: 2.7',
             'Programming Language :: Python :: 3.4',
             'Programming Language :: Python :: 3.5',
             'Programming Language :: Python :: 3.6']


# now call setup
setup(name='osmnx',
      version='0.2.1',
      description='Retrieve, construct, analyze, and visualize street networks from OpenStreetMap',
      long_description=long_description,
      classifiers=classifiers,
      url='https://github.com/gboeing/osmnx',
      author='Geoff Boeing',
      author_email='gboeing@berkeley.edu',
      license='MIT',
      platforms='any',
      data_files = [('', ['LICENSE.txt'])],
      packages=['osmnx'],
      install_requires=['requests>=2.11',
                        'numpy>=1.11',
                        'pandas>=0.19',
                        'geopandas>=0.2.1',
                        'networkx>=1.11',
                        'matplotlib>=1.5',
                        'Shapely>=1.5',
                        'descartes>=1.0',
                        'geopy>=1.11',
                        'Rtree>=0.8.3'])
