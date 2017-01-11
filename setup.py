# OSMnx
# See full license in LICENSE.txt

from setuptools import setup

setup(name='osmnx',
      version='0.2.1',
      description='Retrieve, construct, analyze, and visualize street networks from OpenStreetMap',
      url='https://github.com/gboeing/osmnx',
      author='Geoff Boeing',
      author_email='gboeing@berkeley.edu',
      license='MIT',
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
