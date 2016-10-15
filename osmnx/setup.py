# OSMnx
# See full license in LICENSE.txt

from setuptools import setup, find_packages

setup(name='OSMnx',
      version='0.1a1dev',
      description='Retrieve and construct spatial geometries and street networks from OpenStreetMap',
      url='https://github.com/gboeing/osmnx',
      author='Geoff Boeing',
      author_email='gboeing@berkeley.edu',
      license='MIT',
      packages=find_packages(),
      install_requires=['requests>=2.11',
                        'numpy>=1.11',
                        'pandas>=0.19',
                        'geopandas>=0.2',
                        'networkx>=1.11',
                        'matplotlib>=1.5',
                        'shapely>=1.5',
                        'descartes>=1.0',
                        'geopy>=1.11'])