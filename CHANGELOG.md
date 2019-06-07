# Change log

## 0.11 (T.B.D.)

  - drop formal python 2 support

## 0.10 (2019-05-08)

  - remove deprecated buildings module
  - filter steps ways out of bike queries
  - convert CRS-handling to proj4 strings
  - save graph to xml-formatted .osm file
  - minor refactoring

## 0.9 (2019-01-28)

  - deprecate buildings module and replace with generalized footprints module
  - improve handling of multipolygon footprints
  - new function to find nearest edge(s), given coordinates
  - add "search," "reverse," and "lookup" nominatim queries
  - use unprojected graphs for figure-ground plotting functions
  - allow non-integer osmid values for custom data
  - improve get_route_edge_attributes function
  - improve color mapping by node/edge attribute value
  - make bidirectional network types explicit
  - networkx compatibility fixes to resolve warnings

## 0.8.2 (2018-09-19)

  - add python 3.7 compatibility
  - add convenience function to plot several routes over the same map
  - optimize graph truncation to bounding box
  - give self-loops a null bearing when calculating edge bearings
  - make accept-language http header explicit and configurable
  - add citation function
  - refactor POI module

## 0.8.1 (2018-05-17)

  - add Gephi compatibility argument for saving GraphML
  - handle square bracket encapsulated strings when loading GraphML

## 0.8 (2018-05-05)

  - add ability to retrieve points of interest
  - improve performance for retrieving huge geographies' street networks
  - fix building footprint retrieval query syntax
  - minor bug fixes

## 0.7.4 (2018-04-05)

  - add fast nearest-nodes search
  - allow custom network query filters
  - allow create_graph to return graph with no edges
  - improve figure_ground joint smoothing
  - fix handling of parallel edges when making multidigraph undirected
  - generalize same-geometry checker
  - improve detection of prior topology simplification
  - custom error types for finer-grained handling

## 0.7.3 (2018-03-12)

  - turn off x- and y-axes to improve plotting appearance
  - make floating-point precision and rounding more sensible
  - improve OS path handling cross-platform
  - replace great-circle distance calculator with haversine
  - add access filter as configurable setting
  - improve performance of inducing subgraphs
  - fix utils.get_largest_component for networkx 2.2 compatibility
  - fix config settings namespacing

## 0.7.2 (2018-02-15)

  - compatibility with networkx 2.1

## 0.7.1 (2018-02-04)

  - fix documentation build
  - ignore ways marked access=no

## 0.7 (2018-02-01)

  - ability to load a graph from a .osm file
  - change datum from NAD83 to WGS84
  - make roundabouts one-way
  - conformal plotting for unprojected graphs
  - fix folium web maps rendering

## 0.6 (2017-10-02)

  - migrate to the networkx 2.0 API

## 0.5.4 (2017-09-16)

  - add optional cleaned intersections count to basic stats
  - allow circuity to be calculated for projected or unprojected networks
  - various code clean-up and refactoring

## 0.5.3 (2017-07-22)

  - add requirements files to distribution

## 0.5.2 (2017-07-22)

  - add ability to download other infrastructures besides just roads/paths (e.g., rail lines, power lines, etc.)
  - calculate graph edges' bearings
  - add ability to get nearest node by great circle or euclidean distance
  - move examples/demo notebooks to new repo: osmnx-examples
  - fix docstrings
  - fix building footprint downloads that require multiple calls for large areas
  - fix missing MultiPolygon import in buildings module

## 0.5.1 (2017-05-12)

  - functionality to clean-up and consolidate complex intersections
  - let save_gdf_shapefile save building footprint GeoDataFrames
  - set node color correctly in figure-ground diagrams

## 0.5 (2017-04-25)

  - add elevation module to get node elevations and street grades
  - new color sequence creation and conversion functions in plot module
  - new function to get a path's edge attribute values
  - gracefully handle subpolygons that are invalid or have zero area
  - make truncate_graph_polygon work on projected graphs
  - plot_shape accepts a color or a list of colors
  - make all requests to Overpass API set custom user-agent and referer
  - rewrite algorithms to convert multidigraphs to multigraphs

## 0.4.1 (2017-04-01)

  - fix load_graphml so we can save a graph again after loading it
  - fix load_graphml so edge oneway attribute is not always set to True
  - buildings module gets buildings stored in OSM as relations as well as ways
  - fix figure-ground diagram saving to make perfect square and smooth joints
  - add optional graph argument to plot_figure_ground
  - suppress jupyter notebook deprecation warnings

## 0.4 (2017-03-01)

  - plot entire networks with folium
  - plot routes on top of networks with folium
  - vectorize all great circle calculations
  - new geocode function in utils
  - remove geopy dependency
  - refactor modules
  - simplify before truncating by distance when getting graph by point and network distance
  - project geometries, GeoDataFrames, and graphs to a passed-in CRS

## 0.3.1 (2017-02-15)

  - clean up docstrings throughout
  - remove network code vestiges from buildings.py

## 0.3 (2016-01-29)

  - add route plotting with folium
  - add downloading and visualization of building footprints
  - updates for compatibility with matplotlib 2.0

## 0.2.2 (2017-01-20)

  - fixes for compatibility with networkx 2.0's new API
  - make png default image save format
  - figure-ground plots collect street network from a wider area

## 0.2.1 (2017-01-11)

  - add license file to dist package

## 0.2 (2017-01-10)

  - refactor modules
  - add graph to GDF and GDF to graph functions
  - add encoding argument to save_graph_shapefile
  - add unit tests and continuous integration

## 0.1 (2016-12-19)

  - add street width attribute for ways from OSM

## 0.1b2 (2016-11-29)

  - make simplification error messages explicit

## 0.1b1 (2016-11-28)

  - process land use and area tags from OSM
  - make intersection error messages clear

## 0.1a1 (2016-11-07)

  - first pre-release
