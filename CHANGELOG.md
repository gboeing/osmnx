# Change log

## 1.2.2 (2022-08-05)

  - fix compatibility with rasterio 1.3
  - fix API version when saving OSM XML
  - resolve shapely deprecation warning

## 1.2.1 (2022-06-16)

  - fix rate limit checking and pausing on newest versions of Overpass API
  - allow add_edge_lengths function to be run on a subset of edges
  - resolve pandas deprecation warning

## 1.2.0 (2022-05-23)

  - add ability to load GraphML string data to the load_graphml function
  - add "reversed" edge attribute to support node-order-dependent edge attributes
  - add new edge_color and edge_linewidth arguments to plot_footprints function
  - fix nearest_edges function selecting arbitrary edge when bounding boxes overlap
  - fix get_digraph function's parallel edge handling
  - fix pandas and geopandas version compatibility
  - fix log output appearing in Jupyter notebooks on Unix-like systems
  - remove old functions and arguments previously deprecated in v1.1
  - deprecate utils.config function in favor of using settings module directly

## 1.1.2 (2021-11-17)
  
  - fix geocoding when no geojson is returned
  - fix graph simplification to properly handle travel_time edge attributes
  - fix streets per node not being calculated when clean_periphery=False
  - allow user-defined aggregation function when imputing missing edge speeds
  - allow user to configure requests package keyword arguments when connecting to APIs
  - faster graph projection by calculating UTM zone number with a computationally cheaper method
  - improve efficiency of quadrat-based geometry cutting
  - fall back on google dns resolution when necessary if using a proxy
  - move count_streets_per_node function to stats module
  - resolve shapely and geopandas deprecation warnings

## 1.1.1 (2021-05-19)

  - fix overpass status endpoint checks with explicit IP address resolution
  - fix slot management on local overpass instances by optionally disabling rate limiting
  - parallelize shortest_path calculation for multiple origins/destinations

## 1.1.0 (2021-05-01)

  - add graph-constrained spatial sampling function
  - add add_node_elevations_raster function to add node elevations from local raster file(s)
  - add add_node_elevations_google function and deprecate old add_node_elevations function
  - add faster streamlined nearest_nodes and nearest_edges functions to distance module
  - deprecate old get_nearest_node, get_nearest_nodes, get_nearest_edge, and get_nearest_edges
  - add utils_geo.interpolate_points function and deprecate redistribute_vertices in favor of it
  - add vectorized calculate_bearing function and deprecate get_bearing in favor of it
  - expose individual street network stats functions in stats module
  - deprecate the extended_stats function in stats module
  - add network orientation and entropy stats functions to bearing module
  - add plot_orientation function to bearing module to polar histograms of graph edge bearings
  - add route_linewidths parameter to plot_graph_routes function
  - handle relations of type "boundary" in geometries module
  - multi-index GeoDataFrames returned from geometries module by element type and osmid
  - ensure all nodes have integer IDs after graph intersection consolidation
  - vectorize add_edge_lengths, add_edge_grades, and add_edge_bearings functions
  - improve save_graph_xml speed
  - improve geocoder module error messages
  - improve handling of node geometry when converting graph to/from GeoDataFrames
  - fix network_type filters allowing ways tagged "bus_guideway"
  - fix handling of boolean type conversion in load_graphml
  - fix truncate_graph_dist retaining unreachable nodes
  - fix bug in consolidate_intersections when pygeos is installed
  - move add_edge_lengths function from utils_graph to distance module
  - remove descartes dependency in line with geopandas

## 1.0.1 (2021-01-13)

  - fix network_type filters allowing ways tagged "planned"
  - fix "drive" network_type allowing some alleys
  - fix intersection consolidation for compatibility with v1.0 node ids/indexing
  - fix python 3.6 compatibility
  - deprecate folium polyline styling arguments

## 1.0.0 (2021-01-01)

  - set use_cache=True by default
  - add ability to query a place by OSM ID in geocoder.geocode_to_gdf function
  - add optional setting for download/cache-only mode
  - replace md5 with sha1 for cache filename hashing
  - replace streets_per_node graph attribute with equivalent street_count node attribute
  - remove redundant osmid node attribute
  - make graph_to_gdfs multi-index the edges GeoDataFrame by u, v, key
  - refactor consolidate_intersections function for better speed and efficiency
  - refactor count_streets_per_node function for better speed and efficiency
  - refactor folium module for better speed and efficiency
  - refactor get_undirected functionality for better speed and efficiency
  - extract all private/internal .osm XML functionality into new osm_xml module
  - deprecate io.save_graph_xml with warning (function moved to osm_xml module)
  - remove internal _is_simplified function
  - remove deprecated pois module
  - remove deprecated footprints module
  - remove deprecated utils_graph.induce_subgraph function
  - remove deprecated node_type parameter from io.load_graphml function

## 0.16.2 (2020-11-17)

  - improve graph_from_gdfs speed and efficiency
  - improve plot_route_folium speed and efficiency
  - fix remove_isolated_nodes function mutating the passed-in graph
  - fix gephi compatibility in save_graphml
  - add customizable node/edge attribute data type arguments to load_graphml
  - deprecate old node_type argument in load_graphml
  - expose bidirectional_network_types via config function

## 0.16.1 (2020-10-05)

  - fix handling graphs with no intersections in consolidate_intersections
  - fix consolidate_intersections returning GeoSeries without CRS attribute
  - fix response caching to save only when status code is 200
  - fix elevation module's grade absolute value calculation when grade is null
  - move shortest path functions from utils_graph module to distance module

## 0.16.0 (2020-09-07)

  - new geometries module for creating GeoDataFrames from tag/value queries
  - deprecate old pois and fooptrints modules (replaced by geometries module)
  - auto-select first Polygon/MultiPolygon when geocoding with which_result=None
  - new k_shortest_paths function to solve *k* shortest paths from origin to destination
  - new shortest_path convenience function
  - new get_digraph function to correctly convert MultiDiGraph to DiGraph
  - miscellaneous performance improvements and optimizations
  - deprecate induce_subgraph function
  - remove deprecated boundaries module (replaced by geocoder module in v0.15.0)
  - remove deprecated utils_geo.geocode function (replaced by geocoder.geocode function in v0.15.0)

## 0.15.1 (2020-07-03)

  - fix geopandas future warnings

## 0.15.0 (2020-06-30)

  - improve plotting defaults and streamline plot module speed and efficiency
  - improve color handling in plot module
  - improve route plotting
  - plot_graph_routes function now accepts multiple route colors
  - allow multiple elevation API providers
  - consolidate_intersections replaces update_edge_lengths param with reconnect_edges param
  - fix geopackage file saving after consolidating intersections
  - add new geocoder module and move utils_geo.geocode function into it
  - replace gdf_from_place/s functions with geocoder.geocode_to_gdf
  - deprecate boundaries module
  - remove deprecated timeout, memory, custom_settings, and max_query_area_size function params
  - remove deprecated plotting params and plot_shape function

## 0.14.1 (2020-06-09)

  - fix simplification of graphs with long rural roads
  - reduce memory footprint of graph simplification
  - remove disconnected self-contained rings from graph by default when simplifying
  - improve speed and efficiency of project_graph, graph_to_gdfs, and graph_from_gdfs
  - improve attribute value conversion in load_graphml
  - expose precision parameter for adding bearings, elevations, speeds, and travel times
  - fix config function clobber behavior
  - fix graph periphery cleaning when clean_periphery=True but simplify=False
  - rename settings useful_tags_path to the more appropriate useful_tags_way
  - deprecate the timeout, memory, custom_settings, and max_query_area_size function params
  - the params above are now accessible via config function and settings module
  - deprecate old plot params and plot_shape function
  - remove previously deprecated infrastructure parameter in favor of custom_filter

## 0.14.0 (2020-06-03)

  - better geometry subdividing for huge OSM queries
  - better handling of maxspeed list values for simplified graphs
  - downloader only retrieves url response from cache if no server remark
  - deprecate graph creation infrastructure parameter in favor of flexible custom_filter
  - remove deprecated functions: graph_from_file, clean_intersections, gdfs_to_graph

## 0.13.0 (2020-05-25)

  - major refactor of entire package
  - clean up API and namespace
  - new consolidate_intersections function with topological option
  - new speed module to calculate graph edge speeds and travel times
  - generalize POIs module to query with a flexible tags dict
  - allow folium functions to accept FeatureGroup and kwargs
  - all graph saving functions now take a filepath argument instead of folder/filename
  - save shapefiles in single folder containing both nodes and edges
  - optionally return distance and/or geometry in nearest edge search
  - expose timeout and memory parameters in pois and footprints modules
  - define default crs via epsg code instead of proj4 string
  - update and simplify logging with timestamps
  - graph metadata: add creation date and version, remove name
  - replace inconsistent distance parameters with consistent dist parameters
  - deprecate old clean_intersections function in favor of new consolidate_intersections
  - deprecate old gdfs_to_graph function in favor of graph_from_gdfs
  - deprecate old graph_from_file function in favor of graph_from_xml
  - rename save_as_osm function -> save_graph_xml for consistency
  - rename save_load module -> io
  - remove old save_gdf_shapefile function
  - drop support for python 3.5 and lower

## 0.12.1 (2020-05-01)

  - fix handling relations with missing type tag
  - fix save_graph_geopackage handling numeric attributes
  - fix load_graphml handling elevation and grade attributes
  - improve edge finding algorithms to return edge key
  - more informative graph_from_file data load error message
  - refactor url-in-cache checking
  - add timestamp helper function
  - documentation improvements

## 0.12 (2020-04-10)

  - add ability to save graph as geopackage file
  - add truncate_by_edge implementation in truncate_graph_polygon
  - allow flexible overpass settings (e.g., to query by date)
  - better handling of invalid footprint geometries
  - geocode function now uses nominatim_request function
  - improve .osm xml output
  - improve one-way street handling
  - fix graph projection overwriting original lat/lng
  - fix redistribute_vertices function for MultiLineStrings

## 0.11.4 (2020-01-31)

  - fix .osm xml output
  - fix for pandas 1.0

## 0.11.3 (2020-01-09)

  - fix errant print statement

## 0.11.2 (2020-01-07)

  - fix .osm xml output
  - fix geopandas future compatibility

## 0.11.1 (2020-01-01)

  - fix get_nearest_edges search when not using a spatial index

## 0.11 (2019-12-04)

  - drop formal python 2 support
  - refactor all modules for cleaner package organization
  - make stats betweenness centrality compatible with networkx>=2.4
  - allow configurable overpass and nominatim endpoints
  - allow gdf_from_places to take a which_result list argument
  - handle zero-division in street grade calculation
  - better footprint relation handling
  - improve network type queries for better filtering
  - fix pois_from_polygon returning points outside polygon

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

## 0.3 (2017-01-29)

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
