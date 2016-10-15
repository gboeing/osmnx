# OSMnx

Retrieve and construct spatial geometries and street networks from OpenStreetMap

OSMnx is a Python library that lets you download spatial geometries and construct, project, visualize, 
and analyze street networks from OpenStreetMap's API. Users can construct walkable, drivable, or bikable 
urban networks with a single line of Python code. [Overview demo](examples/01-overview-osmnx.ipynb).

### Use OSMnx to create place boundary shapefiles from OpenStreetMap

OSMnx lets you download spatial "place boundary" geometries from OpenStreetMap, save them to shapefiles, 
project them, and plot them. For a more in-depth demonstration of creating these shapefiles, 
see [this notebook](examples/02-example-osm-to-shapefile.ipynb).

### Use OSMnx to download and construct street networks

OSMnx lets you download street network data and build topologically-corrected street networks, project to UTM and plot the 
networks, and save the street network as SVGs, GraphML files, or shapefiles for later use. The street networks are 
directed and preserve one-way directionality. For a more in-depth demonstration of creating street networks, 
see [this notebook](examples/03-example-osm-place-network.ipynb).

You can download a street network by providing OSMnx any of the following (demonstrated in the examples below):
  - a bounding box
  - a lat-long point plus a distance (either distance along the network, or cardinal)
  - an address plus a distance (either distance along the network, or cardinal)
  - a place name or list of place names (for OSMnx to automatically geocode and get the boundary of)
  - a polygon of the desired street network's boundaries

You can also specify several different network types:
  - `drive` - get drivable public streets (but not service roads)
  - `drive_service` - get drivable streets, including service roads
  - `walk` - get all streets and paths that pedestrians can use (this network type ignores one-way directionality)
  - `bike` - get all streets and paths that cyclists can use
  - `all` - download all non-private OSM streets and paths (this is the default network type unless you specify a different one)
  - `all_private` - download all OSM streets and paths, including private-access ones

### Use OSMnx to correct and simplify street network topology

Simplification is normally done by OSMnx automatically under the hood, but we can break it out to see how it works. 
OpenStreetMap nodes include intersections, but they also include all the points along a single block where 
the street curves. The latter are not nodes in the graph theory sense, so we remove them algorithmically and consolidate the 
set of edges between "true" network nodes into a single edge, but retain the actual spatial geometry. There are two 
simplification modes, strict and non-strict. The main difference is that unlike strict mode, non-strict mode allows 
simplification to an "expansion graph". For a more in-depth demonstration of topological simplification with OSMnx, 
see [this notebook](examples/04-example-simplify-network.ipynb).

### Use OSMnx to save street networks to disk

OSMnx allows users to save street networks to disk as shapefiles to work with in GIS software, as GraphML files
to work with in Gephi or NetworkX, and as SVG files to work with in Illustrator. It also allows you to save place
boundary geometries as shapefiles. For more examples of saving and loading networks to/from disk, 
see [this notebook](examples/05-example-save-load-networks-shapes.ipynb).

### Use OSMnx to analyze street networks

OSMnx allows you to calculate origin-destination routes along the network and quickly visualize them. You can easily
visualize one-way streets, cul de sacs, high/low connectivity intersections, etc. OSMnx provides built-in capabilities
to quickly calculate spatial network metrics like intersection density, average intersection degree, edge density, 
average street segment length, etc. For more examples of analyzing street networks with OSMnx, 
see [this notebook](examples/06-example-osmnx-networkx.ipynb)