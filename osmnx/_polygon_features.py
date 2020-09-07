"""Dictionary of tags to determine if closed ways should become polygons."""

# Note: The dictionary is based on the JSON linked to from the following page:
# https://wiki.openstreetmap.org/wiki/Overpass_turbo/Polygon_Features
_polygon_features = {
    "building": {"polygon": "all"},
    "highway": {"polygon": "passlist", "values": ["services", "rest_area", "escape", "elevator"]},
    "natural": {
        "polygon": "blocklist",
        "values": ["coastline", "cliff", "ridge", "arete", "tree_row"],
    },
    "landuse": {"polygon": "all"},
    "waterway": {"polygon": "passlist", "values": ["riverbank", "dock", "boatyard", "dam"]},
    "amenity": {"polygon": "all"},
    "leisure": {"polygon": "all"},
    "barrier": {
        "polygon": "passlist",
        "values": ["city_wall", "ditch", "hedge", "retaining_wall", "spikes"],
    },
    "railway": {
        "polygon": "passlist",
        "values": ["station", "turntable", "roundhouse", "platform"],
    },
    "area": {"polygon": "all"},
    "boundary": {"polygon": "all"},
    "man_made": {"polygon": "blocklist", "values": ["cutline", "embankment", "pipeline"]},
    "power": {"polygon": "passlist", "values": ["plant", "substation", "generator", "transformer"]},
    "place": {"polygon": "all"},
    "shop": {"polygon": "all"},
    "aeroway": {"polygon": "blocklist", "values": ["taxiway"]},
    "tourism": {"polygon": "all"},
    "historic": {"polygon": "all"},
    "public_transport": {"polygon": "all"},
    "office": {"polygon": "all"},
    "building:part": {"polygon": "all"},
    "military": {"polygon": "all"},
    "ruins": {"polygon": "all"},
    "area:highway": {"polygon": "all"},
    "craft": {"polygon": "all"},
    "golf": {"polygon": "all"},
    "indoor": {"polygon": "all"},
}
