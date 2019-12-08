import xml


class OSMContentHandler(xml.sax.handler.ContentHandler):
    """
    SAX content handler for OSM XML.

    Used to build an Overpass-like response JSON object in self.object. For format
    notes, see http://wiki.openstreetmap.org/wiki/OSM_XML#OSM_XML_file_format_notes
    and http://overpass-api.de/output_formats.html#json
    """

    def __init__(self):
        self._element = None
        self.object = {'elements': []}

    def startElement(self, name, attrs):
        if name == 'osm':
            self.object.update({k: attrs[k] for k in attrs.keys()
                                if k in ('version', 'generator')})

        elif name in ('node', 'way'):
            self._element = dict(type=name, tags={}, nodes=[], **attrs)
            self._element.update({k: float(attrs[k]) for k in attrs.keys()
                                  if k in ('lat', 'lon')})
            self._element.update({k: int(attrs[k]) for k in attrs.keys()
                                  if k in ('id', 'uid', 'version', 'changeset')})

        elif name == 'tag':
            self._element['tags'].update({attrs['k']: attrs['v']})

        elif name == 'nd':
            self._element['nodes'].append(int(attrs['ref']))

        elif name == 'relation':
            # Placeholder for future relation support.
            # Look for nested members and tags.
            pass

    def endElement(self, name):
        if name in ('node', 'way'):
            self.object['elements'].append(self._element)
