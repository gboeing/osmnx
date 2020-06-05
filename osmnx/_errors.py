"""Custom errors."""


class EmptyOverpassResponse(ValueError):  # pragma: no cover
    """Exception for empty overpass response."""

    def __init__(self, *args, **kwargs):
        """Create exception."""
        Exception.__init__(self, *args, **kwargs)
