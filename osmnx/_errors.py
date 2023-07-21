"""Define custom errors and exceptions."""


class ResponseStatusCodeError(ValueError):  # pragma: no cover
    """Exception for an unhandled server response status code."""

    def __init__(self, *args, **kwargs):
        """Create exception."""
        Exception.__init__(self, *args, **kwargs)


class CacheOnlyModeInterrupt(InterruptedError):  # pragma: no cover
    """Exception for settings.cache_only_mode=True interruption."""

    def __init__(self, *args, **kwargs):
        """Create exception."""
        Exception.__init__(self, *args, **kwargs)


class EmptyOverpassResponse(ValueError):  # pragma: no cover
    """Exception for empty Overpass API response."""

    def __init__(self, *args, **kwargs):
        """Create exception."""
        Exception.__init__(self, *args, **kwargs)


class GraphSimplificationError(ValueError):  # pragma: no cover
    """Exception for a problem with graph simplification."""

    def __init__(self, *args, **kwargs):
        """Create exception."""
        Exception.__init__(self, *args, **kwargs)
