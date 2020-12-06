"""Custom errors."""


class CacheOnlyModeInterrupt(InterruptedError):  # pragma: no cover
    """Exception for cache_only_mode interruption."""

    def __init__(self, *args, **kwargs):
        """Create exception."""
        Exception.__init__(self, *args, **kwargs)


class EmptyOverpassResponse(ValueError):  # pragma: no cover
    """Exception for empty overpass response."""

    def __init__(self, *args, **kwargs):
        """Create exception."""
        Exception.__init__(self, *args, **kwargs)
