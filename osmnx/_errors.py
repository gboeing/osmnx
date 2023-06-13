"""Define custom errors and exceptions."""


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
