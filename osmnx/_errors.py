"""Custom errors."""


class EmptyOverpassResponse(ValueError):  # pragma: no cover
    """Exception for empty overpass response."""

    def __init__(self, *args, **kwargs):
        """Create exception."""
        Exception.__init__(self, *args, **kwargs)


class InsufficientNetworkQueryArguments(ValueError):  # pragma: no cover
    """Exception for insufficient network query args."""

    def __init__(self, *args, **kwargs):
        """Create exception."""
        Exception.__init__(self, *args, **kwargs)


class InvalidDistanceType(ValueError):  # pragma: no cover
    """Exception for invalid distance type."""

    def __init__(self, *args, **kwargs):
        """Create exception."""
        Exception.__init__(self, *args, **kwargs)


class UnknownNetworkType(ValueError):  # pragma: no cover
    """Exception for unknown network type."""

    def __init__(self, *args, **kwargs):
        """Create exception."""
        Exception.__init__(self, *args, **kwargs)
