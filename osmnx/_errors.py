"""Define custom errors and exceptions."""


class ResponseStatusCodeError(ValueError):  # pragma: no cover
    """Exception for an unhandled server response status code."""

    def __init__(self, *args, **kwargs):
        """Create exception."""
        Exception.__init__(self, *args, **kwargs)


class CacheOnlyInterruptError(InterruptedError):  # pragma: no cover
    """Exception for settings.cache_only_mode=True interruption."""

    def __init__(self, *args, **kwargs):
        """Create exception."""
        Exception.__init__(self, *args, **kwargs)


class InsufficientResponseError(ValueError):  # pragma: no cover
    """Exception for empty or too few results in server response."""

    def __init__(self, *args, **kwargs):
        """Create exception."""
        Exception.__init__(self, *args, **kwargs)


class GraphSimplificationError(ValueError):  # pragma: no cover
    """Exception for a problem with graph simplification."""

    def __init__(self, *args, **kwargs):
        """Create exception."""
        Exception.__init__(self, *args, **kwargs)
