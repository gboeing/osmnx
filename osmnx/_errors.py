"""Define custom errors and exceptions."""


class CacheOnlyInterruptError(InterruptedError):
    """Exception for `settings.cache_only_mode=True` interruption."""


class GraphSimplificationError(ValueError):
    """Exception for a problem with graph simplification."""


class ValidationError(ValueError):
    """Exception for failed graph or GeoDataFrame validation."""


class InsufficientResponseError(ValueError):
    """Exception for empty or too few results in server response."""


class ResponseStatusCodeError(ValueError):
    """Exception for an unhandled server response status code."""
