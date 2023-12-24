"""Define custom errors and exceptions."""


class ResponseStatusCodeError(ValueError):
    """Exception for an unhandled server response status code."""


class CacheOnlyInterruptError(InterruptedError):
    """Exception for settings.cache_only_mode=True interruption."""


class InsufficientResponseError(ValueError):
    """Exception for empty or too few results in server response."""


class GraphSimplificationError(ValueError):
    """Exception for a problem with graph simplification."""
