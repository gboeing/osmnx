"""Calculate graph edge speeds and travel times."""

from warnings import warn

import numpy as np

from . import routing


def add_edge_speeds(G, hwy_speeds=None, fallback=None, precision=None, agg=np.mean):
    """
    Do not use: deprecated.

    Use the `routing.add_edge_speeds` function instead.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        deprecated, do not use
    hwy_speeds : dict
        deprecated, do not use
    fallback : numeric
        deprecated, do not use
    precision : int
        deprecated, do not use
    agg : function
        deprecated, do not use

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    msg = (
        "The `add_edge_speeds` function has moved to the `routing` module. Calling "
        "`speed.add_edge_speeds` is deprecated and will be removed in the "
        "v2.0.0 release. Call it via `routing.add_edge_speeds` instead. "
        "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
    )
    warn(msg, FutureWarning, stacklevel=2)
    return routing.add_edge_speeds(G, hwy_speeds, fallback, precision, agg)


def add_edge_travel_times(G, precision=None):
    """
    Do not use: deprecated.

    Use the `routing.add_edge_travel_times` function instead.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        deprecated, do not use
    precision : int
        deprecated, do not use

    Returns
    -------
    G : networkx.MultiDiGraph
    """
    msg = (
        "The `add_edge_travel_times` function has moved to the `routing` module. Calling "
        "`speed.add_edge_travel_times` is deprecated and will be removed in the "
        "v2.0.0 release. Call it via `routing.add_edge_travel_times` instead. "
        "See the OSMnx v2 migration guide: https://github.com/gboeing/osmnx/issues/1123"
    )
    warn(msg, FutureWarning, stacklevel=2)
    return routing.add_edge_travel_times(G, precision)
