"""General utility functions."""

from __future__ import annotations

import datetime as dt
import logging as lg
import os
import sys
import unicodedata as ud
from contextlib import redirect_stdout
from pathlib import Path

from . import settings


def citation(style: str = "bibtex") -> None:
    """
    Print the OSMnx package's citation information.

    Boeing, G. (2024). Modeling and Analyzing Urban Networks and Amenities with
    OSMnx. Working paper. https://geoffboeing.com/publications/osmnx-paper/

    Parameters
    ----------
    style
        {"apa", "bibtex", "ieee"}
        The citation format, either APA or BibTeX or IEEE.

    Returns
    -------
    None
    """
    if style == "apa":
        msg = (
            "Boeing, G. (2024). Modeling and Analyzing Urban Networks and Amenities "
            "with OSMnx. Working paper. https://geoffboeing.com/publications/osmnx-paper/"
        )
    elif style == "bibtex":
        msg = (
            "@techreport{boeing_osmnx_2024,\n"
            "    author = {Boeing, Geoff},\n"
            "    title = {{Modeling and Analyzing Urban Networks and Amenities with OSMnx}},\n"
            "    type = {Working paper},\n"
            "    url = {https://geoffboeing.com/publications/osmnx-paper/},\n"
            "    year = {2024}\n"
            "}"
        )
    elif style == "ieee":
        msg = (
            'G. Boeing, "Modeling and Analyzing Urban Networks and Amenities with OSMnx," '
            "Working paper, https://geoffboeing.com/publications/osmnx-paper/"
        )
    else:  # pragma: no cover
        err_msg = f"Invalid citation style {style!r}."
        raise ValueError(err_msg)

    print(msg)  # noqa: T201


def ts(style: str = "datetime", template: str | None = None) -> str:
    """
    Return current local timestamp as a string.

    Parameters
    ----------
    style
        {"datetime", "iso8601", "date", "time"}
        Format the timestamp with this built-in style.
    template
        If not None, format the timestamp with this format string instead of
        one of the built-in styles.

    Returns
    -------
    timestamp
    """
    if template is None:
        if style == "datetime":
            template = "{:%Y-%m-%d %H:%M:%S}"
        elif style == "iso8601":
            template = "{:%Y-%m-%dT%H:%M:%SZ}"
        elif style == "date":
            template = "{:%Y-%m-%d}"
        elif style == "time":
            template = "{:%H:%M:%S}"
        else:  # pragma: no cover
            msg = f"Invalid timestamp style {style!r}."
            raise ValueError(msg)

    return template.format(dt.datetime.now().astimezone())


def log(
    message: str,
    level: int | None = None,
    name: str | None = None,
    filename: str | None = None,
) -> None:
    """
    Write a message to the logger.

    This logs to file and/or prints to the console (terminal), depending on
    the current configuration of `settings.log_file` and
    `settings.log_console`.

    Parameters
    ----------
    message
        The message to log.
    level
        One of the Python `logger.level` constants. If None, set to
        `settings.log_level` value.
    name
        The name of the logger. If None, set to `settings.log_name` value.
    filename
        The name of the log file, without file extension. If None, set to
        `settings.log_filename` value.

    Returns
    -------
    None
    """
    if level is None:
        level = settings.log_level
    if name is None:
        name = settings.log_name
    if filename is None:
        filename = settings.log_filename

    # if logging to file is turned on
    if settings.log_file:
        # get the current logger (or create a new one, if none), then log
        # message at requested level
        logger = _get_logger(name=name, filename=filename)
        if level == lg.DEBUG:
            logger.debug(message)
        elif level == lg.INFO:
            logger.info(message)
        elif level == lg.WARNING:
            logger.warning(message)
        elif level == lg.ERROR:
            logger.error(message)

    # if logging to console (terminal window) is turned on
    if settings.log_console:
        # prepend timestamp then convert to ASCII for Windows command prompts
        message = f"{ts()} {message}"
        message = ud.normalize("NFKD", message).encode("ascii", errors="replace").decode()

        try:
            # print explicitly to terminal in case Jupyter has captured stdout
            if getattr(sys.stdout, "_original_stdstream_copy", None) is not None:
                # redirect the Jupyter-captured pipe back to original
                os.dup2(sys.stdout._original_stdstream_copy, sys.__stdout__.fileno())  # type: ignore[union-attr]
                sys.stdout._original_stdstream_copy = None  # type: ignore[union-attr]
            with redirect_stdout(sys.__stdout__):
                print(message, file=sys.__stdout__, flush=True)
        except OSError:
            # handle pytest on Windows raising OSError from sys.__stdout__
            print(message, flush=True)  # noqa: T201


def _get_logger(name: str, filename: str) -> lg.Logger:
    """
    Create a logger or return the current one if already instantiated.

    Parameters
    ----------
    name
        Name of the logger.
    filename
        Name of the log file, without file extension.

    Returns
    -------
    logger
    """
    logger = lg.getLogger(name)

    # if a logger with this name is not already set up with a handler
    if len(logger.handlers) == 0:
        # make log filepath and create parent folder if it doesn't exist
        filepath = Path(settings.logs_folder) / f"{filename}_{ts(style='date')}.log"
        filepath.parent.mkdir(parents=True, exist_ok=True)

        # create file handler and log formatter and set them up
        handler = lg.FileHandler(filepath, encoding="utf-8")
        handler.setLevel(lg.DEBUG)
        handler.setFormatter(lg.Formatter("%(asctime)s %(levelname)s %(name)s %(message)s"))
        logger.addHandler(handler)
        logger.setLevel(lg.DEBUG)

    return logger
