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

    Boeing, G. (2017). OSMnx: New Methods for Acquiring, Constructing,
    Analyzing, and Visualizing Complex Street Networks. Computers, Environment
    and Urban Systems, 65, 126-139.
    https://doi.org/10.1016/j.compenvurbsys.2017.05.004

    Parameters
    ----------
    style : string {"apa", "bibtex", "ieee"}
        citation format, either APA or BibTeX or IEEE

    Returns
    -------
    None
    """
    if style == "apa":
        msg = (
            "Boeing, G. (2017). OSMnx: New Methods for Acquiring, Constructing, "
            "Analyzing, and Visualizing Complex Street Networks. Computers, "
            "Environment and Urban Systems, 65, 126-139. "
            "https://doi.org/10.1016/j.compenvurbsys.2017.05.004"
        )
    elif style == "bibtex":
        msg = (
            "@article{boeing_osmnx_2017,\n"
            "    title = {{OSMnx}: {New} {Methods} for {Acquiring}, "
            "{Constructing}, {Analyzing}, and {Visualizing} {Complex} "
            "{Street} {Networks}},\n"
            "    volume = {65},\n"
            "    doi = {10.1016/j.compenvurbsys.2017.05.004},\n"
            "    number = {126-139},\n"
            "    journal = {Computers, Environment and Urban Systems},\n"
            "    author = {Boeing, Geoff},\n"
            "    year = {2017},\n"
            "    pages = {126--139}\n"
            "}"
        )
    elif style == "ieee":
        msg = (
            'G. Boeing, "OSMnx: New Methods for Acquiring, Constructing, '
            'Analyzing, and Visualizing Complex Street Networks," Computers, '
            "Environment and Urban Systems, vol. 65, pp. 126-139, 2017, "
            "doi: 10.1016/j.compenvurbsys.2017.05.004."
        )
    else:  # pragma: no cover
        err_msg = f"unrecognized citation style {style!r}"
        raise ValueError(err_msg)

    print(msg)  # noqa: T201


def ts(style: str = "datetime", template: str | None = None) -> str:
    """
    Return current local timestamp as a string.

    Parameters
    ----------
    style : string {"datetime", "date", "time"}
        format the timestamp with this built-in style
    template : string
        if not None, format the timestamp with this format string instead of
        one of the built-in styles

    Returns
    -------
    ts : string
        local timestamp string
    """
    if template is None:
        if style == "datetime":
            template = "{:%Y-%m-%d %H:%M:%S}"
        elif style == "date":
            template = "{:%Y-%m-%d}"
        elif style == "time":
            template = "{:%H:%M:%S}"
        else:  # pragma: no cover
            msg = f"unrecognized timestamp style {style!r}"
            raise ValueError(msg)

    return template.format(dt.datetime.now().astimezone())


def log(
    message: str, level: int | None = None, name: str | None = None, filename: str | None = None
) -> None:
    """
    Write a message to the logger.

    This logs to file and/or prints to the console (terminal), depending on
    the current configuration of settings.log_file and settings.log_console.

    Parameters
    ----------
    message : string
        the message to log
    level : int
        one of Python's logger.level constants
    name : string
        name of the logger
    filename : string
        name of the log file, without file extension

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
                os.dup2(sys.stdout._original_stdstream_copy, sys.__stdout__.fileno())  # type: ignore[attr-defined]
                sys.stdout._original_stdstream_copy = None  # type: ignore[attr-defined]
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
    name : string
        name of the logger
    filename : string
        name of the log file, without file extension

    Returns
    -------
    logger : logging.Logger
    """
    logger = lg.getLogger(name)

    # if a logger with this name is not already set up with a handler
    if not len(logger.handlers) > 0:
        # make log filepath and create parent folder if it doesn't exist
        filepath = Path(settings.logs_folder) / f'{filename}_{ts(style="date")}.log'
        filepath.parent.mkdir(parents=True, exist_ok=True)

        # create file handler and log formatter and set them up
        handler = lg.FileHandler(filepath, encoding="utf-8")
        handler.setLevel(lg.DEBUG)
        handler.setFormatter(lg.Formatter("%(asctime)s %(levelname)s %(name)s %(message)s"))
        logger.addHandler(handler)
        logger.setLevel(lg.DEBUG)

    return logger
