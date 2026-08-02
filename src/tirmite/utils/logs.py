"""
Logging configuration for the tirmite package.

Initialises the ``tirmite`` logger with rich console formatting and an
optional log file. Every module in the package logs through a named child
logger (``logging.getLogger(__name__)``), so configuring the one parent here
configures all of them.

Notes
-----
This deliberately configures the ``tirmite`` logger rather than the root
logger. Configuring root would clear and replace the handlers of whatever
application imported TIRmite, which is hostile to library use: merely calling
a ``main()`` function would silently redirect the host's own logging. Because
the package's loggers propagate by default, records still reach any handler
the host has installed on root.
"""

import logging
from pathlib import Path
from typing import Optional, Union

from rich.console import Console
from rich.logging import RichHandler

# The parent logger for the whole package. Every module logger is a child of
# this one by virtue of being named tirmite.<something>.
PACKAGE_LOGGER_NAME = 'tirmite'

logger = logging.getLogger(__name__)


def init_logging(
    loglevel: str = 'DEBUG', logfile: Optional[Union[str, Path]] = None
) -> None:
    """
    Configure the ``tirmite`` logger with the given level and rich formatting.

    Parameters
    ----------
    loglevel : str, optional
        The log level to use (e.g., "DEBUG", "INFO", "WARNING", "ERROR",
        "CRITICAL"), by default "DEBUG".
    logfile : str or Path, optional
        Path to log file. If provided, logs are written to both console and
        file.

    Returns
    -------
    None
        This function configures the package logger and returns nothing.

    Raises
    ------
    ValueError
        If the provided log level is invalid.

    Notes
    -----
    Safe to call more than once: existing handlers on the package logger are
    removed first, so repeated calls do not duplicate output. Handlers
    installed by the host application on the root logger are left untouched.
    """
    # Convert log level string to numeric value
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {loglevel}')

    package_logger = logging.getLogger(PACKAGE_LOGGER_NAME)

    # Clear only OUR handlers, so calling init_logging twice does not double
    # every message. The root logger's handlers are none of our business.
    for handler in list(package_logger.handlers):
        package_logger.removeHandler(handler)
        handler.close()

    package_logger.setLevel(numeric_level)

    # Propagation stays enabled so that a host application's root handlers,
    # and pytest's caplog, still see TIRmite's records.
    package_logger.propagate = True

    # Rich console output goes to stderr, keeping stdout free for data.
    console_handler = RichHandler(console=Console(stderr=True))
    console_handler.setLevel(numeric_level)
    package_logger.addHandler(console_handler)

    if logfile:
        try:
            logfile_path = Path(logfile)
            logfile_path.parent.mkdir(parents=True, exist_ok=True)

            file_handler = logging.FileHandler(logfile_path, mode='w', encoding='utf-8')
            file_handler.setLevel(numeric_level)
            file_handler.setFormatter(
                logging.Formatter(
                    fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                )
            )
            package_logger.addHandler(file_handler)

            logger.info(f'Logging initialized with level: {loglevel}')
            logger.info(f'Log file: {logfile_path.absolute()}')

        except OSError as e:
            # A failed log file must not abort the run; carry on with console.
            logger.warning(f'Failed to create log file {logfile}: {e}')
            logger.warning('Continuing with console-only logging')
    else:
        logger.info(f'Logging initialized with level: {loglevel}')
