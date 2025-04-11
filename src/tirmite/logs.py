"""
Logging configuration for the tirmite package.

This module provides functionality to initialize and configure logging with rich
formatting for better readability in terminal output. It uses the 'rich' library
to create visually enhanced log messages.
"""

import logging

from rich.console import Console
from rich.logging import RichHandler


def init_logging(loglevel: str = 'DEBUG') -> None:
    """
    Initialize root logger with specified log level and rich formatting.

    Configures the global logging system with rich formatting for all modules
    that use the standard logging calls.

    Parameters
    ----------
    loglevel : str, optional
        The log level to use (e.g., "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
        by default "DEBUG".

    Returns
    -------
    None
        This function configures the global logging system and doesn't return a value.

    Raises
    ------
    ValueError
        If the provided log level is invalid.
    """
    # Convert log level string to numeric value
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {loglevel}')

    # Get the root logger
    root_logger = logging.getLogger()

    # Clear existing handlers if any are present
    if root_logger.hasHandlers():
        root_logger.handlers.clear()

    # Add rich handler that outputs to stderr for better visibility in scripts
    root_logger.addHandler(RichHandler(console=Console(stderr=True)))

    # Set the logger's level according to specified loglevel
    root_logger.setLevel(numeric_level)
