"""Core analysis logic, independent of any command-line interface.

Modules here hold the algorithms that the ``tirmite`` subcommands orchestrate:
hit parsing and filtering, sequence extraction, flank and target-site
reconstruction, and the terminus pairing engine.

Nothing in this subpackage may import from :mod:`tirmite.cli`. The dependency
runs one way only -- ``cli`` builds on ``core``, ``core`` builds on ``utils``
and ``runners`` -- which is what keeps the import graph acyclic.

This module deliberately re-exports nothing. Adding re-exports here would make
``import tirmite.core.anything`` execute every core module, which is the most
likely way to introduce a circular import.
"""
