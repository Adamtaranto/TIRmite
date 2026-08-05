"""Command-line interface for TIRmite.

Each module in this subpackage implements one ``tirmite`` subcommand and
exposes three things: an ``add_<name>_parser`` function that registers the
subcommand on the top-level parser, a ``create_<name>_parser`` function that
builds a standalone parser for direct invocation, and a ``main(args)`` entry
point.  :mod:`tirmite.cli.cli` wires them together.

This module deliberately re-exports nothing.  The subcommand modules are
imported lazily inside :func:`tirmite.cli.cli.create_parser` so that starting
the CLI does not pay the import cost of every subcommand's dependencies.
"""
