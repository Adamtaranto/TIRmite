"""Shared utilities: sequence access, logging setup, and file handling.

- :mod:`tirmite.utils.extract` abstracts over the two ways TIRmite can read
  sequence (an indexed FASTA via ``pyfaidx``, or a BLAST database via
  ``blastdbcmd``) behind a common interface.
- :mod:`tirmite.utils.logs` configures logging for the CLI entry points.
- :mod:`tirmite.utils.utils` holds temporary-directory management, input
  validation and gzip handling.

This module deliberately re-exports nothing, to keep the import graph flat.
"""
