"""Wrappers around the external tools TIRmite shells out to.

Modules here build and execute command lines for HMMER (``hmmbuild``,
``hmmpress``, ``nhmmer``, ``hmmalign``) and BLAST+ (``makeblastdb``,
``blastn``, ``blastdbcmd``).  Command construction is kept separate from
execution so that the argument lists can be unit tested without the binaries
being installed.

This module deliberately re-exports nothing, to keep the import graph flat.
"""
