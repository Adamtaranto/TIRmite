"""MAFFT multiple-sequence-alignment wrappers.

TIRmite aligns sequences in two quite different situations, and they want
different things from MAFFT:

- :func:`align_to_file` is used when building an HMM. It aligns many extracted
  element copies, uppercases the result (HMMER treats lowercase as masked), and
  persists the alignment because it is both an input to ``hmmbuild`` and a
  user-facing output.
- :func:`align_in_memory` is used when validating a junction. It aligns a
  handful of sequences purely to inspect the gap pattern, and never needs the
  alignment on disk.

Both previously existed as separate ``run_mafft_alignment`` functions in
:mod:`tirmite.cli.hmm_build` and :mod:`tirmite.cli.validate`, with the same
name but different signatures. They now share :func:`_run_mafft` for process
handling while keeping the two distinct result contracts.
"""

import io
import logging
import os
from pathlib import Path
import shutil
import subprocess
from typing import List, Optional, Sequence, Union

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Library modules acquire a named logger and attach a NullHandler, so that
# importing TIRmite as a library emits nothing until the host application
# configures logging. Handler setup belongs to the CLI, not here.
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class MafftError(Exception):
    """Raised when MAFFT is unavailable, fails, or emits unparseable output."""


def mafft_available() -> bool:
    """
    Report whether the ``mafft`` executable is on ``PATH``.

    Returns
    -------
    bool
        True if ``mafft`` can be located.
    """
    return shutil.which('mafft') is not None


def _run_mafft(
    sequences: Sequence[SeqRecord],
    input_path: Union[str, Path],
    extra_args: Sequence[str],
) -> str:
    """
    Write sequences to a FASTA file and run MAFFT over it.

    Parameters
    ----------
    sequences : sequence of Bio.SeqRecord.SeqRecord
        Sequences to align. At least two are required.
    input_path : str or Path
        Path to write the MAFFT input FASTA to. The caller owns this file and
        is responsible for removing it.
    extra_args : sequence of str
        MAFFT flags to insert before the input filename, e.g.
        ``['--thread', '4']`` or ``['--auto', '--quiet']``.

    Returns
    -------
    str
        MAFFT's standard output, which is the alignment in FASTA format.

    Raises
    ------
    MafftError
        If MAFFT is not installed, fewer than two sequences were supplied, or
        MAFFT exits non-zero.

    Notes
    -----
    MAFFT writes the alignment to stdout, so output is captured rather than
    redirected to a file. That keeps this function agnostic about whether the
    caller ultimately wants a file or in-memory records.
    """
    if not mafft_available():
        raise MafftError('mafft not found in PATH. Please install MAFFT.')

    if len(sequences) < 2:
        raise MafftError(
            f'Need at least 2 sequences for alignment, got {len(sequences)}'
        )

    input_path = Path(input_path)
    with open(input_path, 'w') as handle:
        SeqIO.write(sequences, handle, 'fasta')

    cmd = ['mafft', *extra_args, str(input_path)]
    logger.info(f'Running MAFFT alignment on {len(sequences)} sequences')
    logger.debug(f'MAFFT command: {" ".join(cmd)}')

    result = subprocess.run(cmd, capture_output=True, text=True, check=False)

    if result.returncode != 0:
        raise MafftError(f'MAFFT failed (exit {result.returncode}): {result.stderr}')

    return result.stdout


def align_to_file(
    sequences: List[SeqRecord], output_file: Path, threads: int = 1
) -> Path:
    """
    Align sequences with MAFFT and write the uppercased result to a file.

    Parameters
    ----------
    sequences : list of Bio.SeqRecord.SeqRecord
        Sequences to align.
    output_file : Path
        Path to write the alignment to, in FASTA format.
    threads : int, default 1
        Number of CPU threads to pass to MAFFT.

    Returns
    -------
    Path
        ``output_file``, for convenient chaining.

    Raises
    ------
    MafftError
        If MAFFT is unavailable, fewer than two sequences were supplied, or
        MAFFT fails.

    Notes
    -----
    The alignment is uppercased before being written. MAFFT preserves the case
    of its input, and ``hmmbuild`` interprets lowercase residues as masked, so
    a soft-masked genomic region would otherwise silently drop out of the
    resulting model.
    """
    # MAFFT input lands beside the output so both share the caller's chosen
    # directory, which is typically a managed temporary directory.
    temp_fasta = output_file.parent / f'{output_file.stem}_input.fasta'

    try:
        stdout = _run_mafft(sequences, temp_fasta, ['--thread', str(threads)])

        alignment_records = [
            SeqRecord(
                Seq(str(record.seq).upper()),
                id=record.id,
                description=record.description,
            )
            for record in SeqIO.parse(io.StringIO(stdout), 'fasta')
        ]

        if not alignment_records:
            raise MafftError('MAFFT produced no alignment records')

        with open(output_file, 'w') as outfile:
            SeqIO.write(alignment_records, outfile, 'fasta')

        logger.info(
            f'Alignment of {len(alignment_records)} sequences written to {output_file}'
        )
        return output_file

    finally:
        # Remove the MAFFT input whether or not the alignment succeeded.
        if temp_fasta.exists():
            temp_fasta.unlink()


def align_in_memory(
    sequences: List[SeqRecord],
    tmpdir: str,
) -> Optional[List[SeqRecord]]:
    """
    Align sequences with MAFFT and return the aligned records.

    Parameters
    ----------
    sequences : list of Bio.SeqRecord.SeqRecord
        Sequences to align. The first is conventionally the query.
    tmpdir : str
        Directory for MAFFT's intermediate input and output files.

    Returns
    -------
    list of Bio.SeqRecord.SeqRecord or None
        The aligned sequences, or None if alignment failed for any reason.

    Notes
    -----
    Returns None rather than raising, because the validation workflow treats a
    failed alignment as an inconclusive result for one candidate junction and
    carries on with the rest rather than aborting the run.
    """
    input_file = os.path.join(tmpdir, 'mafft_input.fasta')
    output_file = os.path.join(tmpdir, 'mafft_output.fasta')

    try:
        # --auto lets MAFFT pick a strategy suited to the (small) input size;
        # --quiet suppresses its progress chatter on stderr.
        stdout = _run_mafft(sequences, input_file, ['--auto', '--quiet'])
    except MafftError as e:
        logger.warning(f'MAFFT alignment failed: {e}')
        return None

    with open(output_file, 'w') as out_handle:
        out_handle.write(stdout)

    try:
        return list(AlignIO.read(output_file, 'fasta'))
    except Exception as e:
        logger.warning(f'Failed to parse MAFFT output: {e}')
        return None
