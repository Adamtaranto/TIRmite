"""
BLAST wrappers for sequence alignment in TIR identification.

Provides subprocess-based BLAST execution with:
- Safety-focused command construction (avoids shell=True)
- Configurable or unlimited timeouts
- Periodic progress logging for long-running searches
- Batch processing capabilities
- Self-alignment support for TIR detection

All functions use Path objects and avoid shell injection vulnerabilities.
"""

import logging
import os
from pathlib import Path
import subprocess
import time
from typing import List, Optional, Union

# Library modules acquire a named logger and attach a NullHandler, so that
# importing TIRmite as a library emits nothing until the host application
# configures logging. Handler setup belongs to the CLI, not here.
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class BlastError(Exception):
    """
    Custom exception for BLAST-related errors.
    """

    pass


def blast_db_exists(db_prefix: Union[str, Path]) -> bool:
    """
    Check if BLAST database files exist for the given prefix.

    Parameters
    ----------
    db_prefix : str or Path
        Path prefix for BLAST database files (without extension).

    Returns
    -------
    bool
        True if any of the expected BLAST database files exist, False otherwise.

    Notes
    -----
    Checks for BLAST+ nucleotide database extensions: .nhr, .nin, .nsq,
    .ndb, .not, .ntf, .nto.
    """
    db_extensions = ['.nhr', '.nin', '.nsq', '.ndb', '.not', '.ntf', '.nto']
    return any(Path(f'{db_prefix}{ext}').exists() for ext in db_extensions)


def check_blast_available() -> bool:
    """
    Check if blastn executable is available in system PATH.

    Returns
    -------
    bool
        True if blastn is found and executable, False otherwise.

    Notes
    -----
    Tests blastn availability by running 'blastn -version' with a 10-second timeout.
    """
    try:
        result = subprocess.run(
            ['blastn', '-version'], capture_output=True, text=True, timeout=10
        )
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False


def run_blastn(
    query: Union[str, Path],
    subject: Union[str, Path],
    output: Union[str, Path],
    word_size: int = 4,
    perc_identity: float = 60.0,
    outfmt: str = '6 qstart qend sstart send length positive pident qlen slen qframe sframe qseqid sseqid',
    additional_args: Optional[List[str]] = None,
    verbose: bool = False,
    num_threads: int = 1,
    timeout: Optional[int] = None,
) -> subprocess.CompletedProcess:
    """
    Execute blastn with specified parameters using subprocess.

    Parameters
    ----------
    query : str or Path
        Path to query sequence file (FASTA format).
    subject : str or Path
        Path to subject sequence file (FASTA format) or BLAST database prefix.
        If a pre-built BLAST database is provided (prefix without extension),
        the ``-db`` flag is used instead of ``-subject``.
    output : str or Path
        Path to output file for BLAST results.
    word_size : int, default 4
        Word size for initial matches (smaller values increase sensitivity).
    perc_identity : float, default 60.0
        Minimum percent identity threshold for reporting alignments.
    outfmt : str, default '6 qstart qend sstart send length positive pident qlen slen qframe sframe qseqid sseqid'
        Output format string for tabular results.
    additional_args : list of str, optional
        Additional command-line arguments to pass to blastn.
    verbose : bool, default False
        If True, logs extra detail including full output to console.
    num_threads : int, default 1
        Number of CPU threads for BLAST to use.
    timeout : int or None, default None
        Maximum number of seconds to wait for BLAST to complete.
        If None, BLAST is allowed to run indefinitely.

    Returns
    -------
    subprocess.CompletedProcess
        Result object containing return code, stdout, and stderr.

    Raises
    ------
    FileNotFoundError
        If query file doesn't exist, or subject is neither an existing file
        nor a valid BLAST database prefix.
    BlastError
        If blastn is not available, execution fails, or output file not created.

    Notes
    -----
    Logs the BLAST command and thread usage at INFO level. Emits a progress
    message every 60 seconds while the search is running so that long jobs are
    distinguishable from frozen processes. Creates the output directory if needed.
    """
    # Validate inputs
    query_path = Path(query)
    subject_path = Path(subject)
    output_path = Path(output)

    if not query_path.exists():
        raise FileNotFoundError(f'Query file not found: {query_path}')

    # Determine whether subject is a FASTA file or a BLAST database prefix
    subject_is_db = blast_db_exists(subject_path)
    if not subject_path.exists() and not subject_is_db:
        raise FileNotFoundError(
            f'Subject file or BLAST database not found: {subject_path}'
        )

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Check if blastn is available
    if not check_blast_available():
        raise BlastError('blastn not found in PATH. Please install BLAST+ tools.')

    # Build command using list (safer than shell=True)
    cmd = [
        'blastn',
        '-word_size',
        str(word_size),
        '-outfmt',
        outfmt,
        '-query',
        str(query_path),
    ]

    # Use -db for BLAST database prefix, -subject for FASTA file
    if subject_is_db:
        cmd.extend(['-db', str(subject_path)])
    else:
        cmd.extend(['-subject', str(subject_path)])

    cmd.extend(
        [
            '-out',
            str(output_path),
            '-perc_identity',
            str(perc_identity),
            '-num_threads',
            str(num_threads),
        ]
    )

    # Add any additional arguments
    if additional_args:
        cmd.extend(additional_args)

    # Always log the command and thread usage
    available_cpus = os.cpu_count() or 1
    logger.info('Running blastn with the following parameters:')
    logger.info(f'  Command: {" ".join(cmd)}')
    logger.info(f'  Query: {query_path}')
    if subject_is_db:
        logger.info(f'  BLAST database: {subject_path}')
    else:
        logger.info(f'  Subject: {subject_path}')
    logger.info(f'  Output: {output_path}')
    logger.info(f'  Word size: {word_size}')
    logger.info(f'  Percent identity: {perc_identity}%')
    logger.info(f'  Threads: {num_threads} requested / {available_cpus} available')
    if additional_args:
        logger.info(f'  Additional args: {" ".join(additional_args)}')
    if timeout is None:
        logger.info('  Timeout: none (will run until complete)')
    else:
        logger.info(f'  Timeout: {timeout}s')

    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        # Poll every 5 seconds; emit a progress log message every 60 seconds
        # so long-running searches don't look frozen.
        _SLEEP_INTERVAL = 5  # seconds between polls
        _LOG_INTERVAL = 60  # seconds between progress messages
        start_time = time.time()
        last_log_time = start_time

        while proc.poll() is None:
            time.sleep(_SLEEP_INTERVAL)
            now = time.time()
            if now - last_log_time >= _LOG_INTERVAL:
                elapsed = int(now - start_time)
                logger.info(f'BLAST search is still running... ({elapsed}s elapsed)')
                last_log_time = now
            if timeout is not None and (now - start_time) >= timeout:
                proc.kill()
                proc.wait()
                unit = 'second' if timeout == 1 else 'seconds'
                raise BlastError(f'BLAST command timed out after {timeout} {unit}')

        stdout, stderr = proc.communicate()
        result = subprocess.CompletedProcess(
            args=cmd,
            returncode=proc.returncode,
            stdout=stdout,
            stderr=stderr,
        )

        # Check for errors
        if result.returncode != 0:
            error_msg = f'blastn failed with exit code {result.returncode}'
            if result.stderr:
                error_msg += f'\nSTDERR: {result.stderr}'
            if result.stdout:
                error_msg += f'\nSTDOUT: {result.stdout}'
            raise BlastError(error_msg)

        if verbose and result.stdout:
            logger.info(f'BLAST output: {result.stdout}')

        # Verify output file was created
        if not output_path.exists():
            raise BlastError(
                f'BLAST completed but output file not created: {output_path}'
            )

        return result

    except BlastError:
        raise
    except Exception as e:
        raise BlastError(f'Error running BLAST: {str(e)}') from e


def run_self_blast(
    sequence_file: Union[str, Path],
    output_file: Union[str, Path],
    perc_identity: float = 60.0,
    verbose: bool = False,
    num_threads: int = 1,  # Add threading support
    timeout: Optional[int] = None,
) -> subprocess.CompletedProcess:
    """
    Perform self-alignment by running blastn with sequence as both query and subject.

    Parameters
    ----------
    sequence_file : str or Path
        Path to sequence file for self-alignment (FASTA format).
    output_file : str or Path
        Path to output file for BLAST results.
    perc_identity : float, default 60.0
        Minimum percent identity threshold for reporting alignments.
    verbose : bool, default False
        If True, logs command and output to console.
    num_threads : int, default 1
        Number of CPU threads for BLAST to use.
    timeout : int or None, default None
        Maximum number of seconds to wait for BLAST to complete.
        If None, BLAST is allowed to run indefinitely.

    Returns
    -------
    subprocess.CompletedProcess
        Result object from blastn execution.

    Notes
    -----
    Convenience wrapper for run_blastn() with query and subject set to same file.
    Commonly used for TIR identification through self-complementarity detection.
    """
    return run_blastn(
        query=sequence_file,
        subject=sequence_file,
        output=output_file,
        perc_identity=perc_identity,
        verbose=verbose,
        num_threads=num_threads,
        timeout=timeout,
    )
