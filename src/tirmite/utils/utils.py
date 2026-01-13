"""
Utility functions for TIRmite workflow management.

Provides tools for:
- Temporary directory management
- Input/output file validation
- Genome indexing and description extraction
- Directory setup and cleanup
- Legacy compatibility functions
"""

from collections import Counter
from contextlib import contextmanager
import logging
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
from typing import Any, Generator, Optional, Tuple, Union

from Bio import SeqIO  # type: ignore[import-not-found]
from pyfaidx import Fasta  # type: ignore[import-not-found]


@contextmanager
def temporary_directory(
    prefix: str = 'tirmite_',
    suffix: Optional[str] = None,
    base_dir: Optional[Union[str, Path]] = None,
    cleanup: bool = True,
) -> 'Generator[Path, None, None]':
    """
    Context manager for creating and managing temporary directories.

    Parameters
    ----------
    prefix : str, default 'tirmite_'
        Prefix for temporary directory name.
    suffix : str, optional
        Suffix for temporary directory name.
    base_dir : str or Path, optional
        Base directory for temp dir creation. Uses system temp if None.
    cleanup : bool, default True
        Whether to automatically clean up the directory on exit.

    Yields
    ------
    Path
        Path to the created temporary directory.

    Examples
    --------
    >>> with temporary_directory(prefix='tirmite_', cleanup=True) as temp_dir:
    ...     # Use temp_dir for operations
    ...     pass
    # Directory is automatically cleaned up
    """
    temp_dir = None
    temp_path = None

    try:
        # Create temporary directory using secure methods
        temp_dir = tempfile.mkdtemp(
            prefix=prefix, suffix=suffix, dir=str(base_dir) if base_dir else None
        )
        temp_path = Path(temp_dir)
        yield temp_path

    except Exception as e:
        # Re-raise the original exception, cleanup happens in finally block
        raise e

    finally:
        # Clean up if requested and directory exists
        if cleanup and temp_dir is not None and Path(temp_dir).exists():
            try:
                shutil.rmtree(temp_dir)
            except OSError as e:
                # Log warning but don't fail the entire operation
                import logging

                logging.warning(
                    f'Failed to clean up temporary directory {temp_dir}: {e}'
                )


def create_output_directory(output_path: Optional[Union[str, Path]] = None) -> Path:
    """
    Create output directory with proper error handling.

    Parameters
    ----------
    output_path : str or Path, optional
        Desired output directory path. Uses current directory if None.

    Returns
    -------
    Path
        Absolute path to the created/validated output directory.

    Raises
    ------
    OSError
        If directory creation fails or path is not writable.
    """
    if output_path:
        out_path = Path(output_path).resolve()
        try:
            out_path.mkdir(parents=True, exist_ok=True)

            # Test if directory is writable
            test_file = out_path / '.tirmite_write_test'
            try:
                test_file.touch()
                test_file.unlink()
            except OSError:
                raise OSError(f'Output directory is not writable: {out_path}') from None

        except OSError as e:
            raise OSError(f'Failed to create output directory {out_path}: {e}') from e
    else:
        out_path = Path.cwd()

    return out_path


def validate_input_files(args: Any) -> None:
    """
    Validate that all required input files exist and are readable.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments containing file paths.

    Raises
    ------
    FileNotFoundError
        If any required input file doesn't exist or isn't readable.
    """
    required_files = []
    optional_files = []

    # Always required
    required_files.append(('genome', args.genome))

    # Conditionally required
    if hasattr(args, 'hmmFile') and args.hmmFile:
        required_files.append(('HMM file', args.hmmFile))
    if hasattr(args, 'alnFile') and args.alnFile:
        required_files.append(('alignment file', args.alnFile))
    if hasattr(args, 'pairbed') and args.pairbed:
        required_files.append(('BED file', args.pairbed))
    if hasattr(args, 'matrix') and args.matrix:
        optional_files.append(('matrix file', args.matrix))

    # Check required files
    for file_type, file_path in required_files:
        if not file_path:
            continue

        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f'{file_type} not found: {file_path}')
        if not path.is_file():
            raise FileNotFoundError(f'{file_type} is not a regular file: {file_path}')
        if not os.access(path, os.R_OK):
            raise PermissionError(f'{file_type} is not readable: {file_path}')

    # Check optional files (warn instead of error)
    for file_type, file_path in optional_files:
        if not file_path:
            continue

        path = Path(file_path)
        if not path.exists():
            import logging

            logging.warning(f'Optional {file_type} not found: {file_path}')


def setup_directories(args: Any) -> Tuple[Path, Path]:
    """
    Set up output and temporary directories with proper error handling.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments.

    Returns
    -------
    tuple[Path, Path]
        - Absolute path to output directory
        - Absolute path to temporary directory

    Raises
    ------
    OSError
        If directory creation fails.
    FileNotFoundError
        If input files don't exist.
    """
    # Validate input files first
    validate_input_files(args)

    # Create output directory
    output_dir = create_output_directory(
        args.outdir if hasattr(args, 'outdir') else None
    )

    # Create temporary directory with proper parent directory handling
    base_temp_dir = None
    if hasattr(args, 'tempdir') and args.tempdir:
        base_temp_dir = Path(args.tempdir)

        # Create the parent directory structure if it doesn't exist
        try:
            base_temp_dir.mkdir(parents=True, exist_ok=True)

            # Test if directory is writable
            test_file = base_temp_dir / '.tirmite_write_test'
            try:
                test_file.touch()
                test_file.unlink()
            except OSError:
                raise OSError(
                    f'Temporary base directory is not writable: {base_temp_dir}'
                ) from None

        except OSError as e:
            raise OSError(
                f'Failed to create temporary base directory {base_temp_dir}: {e}'
            ) from e

    # Create the actual temporary directory
    try:
        temp_dir = Path(
            tempfile.mkdtemp(
                prefix='tirmite_', dir=str(base_temp_dir) if base_temp_dir else None
            )
        )
    except OSError as e:
        raise OSError(f'Failed to create temporary directory: {e}') from e

    return output_dir, temp_dir


def cleanup_temp_directory(temp_dir: Union[str, Path], keep_temp: bool = False) -> None:
    """
    Safely clean up temporary directory.

    Parameters
    ----------
    temp_dir : str or Path
        Path to temporary directory to clean up.
    keep_temp : bool, default False
        If True, skip cleanup and log directory location.
    """
    temp_path = Path(temp_dir)

    if keep_temp:
        import logging

        logging.info(f'Temporary directory preserved: {temp_path}')
        return

    if temp_path.exists() and temp_path.is_dir():
        try:
            shutil.rmtree(temp_path)
            import logging

            logging.debug(f'Cleaned up temporary directory: {temp_path}')
        except OSError as e:
            import logging

            logging.warning(f'Failed to clean up temporary directory {temp_path}: {e}')


# Legacy function for backwards compatibility
def dochecks(args: Any) -> Tuple[Path, Path]:
    """
    DEPRECATED: Use setup_directories() instead.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments.

    Returns
    -------
    tuple
        (outdir, tempdir) paths as strings.
    """
    import warnings

    warnings.warn(
        'dochecks() is deprecated. Use setup_directories() instead.',
        DeprecationWarning,
        stacklevel=2,
    )

    return setup_directories(args)


def isfile(path: str) -> None:
    """
    DEPRECATED: Use validate_input_files() instead.

    Parameters
    ----------
    path : str
        File path to check.

    Returns
    -------
    None
        No return value. Raises SystemExit if file not found.
    """
    import warnings

    warnings.warn(
        'isfile() is deprecated. Use validate_input_files() instead.',
        DeprecationWarning,
        stacklevel=2,
    )

    if not os.path.isfile(path):
        print('Input file not found: %s' % path)
        sys.exit(1)


def cleanID(s: str) -> str:
    """
    Remove non-alphanumeric characters and normalize whitespace in string.

    Parameters
    ----------
    s : str
        Input string to clean.

    Returns
    -------
    str
        Cleaned string with only alphanumeric characters and underscores.
        Whitespace sequences replaced with single underscores.

    Examples
    --------
    >>> cleanID("My-Model Name_v1!")
    'My_Model_Name_v1'
    """
    s = re.sub(r'[^\w\s]', '', s)
    s = re.sub(r'\s+', '_', s)
    return s


def manageTemp(
    record: Any = None, tempPath: Optional[str] = None, scrub: bool = False
) -> None:
    """
    Write sequence record to temporary file or delete temporary file.

    Parameters
    ----------
    record : Bio.SeqRecord.SeqRecord, optional
        Sequence record to write to file. Required if scrub is False.
    tempPath : str, optional
        Path to temporary file for writing or deletion.
    scrub : bool, default False
        If True, deletes file at tempPath. If False, writes record to tempPath.

    Returns
    -------
    None
        No return value.

    Notes
    -----
    This is a legacy function. Consider using temporary_directory()
    context manager for more robust temporary file handling.
    """
    if scrub and tempPath:
        try:
            os.remove(tempPath)
        except OSError:
            pass
    elif tempPath:
        with open(tempPath, 'w') as f:
            SeqIO.write(record, f, 'fasta')


def checkUniqueID(records: Any) -> None:
    """
    Verify that all sequence IDs in record list are unique.

    Parameters
    ----------
    records : list of Bio.SeqRecord.SeqRecord
        List of sequence records to check.

    Returns
    -------
    None
        No return value.

    Raises
    ------
    SystemExit
        If duplicate IDs are found, prints duplicate IDs and exits with code 1.
    """
    seqIDs = [records[x].id for x in range(len(records))]
    IDcounts = Counter(seqIDs)
    duplicates = [k for k, v in IDcounts.items() if v > 1]
    if duplicates:
        print('Input sequence IDs not unique. Quiting.')
        print(duplicates)
        sys.exit(1)
    else:
        pass


def indexGenome(genomePath: Union[str, Path]) -> Tuple[Fasta, dict]:
    """
    Index genome FASTA file and extract sequence descriptions.

    Parameters
    ----------
    genomePath : str or Path
        Path to genome FASTA file to index.

    Returns
    -------
    genome : pyfaidx.Fasta
        Indexed genome object allowing efficient random access to sequences.
    descriptions : dict
        Dictionary mapping sequence IDs to their description strings.

    Raises
    ------
    FileNotFoundError
        If genome file doesn't exist at specified path.

    Notes
    -----
    Creates a .fai index file alongside the genome FASTA for rapid sequence access.
    Descriptions are parsed from FASTA headers (text after first whitespace).
    """
    genome_path = Path(genomePath)

    if not genome_path.exists():
        raise FileNotFoundError(f'Genome file not found: {genome_path}')

    # Index with pyfaidx
    genome = Fasta(str(genome_path))

    # Extract descriptions
    descriptions = extract_genome_descriptions(genome_path)

    logging.info(f'Indexed genome with {len(genome.keys())} sequences')
    logging.debug(f'Extracted descriptions for {len(descriptions)} sequences')

    return genome, descriptions


def extract_genome_descriptions(genome_path: Union[str, Path]) -> dict:
    """
    Parse sequence descriptions from genome FASTA file headers.

    Parameters
    ----------
    genome_path : str or Path
        Path to genome FASTA file.

    Returns
    -------
    dict
        Dictionary mapping sequence IDs to description strings.
        Returns empty dict if parsing fails.

    Notes
    -----
    Extracts text following the sequence ID in FASTA headers.
    Header format: >sequence_id description text
    If no description present, maps to empty string.
    """
    descriptions = {}

    try:
        with open(genome_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Parse header: >ID description...
                    header = line[1:].strip()
                    parts = header.split(None, 1)  # Split on first whitespace
                    seq_id = parts[0]
                    description = parts[1] if len(parts) > 1 else ''
                    descriptions[seq_id] = description

    except Exception as e:
        logging.warning(f'Could not extract genome descriptions: {e}')

    return descriptions
