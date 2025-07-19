from collections import Counter
from contextlib import contextmanager
from datetime import datetime
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
from typing import Optional, Tuple, Union

from Bio import SeqIO
from pyfaidx import Fasta


@contextmanager
def temporary_directory(
    prefix: str = 'tirmite_',
    suffix: Optional[str] = None,
    base_dir: Optional[Union[str, Path]] = None,
    cleanup: bool = True,
):
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


def validate_input_files(args) -> None:
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


def setup_directories(args) -> Tuple[Path, Path]:
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

    # Create temporary directory
    # Use a more secure location by default, but allow override
    base_temp_dir = None
    if hasattr(args, 'tempdir') and args.tempdir:
        base_temp_dir = args.tempdir

    temp_dir = Path(
        tempfile.mkdtemp(
            prefix='tirmite_', dir=str(base_temp_dir) if base_temp_dir else None
        )
    )

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
def dochecks(args):
    """
    DEPRECATED: Use setup_directories() instead.
    Legacy function maintained for backwards compatibility.
    """
    import warnings

    warnings.warn(
        'dochecks() is deprecated. Use setup_directories() instead.',
        DeprecationWarning,
        stacklevel=2,
    )

    return setup_directories(args)


def getTimestring():
    """
    Return int only string of current datetime with milliseconds.

    Note: This function is deprecated for security reasons.
    Use tempfile module for temporary file/directory naming instead.
    """
    import warnings

    warnings.warn(
        'getTimestring() for temp naming is deprecated. Use tempfile module instead.',
        DeprecationWarning,
        stacklevel=2,
    )

    (dt, micro) = datetime.utcnow().strftime('%Y%m%d%H%M%S.%f').split('.')
    dt = '%s%03d' % (dt, int(micro) / 1000)
    return dt


def isfile(path):
    """
    DEPRECATED: Use validate_input_files() instead.
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


def cleanID(s):
    """
    Remove non alphanumeric characters from string.
    Replace whitespace with underscores.
    """
    s = re.sub(r'[^\w\s]', '', s)
    s = re.sub(r'\s+', '_', s)
    return s


def manageTemp(record=None, tempPath=None, scrub=False):
    """
    Create single sequence fasta files or scrub temp files.

    Note: Consider using temporary_directory() context manager instead.
    """
    if scrub and tempPath:
        try:
            os.remove(tempPath)
        except OSError:
            pass
    else:
        with open(tempPath, 'w') as f:
            SeqIO.write(record, f, 'fasta')


def checkUniqueID(records):
    """
    Check that IDs for input elements are unique.
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


def indexGenome(genome_path):
    """
    Create or load pyfaidx index for genome file.
    Returns Fasta object for efficient sequence access.

    Args:
        genome_path: Path to genome FASTA file

    Returns:
        pyfaidx.Fasta: Indexed genome object
    """
    # pyfaidx automatically creates .fai index if it doesn't exist
    # and uses existing index if present
    try:
        genome_index = Fasta(genome_path)
        return genome_index
    except Exception as e:
        print(f'Error indexing genome file {genome_path}: {e}')
        sys.exit(1)
