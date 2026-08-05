"""
Utility functions for TIRmite workflow management.

Provides tools for:
- Temporary directory management
- Input/output file validation
- Genome indexing and description extraction
- Directory setup and cleanup
- Legacy compatibility functions
"""

from contextlib import contextmanager
import gzip
import logging
import os
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any, Generator, Optional, Tuple, Union

from pyfaidx import Fasta

# Library modules acquire a named logger and attach a NullHandler, so that
# importing TIRmite as a library emits nothing until the host application
# configures logging. Handler setup belongs to the CLI, not here.
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


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

                logger.warning(
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
            logger.warning(f'Optional {file_type} not found: {file_path}')


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
        logger.info(f'Temporary directory preserved: {temp_path}')
        return

    if temp_path.exists() and temp_path.is_dir():
        try:
            shutil.rmtree(temp_path)

            logger.debug(f'Cleaned up temporary directory: {temp_path}')
        except OSError as e:
            logger.warning(f'Failed to clean up temporary directory {temp_path}: {e}')


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


def extract_model_name_from_path(model_path: Optional[str]) -> Optional[str]:
    """
    Extract HMM model name from file path by parsing the HMM file header.

    Parameters
    ----------
    model_path : str or Path, optional
        Path to an HMM file. May be None.

    Returns
    -------
    str or None
        Model name taken from the ``NAME`` field in the HMM file. Falls back
        to the filename stem if the field is absent or the file cannot be
        read. Returns None if ``model_path`` is None.

    Notes
    -----
    An HMMER3 model file declares its name on a line beginning with ``NAME``
    followed by two spaces, e.g. ``NAME  MY_TIR``. The name recorded inside
    the file is authoritative: it, not the filename, is what nhmmer reports in
    its output, so hit tables must be keyed on it for the pairing and anchor
    filters to match models to hits.

    Examples
    --------
    >>> extract_model_name_from_path(None) is None
    True
    """
    if not model_path:
        return None

    try:
        with open(model_path, 'r') as f:
            for line in f:
                if line.startswith('NAME  '):
                    return line.split()[1].strip()
    except (FileNotFoundError, IOError):
        # An unreadable model file is not fatal here; the caller only needs a
        # label, and the filename stem is a reasonable approximation.
        return Path(model_path).stem

    # File was readable but carried no NAME field.
    return Path(model_path).stem


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

    logger.info(f'Indexed genome with {len(genome.keys())} sequences')
    logger.debug(f'Extracted descriptions for {len(descriptions)} sequences')

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
        logger.warning(f'Could not extract genome descriptions: {e}')

    return descriptions


def is_gzipped_file(file_path: Union[str, Path]) -> bool:
    """
    Check if a file is gzip-compressed.

    Parameters
    ----------
    file_path : str or Path
        Path to file to check.

    Returns
    -------
    bool
        True if file is gzip-compressed, False otherwise.

    Notes
    -----
    Checks both file extension (.gz) and file magic bytes for gzip format.
    """
    file_path = Path(file_path)

    # First check extension
    if file_path.suffix.lower() == '.gz':
        return True

    # Also check magic bytes if file exists
    if file_path.exists():
        try:
            with open(file_path, 'rb') as f:
                # Gzip files start with 1f 8b magic bytes
                return f.read(2) == b'\x1f\x8b'
        except Exception:
            pass

    return False


def decompress_genome(
    genome_path: Union[str, Path],
    output_dir: Union[str, Path],
) -> Path:
    """
    Decompress a gzip-compressed genome file to output directory.

    Parameters
    ----------
    genome_path : str or Path
        Path to gzip-compressed genome file.
    output_dir : str or Path
        Directory where decompressed file will be written.

    Returns
    -------
    Path
        Path to decompressed genome file.

    Raises
    ------
    FileNotFoundError
        If input genome file doesn't exist.
    OSError
        If decompression fails.

    Notes
    -----
    Creates a decompressed copy with .gz extension removed.
    Original file is not modified.
    """
    genome_path = Path(genome_path)
    output_dir = Path(output_dir)

    if not genome_path.exists():
        raise FileNotFoundError(f'Genome file not found: {genome_path}')

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine output filename (remove .gz extension)
    if genome_path.suffix.lower() == '.gz':
        output_name = genome_path.stem
    else:
        output_name = genome_path.name + '.decompressed'

    output_path = output_dir / output_name

    logger.info(f'Decompressing {genome_path.name} to {output_path}')

    try:
        with gzip.open(genome_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        logger.debug(f'Decompression complete: {output_path}')
        return output_path

    except Exception as e:
        raise OSError(f'Failed to decompress {genome_path}: {e}') from e


def prepare_genome_file(
    genome_path: Union[str, Path],
    temp_dir: Union[str, Path],
) -> Path:
    """
    Prepare genome file for use, decompressing if necessary.

    Parameters
    ----------
    genome_path : str or Path
        Path to genome file (may be gzip-compressed).
    temp_dir : str or Path
        Temporary directory for decompressed files.

    Returns
    -------
    Path
        Path to prepared (decompressed if necessary) genome file.

    Raises
    ------
    FileNotFoundError
        If genome file doesn't exist.
    OSError
        If decompression fails.

    Notes
    -----
    If genome is gzip-compressed, decompresses to temp_dir.
    Otherwise returns original path unchanged.
    """
    genome_path = Path(genome_path)

    if not genome_path.exists():
        raise FileNotFoundError(f'Genome file not found: {genome_path}')

    # Check if file is gzipped
    if is_gzipped_file(genome_path):
        logger.info(f'Detected gzip-compressed genome: {genome_path.name}')
        return decompress_genome(genome_path, temp_dir)
    else:
        return genome_path
