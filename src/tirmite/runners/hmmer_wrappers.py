import logging
from pathlib import Path
import re
import shutil
from typing import List, Optional, Tuple, Union

from tirmite.runners.wrapping import run_commands_sequential


def cleanID(sequence_id: str) -> str:
    """
    Remove non-alphanumeric characters from string and replace whitespace with underscores.

    Parameters
    ----------
    sequence_id : str
        Input string to clean.

    Returns
    -------
    str
        Cleaned string with only alphanumeric characters and underscores.

    Examples
    --------
    >>> cleanID("My Model-Name 1")
    'My_Model_Name_1'
    """
    # Remove non-alphanumeric characters except whitespace
    cleaned = re.sub(r'[^\w\s]', '', sequence_id)
    # Replace whitespace with underscores
    cleaned = re.sub(r'\s+', '_', cleaned)
    return cleaned


def build_hmmbuild_command(
    model_name: str,
    input_alignment: Union[str, Path],
    output_dir: Union[str, Path],
    executable_path: str = 'hmmbuild',
    cores: Optional[int] = None,
) -> Tuple[List[str], Path]:
    """
    Construct hmmbuild command for creating HMM from multiple sequence alignment.

    Parameters
    ----------
    model_name : str
        Name for the HMM model (will be cleaned for filesystem safety).
    input_alignment : str or Path
        Path to input multiple sequence alignment file.
    output_dir : str or Path
        Directory where HMM file will be created.
    executable_path : str, default 'hmmbuild'
        Path or name of hmmbuild executable.
    cores : int, optional
        Number of CPU cores to use for building.

    Returns
    -------
    tuple[list[str], Path]
        - Command as list of arguments for subprocess
        - Path to output HMM file

    Raises
    ------
    FileNotFoundError
        If input alignment file doesn't exist.
    """
    # Clean model name for filesystem safety
    clean_model_name = cleanID(model_name)

    # Convert paths to Path objects
    input_path = Path(input_alignment)
    output_path = Path(output_dir)

    # Validate input file exists
    if not input_path.exists():
        raise FileNotFoundError(f'Input alignment file not found: {input_path}')

    # Create HMM database directory
    hmm_db_dir = output_path / 'hmmDB'
    hmm_db_dir.mkdir(parents=True, exist_ok=True)

    # Construct output file path
    output_hmm = hmm_db_dir / f'{clean_model_name}.hmm'

    # Build command as list (safer than shell strings)
    command = [
        executable_path,
        '--dna',  # DNA alphabet
        '-n',
        clean_model_name,  # Model name
    ]

    # Add optional parameters
    if cores:
        command.extend(['--cpu', str(cores)])

    # Add output file and input alignment
    command.extend([str(output_hmm), str(input_path)])

    return command, output_hmm


def build_hmmpress_command(
    hmm_file: Union[str, Path],
    executable_path: str = 'hmmpress',
) -> List[str]:
    """
    Construct hmmpress command for indexing HMM database.

    Parameters
    ----------
    hmm_file : str or Path
        Path to HMM file to index.
    executable_path : str, default 'hmmpress'
        Path or name of hmmpress executable.

    Returns
    -------
    list[str]
        Command as list of arguments for subprocess.

    Raises
    ------
    FileNotFoundError
        If HMM file doesn't exist.
    """
    hmm_path = Path(hmm_file)

    # Validate HMM file exists
    if not hmm_path.exists():
        raise FileNotFoundError(f'HMM file not found: {hmm_path}')

    # Build command - force overwrite existing index files
    command = [
        executable_path,
        '-f',  # Force overwrite
        str(hmm_path),
    ]

    return command


def build_nhmmer_command(
    model_path: Union[str, Path],
    genome_path: Union[str, Path],
    output_dir: Union[str, Path],
    executable_path: str = 'nhmmer',
    evalue: Optional[float] = None,
    cores: Optional[int] = None,
    nobias: bool = False,
    matrix_file: Optional[Union[str, Path]] = None,
) -> Tuple[List[str], Path]:
    """
    Construct nhmmer command for searching HMM against genome sequence.

    Parameters
    ----------
    model_path : str or Path
        Path to HMM model file.
    genome_path : str or Path
        Path to target genome sequence file.
    output_dir : str or Path
        Directory for output files.
    executable_path : str, default 'nhmmer'
        Path or name of nhmmer executable.
    evalue : float, optional
        E-value threshold for reporting hits.
    cores : int, optional
        Number of CPU cores to use.
    nobias : bool, default False
        Turn off composition bias correction.
    matrix_file : str or Path, optional
        Path to custom scoring matrix file.

    Returns
    -------
    tuple[list[str], Path]
        - Command as list of arguments for subprocess
        - Path to output directory

    Raises
    ------
    FileNotFoundError
        If model or genome files don't exist.
    """
    # Convert to Path objects
    model_path = Path(model_path)
    genome_path = Path(genome_path)
    output_path = Path(output_dir)

    # Validate input files exist
    if not model_path.exists():
        raise FileNotFoundError(f'HMM model file not found: {model_path}')
    if not genome_path.exists():
        raise FileNotFoundError(f'Genome file not found: {genome_path}')

    # Get model basename for output file naming
    model_base = model_path.stem

    # Create nhmmer results directory
    results_dir = output_path / 'nhmmer_results'
    results_dir.mkdir(parents=True, exist_ok=True)

    # Construct output file path
    output_file = results_dir / f'{model_base}.out'

    # Build base command
    command = [
        executable_path,
        '--tblout',
        str(output_file),  # Tabular output format
        '--noali',  # Don't show alignments
        '--notextw',  # Unlimited line width
        '--dna',  # DNA alphabet
        '--max',  # Turn off all acceleration heuristics
    ]

    # Add optional parameters
    if cores:
        command.extend(['--cpu', str(cores)])
    if evalue is not None:
        command.extend(['-E', str(evalue)])
    if nobias:
        command.append('--nobias')
    if matrix_file:
        matrix_path = Path(matrix_file)
        if not matrix_path.exists():
            raise FileNotFoundError(f'Matrix file not found: {matrix_path}')
        command.extend(['--mxfile', str(matrix_path)])

    # Add model and target sequence files
    command.extend([str(model_path), str(genome_path)])

    return command, results_dir


def process_hmmer_workflow(
    hmm_dir: Optional[Union[str, Path]] = None,
    hmm_file: Optional[Union[str, Path]] = None,
    alignment_dir: Optional[Union[str, Path]] = None,
    left_model: Optional[Union[str, Path]] = None,  # Add this parameter
    right_model: Optional[Union[str, Path]] = None,  # Add this parameter
    temp_dir: Optional[Union[str, Path]] = None,
    genome_path: Union[str, Path] = None,
    executable_paths: Optional[dict] = None,
    search_params: Optional[dict] = None,
    verbose: bool = False,
) -> Tuple[Path, Path]:
    """
    Execute complete HMMER workflow: build models, press databases, and search genome.

    Parameters
    ----------
    hmm_dir : str or Path, optional
        Directory containing existing HMM files to copy.
    hmm_file : str or Path, optional
        Single HMM file to copy (mutually exclusive with hmm_dir).
    alignment_dir : str or Path, optional
        Directory containing alignment files to convert to HMMs.
    left_model : str or Path, optional
        Path to left HMM model file.
    right_model : str or Path, optional
        Path to right HMM model file.
    temp_dir : str or Path, optional
        Working directory for temporary files. Defaults to current directory.
    genome_path : str or Path
        Path to target genome sequence file.
    executable_paths : dict, optional
        Dictionary mapping tool names to executable paths.
        Keys: 'hmmbuild', 'hmmpress', 'nhmmer'
    search_params : dict, optional
        Dictionary of search parameters for nhmmer.
        Keys: 'evalue', 'cores', 'nobias', 'matrix'
    verbose : bool, default False
        Enable verbose logging of commands.

    Returns
    -------
    tuple[Path, Path]
        - Path to nhmmer results directory
        - Path to HMM database directory

    Raises
    ------
    ValueError
        If no HMM sources are provided or invalid parameter combinations.
    FileNotFoundError
        If required input files don't exist.

    Examples
    --------
    >>> results_dir, hmm_db = process_hmmer_workflow(
    ...     alignment_dir='alignments/',
    ...     genome_path='genome.fasta',
    ...     search_params={'evalue': 1e-3, 'cores': 4}
    ... )
    """
    # Validate inputs - now includes left_model and right_model
    if not any([hmm_dir, hmm_file, alignment_dir, left_model, right_model]):
        raise ValueError(
            'Must provide at least one of: hmm_dir, hmm_file, alignment_dir, left_model, or right_model'
        )

    if hmm_dir and hmm_file:
        raise ValueError('hmm_dir and hmm_file are mutually exclusive')

    # Set up paths
    if temp_dir:
        temp_path = Path(temp_dir)
        temp_path.mkdir(parents=True, exist_ok=True)
    else:
        temp_path = Path.cwd()

    # Set default executable paths
    default_executables = {
        'hmmbuild': 'hmmbuild',
        'hmmpress': 'hmmpress',
        'nhmmer': 'nhmmer',
    }
    if executable_paths:
        default_executables.update(executable_paths)

    # Set default search parameters
    default_search_params = {
        'evalue': None,
        'cores': None,
        'nobias': False,
        'matrix': None,
    }
    if search_params:
        default_search_params.update(search_params)

    # Create HMM database directory
    hmm_db_path = temp_path / 'hmmDB'
    hmm_db_path.mkdir(parents=True, exist_ok=True)

    # Copy existing HMM files
    if hmm_dir:
        hmm_dir_path = Path(hmm_dir)
        if not hmm_dir_path.exists():
            raise FileNotFoundError(f'HMM directory not found: {hmm_dir_path}')

        for hmm_file_path in hmm_dir_path.glob('*.hmm'):
            shutil.copy2(hmm_file_path, hmm_db_path)
            if verbose:
                logging.info(f'Copied HMM file: {hmm_file_path.name}')

    if hmm_file:
        hmm_file_path = Path(hmm_file)
        if not hmm_file_path.exists():
            raise FileNotFoundError(f'HMM file not found: {hmm_file_path}')

        shutil.copy2(hmm_file_path, hmm_db_path)
        if verbose:
            logging.info(f'Copied HMM file: {hmm_file_path.name}')

    # Copy left and right model files
    if left_model:
        left_model_path = Path(left_model)
        if not left_model_path.exists():
            raise FileNotFoundError(f'Left model file not found: {left_model_path}')

        shutil.copy2(left_model_path, hmm_db_path)
        if verbose:
            logging.info(f'Copied left model file: {left_model_path.name}')

    if right_model:
        right_model_path = Path(right_model)
        if not right_model_path.exists():
            raise FileNotFoundError(f'Right model file not found: {right_model_path}')

        shutil.copy2(right_model_path, hmm_db_path)
        if verbose:
            logging.info(f'Copied right model file: {right_model_path.name}')

    # Build HMMs from alignments
    if alignment_dir:
        alignment_path = Path(alignment_dir)
        if not alignment_path.exists():
            raise FileNotFoundError(f'Alignment directory not found: {alignment_path}')

        build_commands = []

        # Find all alignment files (common formats)
        alignment_patterns = ['*.fasta', '*.fas', '*.fa', '*.aln', '*.sto']
        alignment_files = []
        for pattern in alignment_patterns:
            alignment_files.extend(alignment_path.glob(pattern))

        if not alignment_files:
            logging.warning(f'No alignment files found in {alignment_path}')

        for alignment_file in alignment_files:
            model_name = alignment_file.stem

            try:
                command, output_path = build_hmmbuild_command(
                    model_name=model_name,
                    input_alignment=alignment_file,
                    output_dir=temp_path,
                    executable_path=default_executables['hmmbuild'],
                    cores=default_search_params['cores'],
                )
                build_commands.append(command)

                if verbose:
                    logging.info(f'Prepared hmmbuild command for {model_name}')

            except Exception as e:
                logging.error(f'Failed to prepare hmmbuild for {alignment_file}: {e}')
                continue

        # Execute all build commands
        if build_commands:
            logging.info(f'Building {len(build_commands)} HMM models...')
            try:
                run_commands_sequential(
                    cmds=build_commands, verbose=verbose, stop_on_error=True
                )
                logging.info('HMM building completed successfully')
            except Exception as e:
                logging.error(f'HMM building failed: {e}')
                raise

    # Process all HMM files: press and search
    hmm_files = list(hmm_db_path.glob('*.hmm'))
    if not hmm_files:
        raise ValueError('No HMM files found for processing')

    search_commands = []
    results_dir = None

    for hmm_file_path in hmm_files:
        try:
            # Build hmmpress command
            press_command = build_hmmpress_command(
                hmm_file=hmm_file_path, executable_path=default_executables['hmmpress']
            )

            # Build nhmmer command
            nhmmer_command, current_results_dir = build_nhmmer_command(
                model_path=hmm_file_path,
                genome_path=genome_path,
                output_dir=temp_path,
                executable_path=default_executables['nhmmer'],
                evalue=default_search_params['evalue'],
                cores=default_search_params['cores'],
                nobias=default_search_params['nobias'],
                matrix_file=default_search_params['matrix'],
            )

            # Store commands for execution
            search_commands.extend([press_command, nhmmer_command])

            # All models use same results directory
            if results_dir is None:
                results_dir = current_results_dir

            if verbose:
                logging.info(f'Prepared search commands for {hmm_file_path.name}')

        except Exception as e:
            logging.error(f'Failed to prepare commands for {hmm_file_path}: {e}')
            continue

    # Execute all search commands
    if search_commands:
        logging.info(f'Executing {len(search_commands)} HMMER commands...')
        try:
            run_commands_sequential(
                cmds=search_commands, verbose=verbose, stop_on_error=True
            )
            logging.info('HMMER search completed successfully')
        except Exception as e:
            logging.error(f'HMMER search failed: {e}')
            raise

    return results_dir, hmm_db_path
