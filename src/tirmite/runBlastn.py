import logging
from pathlib import Path
import subprocess
from typing import List, Optional, Union


class BlastError(Exception):
    """Custom exception for BLAST-related errors."""

    pass


def check_blast_available() -> bool:
    """
    Check if blastn is available in PATH.

    Returns:
        bool: True if blastn is found, False otherwise
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
) -> subprocess.CompletedProcess:
    """
    Run blastn with specified parameters using subprocess best practices.

    Args:
        query: Path to query sequence file
        subject: Path to subject sequence file
        output: Path to output file
        word_size: Word size for blastn (default: 4)
        perc_identity: Minimum percent identity (default: 60.0)
        outfmt: Output format string (default: tabular with specific fields)
        additional_args: Additional command line arguments
        verbose: Print command and output if True
        num_threads: Number of CPU threads to use (default: 1)

    Returns:
        subprocess.CompletedProcess: Result of the blastn command

    Raises:
        BlastError: If blastn fails or is not available
        FileNotFoundError: If input files don't exist
    """
    # Validate inputs
    query_path = Path(query)
    subject_path = Path(subject)
    output_path = Path(output)

    if not query_path.exists():
        raise FileNotFoundError(f'Query file not found: {query_path}')
    if not subject_path.exists():
        raise FileNotFoundError(f'Subject file not found: {subject_path}')

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
        '-subject',
        str(subject_path),
        '-out',
        str(output_path),
        '-perc_identity',
        str(perc_identity),
        '-num_threads',
        str(num_threads),
    ]

    # Add any additional arguments
    if additional_args:
        cmd.extend(additional_args)

    if verbose:
        logging.info(
            f'Running BLAST command with {num_threads} threads: {" ".join(cmd)}'
        )

    try:
        # Use subprocess.run with proper error handling
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout
            check=False,  # Don't raise on non-zero exit
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
            logging.info(f'BLAST output: {result.stdout}')

        # Verify output file was created
        if not output_path.exists():
            raise BlastError(
                f'BLAST completed but output file not created: {output_path}'
            )

        return result

    except subprocess.TimeoutExpired as err:
        raise BlastError('BLAST command timed out after 5 minutes') from err
    except Exception as e:
        raise BlastError(f'Error running BLAST: {str(e)}') from e


def run_self_blast(
    sequence_file: Union[str, Path],
    output_file: Union[str, Path],
    perc_identity: float = 60.0,
    verbose: bool = False,
    num_threads: int = 1,  # Add threading support
) -> subprocess.CompletedProcess:
    """
    Run blastn with sequence file as both query and subject (self-alignment).

    This is a convenience wrapper for the common case of self-alignment
    used in TIR identification.

    Args:
        sequence_file: Path to sequence file for self-alignment
        output_file: Path to output file
        perc_identity: Minimum percent identity (default: 60.0)
        verbose: Print command and output if True
        num_threads: Number of CPU threads to use (default: 1)

    Returns:
        subprocess.CompletedProcess: Result of the blastn command
    """
    return run_blastn(
        query=sequence_file,
        subject=sequence_file,
        output=output_file,
        perc_identity=perc_identity,
        verbose=verbose,
        num_threads=num_threads,
    )


def run_blast_batch(
    sequence_files: List[Union[str, Path]],
    output_dir: Union[str, Path],
    perc_identity: float = 60.0,
    verbose: bool = False,
    max_workers: int = 1,
) -> List[subprocess.CompletedProcess]:
    """
    Run BLAST on multiple sequence files.

    Args:
        sequence_files: List of sequence files for self-alignment
        output_dir: Directory to write output files
        perc_identity: Minimum percent identity (default: 60.0)
        verbose: Print commands and output if True
        max_workers: Number of parallel processes (default: 1 for sequential)

    Returns:
        List of subprocess.CompletedProcess results
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    results = []

    if max_workers == 1:
        # Sequential execution
        for seq_file in sequence_files:
            seq_path = Path(seq_file)
            output_file = output_path / f'{seq_path.stem}_blast.tab'

            try:
                result = run_self_blast(
                    sequence_file=seq_path,
                    output_file=output_file,
                    perc_identity=perc_identity,
                    verbose=verbose,
                )
                results.append(result)

                if verbose:
                    logging.info(f'Completed BLAST for {seq_path.name}')

            except Exception as e:
                logging.error(f'BLAST failed for {seq_path.name}: {str(e)}')
                # Continue with remaining files
                continue
    else:
        # Parallel execution using ThreadPoolExecutor
        import concurrent.futures

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all jobs
            future_to_file = {}
            for seq_file in sequence_files:
                seq_path = Path(seq_file)
                output_file = output_path / f'{seq_path.stem}_blast.tab'

                future = executor.submit(
                    run_self_blast,
                    sequence_file=seq_path,
                    output_file=output_file,
                    perc_identity=perc_identity,
                    verbose=verbose,
                )
                future_to_file[future] = seq_path

            # Collect results
            for future in concurrent.futures.as_completed(future_to_file):
                seq_path = future_to_file[future]
                try:
                    result = future.result()
                    results.append(result)
                    if verbose:
                        logging.info(f'Completed BLAST for {seq_path.name}')
                except Exception as e:
                    logging.error(f'BLAST failed for {seq_path.name}: {str(e)}')
                    continue

    return results


# Deprecated functions kept for backwards compatibility
def makeBlast(seq=None, outfile=None, pid=60):
    """
    DEPRECATED: Use run_self_blast() instead.
    Legacy function maintained for backwards compatibility.
    """
    import warnings

    warnings.warn(
        'makeBlast() is deprecated. Use run_self_blast() instead.',
        DeprecationWarning,
        stacklevel=2,
    )

    try:
        run_self_blast(sequence_file=seq, output_file=outfile, perc_identity=pid)
    except Exception as e:
        raise Error(f'BLAST failed: {str(e)}') from e


def run_blast(cmds, verbose=False):
    """
    DEPRECATED: Use run_blastn() or run_self_blast() instead.
    Legacy function maintained for backwards compatibility.
    """
    import warnings

    warnings.warn(
        'run_blast() is deprecated. Use run_blastn() or run_self_blast() instead.',
        DeprecationWarning,
        stacklevel=2,
    )

    # This is a simplified implementation for backwards compatibility
    # The original function used shell scripts which is not recommended
    for cmd in cmds:
        try:
            result = subprocess.run(
                cmd,
                shell=True,  # Kept for compatibility but not recommended
                capture_output=True,
                text=True,
                timeout=300,
            )
            if result.returncode != 0:
                raise Error(f'Command failed: {cmd}\nError: {result.stderr}')
            if verbose:
                print(result.stdout)
        except subprocess.TimeoutExpired as err:
            raise Error(f'Command timed out: {cmd}') from err
        except Exception as e:
            raise Error(f'Error running command: {cmd}\nError: {str(e)}') from e


# Keep original Error class for backwards compatibility
class Error(Exception):
    pass


"""
Recreate pymummer-like coords file output from blast.

#Blast field, coords field, Description

qstart			[S1] 	Start of the alignment region in the reference sequence
qend			[E1] 	End of the alignment region in the reference sequence
sstart			[S2] 	Start of the alignment region in the query sequence
send 			[E2] 	End of the alignment region in the query sequence
length			[LEN 1] Length of the alignment region in the reference sequence
positive		[LEN 2] Length of the alignment region in the query sequence (#"positive" is just a filler as blast won't all repeated fields)
pident			[% IDY] Percent identity of the alignment
qlen			[LEN R] Length of the reference sequence
slen			[LEN Q] Length of the query sequence
qframe sframe	[FRM] 	Reading frame for the reference AND query sequence alignments respectively
qseqid sseqid	[TAGS] 	The reference AND query FastA IDs respectively. All output coordinates and lengths are relative to the forward strand of the reference DNA sequence.
"""
