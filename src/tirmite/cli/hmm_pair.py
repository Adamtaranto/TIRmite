#!/usr/bin/env python3
"""
TIRmite-pair: Pair precomputed nhmmer hits for TIR detection.

This module processes precomputed nhmmer search results:
1. Imports nhmmer hits from tabular files
2. Filters hits by model coverage and e-value thresholds
3. Applies pairing algorithms (symmetric or asymmetric)
4. Outputs paired elements, TIR hits, and GFF3 annotations

Supports both canonical (F,R) and custom strand orientations
for flexible transposon architecture detection.
"""

import argparse
import logging
import os
from pathlib import Path
import sys
from typing import Any, Dict, Optional, cast

from tirmite._version import __version__  # type: ignore[import-not-found]
import tirmite.tirmitetools as tirmite
from tirmite.utils.logs import init_logging
from tirmite.utils.utils import (
    cleanup_temp_directory,
    indexGenome,
    setup_directories,
)


def get_hmm_model_length(hmm_file_path: str) -> Dict[str, int]:
    """
    Extract HMM model lengths from HMM file by parsing LENG fields.

    Parameters
    ----------
    hmm_file_path : str or Path
        Path to HMM file containing one or more models.

    Returns
    -------
    dict
        Dictionary mapping model names to their lengths (in alignment columns).

    Notes
    -----
    Parses HMM file format looking for NAME and LENG lines.
    Handles multi-model HMM files, extracting length for each named model.
    """
    model_lengths = {}

    try:
        with open(hmm_file_path, 'r') as f:
            current_model = None
            for line in f:
                line = line.strip()
                if line.startswith('NAME  '):
                    current_model = line.split()[1]
                elif line.startswith('LENG  '):
                    if current_model:
                        length = int(line.split()[1])
                        model_lengths[current_model] = length
                elif line == '//':
                    current_model = None

    except Exception as e:
        logging.error(f'Error reading HMM file {hmm_file_path}: {e}')

    return model_lengths


def load_model_lengths_file(lengths_file: str) -> Dict[str, int]:
    """
    Load model lengths from tab-delimited text file.

    Parameters
    ----------
    lengths_file : str or Path
        Path to tab-delimited file with format: model_name<TAB>model_length.

    Returns
    -------
    dict
        Dictionary mapping model names to integer lengths.

    Notes
    -----
    Skips comment lines (starting with #) and blank lines.
    Logs warnings for malformed lines but continues parsing.
    """
    model_lengths = {}

    try:
        with open(lengths_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) != 2:
                    logging.warning(
                        f'Skipping malformed line {line_num} in {lengths_file}: {line}'
                    )
                    continue

                model_name, model_length_str = parts
                try:
                    model_length = int(model_length_str)
                    model_lengths[model_name] = model_length
                except ValueError:
                    logging.warning(
                        f'Invalid model length on line {line_num}: {model_length_str}'
                    )

    except Exception as e:
        logging.error(f'Error reading model lengths file {lengths_file}: {e}')

    return model_lengths


def calculate_hit_coverage(hitTable: Any, model_lengths: Dict[str, int]) -> Any:
    """
    Calculate coverage for hits based on model lengths.

    Parameters
    ----------
    hitTable : pandas.DataFrame
        DataFrame with hits containing model, hitStart, hitEnd columns.
    model_lengths : dict
        Dictionary mapping model names to their lengths.

    Returns
    -------
    pandas.DataFrame
        DataFrame with coverage column added.
    """
    hitTable = hitTable.copy()
    coverage_values = []

    for _idx, row in hitTable.iterrows():
        model = row['model']
        hit_length = abs(int(row['hitEnd']) - int(row['hitStart'])) + 1

        if model in model_lengths:
            model_length = model_lengths[model]
            coverage = hit_length / model_length
        else:
            logging.warning(f'Model length not found for {model}, using coverage = 0')
            coverage = 0.0

        coverage_values.append(coverage)

    hitTable['coverage'] = coverage_values
    return hitTable


def filter_hits_coverage(hitTable: Any, mincov: float) -> Any:
    """
    Filter hits by coverage threshold.

    Parameters
    ----------
    hitTable : pandas.DataFrame
        DataFrame with coverage column.
    mincov : float
        Minimum coverage threshold (0.0 to 1.0).

    Returns
    -------
    pandas.DataFrame
        Filtered DataFrame containing only hits with coverage >= mincov.
    """
    return hitTable[hitTable['coverage'] >= mincov]


def extract_model_name_from_path(model_path: Optional[str]) -> Optional[str]:
    """
    Extract model name from HMM file path by reading the HMM file.

    Parameters
    ----------
    model_path : str or Path
        Path to HMM file.

    Returns
    -------
    str or None
        Model name from NAME field, filename stem if not found, or None if no path.
    """
    if not model_path:
        return None

    try:
        with open(model_path, 'r') as f:
            for line in f:
                if line.startswith('NAME  '):
                    return line.split()[1].strip()
    except (FileNotFoundError, IOError):
        # Fallback to basename without extension
        return Path(model_path).stem

    return Path(model_path).stem


def load_pairing_map(pairing_map_file: str) -> list[tuple[str, str]]:
    """
    Load pairing map from tab-delimited file.

    Parameters
    ----------
    pairing_map_file : str
        Path to tab-delimited file with left and right feature names.

    Returns
    -------
    list of tuple
        List of (left_feature, right_feature) tuples.

    Raises
    ------
    ValueError
        If file format is invalid.
    FileNotFoundError
        If file doesn't exist.

    Notes
    -----
    File format: left_feature<TAB>right_feature
    For symmetric pairing, both columns have the same value.
    Skips comment lines (starting with #) and blank lines.
    """
    pairings = []
    seen_features = {}  # Track occurrences of each feature

    try:
        with open(pairing_map_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                # Skip comments and blank lines
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) != 2:
                    raise ValueError(
                        f'Invalid format on line {line_num}: expected 2 tab-delimited columns, '
                        f'got {len(parts)}'
                    )

                left_feature, right_feature = parts[0].strip(), parts[1].strip()

                if not left_feature or not right_feature:
                    raise ValueError(f'Empty feature name on line {line_num}')

                pairings.append((left_feature, right_feature))

                # Track feature occurrences for warning
                for feature in [left_feature, right_feature]:
                    if feature not in seen_features:
                        seen_features[feature] = []
                    seen_features[feature].append(line_num)

    except FileNotFoundError:
        raise FileNotFoundError(
            f'Pairing map file not found: {pairing_map_file}'
        ) from None
    except Exception as e:
        raise ValueError(f'Error reading pairing map file: {e}') from e

    if not pairings:
        raise ValueError(f'No valid pairings found in {pairing_map_file}')

    # Warn about features appearing in multiple pairings
    for feature, line_nums in seen_features.items():
        if len(line_nums) > 1:
            logging.warning(
                f'Feature "{feature}" appears in multiple pairing combinations '
                f'(lines: {", ".join(map(str, line_nums))})'
            )

    logging.info(f'Loaded {len(pairings)} pairing combinations from {pairing_map_file}')
    return pairings


def check_multiple_models(hitTable: Any) -> list[str]:
    """
    Check if hitTable contains hits from multiple models/queries.

    Parameters
    ----------
    hitTable : pandas.DataFrame
        DataFrame containing hit data with 'model' column.

    Returns
    -------
    list of str
        List of unique model names found in hitTable.
    """
    unique_models = hitTable['model'].unique().tolist()
    return unique_models


def check_overlapping_hits(left_hitTable: Any, right_hitTable: Any) -> tuple[int, int]:
    """
    Check for overlapping hits between left and right hit tables.

    Two hits overlap if they are on the same target sequence and their
    genomic coordinates overlap (considering both forward and reverse strands).

    Parameters
    ----------
    left_hitTable : pandas.DataFrame
        DataFrame containing hits from left file.
    right_hitTable : pandas.DataFrame
        DataFrame containing hits from right file.

    Returns
    -------
    tuple of (int, int)
        Number of left hits that overlap with right hits,
        Number of right hits that overlap with left hits.

    Notes
    -----
    This function uses an efficient algorithm that groups hits by target
    and checks for coordinate overlaps only within each target group.
    """
    if left_hitTable.empty or right_hitTable.empty:
        return 0, 0

    # Convert hitStart and hitEnd to integers for comparison
    left_hitTable = left_hitTable.copy()
    right_hitTable = right_hitTable.copy()
    left_hitTable['hitStart'] = left_hitTable['hitStart'].astype(int)
    left_hitTable['hitEnd'] = left_hitTable['hitEnd'].astype(int)
    right_hitTable['hitStart'] = right_hitTable['hitStart'].astype(int)
    right_hitTable['hitEnd'] = right_hitTable['hitEnd'].astype(int)

    # Group hits by target sequence for efficient lookup
    left_by_target = left_hitTable.groupby('target')
    right_by_target = right_hitTable.groupby('target')

    # Find common targets
    left_targets = set(left_hitTable['target'].unique())
    right_targets = set(right_hitTable['target'].unique())
    common_targets = left_targets & right_targets

    if not common_targets:
        return 0, 0

    left_overlapping_indices = set()
    right_overlapping_indices = set()

    # Check for overlaps only on common targets
    for target in common_targets:
        left_hits = left_by_target.get_group(target)
        right_hits = right_by_target.get_group(target)

        # Check each left hit against each right hit on this target
        for left_idx, left_row in left_hits.iterrows():
            left_start = left_row['hitStart']
            left_end = left_row['hitEnd']

            for right_idx, right_row in right_hits.iterrows():
                right_start = right_row['hitStart']
                right_end = right_row['hitEnd']

                # Check if coordinates overlap
                # Two intervals [a,b] and [c,d] overlap if max(a,c) <= min(b,d)
                if max(left_start, right_start) <= min(left_end, right_end):
                    left_overlapping_indices.add(left_idx)
                    right_overlapping_indices.add(right_idx)

    return len(left_overlapping_indices), len(right_overlapping_indices)


def create_pair_parser() -> argparse.ArgumentParser:
    """
    Create standalone argument parser for pair command.

    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser for pair workflow options.
    """
    parser = argparse.ArgumentParser(
        prog='tirmite-pair',
        description='Pair precomputed nhmmer hits for transposon detection',
    )
    _configure_pair_parser(parser)
    return parser


def add_pair_parser(subparsers: Any) -> argparse.ArgumentParser:
    """
    Add pair subcommand parser.

    Parameters
    ----------
    subparsers : argparse._SubParsersAction
        Subparser object to add pair command to.

    Returns
    -------
    argparse.ArgumentParser
        The configured pair subcommand parser.
    """
    parser = cast(
        argparse.ArgumentParser,
        subparsers.add_parser(
            'pair',
            help='Pair precomputed nhmmer hits',
            description='Pair precomputed nhmmer hits for transposon detection',
        ),
    )
    _configure_pair_parser(parser)
    return parser


def _configure_pair_parser(parser: argparse.ArgumentParser) -> None:
    """
    Configure parser with pair command arguments.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser to configure.

    Returns
    -------
    None
        Modifies parser in place.
    """

    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__),
    )

    parser.add_argument(
        '--loglevel',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        help="Set logging level. Default: 'INFO'",
    )

    parser.add_argument(
        '--logfile',
        action='store_true',
        default=False,
        help='Write log messages to file in output directory.',
    )

    # Required inputs
    parser.add_argument(
        '--genome',
        type=str,
        required=False,
        help='Path to target genome FASTA file.',
    )

    parser.add_argument(
        '--blastdb',
        type=str,
        required=False,
        help='Path to BLAST database (alternative to --genome for sequence extraction).',
    )

    # Search result input files (mutually exclusive groups)
    search_group = parser.add_mutually_exclusive_group(required=True)
    search_group.add_argument(
        '--nhmmerFile',
        type=str,
        help='Path to single nhmmer output file (requires --hmmFile or --lengthsFile).',
    )
    search_group.add_argument(
        '--leftNhmmer',
        type=str,
        help='Path to nhmmer output for left model (use with --rightNhmmer).',
    )
    search_group.add_argument(
        '--blastFile',
        type=str,
        help='Path to single BLAST tabular output file (requires --queryLen or --lengthsFile).',
    )
    search_group.add_argument(
        '--leftBlast',
        type=str,
        help='Path to BLAST output for left query (use with --rightBlast).',
    )

    parser.add_argument(
        '--rightNhmmer',
        type=str,
        help='Path to nhmmer output for right model (use with --leftNhmmer).',
    )

    parser.add_argument(
        '--rightBlast',
        type=str,
        help='Path to BLAST output for right query (use with --leftBlast).',
    )

    # Model/Query length sources (mutually exclusive)
    length_group = parser.add_mutually_exclusive_group()
    length_group.add_argument(
        '--hmmFile',
        type=str,
        help='Path to HMM file for extracting model lengths (for single model pairing).',
    )
    length_group.add_argument(
        '--leftModel',
        type=str,
        help='Path to left HMM model file (for asymmetric pairing).',
    )
    length_group.add_argument(
        '--lengthsFile',
        type=str,
        help='Path to tab-delimited file with model_name and model_length columns.',
    )
    length_group.add_argument(
        '--queryLen',
        type=int,
        help='Length of BLAST query sequence (for single query pairing).',
    )

    parser.add_argument(
        '--rightModel',
        type=str,
        help='Path to right HMM model file (for asymmetric pairing).',
    )

    # Filtering parameters
    parser.add_argument(
        '--maxeval',
        type=float,
        default=0.001,
        help='Maximum e-value allowed for valid hit. Default: 0.001',
    )

    parser.add_argument(
        '--mincov',
        type=float,
        default=0.5,
        help='Minimum hit coverage as proportion of model length. Default: 0.5',
    )

    parser.add_argument(
        '--maxdist',
        type=int,
        default=None,
        help='Maximum distance allowed between termini for pairing.',
    )

    # Pairing configuration
    parser.add_argument(
        '--orientation',
        type=str,
        default='F,R',
        help='Orientation pattern for pairing. F=Forward(+), R=Reverse(-). Options: F,R (TIR), F,F (LTR), R,R, R,F. Default: F,R',
    )

    parser.add_argument(
        '--stableReps',
        type=int,
        default=0,
        help='Number of iterations when no new pairs found. Default: 0',
    )

    parser.add_argument(
        '--pairing_map',
        type=str,
        default=None,
        help='Tab-delimited file mapping left to right feature names for pairing. Required when input contains multiple models/queries.',
    )

    # Output options
    parser.add_argument(
        '--outdir',
        type=str,
        default=None,
        help='Output directory. Default: current directory',
    )

    parser.add_argument(
        '--prefix',
        type=str,
        default=None,
        help='Prefix for output files.',
    )

    parser.add_argument(
        '--nopairing',
        action='store_true',
        default=False,
        help='Only report individual hits, skip pairing.',
    )

    parser.add_argument(
        '--gffOut',
        action='store_true',
        default=False,
        help='Generate GFF3 output file.',
    )

    parser.add_argument(
        '--report',
        default='all',
        choices=['all', 'paired', 'unpaired'],
        help='Types of hits to include in GFF output.',
    )

    parser.add_argument(
        '--padlen',
        type=int,
        default=None,
        help='Extract N bases flanking each hit in FASTA output.',
    )

    # Utility options
    parser.add_argument(
        '--tempdir',
        type=str,
        default=None,
        help='Base directory for temporary files.',
    )

    parser.add_argument(
        '--keep-temp',
        action='store_true',
        default=False,
        help='Preserve temporary directory.',
    )


def validate_arguments(args: Any) -> None:
    """
    Validate argument combinations and file existence.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    None
        No return value.

    Raises
    ------
    ValueError
        If argument combinations are invalid.
    FileNotFoundError
        If required files don't exist.
    """
    # Check that either genome or blastdb is provided
    if not args.genome and not args.blastdb:
        raise ValueError('Either --genome or --blastdb must be provided')

    # Check asymmetric pairing requirements
    if args.leftNhmmer and not args.rightNhmmer:
        raise ValueError('--leftNhmmer requires --rightNhmmer')
    if args.rightNhmmer and not args.leftNhmmer:
        raise ValueError('--rightNhmmer requires --leftNhmmer')
    if args.leftBlast and not args.rightBlast:
        raise ValueError('--leftBlast requires --rightBlast')
    if args.rightBlast and not args.leftBlast:
        raise ValueError('--rightBlast requires --leftBlast')
    if args.leftModel and not args.rightModel:
        raise ValueError('--leftModel requires --rightModel')
    if args.rightModel and not args.leftModel:
        raise ValueError('--rightModel requires --leftModel')

    # Check model/query length source requirements
    if args.nhmmerFile:
        if not (args.hmmFile or args.lengthsFile):
            raise ValueError('--nhmmerFile requires either --hmmFile or --lengthsFile')

    if args.blastFile:
        if not (args.queryLen or args.lengthsFile):
            raise ValueError('--blastFile requires either --queryLen or --lengthsFile')

    if args.leftNhmmer and args.rightNhmmer:
        if not (args.leftModel and args.rightModel) and not args.lengthsFile:
            raise ValueError(
                'Asymmetric pairing requires --leftModel/--rightModel or --lengthsFile'
            )

    if args.leftBlast and args.rightBlast:
        if not args.lengthsFile:
            raise ValueError(
                'Asymmetric BLAST pairing requires --lengthsFile with query lengths'
            )

    # Check file existence
    required_files = []
    if args.genome:
        required_files.append(args.genome)

    if args.nhmmerFile:
        required_files.append(args.nhmmerFile)
    if args.leftNhmmer:
        required_files.extend([args.leftNhmmer, args.rightNhmmer])
    if args.blastFile:
        required_files.append(args.blastFile)
    if args.leftBlast:
        required_files.extend([args.leftBlast, args.rightBlast])
    if args.hmmFile:
        required_files.append(args.hmmFile)
    if args.leftModel:
        required_files.extend([args.leftModel, args.rightModel])
    if args.lengthsFile:
        required_files.append(args.lengthsFile)
    if args.pairing_map:
        required_files.append(args.pairing_map)

    for file_path in required_files:
        if not Path(file_path).exists():
            raise FileNotFoundError(f'Required file not found: {file_path}')

    # Validate blastdb existence if provided
    if args.blastdb:
        # Check if any of the expected BLAST DB files exist
        db_extensions = ['.nhr', '.nin', '.nsq', '.ndb', '.not', '.ntf', '.nto']
        db_exists = any(Path(f'{args.blastdb}{ext}').exists() for ext in db_extensions)
        if not db_exists:
            raise FileNotFoundError(
                f'BLAST database files not found for: {args.blastdb}'
            )


def main(args: Optional[argparse.Namespace] = None) -> int:
    """
    Main entry point for tirmite-pair.

    Parameters
    ----------
    args : argparse.Namespace, optional
        Parsed command-line arguments. If None, parses from sys.argv.

    Returns
    -------
    int
        Exit code (0 for success, 1 for error).
    """
    # Parse arguments if not provided
    if args is None:
        parser = create_pair_parser()
        args = parser.parse_args()

    # Mypy assertion: args is guaranteed non-None after parsing
    assert args is not None

    try:
        # Validate arguments
        try:
            validate_arguments(args)
        except (ValueError, FileNotFoundError) as e:
            logging.error(f'Argument validation failed: {e}')
            sys.exit(1)

        # Set up directories first (we need output dir for default logfile location)
        try:
            outDir, tempDir = setup_directories(args)
        except (OSError, FileNotFoundError, PermissionError) as e:
            # If directory setup fails, we can't even set up logging properly
            print(f'Directory setup failed: {e}', file=sys.stderr)
            sys.exit(1)

        # Determine logfile path
        logfile_path = None
        if args.logfile:
            # Create logfile in output directory
            if args.prefix:
                logfile_name = f'{args.prefix}_tirmite_pair.log'
            else:
                logfile_name = 'tirmite_pair.log'
            logfile_path = outDir / logfile_name

        # Set up logging with or without file output
        init_logging(loglevel=args.loglevel, logfile=logfile_path)

        logging.info(f'TIRmite-pair version {__version__}')
        logging.info(f'Output directory: {outDir}')
        logging.debug(f'Temporary directory: {tempDir}')

        # Index genome if provided
        genome = None
        genome_descriptions = None
        if args.genome:
            logging.info(f'Indexing genome: {args.genome}')
            genome, genome_descriptions = indexGenome(args.genome)
        elif args.blastdb:
            logging.info(f'Using BLAST database: {args.blastdb}')
            # Note: genome will remain None, sequence extraction will use blastdbcmd

        # Load model/query lengths
        logging.info('Loading model/query lengths...')
        model_lengths = {}

        if args.lengthsFile:
            model_lengths = load_model_lengths_file(args.lengthsFile)
        elif args.hmmFile:
            model_lengths = get_hmm_model_length(args.hmmFile)
        elif args.queryLen:
            # For single BLAST query, we'll assign the length after importing hits
            # to get the query name
            pass
        elif args.leftModel and args.rightModel:
            left_lengths = get_hmm_model_length(args.leftModel)
            right_lengths = get_hmm_model_length(args.rightModel)
            model_lengths = {**left_lengths, **right_lengths}

        # Import search hits
        hitTable = None
        input_format = None

        if args.nhmmerFile:
            # Single nhmmer file mode
            logging.info('Importing nhmmer hits...')
            input_format = 'nhmmer'
            hitTable = tirmite.import_nhmmer(
                infile=args.nhmmerFile, hitTable=None, prefix=args.prefix
            )
        elif args.leftNhmmer:
            # Asymmetric nhmmer mode - import from both files
            logging.info('Importing nhmmer hits from left and right models...')
            input_format = 'nhmmer'

            # Check if left and right files are the same
            if args.leftNhmmer == args.rightNhmmer:
                raise ValueError(
                    f'Left and right nhmmer files cannot be the same: {args.leftNhmmer}'
                )

            left_model_name = extract_model_name_from_path(args.leftModel)
            right_model_name = extract_model_name_from_path(args.rightModel)

            # Import left file
            left_hitTable = tirmite.import_nhmmer(
                infile=args.leftNhmmer,
                hitTable=None,
                prefix=args.prefix,
            )
            left_models = check_multiple_models(left_hitTable)
            logging.info(
                f'Left nhmmer file: {len(left_hitTable)} hits, {len(left_models)} unique query/model name(s)'
            )

            # Import right file
            right_hitTable = tirmite.import_nhmmer(
                infile=args.rightNhmmer,
                hitTable=None,
                prefix=args.prefix,
            )
            right_models = check_multiple_models(right_hitTable)
            logging.info(
                f'Right nhmmer file: {len(right_hitTable)} hits, {len(right_models)} unique query/model name(s)'
            )

            # Check for overlapping query names
            overlapping_models = set(left_models) & set(right_models)
            if overlapping_models:
                logging.warning(
                    f'Query/model names appear in both left and right files: {", ".join(overlapping_models)}'
                )

            # Check for overlapping hits (same target, overlapping coordinates)
            left_overlap_count, right_overlap_count = check_overlapping_hits(
                left_hitTable, right_hitTable
            )
            if left_overlap_count > 0 or right_overlap_count > 0:
                logging.warning(
                    f'Found {left_overlap_count} left hit(s) and {right_overlap_count} right hit(s) '
                    'with overlapping genomic coordinates'
                )

            # Validate single query per file or require pairing_map
            if len(left_models) > 1 or len(right_models) > 1:
                if not args.pairing_map:
                    raise ValueError(
                        f'Left file contains {len(left_models)} query/model name(s), '
                        f'right file contains {len(right_models)} query/model name(s). '
                        'When either file contains multiple queries, --pairing_map is required.'
                    )

            # Combine hit tables
            hitTable = tirmite.import_nhmmer(
                infile=args.leftNhmmer,
                hitTable=None,
                prefix=args.prefix,
            )
            hitTable = tirmite.import_nhmmer(
                infile=args.rightNhmmer,
                hitTable=hitTable,
                prefix=args.prefix,
            )
        elif args.blastFile:
            # Single BLAST file mode
            logging.info('Importing BLAST hits...')
            input_format = 'blast'

            # Detect format and warn if mismatch
            detected_format = tirmite.detect_input_format(args.blastFile)
            if detected_format != 'blast' and detected_format != 'unknown':
                logging.warning(
                    f'File format appears to be {detected_format}, but --blastFile was specified. '
                    'Consider using --nhmmerFile instead.'
                )

            hitTable = tirmite.import_blast(
                infile=args.blastFile, hitTable=None, prefix=args.prefix
            )

        elif args.leftBlast:
            # Asymmetric BLAST mode - import from both files
            logging.info('Importing BLAST hits from left and right queries...')
            input_format = 'blast'

            # Check if left and right files are the same
            if args.leftBlast == args.rightBlast:
                raise ValueError(
                    f'Left and right BLAST files cannot be the same: {args.leftBlast}'
                )

            # Detect format for both files
            detected_left = tirmite.detect_input_format(args.leftBlast)
            detected_right = tirmite.detect_input_format(args.rightBlast)

            if detected_left != 'blast' and detected_left != 'unknown':
                logging.warning(
                    f'Left file format appears to be {detected_left}, but --leftBlast was specified.'
                )
            if detected_right != 'blast' and detected_right != 'unknown':
                logging.warning(
                    f'Right file format appears to be {detected_right}, but --rightBlast was specified.'
                )

            # Import left file
            left_hitTable = tirmite.import_blast(
                infile=args.leftBlast,
                hitTable=None,
                prefix=args.prefix,
            )
            left_models = check_multiple_models(left_hitTable)
            logging.info(
                f'Left BLAST file: {len(left_hitTable)} hits, {len(left_models)} unique query/model name(s)'
            )

            # Import right file
            right_hitTable = tirmite.import_blast(
                infile=args.rightBlast,
                hitTable=None,
                prefix=args.prefix,
            )
            right_models = check_multiple_models(right_hitTable)
            logging.info(
                f'Right BLAST file: {len(right_hitTable)} hits, {len(right_models)} unique query/model name(s)'
            )

            # Check for overlapping query names
            overlapping_models = set(left_models) & set(right_models)
            if overlapping_models:
                logging.warning(
                    f'Query/model names appear in both left and right files: {", ".join(overlapping_models)}'
                )

            # Check for overlapping hits (same target, overlapping coordinates)
            left_overlap_count, right_overlap_count = check_overlapping_hits(
                left_hitTable, right_hitTable
            )
            if left_overlap_count > 0 or right_overlap_count > 0:
                logging.warning(
                    f'Found {left_overlap_count} left hit(s) and {right_overlap_count} right hit(s) '
                    'with overlapping genomic coordinates'
                )

            # Validate single query per file or require pairing_map
            if len(left_models) > 1 or len(right_models) > 1:
                if not args.pairing_map:
                    raise ValueError(
                        f'Left file contains {len(left_models)} query/model name(s), '
                        f'right file contains {len(right_models)} query/model name(s). '
                        'When either file contains multiple queries, --pairing_map is required.'
                    )

            # Combine hit tables
            hitTable = tirmite.import_blast(
                infile=args.leftBlast,
                hitTable=None,
                prefix=args.prefix,
            )
            hitTable = tirmite.import_blast(
                infile=args.rightBlast,
                hitTable=hitTable,
                prefix=args.prefix,
            )

        if hitTable is None or len(hitTable) == 0:
            logging.error('No hits were imported from input files')
            cleanup_temp_directory(tempDir, args.keep_temp)
            sys.exit(1)

        # Log import statistics
        logging.info(f'Imported {len(hitTable)} total hits')

        # Get unique models and log statistics
        unique_models = check_multiple_models(hitTable)
        logging.info(f'Found {len(unique_models)} unique query/model name(s)')

        # Log per-query hit counts at debug level
        for model in unique_models:
            hit_count = len(hitTable[hitTable['model'] == model])
            logging.debug(f'Query/model "{model}": {hit_count} hits')

        # If queryLen was provided for BLAST input, assign it to ALL queries
        if args.blastFile and args.queryLen:
            for query_name in unique_models:
                model_lengths[query_name] = args.queryLen
                logging.debug(f'Set length for query {query_name}: {args.queryLen}')
            logging.info(
                f'Applied query length {args.queryLen} to {len(unique_models)} query name(s)'
            )

        # Validate that we have model lengths for all models in hitTable
        if not model_lengths:
            logging.error('No model/query lengths could be loaded')
            cleanup_temp_directory(tempDir, args.keep_temp)
            sys.exit(1)

        logging.info(f'Loaded lengths for models: {", ".join(model_lengths.keys())}')

        # Calculate hit coverage
        logging.info('Calculating hit coverage...')
        hitTable = calculate_hit_coverage(hitTable, model_lengths)

        # Apply coverage filtering
        logging.info(f'Filtering hits with coverage < {args.mincov}')
        hitCount = len(hitTable)

        # Log pre-filtering counts
        model_counts_before = hitTable['model'].value_counts().to_dict()
        for model, count in model_counts_before.items():
            logging.info(f'Model {model}: {count} hits before coverage filtering')

        hitTable = filter_hits_coverage(hitTable, args.mincov)

        # Log post-filtering counts
        model_counts_after = hitTable['model'].value_counts().to_dict()
        total_excluded = hitCount - len(hitTable)
        logging.info(f'Excluded {total_excluded} hits on coverage criteria')

        for model in model_counts_before:
            before = model_counts_before[model]
            after = model_counts_after.get(model, 0)
            excluded = before - after
            logging.info(f'Model {model}: {excluded} excluded, {after} remaining')

        # Apply e-value filtering
        logging.info(f'Filtering hits with e-value > {args.maxeval}')
        hitCount = len(hitTable)

        model_counts_before = hitTable['model'].value_counts().to_dict()
        hitTable = tirmite.filterHitsEval(maxeval=args.maxeval, hitTable=hitTable)

        model_counts_after = hitTable['model'].value_counts().to_dict()
        total_excluded = hitCount - len(hitTable)
        logging.info(f'Excluded {total_excluded} hits on e-value criteria')

        for model in model_counts_before:
            before = model_counts_before[model]
            after = model_counts_after.get(model, 0)
            excluded = before - after
            logging.info(f'Model {model}: {excluded} excluded, {after} remaining')

        logging.info(f'Total remaining hits: {len(hitTable)}')

        # Convert to dict structure
        hitsDict, hitIndex = tirmite.table2dict(hitTable)

        # Check for multiple models and validate pairing map requirement
        # Note: unique_models was already determined and logged after import
        # For asymmetric modes, validation was already done per-file
        is_asymmetric = (args.leftNhmmer and args.rightNhmmer) or (
            args.leftBlast and args.rightBlast
        )

        # Load pairing map if provided
        pairing_map = None
        if args.pairing_map:
            pairing_map = load_pairing_map(args.pairing_map)
        elif len(unique_models) > 1 and not is_asymmetric:
            # Multiple models in single file without pairing map - raise error
            # (Asymmetric mode already validated per-file)
            raise ValueError(
                f'Input contains {len(unique_models)} distinct models/queries: {", ".join(unique_models)}. '
                'Multiple models require --pairing_map to specify which features should be paired together.'
            )

        # Create pairing configuration
        if pairing_map:
            # Use pairing map - will create configs for each pairing later
            # Pairing map workflow handles both single and multiple pairings
            logging.info(
                f'Will execute {len(pairing_map)} independent pairing procedure(s) based on pairing map'
            )
            config = None
        elif args.leftNhmmer and args.rightNhmmer:
            # Asymmetric nhmmer pairing
            left_model_name = extract_model_name_from_path(args.leftModel)
            right_model_name = extract_model_name_from_path(args.rightModel)

            config = tirmite.PairingConfig(
                orientation=args.orientation,
                left_model=left_model_name,
                right_model=right_model_name,
            )
        elif args.leftBlast and args.rightBlast:
            # Asymmetric BLAST pairing - extract model names from hitTable
            blast_models = hitTable['model'].unique()
            if len(blast_models) < 2:
                logging.error(
                    f'Expected 2 models for asymmetric pairing, found {len(blast_models)}'
                )
                cleanup_temp_directory(tempDir, args.keep_temp)
                sys.exit(1)

            # Note: Model assignment is based on order of appearance in combined hit table
            # The first model encountered becomes 'left', second becomes 'right'
            # For deterministic results, ensure leftBlast file contains only left query hits
            # and rightBlast file contains only right query hits
            left_model_name = blast_models[0]
            right_model_name = blast_models[1]
            logging.info(
                f'Assigning models for asymmetric pairing: '
                f'left={left_model_name}, right={right_model_name}'
            )
            logging.info(
                'Note: First unique model becomes "left", second becomes "right" '
                'based on order in input files'
            )

            config = tirmite.PairingConfig(
                orientation=args.orientation,
                left_model=left_model_name,
                right_model=right_model_name,
            )
        else:
            # Symmetric pairing
            available_models = list(hitsDict.keys())
            if not available_models:
                logging.error('No models found in hits')
                cleanup_temp_directory(tempDir, args.keep_temp)
                sys.exit(1)

            single_model = available_models[0]
            config = tirmite.PairingConfig(
                orientation=args.orientation, single_model=single_model
            )

        # Write individual hits
        logging.info('Writing individual hits to FASTA...')
        tirmite.writeTIRs(
            outDir=outDir,
            hitTable=hitTable,
            maxeval=args.maxeval,
            genome=genome,
            prefix=args.prefix,
            padlen=args.padlen,
            genome_descriptions=genome_descriptions,
            blastdb=args.blastdb if args.blastdb else None,
        )

        # Skip pairing if requested
        if args.nopairing:
            logging.info('Pairing disabled. Analysis complete.')
            cleanup_temp_directory(tempDir, args.keep_temp)
            return 0

        # Run pairing - handle single or multiple pairing procedures
        if pairing_map:
            # Multiple pairing procedures based on pairing map
            logging.info(
                f'Running {len(pairing_map)} independent pairing procedures based on pairing map'
            )

            all_paired = {}  # Accumulate all paired results by model
            all_paired_hits = set()  # Track which hit indices have been paired
            original_hitIndex = hitIndex  # Preserve original for unpaired hit tracking

            for pair_idx, (left_feature, right_feature) in enumerate(pairing_map, 1):
                logging.info(
                    f'Pairing procedure {pair_idx}/{len(pairing_map)}: {left_feature} <-> {right_feature}'
                )

                # Check if both features exist in hitsDict
                if left_feature not in hitsDict:
                    logging.warning(
                        f'Feature "{left_feature}" not found in hits, skipping this pairing'
                    )
                    continue
                if right_feature not in hitsDict:
                    logging.warning(
                        f'Feature "{right_feature}" not found in hits, skipping this pairing'
                    )
                    continue

                # Create config for this pairing
                if left_feature == right_feature:
                    # Symmetric pairing
                    pair_config = tirmite.PairingConfig(
                        orientation=args.orientation, single_model=left_feature
                    )
                else:
                    # Asymmetric pairing
                    pair_config = tirmite.PairingConfig(
                        orientation=args.orientation,
                        left_model=left_feature,
                        right_model=right_feature,
                    )

                # Run pairing for this combination (use fresh copy of hitIndex for each)
                logging.info(
                    f'Searching for candidate pairings: {left_feature} <-> {right_feature}'
                )
                pair_hitIndex = tirmite.parseHitsGeneral(
                    hitsDict=hitsDict,
                    hitIndex=original_hitIndex,
                    maxDist=args.maxdist,
                    config=pair_config,
                )

                logging.info('Performing iterative pairing...')
                if pair_config.is_asymmetric:
                    logging.info(
                        f'Using asymmetric pairing: {pair_config.left_model} + {pair_config.right_model}'
                    )
                    pair_hitIndex, pair_paired, pair_unpaired = (
                        tirmite.iterateGetPairsAsymmetric(
                            pair_hitIndex, pair_config, stableReps=args.stableReps
                        )
                    )
                else:
                    logging.info(
                        f'Using symmetric pairing with orientation {pair_config.orientation}'
                    )
                    pair_hitIndex, pair_paired, pair_unpaired = (
                        tirmite.iterateGetPairsCustom(
                            pair_hitIndex, pair_config, stableReps=args.stableReps
                        )
                    )

                # Accumulate paired results
                for model, pairs in pair_paired.items():
                    if model not in all_paired:
                        all_paired[model] = []
                    all_paired[model].extend(pairs)
                    # Track paired hit indices
                    for pair_set in pairs:
                        all_paired_hits.update(pair_set)

                # Log this pairing's results
                total_pairs = sum(len(pairs) for pairs in pair_paired.values())
                logging.info(
                    f'Pairing procedure {pair_idx} completed: {total_pairs} pairs found'
                )

            # Final pairing results
            paired = all_paired
            hitIndex = original_hitIndex  # Use original hitIndex for output

            # Collect truly unpaired hits (not paired in any procedure)
            unpaired = []
            for model in hitsDict.keys():
                if model in original_hitIndex:
                    for hit_id in original_hitIndex[model].keys():
                        if hit_id not in all_paired_hits:
                            unpaired.append(hit_id)

            # Log final pairing results
            total_pairs = sum(len(pairs) for pairs in paired.values())
            total_unpaired = len(unpaired)
            logging.info(
                f'All pairing procedures completed: {total_pairs} total pairs, {total_unpaired} unpaired hits'
            )

        else:
            # Single pairing procedure (original logic)
            logging.info('Searching for candidate pairings...')
            hitIndex = tirmite.parseHitsGeneral(
                hitsDict=hitsDict,
                hitIndex=hitIndex,
                maxDist=args.maxdist,
                config=config,
            )

            logging.info('Performing iterative pairing...')
            if config.is_asymmetric:
                logging.info(
                    f'Using asymmetric pairing: {config.left_model} + {config.right_model}'
                )
                hitIndex, paired, unpaired = tirmite.iterateGetPairsAsymmetric(
                    hitIndex, config, stableReps=args.stableReps
                )
            else:
                logging.info(
                    f'Using symmetric pairing with orientation {config.orientation}'
                )
                hitIndex, paired, unpaired = tirmite.iterateGetPairsCustom(
                    hitIndex, config, stableReps=args.stableReps
                )

            # Log pairing results
            total_pairs = sum(len(pairs) for pairs in paired.values())
            total_unpaired = len(unpaired)
            logging.info(
                f'Pairing completed: {total_pairs} pairs, {total_unpaired} unpaired'
            )

        # Write paired TIRs
        if args.report in ['all', 'paired']:
            logging.info('Writing paired termini to FASTA...')
            tirmite.writePairedTIRs(
                outDir=outDir,
                paired=paired,
                hitIndex=hitIndex,
                genome=genome,
                prefix=args.prefix,
                padlen=args.padlen,
                genome_descriptions=genome_descriptions,
                blastdb=args.blastdb if args.blastdb else None,
            )

        # Extract and write elements
        logging.info('Extracting paired elements...')
        pairedEles = tirmite.fetchElements(
            paired=paired,
            hitIndex=hitIndex,
            genome=genome,
            genome_descriptions=genome_descriptions,
            blastdb=args.blastdb if args.blastdb else None,
        )

        logging.info('Writing paired elements to FASTA...')
        tirmite.writeElements(outDir, eleDict=pairedEles, prefix=args.prefix)

        # Write GFF if requested
        if args.gffOut:
            # Get unpaired TIRs if needed
            unpairedTIRs = None
            if args.report in ['all', 'unpaired']:
                unpairedTIRs = tirmite.fetchUnpaired(hitIndex=hitIndex)
                logging.info(f'Found {len(unpairedTIRs)} unpaired termini')

            # Set GFF output path
            if args.prefix:
                gffOutPath = os.path.join(outDir, f'{args.prefix}.gff3')
            else:
                gffOutPath = os.path.join(outDir, 'tirmite_pair_report.gff3')

            logging.info(f'Writing GFF3 output: {gffOutPath}')
            tirmite.gffWrite(
                outpath=gffOutPath,
                featureList=pairedEles,
                writeTIRs=args.report,
                unpaired=unpairedTIRs,
                prefix=args.prefix,
            )

        logging.info('TIRmite-pair analysis completed successfully')

    except KeyboardInterrupt:
        logging.info('Analysis interrupted by user')
        sys.exit(130)
    except Exception as e:
        logging.error(f'Unexpected error: {e}')
        logging.exception('Full traceback:')  # This will log the full stack trace
        sys.exit(1)
    finally:
        if 'tempDir' in locals():
            cleanup_temp_directory(tempDir, args.keep_temp)

        # Log completion message with logfile location if enabled
        if 'logfile_path' in locals() and logfile_path and args.logfile:
            logging.info(f'Log file saved to: {logfile_path}')

    return 0


if __name__ == '__main__':
    main()
