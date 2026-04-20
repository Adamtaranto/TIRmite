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
from typing import Any, Dict, List, Optional, cast

import pandas as pd  # type: ignore[import-untyped]

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


def compute_outer_edge_offset(
    hmm_start: int,
    hmm_end: int,
    model_len: int,
    strand: str,
    terminus_type: str,
) -> int:
    """
    Compute the offset between a hit alignment boundary and the outer edge of the query model.

    The "outer edge" is the external end of the terminus model:
    - For a left terminus the outer edge faces upstream (model position 1 on + strand).
    - For a right terminus the outer edge faces downstream (model position model_len on + strand).

    Parameters
    ----------
    hmm_start : int
        1-based start position of the alignment on the query/HMM model.
    hmm_end : int
        1-based end position of the alignment on the query/HMM model.
    model_len : int
        Total length of the query/HMM model.
    strand : str
        Strand of the hit: '+' or '-'.
    terminus_type : str
        'left' or 'right' terminus.

    Returns
    -------
    int
        Number of unaligned model positions between the hit and the outer edge.
        Zero means the hit reaches the outer edge exactly.
    """
    if terminus_type == 'left':
        if strand == '+':
            return hmm_start - 1
        else:  # '-'
            return model_len - hmm_end
    else:  # 'right'
        if strand == '+':
            return model_len - hmm_end
        else:  # '-'
            return hmm_start - 1


def filter_hits_by_anchor(
    hit_table: pd.DataFrame,
    model_lengths: Dict[str, int],
    max_offset: int,
    orientation: str = 'F,R',
    pairing_map: Optional[List] = None,
) -> pd.DataFrame:
    """
    Filter hits to only those anchored near the outer edge of the query model.

    For asymmetric pairings (different left/right models) or symmetric pairings
    with different strands (e.g. F,R), the terminus type is determined and only
    the external edge offset is checked.

    For symmetric same-strand pairings (F,F or R,R) without a pairing map, the
    hit must be within ``max_offset`` bases of **both** ends of the query model
    (i.e. both ``hmmStart - 1 <= max_offset`` and ``model_len - hmmEnd <= max_offset``).

    Parameters
    ----------
    hit_table : pandas.DataFrame
        Hit table with columns: model, strand, hmmStart, hmmEnd.
    model_lengths : dict
        Mapping of model name to model length.
    max_offset : int
        Maximum allowed offset from the outer edge.
    orientation : str, default 'F,R'
        Comma-separated orientation codes (F=Forward/+, R=Reverse/-).
    pairing_map : list of tuple, optional
        List of (left_feature, right_feature) tuples for asymmetric pairing.

    Returns
    -------
    pandas.DataFrame
        Filtered hit table.
    """
    if hit_table.empty:
        return hit_table

    # Parse orientation
    orientation_parts = orientation.upper().split(',')
    if len(orientation_parts) != 2:
        logging.warning(
            f'Invalid orientation "{orientation}"; expected two comma-separated codes. '
            'Skipping anchor filter.'
        )
        return hit_table

    left_strand = '+' if orientation_parts[0] == 'F' else '-'
    right_strand = '+' if orientation_parts[1] == 'F' else '-'
    strands_differ = left_strand != right_strand

    # Build model-to-terminus map from pairing map
    model_terminus: Dict[str, str] = {}
    if pairing_map:
        for left_feature, right_feature in pairing_map:
            model_terminus[left_feature] = 'left'
            model_terminus[right_feature] = 'right'

    kept: List[bool] = []
    skipped_no_terminus = 0
    removed = 0
    removed_per_model: Dict[str, int] = {}

    for _, row in hit_table.iterrows():
        model = row['model']
        strand = row['strand']

        model_len = model_lengths.get(model)
        if model_len is None:
            logging.warning(
                f'Model length not found for {model}, keeping hit without anchor check'
            )
            kept.append(True)
            continue

        try:
            hmm_start = int(row['hmmStart'])
            hmm_end = int(row['hmmEnd'])
        except (ValueError, TypeError):
            kept.append(True)
            continue

        # Determine terminus type
        if model in model_terminus:
            # Asymmetric: model name determines terminus type
            terminus_type: Optional[str] = model_terminus[model]
        elif strands_differ:
            # Symmetric with different strands: use strand to distinguish
            if strand == left_strand:
                terminus_type = 'left'
            elif strand == right_strand:
                terminus_type = 'right'
            else:
                terminus_type = None
        else:
            # Same-strand symmetric (F,F or R,R) without pairing map:
            # check both ends of the query model
            offset_start = hmm_start - 1
            offset_end = model_len - hmm_end
            if offset_start <= max_offset and offset_end <= max_offset:
                kept.append(True)
            else:
                kept.append(False)
                removed += 1
                removed_per_model[model] = removed_per_model.get(model, 0) + 1
            continue

        if terminus_type is None:
            skipped_no_terminus += 1
            kept.append(True)
            continue

        offset = compute_outer_edge_offset(
            hmm_start, hmm_end, model_len, strand, terminus_type
        )

        if offset <= max_offset:
            kept.append(True)
        else:
            kept.append(False)
            removed += 1
            removed_per_model[model] = removed_per_model.get(model, 0) + 1

    if skipped_no_terminus:
        logging.debug(
            f'Anchor filter: {skipped_no_terminus} hit(s) kept without checking – '
            'terminus type could not be determined.'
        )

    result = hit_table[kept].copy()

    logging.info(
        f'Anchor filter (max_offset={max_offset}): {len(hit_table)} -> {len(result)} hits '
        f'({removed} removed)'
    )
    if removed_per_model:
        for model_name, count in sorted(removed_per_model.items()):
            logging.info(f'  {model_name}: {count} hit(s) excluded by anchor filter')
    return result


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
        '--nhmmer-file',
        type=str,
        help='Path to single nhmmer output file (requires --hmm-file or --lengths-file).',
    )
    search_group.add_argument(
        '--left-nhmmer',
        type=str,
        help='Path to nhmmer output for left model (use with --right-nhmmer).',
    )
    search_group.add_argument(
        '--blast-file',
        type=str,
        help='Path to single BLAST tabular output file (requires --query-len or --lengths-file).',
    )
    search_group.add_argument(
        '--left-blast',
        type=str,
        help='Path to BLAST output for left query (use with --right-blast).',
    )

    parser.add_argument(
        '--right-nhmmer',
        type=str,
        help='Path to nhmmer output for right model (use with --left-nhmmer).',
    )

    parser.add_argument(
        '--right-blast',
        type=str,
        help='Path to BLAST output for right query (use with --left-blast).',
    )

    # Model/Query length sources (mutually exclusive)
    length_group = parser.add_mutually_exclusive_group()
    length_group.add_argument(
        '--hmm-file',
        type=str,
        help='Path to HMM file for extracting model lengths (for single model pairing).',
    )
    length_group.add_argument(
        '--left-model',
        type=str,
        help='Path to left HMM model file (for asymmetric pairing).',
    )
    length_group.add_argument(
        '--lengths-file',
        type=str,
        help='Path to tab-delimited file with model_name and model_length columns.',
    )
    length_group.add_argument(
        '--query-len',
        type=int,
        help='Length of BLAST query sequence (for single query pairing).',
    )

    parser.add_argument(
        '--right-model',
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

    parser.add_argument(
        '--max-offset',
        type=int,
        default=None,
        dest='max_offset',
        help=(
            'Maximum allowed offset (in bases) between the hit alignment boundary and '
            'the outer edge of the query model. Hits that do not reach within this '
            'many bases of the external edge are removed before pairing. '
            'For asymmetric queries, the outer edge is determined by the terminus type '
            '(left or right model). For symmetric same-strand queries (F,F or R,R), '
            'the hit must be within this offset of both ends of the query model. '
            'Default: no restriction.'
        ),
    )

    # Pairing configuration
    parser.add_argument(
        '--orientation',
        type=str,
        default='F,R',
        help='Orientation pattern for pairing. F=Forward(+), R=Reverse(-). Options: F,R (TIR), F,F (LTR), R,R, R,F. Default: F,R',
    )

    parser.add_argument(
        '--stable-reps',
        type=int,
        default=0,
        help='Number of iterations when no new pairs found. Default: 0',
    )

    parser.add_argument(
        '--pairing-map',
        type=str,
        default=None,
        dest='pairing_map',
        help=(
            'Tab-delimited file mapping left to right feature names for pairing. '
            'Each row: left_model<TAB>right_model. '
            'When provided, independent pairing is performed for each pair; '
            'models not listed in the map are skipped. '
            'Required when either input file contains hits to multiple query/model names.'
        ),
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
        '--gff-out',
        action='store_true',
        default=False,
        help='Generate GFF3 output file.',
    )

    parser.add_argument(
        '--gff-report',
        default='all',
        dest='gff_report',
        choices=['all', 'paired', 'unpaired'],
        help='Types of hits to include in GFF output.',
    )

    parser.add_argument(
        '--padlen',
        type=int,
        default=None,
        help='Extract N bases flanking each hit in FASTA output.',
    )

    parser.add_argument(
        '--flanks',
        action='store_true',
        default=False,
        help=(
            'Enable writing of external flanking sequences for all hits. '
            'Flanks are extracted using the length set by --flank-len (default 50). '
            'For symmetric same-strand orientations (F,F or R,R) both left and right '
            'flanks are written to separate files and a warning is raised advising '
            'use of --flanks-paired instead.'
        ),
    )

    parser.add_argument(
        '--flanks-paired',
        action='store_true',
        default=False,
        dest='flanks_paired',
        help=(
            'Write outer flanking sequences for termini that were assigned to pairs. '
            'If --flanks is also set, all-hit flanks are written first and paired '
            'flanks are written to separate files (reusing already-extracted sequences). '
            'If --flanks is not set, only paired flanks are written.'
        ),
    )

    parser.add_argument(
        '--flank-len',
        type=int,
        default=50,
        dest='flank_len',
        help=(
            'Length of the external flanking region to extract (in bases). '
            'The flank is the genomic sequence immediately outside each terminus hit '
            '(upstream of the left terminus, downstream of the right terminus). '
            'Default: 50.'
        ),
    )

    parser.add_argument(
        '--flank-max-offset',
        type=int,
        default=None,
        dest='flank_max_offset',
        help=(
            'Maximum allowed offset (in bases) between the hit alignment end and '
            'the external end of the query model. When a hit does not reach position 1 '
            'of the query, the flank start is corrected by this offset. '
            'Hits with offset greater than this value are skipped. '
            'Default: no limit.'
        ),
    )

    # Target site reconstruction options
    parser.add_argument(
        '--insertion-site',
        action='store_true',
        default=False,
        dest='insertion_site',
        help=(
            'Enable insertion site reconstruction and reporting. '
            'When set, flanking sequences are used to reconstruct the original '
            'target site for each paired element. '
            'Requires --flanks or --flanks-paired to be set. '
            'Use --tsd-length / --tsd-length-map / --tsd-in-model to configure '
            'TSD handling. '
            'Default: off.'
        ),
    )

    parser.add_argument(
        '--tsd-length',
        type=int,
        default=0,
        dest='tsd_length',
        help=(
            'Length of the Target Site Duplication (TSD) or Direct Repeat (DR) '
            'feature to account for when reconstructing target sites. '
            'Requires --insertion-site and --flank-len to be set. Default: 0 (no TSD).'
        ),
    )

    parser.add_argument(
        '--tsd-length-map',
        type=str,
        default=None,
        dest='tsd_length_map',
        help=(
            'Path to tab-delimited file mapping model pairs to TSD lengths. '
            'Format: left_model<TAB>right_model<TAB>tsd_length. '
            'Use when processing multiple model pairs with different TSD lengths.'
        ),
    )

    parser.add_argument(
        '--tsd-in-model',
        action='store_true',
        default=False,
        dest='tsd_in_model',
        help=(
            'If set, the TSD/DR feature is part of the termini model. '
            'If not set (default), the TSD/DR occurs in the flanking region '
            'immediately adjacent to the terminus hits.'
        ),
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
    if args.left_nhmmer and not args.right_nhmmer:
        raise ValueError('--left-nhmmer requires --right-nhmmer')
    if args.right_nhmmer and not args.left_nhmmer:
        raise ValueError('--right-nhmmer requires --left-nhmmer')
    if args.left_blast and not args.right_blast:
        raise ValueError('--left-blast requires --right-blast')
    if args.right_blast and not args.left_blast:
        raise ValueError('--right-blast requires --left-blast')
    if args.left_model and not args.right_model:
        raise ValueError('--left-model requires --right-model')
    if args.right_model and not args.left_model:
        raise ValueError('--right-model requires --left-model')

    # Check model/query length source requirements
    if args.nhmmer_file:
        if not (args.hmm_file or args.lengths_file):
            raise ValueError(
                '--nhmmer-file requires either --hmm-file or --lengths-file'
            )

    if args.blast_file:
        if not (args.query_len or args.lengths_file):
            raise ValueError(
                '--blast-file requires either --query-len or --lengths-file'
            )

    if args.left_nhmmer and args.right_nhmmer:
        if not (args.left_model and args.right_model) and not args.lengths_file:
            raise ValueError(
                'Asymmetric pairing requires --left-model/--right-model or --lengths-file'
            )

    if args.left_blast and args.right_blast:
        if not args.lengths_file:
            raise ValueError(
                'Asymmetric BLAST pairing requires --lengths-file with query lengths'
            )

    # Validate target site reconstruction options
    tsd_options_set = args.tsd_length > 0 or args.tsd_length_map or args.tsd_in_model
    if tsd_options_set and not args.insertion_site:
        logging.warning(
            'TSD options (--tsd-length, --tsd-length-map, --tsd-in-model) are set '
            'but --insertion-site is not enabled. These options will be ignored. '
            'Add --insertion-site to enable insertion site reconstruction.'
        )
    if args.insertion_site and not (args.flanks or args.flanks_paired):
        raise ValueError(
            '--insertion-site requires --flanks or --flanks-paired to be set'
        )

    # Check file existence
    required_files = []
    if args.genome:
        required_files.append(args.genome)

    if args.nhmmer_file:
        required_files.append(args.nhmmer_file)
    if args.left_nhmmer:
        required_files.extend([args.left_nhmmer, args.right_nhmmer])
    if args.blast_file:
        required_files.append(args.blast_file)
    if args.left_blast:
        required_files.extend([args.left_blast, args.right_blast])
    if args.hmm_file:
        required_files.append(args.hmm_file)
    if args.left_model:
        required_files.extend([args.left_model, args.right_model])
    if args.lengths_file:
        required_files.append(args.lengths_file)
    if args.pairing_map:
        required_files.append(args.pairing_map)
    if args.tsd_length_map:
        required_files.append(args.tsd_length_map)

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


def _write_pair_summary(
    outdir: str,
    prefix: Optional[str],
    left_feature: str,
    right_feature: str,
    hitTable: Any,
    pair_paired: Dict[str, List[Any]],
    pair_hitIndex: Dict[str, Dict[int, Dict[str, Any]]],
    total_pairs: int,
    total_elements: int,
    filter_stats: Optional[Dict[str, Any]] = None,
) -> None:
    """
    Write a summary report for a single model pairing.

    Parameters
    ----------
    outdir : str
        Output directory for the summary file.
    prefix : str or None
        Prefix for output filename.
    left_feature : str
        Name of the left model/feature.
    right_feature : str
        Name of the right model/feature.
    hitTable : pandas.DataFrame
        DataFrame of all hits.
    pair_paired : dict
        Paired hits for this model pair.
    pair_hitIndex : dict
        Hit index for this model pair.
    total_pairs : int
        Total number of pairs found.
    total_elements : int
        Total number of elements extracted.
    filter_stats : dict, optional
        Dictionary with hit-filtering statistics to include in the report.
        Keys used: 'initial_hits', 'mincov', 'coverage_excluded',
        'maxeval', 'evalue_excluded', 'max_offset', 'anchor_excluded',
        'pairing_map_models_ignored', 'pairing_map_hits_ignored'.

    Returns
    -------
    None
        Writes summary text file to disk.
    """
    prefix_str = f'{prefix}_' if prefix else ''
    pair_label = (
        f'{left_feature}_{right_feature}'
        if left_feature != right_feature
        else left_feature
    )
    summary_path = os.path.join(outdir, f'{prefix_str}{pair_label}_summary.txt')

    # Count hits per model in this pairing
    models = set()
    models.add(left_feature)
    models.add(right_feature)

    hits_per_model: Dict[str, int] = {}
    for model in models:
        if model in pair_hitIndex:
            hits_per_model[model] = len(pair_hitIndex[model])
        else:
            hits_per_model[model] = 0

    # Count paired hits
    paired_hit_ids: set = set()
    for _model, pairs in pair_paired.items():
        for pair_set in pairs:
            paired_hit_ids.update(pair_set)

    total_hits = sum(hits_per_model.values())
    total_unpaired = total_hits - len(paired_hit_ids)

    with open(summary_path, 'w') as f:
        f.write('TIRmite Pair Summary Report\n')
        f.write('===========================\n\n')
        f.write(f'Model pair: {left_feature} <-> {right_feature}\n\n')

        # --- Filtering criteria section ---
        if filter_stats:
            f.write('Filtering criteria applied\n')
            f.write('--------------------------\n')
            initial = filter_stats.get('initial_hits')
            if initial is not None:
                f.write(f'  Initial hits imported: {initial}\n')

            # Pairing-map pre-filter
            ignored_models = filter_stats.get('pairing_map_models_ignored')
            if ignored_models:
                f.write(
                    f'  Pairing-map model filter: {len(ignored_models)} model(s) ignored '
                    f'({", ".join(sorted(ignored_models))}), '
                    f'{filter_stats.get("pairing_map_hits_ignored", 0)} hits excluded\n'
                )

            # Coverage filter
            mincov = filter_stats.get('mincov')
            cov_excluded = filter_stats.get('coverage_excluded', 0)
            if mincov is not None:
                f.write(
                    f'  Coverage filter (min coverage >= {mincov}): '
                    f'{cov_excluded} hit(s) excluded\n'
                )

            # E-value filter
            maxeval = filter_stats.get('maxeval')
            eval_excluded = filter_stats.get('evalue_excluded', 0)
            if maxeval is not None:
                f.write(
                    f'  E-value filter (max e-value <= {maxeval}): '
                    f'{eval_excluded} hit(s) excluded\n'
                )

            # Anchor offset filter
            max_offset = filter_stats.get('max_offset')
            anchor_excluded = filter_stats.get('anchor_excluded')
            if max_offset is not None:
                excl_str = (
                    str(anchor_excluded) if anchor_excluded is not None else 'unknown'
                )
                f.write(
                    f'  Anchor offset filter (max offset <= {max_offset}): '
                    f'{excl_str} hit(s) excluded\n'
                )

            after_filtering = filter_stats.get('after_filtering')
            if after_filtering is not None:
                f.write(f'  Hits remaining after all filters: {after_filtering}\n')
            f.write('\n')

        # --- Pairing results section ---
        f.write('Hits per model (after all filters)\n')
        f.write('----------------------------------\n')
        for model, count in sorted(hits_per_model.items()):
            f.write(f'  {model}: {count}\n')
        f.write(f'\nTotal hits for this pair: {total_hits}\n')
        f.write(f'Total pairs found: {total_pairs}\n')
        f.write(f'Total paired elements extracted: {total_elements}\n')
        f.write(f'Total unpaired hits: {total_unpaired}\n')

    logging.info(f'Wrote summary report to {summary_path}')


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

        if args.lengths_file:
            model_lengths = load_model_lengths_file(args.lengths_file)
        elif args.hmm_file:
            model_lengths = get_hmm_model_length(args.hmm_file)
        elif args.query_len:
            # For single BLAST query, we'll assign the length after importing hits
            # to get the query name
            pass
        elif args.left_model and args.right_model:
            left_lengths = get_hmm_model_length(args.left_model)
            right_lengths = get_hmm_model_length(args.right_model)
            model_lengths = {**left_lengths, **right_lengths}

        # Import search hits
        hitTable = None
        input_format = None

        if args.nhmmer_file:
            # Single nhmmer file mode
            logging.info('Importing nhmmer hits...')
            input_format = 'nhmmer'
            hitTable = tirmite.import_nhmmer(
                infile=args.nhmmer_file, hitTable=None, prefix=args.prefix
            )
        elif args.left_nhmmer:
            # Asymmetric nhmmer mode - import from both files
            logging.info('Importing nhmmer hits from left and right models...')
            input_format = 'nhmmer'

            # Check if left and right files are the same
            if args.left_nhmmer == args.right_nhmmer:
                raise ValueError(
                    f'Left and right nhmmer files cannot be the same: {args.left_nhmmer}'
                )

            left_model_name = extract_model_name_from_path(args.left_model)
            right_model_name = extract_model_name_from_path(args.right_model)

            # Import left file
            left_hitTable = tirmite.import_nhmmer(
                infile=args.left_nhmmer,
                hitTable=None,
                prefix=args.prefix,
            )
            left_models = check_multiple_models(left_hitTable)
            logging.info(
                f'Left nhmmer file: {len(left_hitTable)} hits, {len(left_models)} unique query/model name(s)'
            )

            # Import right file
            right_hitTable = tirmite.import_nhmmer(
                infile=args.right_nhmmer,
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

            # Validate single query per file or require --pairing-map
            if len(left_models) > 1 or len(right_models) > 1:
                if len(left_models) > 1:
                    logging.warning(
                        f'Left nhmmer file contains {len(left_models)} query/model names: '
                        + ', '.join(sorted(left_models))
                    )
                if len(right_models) > 1:
                    logging.warning(
                        f'Right nhmmer file contains {len(right_models)} query/model names: '
                        + ', '.join(sorted(right_models))
                    )
                if not args.pairing_map:
                    raise ValueError(
                        f'Left file contains {len(left_models)} query/model name(s), '
                        f'right file contains {len(right_models)} query/model name(s). '
                        'When either file contains multiple queries, --pairing-map is required.'
                    )

            # Combine hit tables
            hitTable = tirmite.import_nhmmer(
                infile=args.left_nhmmer,
                hitTable=None,
                prefix=args.prefix,
            )
            hitTable = tirmite.import_nhmmer(
                infile=args.right_nhmmer,
                hitTable=hitTable,
                prefix=args.prefix,
            )
        elif args.blast_file:
            # Single BLAST file mode
            logging.info('Importing BLAST hits...')
            input_format = 'blast'

            # Detect format and warn if mismatch
            detected_format = tirmite.detect_input_format(args.blast_file)
            if detected_format != 'blast' and detected_format != 'unknown':
                logging.warning(
                    f'File format appears to be {detected_format}, but --blast-file was specified. '
                    'Consider using --nhmmer-file instead.'
                )

            hitTable = tirmite.import_blast(
                infile=args.blast_file, hitTable=None, prefix=args.prefix
            )

        elif args.left_blast:
            # Asymmetric BLAST mode - import from both files
            logging.info('Importing BLAST hits from left and right queries...')
            input_format = 'blast'

            # Check if left and right files are the same
            if args.left_blast == args.right_blast:
                raise ValueError(
                    f'Left and right BLAST files cannot be the same: {args.left_blast}'
                )

            # Detect format for both files
            detected_left = tirmite.detect_input_format(args.left_blast)
            detected_right = tirmite.detect_input_format(args.right_blast)

            if detected_left != 'blast' and detected_left != 'unknown':
                logging.warning(
                    f'Left file format appears to be {detected_left}, but --left-blast was specified.'
                )
            if detected_right != 'blast' and detected_right != 'unknown':
                logging.warning(
                    f'Right file format appears to be {detected_right}, but --right-blast was specified.'
                )

            # Import left file
            left_hitTable = tirmite.import_blast(
                infile=args.left_blast,
                hitTable=None,
                prefix=args.prefix,
            )
            left_models = check_multiple_models(left_hitTable)
            logging.info(
                f'Left BLAST file: {len(left_hitTable)} hits, {len(left_models)} unique query/model name(s)'
            )

            # Import right file
            right_hitTable = tirmite.import_blast(
                infile=args.right_blast,
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

            # Validate single query per file or require --pairing-map
            if len(left_models) > 1 or len(right_models) > 1:
                if len(left_models) > 1:
                    logging.warning(
                        f'Left BLAST file contains {len(left_models)} query/model names: '
                        + ', '.join(sorted(left_models))
                    )
                if len(right_models) > 1:
                    logging.warning(
                        f'Right BLAST file contains {len(right_models)} query/model names: '
                        + ', '.join(sorted(right_models))
                    )
                if not args.pairing_map:
                    raise ValueError(
                        f'Left file contains {len(left_models)} query/model name(s), '
                        f'right file contains {len(right_models)} query/model name(s). '
                        'When either file contains multiple queries, --pairing-map is required.'
                    )

            # Combine hit tables
            hitTable = tirmite.import_blast(
                infile=args.left_blast,
                hitTable=None,
                prefix=args.prefix,
            )
            hitTable = tirmite.import_blast(
                infile=args.right_blast,
                hitTable=hitTable,
                prefix=args.prefix,
            )

        if hitTable is None or len(hitTable) == 0:
            logging.error('No hits were imported from input files')
            cleanup_temp_directory(tempDir, args.keep_temp)
            sys.exit(1)

        # Log import statistics
        initial_hit_count = len(hitTable)
        logging.info(f'Imported {initial_hit_count} total hits')

        # Get unique models and log statistics
        unique_models = check_multiple_models(hitTable)
        logging.info(f'Found {len(unique_models)} unique query/model name(s)')

        # Log per-query hit counts at debug level
        for model in unique_models:
            hit_count = len(hitTable[hitTable['model'] == model])
            logging.debug(f'Query/model "{model}": {hit_count} hits')

        # Filter stats dict – populated incrementally and passed to summary reports
        filter_stats: Dict[str, Any] = {
            'initial_hits': initial_hit_count,
            'mincov': args.mincov,
            'coverage_excluded': 0,
            'maxeval': args.maxeval,
            'evalue_excluded': 0,
            'max_offset': args.max_offset if args.max_offset is not None else None,
            'anchor_excluded': None,
            'pairing_map_models_ignored': None,
            'pairing_map_hits_ignored': 0,
            'after_filtering': initial_hit_count,
        }

        # Filter hitTable to only include models present in the pairing map (if provided)
        # This must happen before any processing to avoid wasting work on unused models
        if args.pairing_map:
            _pre_filter_pairing_map = load_pairing_map(args.pairing_map)
            # Collect all unique model names referenced in the pairing map
            map_models: set = set()
            for _lf, _rf in _pre_filter_pairing_map:
                map_models.add(_lf)
                map_models.add(_rf)

            # Determine which loaded models are not in the pairing map
            models_in_hits = set(unique_models)
            models_to_ignore = models_in_hits - map_models
            models_to_process = models_in_hits & map_models

            if models_to_ignore:
                ignored_hits = len(hitTable[hitTable['model'].isin(models_to_ignore)])
                logging.warning(
                    f'Pairing map references {len(map_models)} model(s). '
                    f'Ignoring {len(models_to_ignore)} model(s) not in pairing map: '
                    + ', '.join(sorted(models_to_ignore))
                )
                logging.info(
                    f'Excluding {ignored_hits} hits for {len(models_to_ignore)} ignored model(s)'
                )
                hitTable = hitTable[hitTable['model'].isin(models_to_process)].copy()
                logging.info(
                    f'{len(hitTable)} hits remaining for {len(models_to_process)} model(s) in pairing map'
                )
                if len(hitTable) == 0:
                    logging.error('No hits remaining after pairing map model filter')
                    cleanup_temp_directory(tempDir, args.keep_temp)
                    sys.exit(1)
                # Update unique_models after filtering
                unique_models = check_multiple_models(hitTable)
                filter_stats['pairing_map_models_ignored'] = sorted(models_to_ignore)
                filter_stats['pairing_map_hits_ignored'] = ignored_hits
            else:
                logging.info(
                    f'All {len(models_in_hits)} model(s) in hits are covered by pairing map'
                )

        # If --query-len was provided for BLAST input, assign it to ALL queries
        if args.blast_file and args.query_len:
            for query_name in unique_models:
                model_lengths[query_name] = args.query_len
                logging.debug(f'Set length for query {query_name}: {args.query_len}')
            logging.info(
                f'Applied query length {args.query_len} to {len(unique_models)} query name(s)'
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
        coverage_excluded = hitCount - len(hitTable)
        logging.info(f'Excluded {coverage_excluded} hits on coverage criteria')
        filter_stats['coverage_excluded'] = coverage_excluded

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
        evalue_excluded = hitCount - len(hitTable)
        logging.info(f'Excluded {evalue_excluded} hits on e-value criteria')
        filter_stats['evalue_excluded'] = evalue_excluded

        for model in model_counts_before:
            before = model_counts_before[model]
            after = model_counts_after.get(model, 0)
            excluded = before - after
            logging.info(f'Model {model}: {excluded} excluded, {after} remaining')

        logging.info(f'Total remaining hits: {len(hitTable)}')

        # Apply anchor filter if --max-offset is set
        if args.max_offset is not None:
            logging.info(
                f'Filtering hits by anchor offset (max_offset={args.max_offset})...'
            )
            hitCount_before_anchor = len(hitTable)

            # Build pairing map for asymmetric model identification
            anchor_pairing_map = None
            if args.pairing_map:
                anchor_pairing_map = load_pairing_map(args.pairing_map)

            hitTable = filter_hits_by_anchor(
                hit_table=hitTable,
                model_lengths=model_lengths,
                max_offset=args.max_offset,
                orientation=args.orientation,
                pairing_map=anchor_pairing_map,
            )

            filter_stats['anchor_excluded'] = hitCount_before_anchor - len(hitTable)

            if len(hitTable) == 0:
                logging.error('No hits remaining after anchor filter')
                cleanup_temp_directory(tempDir, args.keep_temp)
                sys.exit(1)

        # Finalise after_filtering count after all filter steps
        filter_stats['after_filtering'] = len(hitTable)

        # Convert to dict structure
        hitsDict, hitIndex = tirmite.table2dict(hitTable)

        # Check for multiple models and validate pairing map requirement
        # Note: unique_models was already determined and logged after import
        # For asymmetric modes, validation was already done per-file
        is_asymmetric = (args.left_nhmmer and args.right_nhmmer) or (
            args.left_blast and args.right_blast
        )

        # Load pairing map if provided
        pairing_map = None
        if args.pairing_map:
            pairing_map = load_pairing_map(args.pairing_map)
        elif len(unique_models) > 1 and not is_asymmetric:
            # Multiple models in single file without pairing map - warn and raise error
            # (Asymmetric mode already validated per-file)
            logging.warning(
                f'Input contains {len(unique_models)} distinct models/queries: '
                + ', '.join(sorted(unique_models))
            )
            raise ValueError(
                f'Input contains {len(unique_models)} distinct models/queries: {", ".join(sorted(unique_models))}. '
                'Multiple models require --pairing-map to specify which features should be paired together.'
            )

        # Create pairing configuration
        if pairing_map:
            # Use pairing map - will create configs for each pairing later
            # Pairing map workflow handles both single and multiple pairings
            logging.info(
                f'Will execute {len(pairing_map)} independent pairing procedure(s) based on pairing map'
            )
            config = None
        elif args.left_nhmmer and args.right_nhmmer:
            # Asymmetric nhmmer pairing
            left_model_name = extract_model_name_from_path(args.left_model)
            right_model_name = extract_model_name_from_path(args.right_model)

            config = tirmite.PairingConfig(
                orientation=args.orientation,
                left_model=left_model_name,
                right_model=right_model_name,
            )
        elif args.left_blast and args.right_blast:
            # Asymmetric BLAST pairing - extract model names from each file's hitTable
            # Use left_hitTable and right_hitTable (imported earlier) to ensure the
            # left model comes from --left-blast and right model from --right-blast.
            # Do NOT use hitTable['model'].unique() here because the combined table is
            # sorted alphabetically, which would assign models based on query name order
            # rather than which file they came from.
            left_model_names = left_hitTable['model'].unique()
            right_model_names = right_hitTable['model'].unique()

            if len(left_model_names) == 0 or len(right_model_names) == 0:
                logging.error(
                    'Expected at least 1 model in each BLAST file for asymmetric pairing'
                )
                cleanup_temp_directory(tempDir, args.keep_temp)
                sys.exit(1)

            left_model_name = left_model_names[0]
            right_model_name = right_model_names[0]
            logging.info(
                f'Assigning models for asymmetric pairing: '
                f'left={left_model_name} (from --left-blast), '
                f'right={right_model_name} (from --right-blast)'
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
        # In pairing_map mode individual hits are written per pair (below),
        # so only write to the base outDir when no pairing map is used.
        if not pairing_map:
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
                            pair_hitIndex, pair_config, stableReps=args.stable_reps
                        )
                    )
                else:
                    logging.info(
                        f'Using symmetric pairing with orientation {pair_config.orientation}'
                    )
                    pair_hitIndex, pair_paired, pair_unpaired = (
                        tirmite.iterateGetPairsCustom(
                            pair_hitIndex, pair_config, stableReps=args.stable_reps
                        )
                    )

                # Create sub-directory for this pair
                pair_label = (
                    f'{left_feature}_{right_feature}'
                    if left_feature != right_feature
                    else left_feature
                )
                pair_outDir = os.path.join(outDir, pair_label)
                os.makedirs(pair_outDir, exist_ok=True)

                # Write per-pair output to sub-directory
                # Write individual hits for the models in this pair
                pair_hit_models = {left_feature, right_feature}
                pair_hitTable_tirs = hitTable[hitTable['model'].isin(pair_hit_models)]
                logging.info(f'Writing individual hits for pair {pair_label}...')
                tirmite.writeTIRs(
                    outDir=pair_outDir,
                    hitTable=pair_hitTable_tirs,
                    maxeval=args.maxeval,
                    genome=genome,
                    prefix=args.prefix,
                    padlen=args.padlen,
                    genome_descriptions=genome_descriptions,
                    blastdb=args.blastdb if args.blastdb else None,
                )

                # Write paired TIRs
                if args.gff_report in ['all', 'paired']:
                    tirmite.writePairedTIRs(
                        outDir=pair_outDir,
                        paired=pair_paired,
                        hitIndex=pair_hitIndex,
                        genome=genome,
                        prefix=args.prefix,
                        padlen=args.padlen,
                        genome_descriptions=genome_descriptions,
                        blastdb=args.blastdb if args.blastdb else None,
                    )

                # Extract and write elements for this pair
                pair_pairedEles = tirmite.fetchElements(
                    paired=pair_paired,
                    hitIndex=pair_hitIndex,
                    genome=genome,
                    genome_descriptions=genome_descriptions,
                    blastdb=args.blastdb if args.blastdb else None,
                )
                tirmite.writeElements(
                    pair_outDir, eleDict=pair_pairedEles, prefix=args.prefix
                )

                # Extract and write flanks for this pair
                if args.flanks or args.flanks_paired:
                    tirmite.writeFlanks(
                        outDir=pair_outDir,
                        hitTable=hitTable,
                        model_lengths=model_lengths,
                        paired=pair_paired,
                        hitIndex=pair_hitIndex,
                        config=pair_config,
                        genome=genome,
                        prefix=args.prefix,
                        flank_len=args.flank_len,
                        flank_max_offset=args.flank_max_offset,
                        write_all=args.flanks,
                        write_paired=args.flanks_paired,
                        genome_descriptions=genome_descriptions,
                        blastdb=args.blastdb if args.blastdb else None,
                    )

                    # Reconstruct target sites for this pair
                    if args.insertion_site:
                        tsd_length_map_data = None
                        if args.tsd_length_map:
                            tsd_length_map_data = tirmite.load_tsd_length_map(
                                args.tsd_length_map
                            )

                        tirmite.writeTargetSites(
                            outDir=pair_outDir,
                            hitTable=hitTable,
                            model_lengths=model_lengths,
                            paired=pair_paired,
                            hitIndex=pair_hitIndex,
                            config=pair_config,
                            genome=genome,
                            prefix=args.prefix,
                            flank_len=args.flank_len,
                            flank_max_offset=args.flank_max_offset,
                            tsd_length=args.tsd_length,
                            tsd_in_model=args.tsd_in_model,
                            tsd_length_map=tsd_length_map_data,
                            genome_descriptions=genome_descriptions,
                            blastdb=args.blastdb if args.blastdb else None,
                        )

                # Write summary report for this pair
                total_pair_pairs = sum(len(pairs) for pairs in pair_paired.values())
                total_pair_elements = sum(
                    len(eles) for eles in pair_pairedEles.values()
                )
                _write_pair_summary(
                    outdir=pair_outDir,
                    prefix=args.prefix,
                    left_feature=left_feature,
                    right_feature=right_feature,
                    hitTable=hitTable,
                    pair_paired=pair_paired,
                    pair_hitIndex=pair_hitIndex,
                    total_pairs=total_pair_pairs,
                    total_elements=total_pair_elements,
                    filter_stats=filter_stats,
                )

                logging.info(f'Pair {pair_label}: wrote output to {pair_outDir}')

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
                    hitIndex, config, stableReps=args.stable_reps
                )
            else:
                logging.info(
                    f'Using symmetric pairing with orientation {config.orientation}'
                )
                hitIndex, paired, unpaired = tirmite.iterateGetPairsCustom(
                    hitIndex, config, stableReps=args.stable_reps
                )

            # Log pairing results
            total_pairs = sum(len(pairs) for pairs in paired.values())
            total_unpaired = len(unpaired)
            logging.info(
                f'Pairing completed: {total_pairs} pairs, {total_unpaired} unpaired'
            )

        # For pairing-map mode all output has already been written to per-pair
        # subdirectories inside the loop above.  The block below handles the
        # single-pairing (no pairing map) case only.
        if not pairing_map:
            # Write paired TIRs
            if args.gff_report in ['all', 'paired']:
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

            # Write summary report for single-pairing mode
            if config is not None:
                left_feature = config.left_model or ''
                right_feature = config.right_model or ''
                total_ele_count = sum(len(eles) for eles in pairedEles.values())
                _write_pair_summary(
                    outdir=outDir,
                    prefix=args.prefix,
                    left_feature=left_feature,
                    right_feature=right_feature,
                    hitTable=hitTable,
                    pair_paired=paired,
                    pair_hitIndex=hitIndex,
                    total_pairs=sum(len(p) for p in paired.values()),
                    total_elements=total_ele_count,
                    filter_stats=filter_stats,
                )

            # Extract and write external flanks if requested
            if args.flanks or args.flanks_paired:
                logging.info('Extracting external flanking sequences...')
                tirmite.writeFlanks(
                    outDir=outDir,
                    hitTable=hitTable,
                    model_lengths=model_lengths,
                    paired=paired,
                    hitIndex=hitIndex,
                    config=config,
                    genome=genome,
                    prefix=args.prefix,
                    flank_len=args.flank_len,
                    flank_max_offset=args.flank_max_offset,
                    write_all=args.flanks,
                    write_paired=args.flanks_paired,
                    genome_descriptions=genome_descriptions,
                    blastdb=args.blastdb if args.blastdb else None,
                )

                # Reconstruct target sites if explicitly requested
                if args.insertion_site:
                    logging.info('Reconstructing target sites...')

                    # Load TSD length map if provided
                    tsd_length_map = None
                    if args.tsd_length_map:
                        tsd_length_map = tirmite.load_tsd_length_map(
                            args.tsd_length_map
                        )

                    tirmite.writeTargetSites(
                        outDir=outDir,
                        hitTable=hitTable,
                        model_lengths=model_lengths,
                        paired=paired,
                        hitIndex=hitIndex,
                        config=config,
                        genome=genome,
                        prefix=args.prefix,
                        flank_len=args.flank_len,
                        flank_max_offset=args.flank_max_offset,
                        tsd_length=args.tsd_length,
                        tsd_in_model=args.tsd_in_model,
                        tsd_length_map=tsd_length_map,
                        genome_descriptions=genome_descriptions,
                        blastdb=args.blastdb if args.blastdb else None,
                    )

            # Write GFF if requested
            if args.gff_out:
                # Get unpaired TIRs if needed
                unpairedTIRs = None
                if args.gff_report in ['all', 'unpaired']:
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
                    writeTIRs=args.gff_report,
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
