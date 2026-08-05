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
import shutil
import sys
from typing import Any, Dict, List, Optional, Set, cast

from tirmite._version import __version__
from tirmite.core.extraction import (
    fetchElements,
    writeElements,
    writePairedTIRs,
    writeTIRs,
)
from tirmite.core.filters import (
    filterHitsEval,
)
from tirmite.core.flanks import (
    writeFlanks,
)

# Re-exported rather than imported for use only: tests and downstream code
# import these names from this module, and they were defined here before the
# two diverged copies were unified into tirmite.core.hit_filters.
from tirmite.core.hit_filters import (  # noqa: F401
    compute_outer_edge_offset,
    filter_hits_by_anchor,
)
from tirmite.core.output import (
    fetchUnpaired,
    gffWrite,
)
from tirmite.core.pairing import (
    PairingConfig,
    iterateGetPairsAsymmetric,
    iterateGetPairsCustom,
    parseHitsGeneral,
    table2dict,
)
from tirmite.core.parsers import (
    detect_input_format,
    import_blast,
    import_nhmmer,
)
from tirmite.core.tsd import (
    load_tsd_length_map,
    writeTargetSites,
)
from tirmite.report.stats import (
    PairSummary,
    format_pair_summary,
    pair_summary_stats,
)
from tirmite.utils.extract import check_ids, make_source
from tirmite.utils.logs import init_logging
from tirmite.utils.utils import (
    cleanup_temp_directory,
    extract_model_name_from_path,
    indexGenome,
    setup_directories,
)

logger = logging.getLogger(__name__)


def log_startup_info(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    """
    Log startup information: version, package location, full command, and non-default args.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    parser : argparse.ArgumentParser
        The argument parser (used to retrieve defaults for comparison).
    """
    import tirmite as _tirmite_pkg

    pkg_path = Path(_tirmite_pkg.__file__).parent
    logger.info(f'TIRmite package location: {pkg_path}')
    logger.info(f'TIRmite-pair version: {__version__}')
    logger.info(f'Command: {" ".join(sys.argv)}')

    # Collect args that differ from their defaults
    defaults = {action.dest: action.default for action in parser._actions}
    non_defaults = []
    for dest, default_val in defaults.items():
        current_val = getattr(args, dest, default_val)
        if current_val != default_val:
            flag = dest.replace('_', '-')
            non_defaults.append(f'  --{flag} {current_val}')

    if non_defaults:
        logger.info('Non-default arguments:\n' + '\n'.join(non_defaults))
    else:
        logger.info('All arguments are at their default values.')


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
        logger.error(f'Error reading HMM file {hmm_file_path}: {e}')

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
                    logger.warning(
                        f'Skipping malformed line {line_num} in {lengths_file}: {line}'
                    )
                    continue

                model_name, model_length_str = parts
                try:
                    model_length = int(model_length_str)
                    model_lengths[model_name] = model_length
                except ValueError:
                    logger.warning(
                        f'Invalid model length on line {line_num}: {model_length_str}'
                    )

    except Exception as e:
        logger.error(f'Error reading model lengths file {lengths_file}: {e}')

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
            logger.warning(f'Model length not found for {model}, using coverage = 0')
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
    pairings: list[tuple[str, str]] = []
    # Feature name -> line numbers it appeared on, for the duplicate warning.
    seen_features: dict[str, list[int]] = {}

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
            logger.warning(
                f'Feature "{feature}" appears in multiple pairing combinations '
                f'(lines: {", ".join(map(str, line_nums))})'
            )

    logger.info(f'Loaded {len(pairings)} pairing combinations from {pairing_map_file}')
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
    # .tolist() is untyped in pandas, so state the element type explicitly
    # rather than leaking Any out of a function declared to return list[str].
    unique_models: list[str] = list(hitTable['model'].unique())
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
        help=(
            'Maximum distance allowed between termini for pairing, measured '
            'between the facing inner edges of the two terminus hits. This is '
            'the length of the element interior and excludes the termini '
            'themselves, so it does not depend on model length or strand. '
            'Default: no limit.'
        ),
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
        '--no-hits',
        action='store_true',
        default=False,
        dest='no_hits',
        help=(
            'Skip writing individual hit sequences to FASTA. '
            'Useful when extracting sequences from a large genome or BLAST database '
            'would be time-consuming and hit sequences are not required.'
        ),
    )

    parser.add_argument(
        '--no-elements',
        action='store_true',
        default=False,
        dest='no_elements',
        help=(
            'Skip extraction and writing of paired element sequences to FASTA. '
            'Useful when extracting full-length elements from a large genome or '
            'BLAST database would be time-consuming and element sequences are not required.'
        ),
    )

    parser.add_argument(
        '--gff',
        action='store_true',
        default=False,
        dest='gff_out',
        help='Generate GFF3 output file.',
    )

    parser.add_argument(
        '--gff-report',
        default='all',
        dest='gff_report',
        choices=['all', 'paired', 'unpaired'],
        help='Types of hits to include in GFF output.',
    )

    # Report options
    parser.add_argument(
        '--report',
        action='store_true',
        default=False,
        help=(
            'Write a self-contained HTML report of the run: annotation tracks '
            'for every sequence with a hit, terminus alignments, distribution '
            'plots and summary statistics.'
        ),
    )

    parser.add_argument(
        '--report-out',
        type=str,
        default=None,
        dest='report_out',
        help=(
            'Path for the HTML report. Implies --report. '
            'Default: <outdir>/<prefix>tirmite_pair_report.html'
        ),
    )

    parser.add_argument(
        '--report-title',
        type=str,
        default=None,
        dest='report_title',
        help='Heading shown at the top of the HTML report.',
    )

    parser.add_argument(
        '--no-report-sequences',
        action='store_true',
        default=False,
        dest='report_no_sequences',
        help=(
            'Do not embed element sequences in the HTML report. The report '
            'stays much smaller, but its copy-to-clipboard buttons only show '
            'coordinates.'
        ),
    )

    parser.add_argument(
        '--report-max-seq-mb',
        type=float,
        default=20.0,
        dest='report_max_seq_mb',
        help=(
            'Budget in megabytes for element sequences embedded in the HTML '
            'report. Shorter elements are embedded first. Default: 20'
        ),
    )

    parser.add_argument(
        '--report-msa',
        default='auto',
        dest='report_msa',
        choices=['auto', 'mafft', 'anchor', 'off'],
        help=(
            'How to build the terminus alignment panels. "auto" uses MAFFT '
            'when it is on PATH and falls back to placing hits by their model '
            'coordinates; "mafft" requires MAFFT; "anchor" never calls it; '
            '"off" omits the panels. Default: auto'
        ),
    )

    parser.add_argument(
        '--report-msa-max-rows',
        type=int,
        default=500,
        dest='report_msa_max_rows',
        help=(
            'Maximum hits shown per terminus alignment panel, most '
            'significant first. Default: 500'
        ),
    )

    parser.add_argument(
        '--report-max-hits',
        type=int,
        default=200000,
        dest='report_max_hits',
        help=(
            'Maximum hits included in the HTML report. Paired hits and the '
            'most significant unpaired hits are kept. Default: 200000'
        ),
    )

    parser.add_argument(
        '--report-max-rows',
        type=int,
        default=30,
        dest='report_max_rows',
        help=(
            'Maximum stacked annotation rows per sequence in the HTML report. '
            'Default: 30'
        ),
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

    parser.add_argument(
        '--no-pad-flanks',
        action='store_true',
        default=False,
        dest='no_pad_flanks',
        help=(
            'Do not pad flanking and target-site sequences that extend past a '
            'contig boundary. By default such regions are padded with N so that '
            'every flank is --flank-len bases and records remain comparable '
            'position by position; padded records are marked in their FASTA '
            'description. With this flag they are truncated instead, so record '
            'lengths vary near contig ends. Regions that fall entirely outside '
            'a contig are skipped either way.'
        ),
    )

    parser.add_argument(
        '--extend-hits-to-model',
        action='store_true',
        default=False,
        dest='extend_hits_to_model',
        help=(
            'Extend extracted terminus sequences outward by the offset between '
            'the hit alignment and the external end of the model, so that hits '
            'which only partially cover the model are emitted at full model '
            'length and are directly comparable. Affects the hit and paired '
            'terminus FASTA output only, not flanks or pairing. '
            'Default: extract only the aligned region.'
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
        logger.warning(
            'TSD options (--tsd-length, --tsd-length-map, --tsd-in-model) are set '
            'but --insertion-site is not enabled. These options will be ignored. '
            'Add --insertion-site to enable insertion site reconstruction.'
        )
    if args.insertion_site and not (args.flanks or args.flanks_paired):
        raise ValueError(
            '--insertion-site requires --flanks or --flanks-paired to be set'
        )

    _validate_report_arguments(args)

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


def _validate_report_arguments(args: argparse.Namespace) -> None:
    """
    Validate and normalise the HTML report options.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments. Mutated in place: `--report-out` implies `--report`,
        and options that cannot be honoured are downgraded with a warning.

    Returns
    -------
    None
        Updates `args` in place.

    Raises
    ------
    ValueError
        If a report option was given a value that cannot be honoured.

    Notes
    -----
    Conflicts here downgrade rather than abort wherever the run can still
    produce a useful report. The report is a secondary output; refusing to run
    a whole pairing analysis because a visualisation option is redundant would
    be the wrong trade.
    """
    # Tolerate namespaces built before these options existed, as the tests and
    # any library caller may do.
    if not hasattr(args, 'report'):
        return

    if getattr(args, 'report_out', None):
        args.report = True

    if not args.report:
        return

    if args.report_max_seq_mb <= 0:
        raise ValueError('--report-max-seq-mb must be greater than 0')
    if args.report_max_hits < 1:
        raise ValueError('--report-max-hits must be at least 1')
    if args.report_max_rows < 1:
        raise ValueError('--report-max-rows must be at least 1')
    if args.report_msa_max_rows < 1:
        raise ValueError('--report-msa-max-rows must be at least 1')

    if args.report_msa == 'mafft':
        from tirmite.runners.mafft import mafft_available

        if not mafft_available():
            raise ValueError(
                '--report-msa mafft was requested but mafft was not found on '
                'PATH. Install MAFFT, or use --report-msa auto to fall back to '
                'placing hits by their model coordinates.'
            )

    if args.nopairing:
        logger.warning(
            'The HTML report will show individual hits only: --nopairing '
            'disables the pairing step, so it can contain no elements or '
            'terminus pairs.'
        )

    if args.no_elements and not args.report_no_sequences:
        # fetchElements is what produces the sequences the copy buttons use.
        logger.warning(
            '--no-elements disables element extraction, so element sequences '
            'cannot be embedded in the HTML report. Its copy buttons will show '
            'coordinates only.'
        )
        args.report_no_sequences = True


def _make_report_accumulator(
    args: argparse.Namespace,
    hitTable: Any,
    model_lengths: Dict[str, int],
    filter_stats: Dict[str, Any],
    extraction_source: Any,
) -> Optional[Any]:
    """
    Create the HTML report accumulator, or None when no report was requested.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments.
    hitTable : pandas.DataFrame
        The filtered hit table.
    model_lengths : dict
        Model name to declared length.
    filter_stats : dict
        Hit-filtering statistics accumulated so far.
    extraction_source : SequenceSource
        Source used to look up sequence lengths and hit sequences.

    Returns
    -------
    PairReportAccumulator or None
        The accumulator, or None if `--report` was not given or the report
        machinery could not be imported.
    """
    if not getattr(args, 'report', False):
        return None

    try:
        from tirmite.report.collect import PairReportAccumulator
    except Exception as exc:  # noqa: BLE001 - a report must not fail a run
        logger.error(f'HTML report unavailable: {exc}')
        return None

    genome_label = args.genome or args.blastdb or 'sequences'
    return PairReportAccumulator(
        tirmite_version=__version__,
        command=' '.join(sys.argv),
        title=args.report_title or f'TIRmite pair report — {Path(genome_label).name}',
        params={
            'orientation': args.orientation,
            'maxdist': args.maxdist,
            'mincov': args.mincov,
            'maxeval': args.maxeval,
            'max_offset': args.max_offset,
            'stable_reps': args.stable_reps,
        },
        filter_stats=filter_stats,
        hit_table=hitTable,
        model_lengths=model_lengths,
        contig_length=(
            extraction_source.contig_length if extraction_source is not None else None
        ),
        max_hits=args.report_max_hits,
        max_rows=args.report_max_rows,
    )


def _write_report(
    args: argparse.Namespace,
    accumulator: Any,
    unpaired_index: Dict[str, Dict[int, Dict[str, Any]]],
    outDir: str,
    tempDir: str,
    extraction_source: Any,
) -> None:
    """
    Finalise and write the HTML report.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments.
    accumulator : PairReportAccumulator or None
        The accumulator, or None if no report was requested.
    unpaired_index : dict
        The run's full hit index, so hits that no pairing procedure claimed
        still appear on the tracks.
    outDir : str
        Output directory.
    tempDir : str
        The run's temporary directory, used for MAFFT's intermediate files.
    extraction_source : SequenceSource
        Source used to read hit sequences for the alignment panels.

    Returns
    -------
    None
        Writes the report to disk.

    Notes
    -----
    Failures are logged and swallowed. By this point the run's real
    outputs -- GFF3, FASTA and the text summaries -- are already on disk, and a
    visualisation is not worth failing a completed analysis for.
    """
    if accumulator is None:
        return

    try:
        from tirmite.report.render import write_pair_report

        accumulator.add_unpaired(unpaired_index)

        prefix_str = f'{args.prefix}_' if args.prefix else ''
        outpath = Path(
            args.report_out
            or os.path.join(outDir, f'{prefix_str}tirmite_pair_report.html')
        )

        logger.info('Building HTML report...')
        data = accumulator.finalise(
            embed_sequences=not args.report_no_sequences,
            max_seq_bytes=int(args.report_max_seq_mb * 1024 * 1024),
            elements_fasta_path=outDir,
            source=extraction_source,
            tempdir=tempDir,
            msa_mode=args.report_msa,
            msa_max_rows=args.report_msa_max_rows,
        )
        write_pair_report(data, outpath)
    except Exception as exc:  # noqa: BLE001 - see the note above
        logger.error(f'HTML report generation failed: {exc}')
        logger.exception('Full traceback:')


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
) -> PairSummary:
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
    PairSummary
        The counted summary that was written. Returned so callers that also
        build the HTML report use exactly the numbers reported in the text
        file rather than recomputing them.

    Notes
    -----
    Counting and formatting live in :mod:`tirmite.report.stats`; this function
    is the file-writing wrapper around them.
    """
    summary = pair_summary_stats(
        left_feature=left_feature,
        right_feature=right_feature,
        pair_paired=pair_paired,
        pair_hitIndex=pair_hitIndex,
        total_pairs=total_pairs,
        total_elements=total_elements,
        filter_stats=filter_stats,
    )

    prefix_str = f'{prefix}_' if prefix else ''
    summary_path = os.path.join(outdir, f'{prefix_str}{summary.pair_label}_summary.txt')
    with open(summary_path, 'w') as f:
        f.write(format_pair_summary(summary))

    logger.info(f'Wrote summary report to {summary_path}')
    return summary


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
    parser = create_pair_parser()
    if args is None:
        args = parser.parse_args()

    # Mypy assertion: args is guaranteed non-None after parsing
    assert args is not None

    try:
        # Validate arguments
        try:
            validate_arguments(args)
        except (ValueError, FileNotFoundError) as e:
            logger.error(f'Argument validation failed: {e}')
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

        # Log startup information: version, package path, command, non-default args
        log_startup_info(args, parser)
        logger.info(f'Output directory: {outDir}')
        logger.debug(f'Temporary directory: {tempDir}')

        # Index genome if provided
        genome = None
        genome_descriptions = None
        if args.genome:
            logger.info(f'Indexing genome: {args.genome}')
            genome, genome_descriptions = indexGenome(args.genome)
        elif args.blastdb:
            logger.info(f'Using BLAST database: {args.blastdb}')
            blastdbcmd_path = shutil.which('blastdbcmd')
            if blastdbcmd_path:
                logger.info(f'blastdbcmd found: {blastdbcmd_path}')
            else:
                logger.warning(
                    'blastdbcmd not found in PATH; sequence extraction will fail'
                )

        # Load model/query lengths
        logger.info('Loading model/query lengths...')
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
            logger.info('Importing nhmmer hits...')
            input_format = 'nhmmer'
            hitTable = import_nhmmer(
                infile=args.nhmmer_file, hitTable=None, prefix=args.prefix
            )
        elif args.left_nhmmer:
            # Asymmetric nhmmer mode - import from both files
            logger.info('Importing nhmmer hits from left and right models...')
            input_format = 'nhmmer'

            # Check if left and right files are the same
            if args.left_nhmmer == args.right_nhmmer:
                raise ValueError(
                    f'Left and right nhmmer files cannot be the same: {args.left_nhmmer}'
                )

            left_model_name = extract_model_name_from_path(args.left_model)
            right_model_name = extract_model_name_from_path(args.right_model)

            # Import left file
            left_hitTable = import_nhmmer(
                infile=args.left_nhmmer,
                hitTable=None,
                prefix=args.prefix,
            )
            left_models = check_multiple_models(left_hitTable)
            logger.info(
                f'Left nhmmer file: {len(left_hitTable)} hits, {len(left_models)} unique query/model name(s)'
            )

            # Import right file
            right_hitTable = import_nhmmer(
                infile=args.right_nhmmer,
                hitTable=None,
                prefix=args.prefix,
            )
            right_models = check_multiple_models(right_hitTable)
            logger.info(
                f'Right nhmmer file: {len(right_hitTable)} hits, {len(right_models)} unique query/model name(s)'
            )

            # Check for overlapping query names
            overlapping_models = set(left_models) & set(right_models)
            if overlapping_models:
                logger.warning(
                    f'Query/model names appear in both left and right files: {", ".join(overlapping_models)}'
                )

            # Check for overlapping hits (same target, overlapping coordinates)
            left_overlap_count, right_overlap_count = check_overlapping_hits(
                left_hitTable, right_hitTable
            )
            if left_overlap_count > 0 or right_overlap_count > 0:
                logger.warning(
                    f'Found {left_overlap_count} left hit(s) and {right_overlap_count} right hit(s) '
                    'with overlapping genomic coordinates'
                )

            # Validate single query per file or require --pairing-map
            if len(left_models) > 1 or len(right_models) > 1:
                if len(left_models) > 1:
                    logger.warning(
                        f'Left nhmmer file contains {len(left_models)} query/model names: '
                        + ', '.join(sorted(left_models))
                    )
                if len(right_models) > 1:
                    logger.warning(
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
            hitTable = import_nhmmer(
                infile=args.left_nhmmer,
                hitTable=None,
                prefix=args.prefix,
            )
            hitTable = import_nhmmer(
                infile=args.right_nhmmer,
                hitTable=hitTable,
                prefix=args.prefix,
            )
        elif args.blast_file:
            # Single BLAST file mode
            logger.info('Importing BLAST hits...')
            input_format = 'blast'

            # Detect format and warn if mismatch
            detected_format = detect_input_format(args.blast_file)
            if detected_format != 'blast' and detected_format != 'unknown':
                logger.warning(
                    f'File format appears to be {detected_format}, but --blast-file was specified. '
                    'Consider using --nhmmer-file instead.'
                )

            hitTable = import_blast(
                infile=args.blast_file, hitTable=None, prefix=args.prefix
            )

        elif args.left_blast:
            # Asymmetric BLAST mode - import from both files
            logger.info('Importing BLAST hits from left and right queries...')
            input_format = 'blast'

            # Check if left and right files are the same
            if args.left_blast == args.right_blast:
                raise ValueError(
                    f'Left and right BLAST files cannot be the same: {args.left_blast}'
                )

            # Detect format for both files
            detected_left = detect_input_format(args.left_blast)
            detected_right = detect_input_format(args.right_blast)

            if detected_left != 'blast' and detected_left != 'unknown':
                logger.warning(
                    f'Left file format appears to be {detected_left}, but --left-blast was specified.'
                )
            if detected_right != 'blast' and detected_right != 'unknown':
                logger.warning(
                    f'Right file format appears to be {detected_right}, but --right-blast was specified.'
                )

            # Import left file
            left_hitTable = import_blast(
                infile=args.left_blast,
                hitTable=None,
                prefix=args.prefix,
            )
            left_models = check_multiple_models(left_hitTable)
            logger.info(
                f'Left BLAST file: {len(left_hitTable)} hits, {len(left_models)} unique query/model name(s)'
            )

            # Import right file
            right_hitTable = import_blast(
                infile=args.right_blast,
                hitTable=None,
                prefix=args.prefix,
            )
            right_models = check_multiple_models(right_hitTable)
            logger.info(
                f'Right BLAST file: {len(right_hitTable)} hits, {len(right_models)} unique query/model name(s)'
            )

            # Check for overlapping query names
            overlapping_models = set(left_models) & set(right_models)
            if overlapping_models:
                logger.warning(
                    f'Query/model names appear in both left and right files: {", ".join(overlapping_models)}'
                )

            # Check for overlapping hits (same target, overlapping coordinates)
            left_overlap_count, right_overlap_count = check_overlapping_hits(
                left_hitTable, right_hitTable
            )
            if left_overlap_count > 0 or right_overlap_count > 0:
                logger.warning(
                    f'Found {left_overlap_count} left hit(s) and {right_overlap_count} right hit(s) '
                    'with overlapping genomic coordinates'
                )

            # Validate single query per file or require --pairing-map
            if len(left_models) > 1 or len(right_models) > 1:
                if len(left_models) > 1:
                    logger.warning(
                        f'Left BLAST file contains {len(left_models)} query/model names: '
                        + ', '.join(sorted(left_models))
                    )
                if len(right_models) > 1:
                    logger.warning(
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
            hitTable = import_blast(
                infile=args.left_blast,
                hitTable=None,
                prefix=args.prefix,
            )
            hitTable = import_blast(
                infile=args.right_blast,
                hitTable=hitTable,
                prefix=args.prefix,
            )

        if hitTable is None or len(hitTable) == 0:
            logger.error('No hits were imported from input files')
            cleanup_temp_directory(tempDir, args.keep_temp)
            sys.exit(1)

        # Log import statistics
        initial_hit_count = len(hitTable)
        logger.info(f'Imported {initial_hit_count} total hits')

        # Get unique models and log statistics
        unique_models = check_multiple_models(hitTable)
        logger.info(f'Found {len(unique_models)} unique query/model name(s)')

        # Log per-query hit counts at debug level
        for model in unique_models:
            hit_count = len(hitTable[hitTable['model'] == model])
            logger.debug(f'Query/model "{model}": {hit_count} hits')

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
                logger.warning(
                    f'Pairing map references {len(map_models)} model(s). '
                    f'Ignoring {len(models_to_ignore)} model(s) not in pairing map: '
                    + ', '.join(sorted(models_to_ignore))
                )
                logger.info(
                    f'Excluding {ignored_hits} hits for {len(models_to_ignore)} ignored model(s)'
                )
                hitTable = hitTable[hitTable['model'].isin(models_to_process)].copy()
                logger.info(
                    f'{len(hitTable)} hits remaining for {len(models_to_process)} model(s) in pairing map'
                )
                if len(hitTable) == 0:
                    logger.error('No hits remaining after pairing map model filter')
                    cleanup_temp_directory(tempDir, args.keep_temp)
                    sys.exit(1)
                # Update unique_models after filtering
                unique_models = check_multiple_models(hitTable)
                filter_stats['pairing_map_models_ignored'] = sorted(models_to_ignore)
                filter_stats['pairing_map_hits_ignored'] = ignored_hits
            else:
                logger.info(
                    f'All {len(models_in_hits)} model(s) in hits are covered by pairing map'
                )

        # If --query-len was provided for BLAST input, assign it to ALL queries
        if args.blast_file and args.query_len:
            for query_name in unique_models:
                model_lengths[query_name] = args.query_len
                logger.debug(f'Set length for query {query_name}: {args.query_len}')
            logger.info(
                f'Applied query length {args.query_len} to {len(unique_models)} query name(s)'
            )

        # Validate that we have model lengths for all models in hitTable
        if not model_lengths:
            logger.error('No model/query lengths could be loaded')
            cleanup_temp_directory(tempDir, args.keep_temp)
            sys.exit(1)

        logger.info(f'Loaded lengths for models: {", ".join(model_lengths.keys())}')

        # Calculate hit coverage
        logger.info('Calculating hit coverage...')
        hitTable = calculate_hit_coverage(hitTable, model_lengths)

        # Apply coverage filtering
        logger.info(f'Filtering hits with coverage < {args.mincov}')
        hitCount = len(hitTable)

        # Log pre-filtering counts
        model_counts_before = hitTable['model'].value_counts().to_dict()
        for model, count in model_counts_before.items():
            logger.info(f'Model {model}: {count} hits before coverage filtering')

        hitTable = filter_hits_coverage(hitTable, args.mincov)

        # Log post-filtering counts
        model_counts_after = hitTable['model'].value_counts().to_dict()
        coverage_excluded = hitCount - len(hitTable)
        logger.info(f'Excluded {coverage_excluded} hits on coverage criteria')
        filter_stats['coverage_excluded'] = coverage_excluded

        for model in model_counts_before:
            before = model_counts_before[model]
            after = model_counts_after.get(model, 0)
            excluded = before - after
            logger.info(f'Model {model}: {excluded} excluded, {after} remaining')

        # Apply e-value filtering
        logger.info(f'Filtering hits with e-value > {args.maxeval}')
        hitCount = len(hitTable)

        model_counts_before = hitTable['model'].value_counts().to_dict()
        hitTable = filterHitsEval(maxeval=args.maxeval, hitTable=hitTable)

        model_counts_after = hitTable['model'].value_counts().to_dict()
        evalue_excluded = hitCount - len(hitTable)
        logger.info(f'Excluded {evalue_excluded} hits on e-value criteria')
        filter_stats['evalue_excluded'] = evalue_excluded

        for model in model_counts_before:
            before = model_counts_before[model]
            after = model_counts_after.get(model, 0)
            excluded = before - after
            logger.info(f'Model {model}: {excluded} excluded, {after} remaining')

        logger.info(f'Total remaining hits: {len(hitTable)}')

        # Verify that every target sequence can be resolved in the chosen
        # sequence source before extraction begins. On the BLAST side the usual
        # cause of failure is a database built without -parse_seqids, or an
        # accession that differs from the FASTA header token.
        # Built unconditionally so the report can ask it for sequence lengths
        # and hit sequences even when the hit table came back empty.
        extraction_source = make_source(genome=genome, blastdb=args.blastdb)
        if len(hitTable):
            missing_targets = check_ids(extraction_source, hitTable['target'])
            if missing_targets:
                logger.warning(
                    f'{len(missing_targets)} target sequence(s) could not be '
                    'resolved; hits on those sequences will be skipped.'
                )

        # Apply anchor filter if --max-offset is set
        if args.max_offset is not None:
            logger.info(
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
                logger.error('No hits remaining after anchor filter')
                cleanup_temp_directory(tempDir, args.keep_temp)
                sys.exit(1)

        # Finalise after_filtering count after all filter steps
        filter_stats['after_filtering'] = len(hitTable)

        # Convert to dict structure
        hitsDict, hitIndex = table2dict(hitTable)

        # The report accumulates each pairing procedure's outcome as it runs,
        # and is written once at the end. Built here so it sees the same
        # filtered hit table and the same filter statistics the text summaries
        # report.
        report_acc = _make_report_accumulator(
            args, hitTable, model_lengths, filter_stats, extraction_source
        )

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
            logger.warning(
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
            logger.info(
                f'Will execute {len(pairing_map)} independent pairing procedure(s) based on pairing map'
            )
            config = None
        elif args.left_nhmmer and args.right_nhmmer:
            # Asymmetric nhmmer pairing
            left_model_name = extract_model_name_from_path(args.left_model)
            right_model_name = extract_model_name_from_path(args.right_model)

            config = PairingConfig(
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
                logger.error(
                    'Expected at least 1 model in each BLAST file for asymmetric pairing'
                )
                cleanup_temp_directory(tempDir, args.keep_temp)
                sys.exit(1)

            left_model_name = left_model_names[0]
            right_model_name = right_model_names[0]
            logger.info(
                f'Assigning models for asymmetric pairing: '
                f'left={left_model_name} (from --left-blast), '
                f'right={right_model_name} (from --right-blast)'
            )

            config = PairingConfig(
                orientation=args.orientation,
                left_model=left_model_name,
                right_model=right_model_name,
            )
        else:
            # Symmetric pairing
            available_models = list(hitsDict.keys())
            if not available_models:
                logger.error('No models found in hits')
                cleanup_temp_directory(tempDir, args.keep_temp)
                sys.exit(1)

            single_model = available_models[0]
            config = PairingConfig(
                orientation=args.orientation, single_model=single_model
            )

        # Write individual hits
        # In pairing_map mode individual hits are written per pair (below),
        # so only write to the base outDir when no pairing map is used.
        if not pairing_map and not args.no_hits:
            logger.info('Writing individual hits to FASTA...')
            writeTIRs(
                outDir=outDir,
                hitTable=hitTable,
                maxeval=args.maxeval,
                genome=genome,
                prefix=args.prefix,
                padlen=args.padlen,
                genome_descriptions=genome_descriptions,
                blastdb=args.blastdb if args.blastdb else None,
                model_lengths=model_lengths,
                extend_hits_to_model=args.extend_hits_to_model,
                pad=not args.no_pad_flanks,
            )

        # Skip pairing if requested
        if args.nopairing:
            logger.info('Pairing disabled. Analysis complete.')
            # A hits-only report: no groups were added, so it shows the
            # annotation tracks, the alignment panels and the hit statistics.
            _write_report(
                args, report_acc, hitIndex, outDir, tempDir, extraction_source
            )
            cleanup_temp_directory(tempDir, args.keep_temp)
            return 0

        # Run pairing - handle single or multiple pairing procedures
        if pairing_map:
            # Multiple pairing procedures based on pairing map
            logger.info(
                f'Running {len(pairing_map)} independent pairing procedures based on pairing map'
            )

            # Accumulate all paired results by model.
            all_paired: Dict[str, List[Set[int]]] = {}
            # Track which hit indices have been paired.
            all_paired_hits: Set[int] = set()
            original_hitIndex = hitIndex  # Preserve original for unpaired hit tracking

            for pair_idx, (left_feature, right_feature) in enumerate(pairing_map, 1):
                logger.info(
                    f'Pairing procedure {pair_idx}/{len(pairing_map)}: {left_feature} <-> {right_feature}'
                )

                # Check if both features exist in hitsDict
                if left_feature not in hitsDict:
                    logger.warning(
                        f'Feature "{left_feature}" not found in hits, skipping this pairing'
                    )
                    continue
                if right_feature not in hitsDict:
                    logger.warning(
                        f'Feature "{right_feature}" not found in hits, skipping this pairing'
                    )
                    continue

                # Create config for this pairing
                if left_feature == right_feature:
                    # Symmetric pairing
                    pair_config = PairingConfig(
                        orientation=args.orientation, single_model=left_feature
                    )
                else:
                    # Asymmetric pairing
                    pair_config = PairingConfig(
                        orientation=args.orientation,
                        left_model=left_feature,
                        right_model=right_feature,
                    )

                # Run pairing for this combination (use fresh copy of hitIndex for each)
                logger.info(
                    f'Searching for candidate pairings: {left_feature} <-> {right_feature}'
                )
                pair_hitIndex = parseHitsGeneral(
                    hitsDict=hitsDict,
                    hitIndex=original_hitIndex,
                    maxDist=args.maxdist,
                    config=pair_config,
                )

                logger.info('Performing iterative pairing...')
                if pair_config.is_asymmetric:
                    logger.info(
                        f'Using asymmetric pairing: {pair_config.left_model} + {pair_config.right_model}'
                    )
                    pair_hitIndex, pair_paired, pair_unpaired = (
                        iterateGetPairsAsymmetric(
                            pair_hitIndex, pair_config, stableReps=args.stable_reps
                        )
                    )
                else:
                    logger.info(
                        f'Using symmetric pairing with orientation {pair_config.orientation}'
                    )
                    pair_hitIndex, pair_paired, pair_unpaired = iterateGetPairsCustom(
                        pair_hitIndex, pair_config, stableReps=args.stable_reps
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
                if not args.no_hits:
                    logger.info(f'Writing individual hits for pair {pair_label}...')
                    writeTIRs(
                        outDir=pair_outDir,
                        hitTable=pair_hitTable_tirs,
                        maxeval=args.maxeval,
                        genome=genome,
                        prefix=args.prefix,
                        padlen=args.padlen,
                        genome_descriptions=genome_descriptions,
                        blastdb=args.blastdb if args.blastdb else None,
                        model_lengths=model_lengths,
                        extend_hits_to_model=args.extend_hits_to_model,
                        pad=not args.no_pad_flanks,
                    )

                # Write paired TIRs
                if args.gff_report in ['all', 'paired']:
                    writePairedTIRs(
                        outDir=pair_outDir,
                        paired=pair_paired,
                        hitIndex=pair_hitIndex,
                        genome=genome,
                        prefix=args.prefix,
                        padlen=args.padlen,
                        genome_descriptions=genome_descriptions,
                        blastdb=args.blastdb if args.blastdb else None,
                        config=pair_config,
                    )

                # Extract and write elements for this pair
                if not args.no_elements:
                    pair_pairedEles = fetchElements(
                        paired=pair_paired,
                        hitIndex=pair_hitIndex,
                        genome=genome,
                        genome_descriptions=genome_descriptions,
                        blastdb=args.blastdb if args.blastdb else None,
                    )
                    writeElements(
                        pair_outDir, eleDict=pair_pairedEles, prefix=args.prefix
                    )
                else:
                    pair_pairedEles = {}

                # Extract and write flanks for this pair
                if args.flanks or args.flanks_paired:
                    writeFlanks(
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
                        pad_flanks=not args.no_pad_flanks,
                    )

                    # Reconstruct target sites for this pair
                    if args.insertion_site:
                        tsd_length_map_data = None
                        if args.tsd_length_map:
                            tsd_length_map_data = load_tsd_length_map(
                                args.tsd_length_map
                            )

                        writeTargetSites(
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
                            pad_flanks=not args.no_pad_flanks,
                        )

                # Write summary report for this pair
                total_pair_pairs = sum(len(pairs) for pairs in pair_paired.values())
                total_pair_elements = sum(
                    len(eles) for eles in pair_pairedEles.values()
                )
                pair_summary = _write_pair_summary(
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

                # Recorded here rather than after the loop: this loop reuses
                # one hit index across every row, so a later row overwrites the
                # pairing state this one produced. The accumulator snapshots
                # what it is given.
                if report_acc is not None:
                    report_acc.add_group(
                        left_feature=left_feature,
                        right_feature=right_feature,
                        config=pair_config,
                        paired=pair_paired,
                        hit_index=pair_hitIndex,
                        elements=pair_pairedEles,
                        summary=pair_summary,
                    )

                # Write GFF for this pair if requested
                if args.gff_out:
                    pair_unpairedTIRs = None
                    if args.gff_report in ['all', 'unpaired']:
                        pair_unpairedTIRs = fetchUnpaired(hitIndex=pair_hitIndex)
                        logger.info(
                            f'Pair {pair_label}: {len(pair_unpairedTIRs)} unpaired termini'
                        )

                    if args.prefix:
                        pair_gffOutPath = os.path.join(
                            pair_outDir, f'{args.prefix}_{pair_label}.gff3'
                        )
                    else:
                        pair_gffOutPath = os.path.join(
                            pair_outDir, f'{pair_label}.gff3'
                        )

                    logger.info(f'Writing GFF3 output: {pair_gffOutPath}')
                    gffWrite(
                        outpath=pair_gffOutPath,
                        featureList=pair_pairedEles,
                        writeTIRs=args.gff_report,
                        unpaired=pair_unpairedTIRs,
                        prefix=args.prefix,
                    )

                logger.info(f'Pair {pair_label}: wrote output to {pair_outDir}')

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
                logger.info(
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
            logger.info(
                f'All pairing procedures completed: {total_pairs} total pairs, {total_unpaired} unpaired hits'
            )

        else:
            # Single pairing procedure (original logic)
            logger.info('Searching for candidate pairings...')
            hitIndex = parseHitsGeneral(
                hitsDict=hitsDict,
                hitIndex=hitIndex,
                maxDist=args.maxdist,
                config=config,
            )

            logger.info('Performing iterative pairing...')
            if config.is_asymmetric:
                logger.info(
                    f'Using asymmetric pairing: {config.left_model} + {config.right_model}'
                )
                hitIndex, paired, unpaired = iterateGetPairsAsymmetric(
                    hitIndex, config, stableReps=args.stable_reps
                )
            else:
                logger.info(
                    f'Using symmetric pairing with orientation {config.orientation}'
                )
                hitIndex, paired, unpaired = iterateGetPairsCustom(
                    hitIndex, config, stableReps=args.stable_reps
                )

            # Log pairing results
            total_pairs = sum(len(pairs) for pairs in paired.values())
            total_unpaired = len(unpaired)
            logger.info(
                f'Pairing completed: {total_pairs} pairs, {total_unpaired} unpaired'
            )

        # For pairing-map mode all output has already been written to per-pair
        # subdirectories inside the loop above.  The block below handles the
        # single-pairing (no pairing map) case only.
        if not pairing_map:
            # Write paired TIRs
            if args.gff_report in ['all', 'paired']:
                logger.info('Writing paired termini to FASTA...')
                writePairedTIRs(
                    outDir=outDir,
                    paired=paired,
                    hitIndex=hitIndex,
                    genome=genome,
                    prefix=args.prefix,
                    padlen=args.padlen,
                    genome_descriptions=genome_descriptions,
                    blastdb=args.blastdb if args.blastdb else None,
                    config=config,
                )

            # Extract and write elements
            if not args.no_elements:
                logger.info('Extracting paired elements...')
                pairedEles = fetchElements(
                    paired=paired,
                    hitIndex=hitIndex,
                    genome=genome,
                    genome_descriptions=genome_descriptions,
                    blastdb=args.blastdb if args.blastdb else None,
                )

                logger.info('Writing paired elements to FASTA...')
                writeElements(outDir, eleDict=pairedEles, prefix=args.prefix)
            else:
                pairedEles = {}

            # Write summary report for single-pairing mode
            if config is not None:
                left_feature = config.left_model or ''
                right_feature = config.right_model or ''
                total_ele_count = sum(len(eles) for eles in pairedEles.values())
                single_summary = _write_pair_summary(
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
                if report_acc is not None:
                    report_acc.add_group(
                        left_feature=left_feature,
                        right_feature=right_feature,
                        config=config,
                        paired=paired,
                        hit_index=hitIndex,
                        elements=pairedEles,
                        summary=single_summary,
                    )

            # Extract and write external flanks if requested
            if args.flanks or args.flanks_paired:
                logger.info('Extracting external flanking sequences...')
                writeFlanks(
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
                    pad_flanks=not args.no_pad_flanks,
                )

                # Reconstruct target sites if explicitly requested
                if args.insertion_site:
                    logger.info('Reconstructing target sites...')

                    # Load TSD length map if provided
                    tsd_length_map = None
                    if args.tsd_length_map:
                        tsd_length_map = load_tsd_length_map(args.tsd_length_map)

                    writeTargetSites(
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
                        pad_flanks=not args.no_pad_flanks,
                    )

            # Write GFF if requested
            if args.gff_out:
                # Get unpaired TIRs if needed
                unpairedTIRs = None
                if args.gff_report in ['all', 'unpaired']:
                    unpairedTIRs = fetchUnpaired(hitIndex=hitIndex)
                    logger.info(f'Found {len(unpairedTIRs)} unpaired termini')

                # Set GFF output path
                if args.prefix:
                    gffOutPath = os.path.join(outDir, f'{args.prefix}.gff3')
                else:
                    gffOutPath = os.path.join(outDir, 'tirmite_pair_report.gff3')

                logger.info(f'Writing GFF3 output: {gffOutPath}')
                gffWrite(
                    outpath=gffOutPath,
                    featureList=pairedEles,
                    writeTIRs=args.gff_report,
                    unpaired=unpairedTIRs,
                    prefix=args.prefix,
                )

        # Written last, while tempDir is still alive for MAFFT: every other
        # output is already on disk, so a failure here costs only the report.
        _write_report(args, report_acc, hitIndex, outDir, tempDir, extraction_source)

        logger.info('TIRmite-pair analysis completed successfully')

    except KeyboardInterrupt:
        logger.info('Analysis interrupted by user')
        sys.exit(130)
    except Exception as e:
        logger.error(f'Unexpected error: {e}')
        logger.exception('Full traceback:')  # This will log the full stack trace
        sys.exit(1)
    finally:
        if 'tempDir' in locals():
            cleanup_temp_directory(tempDir, args.keep_temp)

        # Log completion message with logfile location if enabled
        if 'logfile_path' in locals() and logfile_path and args.logfile:
            logger.info(f'Log file saved to: {logfile_path}')

    return 0


if __name__ == '__main__':
    main()
