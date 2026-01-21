#!/usr/bin/env python3
"""
Ensemble search workflow for TIRmite.

This module implements an ensemble search strategy for transposon terminal detection:
1. Run BLAST and/or nhmmer searches against target genomes
2. Load and filter hits by quality thresholds
3. Optionally merge overlapping hits from clustered component features
4. Remove nested weak hits between paired terminal models
5. Output merged hits in BLAST-compatible tabular format

Supports both direct search execution and loading of precomputed search results.
"""

import argparse
import logging
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, cast

import pandas as pd  # type: ignore[import-untyped]

from tirmite.runners.hmmer_wrappers import build_nhmmer_command
from tirmite.runners.runBlastn import BlastError, run_blastn
from tirmite.runners.wrapping import run_command
import tirmite.tirmitetools as tirmite
from tirmite.utils.logs import init_logging
from tirmite.utils.utils import temporary_directory


class EnsembleSearchError(Exception):
    """Custom exception for ensemble search errors."""

    pass


# -----------------------------------------------------------------------------
# Cluster Mapping Functions
# -----------------------------------------------------------------------------


def parse_cluster_mapping(mapping_file: Path) -> Dict[str, List[str]]:
    """
    Parse cluster mapping file linking component features to cluster names.

    Parameters
    ----------
    mapping_file : Path
        Path to tab-delimited cluster mapping file.
        Format: cluster_name<TAB>component1<TAB>component2<TAB>...

    Returns
    -------
    dict
        Dictionary mapping cluster names to lists of component feature names.

    Raises
    ------
    EnsembleSearchError
        If file cannot be read or has invalid format.
    """
    cluster_map: Dict[str, List[str]] = {}

    try:
        with open(mapping_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                # Skip comments and empty lines
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) < 2:
                    logging.warning(
                        f'Skipping line {line_num}: insufficient columns (need cluster name + at least one component)'
                    )
                    continue

                cluster_name = parts[0].strip()
                components = [p.strip() for p in parts[1:] if p.strip()]

                if not cluster_name:
                    logging.warning(f'Skipping line {line_num}: empty cluster name')
                    continue

                if not components:
                    logging.warning(
                        f'Skipping line {line_num}: no component features specified for cluster {cluster_name}'
                    )
                    continue

                cluster_map[cluster_name] = components

    except Exception as e:
        raise EnsembleSearchError(
            f'Failed to parse cluster mapping file {mapping_file}: {e}'
        ) from e

    return cluster_map


def validate_cluster_mapping(
    cluster_map: Dict[str, List[str]], available_features: Set[str]
) -> Tuple[bool, List[str]]:
    """
    Validate cluster mapping for uniqueness and completeness.

    Parameters
    ----------
    cluster_map : dict
        Dictionary mapping cluster names to component lists.
    available_features : set
        Set of feature names present in the hit data.

    Returns
    -------
    tuple
        (is_valid, list of warning messages)

    Notes
    -----
    Checks:
    - Component names appear in only one cluster
    - Cluster names are unique (implicit in dict structure)
    - All available features are assigned to clusters (warns if not)
    """
    warnings: List[str] = []
    is_valid = True

    # Check for duplicate component assignments
    seen_components: Dict[str, str] = {}  # component -> cluster
    for cluster_name, components in cluster_map.items():
        for component in components:
            if component in seen_components:
                warnings.append(
                    f"Component '{component}' appears in multiple clusters: "
                    f"'{seen_components[component]}' and '{cluster_name}'"
                )
                is_valid = False
            else:
                seen_components[component] = cluster_name

    # Check for features not assigned to any cluster
    all_components = set(seen_components.keys())
    unassigned = available_features - all_components
    if unassigned:
        warnings.append(
            f'Features not assigned to any cluster will be ignored: {sorted(unassigned)}'
        )

    # Check for cluster components not in available features
    missing = all_components - available_features
    if missing:
        warnings.append(
            f'Cluster components not found in hit data: {sorted(missing)}'
        )

    return is_valid, warnings


def build_component_to_cluster_map(
    cluster_map: Dict[str, List[str]]
) -> Dict[str, str]:
    """
    Build reverse mapping from component names to cluster names.

    Parameters
    ----------
    cluster_map : dict
        Dictionary mapping cluster names to component lists.

    Returns
    -------
    dict
        Dictionary mapping component names to their cluster name.
    """
    component_to_cluster: Dict[str, str] = {}
    for cluster_name, components in cluster_map.items():
        for component in components:
            component_to_cluster[component] = cluster_name
    return component_to_cluster


# -----------------------------------------------------------------------------
# Pairing Map Functions
# -----------------------------------------------------------------------------


def parse_pairing_map(pairing_file: Path) -> Dict[str, str]:
    """
    Parse pairing map file linking left and right terminal features.

    Parameters
    ----------
    pairing_file : Path
        Path to tab-delimited pairing map file.
        Format: left_feature<TAB>right_feature

    Returns
    -------
    dict
        Dictionary mapping left feature names to right feature names.

    Raises
    ------
    EnsembleSearchError
        If file cannot be read or has invalid format.
    """
    pairing_map: Dict[str, str] = {}

    try:
        with open(pairing_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) != 2:
                    logging.warning(
                        f'Skipping line {line_num}: expected 2 columns (left, right), got {len(parts)}'
                    )
                    continue

                left_name = parts[0].strip()
                right_name = parts[1].strip()

                if not left_name or not right_name:
                    logging.warning(f'Skipping line {line_num}: empty feature name')
                    continue

                pairing_map[left_name] = right_name

    except Exception as e:
        raise EnsembleSearchError(
            f'Failed to parse pairing map file {pairing_file}: {e}'
        ) from e

    return pairing_map


# -----------------------------------------------------------------------------
# Hit Loading and Filtering
# -----------------------------------------------------------------------------


def load_hits_from_files(
    blast_files: Optional[List[Path]] = None,
    nhmmer_files: Optional[List[Path]] = None,
) -> pd.DataFrame:
    """
    Load hits from BLAST and/or nhmmer output files.

    Parameters
    ----------
    blast_files : list of Path, optional
        List of BLAST tabular output files to load.
    nhmmer_files : list of Path, optional
        List of nhmmer tabular output files to load.

    Returns
    -------
    pandas.DataFrame
        Combined hit table with standardized columns.

    Raises
    ------
    EnsembleSearchError
        If no input files provided or loading fails.
    """
    if not blast_files and not nhmmer_files:
        raise EnsembleSearchError(
            'Must provide at least one BLAST or nhmmer file to load'
        )

    hit_table: Optional[pd.DataFrame] = None

    # Load BLAST files
    if blast_files:
        for blast_file in blast_files:
            if not blast_file.exists():
                logging.warning(f'BLAST file not found, skipping: {blast_file}')
                continue

            logging.info(f'Loading BLAST hits from: {blast_file}')
            try:
                hit_table = tirmite.import_blast(str(blast_file), hitTable=hit_table)
            except Exception as e:
                logging.error(f'Failed to load BLAST file {blast_file}: {e}')
                continue

    # Load nhmmer files
    if nhmmer_files:
        for nhmmer_file in nhmmer_files:
            if not nhmmer_file.exists():
                logging.warning(f'nhmmer file not found, skipping: {nhmmer_file}')
                continue

            logging.info(f'Loading nhmmer hits from: {nhmmer_file}')
            try:
                hit_table = tirmite.import_nhmmer(str(nhmmer_file), hitTable=hit_table)
            except Exception as e:
                logging.error(f'Failed to load nhmmer file {nhmmer_file}: {e}')
                continue

    if hit_table is None or hit_table.empty:
        logging.warning('No hits loaded from input files')
        # Return empty DataFrame with correct columns
        cols = [
            'model',
            'target',
            'hitStart',
            'hitEnd',
            'strand',
            'evalue',
            'score',
            'bias',
            'hmmStart',
            'hmmEnd',
        ]
        return pd.DataFrame(columns=cols)

    return hit_table


def filter_hits_by_evalue(hit_table: pd.DataFrame, max_evalue: float) -> pd.DataFrame:
    """
    Filter hits by e-value threshold.

    Parameters
    ----------
    hit_table : pandas.DataFrame
        Hit table with 'evalue' column.
    max_evalue : float
        Maximum e-value threshold for keeping hits.

    Returns
    -------
    pandas.DataFrame
        Filtered hit table.
    """
    # Convert evalue column to float for comparison
    hit_table = hit_table.copy()
    hit_table['evalue_float'] = hit_table['evalue'].astype(float)

    filtered = hit_table[hit_table['evalue_float'] <= max_evalue].copy()
    filtered = filtered.drop(columns=['evalue_float'])

    logging.info(
        f'E-value filter (max={max_evalue}): {len(hit_table)} -> {len(filtered)} hits'
    )
    return filtered


def report_hit_statistics(hit_table: pd.DataFrame, stage: str = '') -> None:
    """
    Report statistics about hits in the table.

    Parameters
    ----------
    hit_table : pandas.DataFrame
        Hit table to analyze.
    stage : str, optional
        Label for the processing stage (e.g., 'before filtering', 'after filtering').
    """
    if hit_table.empty:
        logging.info(f'Hit statistics {stage}: 0 total hits')
        return

    total_hits = len(hit_table)
    unique_models = hit_table['model'].nunique()
    unique_targets = hit_table['target'].nunique()

    logging.info(f'Hit statistics {stage}:')
    logging.info(f'  Total hits: {total_hits}')
    logging.info(f'  Unique query features (models): {unique_models}')
    logging.info(f'  Unique target sequences: {unique_targets}')

    # Report hits per model
    hits_per_model = hit_table.groupby('model').size()
    logging.info('  Hits per feature:')
    for model, count in hits_per_model.items():
        logging.info(f'    {model}: {count}')


# -----------------------------------------------------------------------------
# Hit Merging Functions
# -----------------------------------------------------------------------------


def merge_overlapping_cluster_hits(
    hit_table: pd.DataFrame,
    cluster_map: Dict[str, List[str]],
    min_overlap: int = 1,
) -> pd.DataFrame:
    """
    Merge overlapping hits from component features within the same cluster.

    Parameters
    ----------
    hit_table : pandas.DataFrame
        Hit table with 'model', 'target', 'hitStart', 'hitEnd', 'strand', 'score' columns.
    cluster_map : dict
        Dictionary mapping cluster names to lists of component feature names.
    min_overlap : int, default 1
        Minimum overlap in base pairs to trigger merging.

    Returns
    -------
    pandas.DataFrame
        Hit table with overlapping same-cluster hits merged.
        Merged hits inherit properties from highest-scoring component.
    """
    if hit_table.empty:
        return hit_table

    component_to_cluster = build_component_to_cluster_map(cluster_map)

    # Add cluster column to hits
    hit_table = hit_table.copy()
    hit_table['cluster'] = hit_table['model'].map(component_to_cluster)

    # Separate clustered and unclustered hits
    clustered_hits = hit_table[hit_table['cluster'].notna()].copy()
    unclustered_hits = hit_table[hit_table['cluster'].isna()].copy()

    if unclustered_hits.shape[0] > 0:
        logging.warning(
            f'{len(unclustered_hits)} hits from unclustered features will be ignored'
        )

    if clustered_hits.empty:
        logging.warning('No clustered hits to merge')
        return pd.DataFrame(columns=hit_table.columns.drop('cluster'))

    # Group by target, strand, and cluster for merging
    merged_records: List[Dict[str, Any]] = []

    for (_target, _strand, cluster), group in clustered_hits.groupby(
        ['target', 'strand', 'cluster']
    ):
        # Sort by start position
        group = group.copy()
        group['hitStart_int'] = group['hitStart'].astype(int)
        group['hitEnd_int'] = group['hitEnd'].astype(int)
        group = group.sort_values('hitStart_int')

        # Convert score to float for comparison
        group['score_float'] = pd.to_numeric(group['score'], errors='coerce').fillna(0)

        # Merge overlapping hits
        current_hits: List[pd.Series] = []
        current_start: Optional[int] = None
        current_end: Optional[int] = None

        for _, hit in group.iterrows():
            hit_start = int(hit['hitStart_int'])
            hit_end = int(hit['hitEnd_int'])

            if current_start is None:
                # First hit
                current_start = hit_start
                current_end = hit_end
                current_hits = [hit]
            else:
                # current_end is guaranteed to be set when current_start is set
                assert current_end is not None
                if hit_start <= current_end + min_overlap:
                    # Overlapping or adjacent - extend region
                    current_end = max(current_end, hit_end)
                    current_hits.append(hit)
                else:
                    # Non-overlapping - save current merged hit and start new one
                    merged_records.append(
                        _create_merged_hit(current_hits, current_start, current_end, cluster)
                    )
                    current_start = hit_start
                    current_end = hit_end
                    current_hits = [hit]

        # Don't forget the last merged region
        if current_hits and current_start is not None and current_end is not None:
            merged_records.append(
                _create_merged_hit(current_hits, current_start, current_end, cluster)
            )

    if not merged_records:
        cols = [
            'model',
            'target',
            'hitStart',
            'hitEnd',
            'strand',
            'evalue',
            'score',
            'bias',
            'hmmStart',
            'hmmEnd',
        ]
        return pd.DataFrame(columns=cols)

    # Create merged DataFrame
    merged_df = pd.DataFrame(merged_records)

    # Sort by model, target, location
    merged_df = merged_df.sort_values(
        ['model', 'target', 'hitStart', 'hitEnd', 'strand'],
        ascending=[True, True, True, True, True],
    )
    merged_df = merged_df.reset_index(drop=True)

    logging.info(
        f'Merged overlapping cluster hits: {len(clustered_hits)} -> {len(merged_df)} hits'
    )

    return merged_df


def _create_merged_hit(
    hits: List[pd.Series],
    merged_start: int,
    merged_end: int,
    cluster_name: str,
) -> Dict[str, Any]:
    """
    Create a merged hit record from overlapping component hits.

    Parameters
    ----------
    hits : list of pandas.Series
        List of overlapping hit records to merge.
    merged_start : int
        Start position of merged region.
    merged_end : int
        End position of merged region.
    cluster_name : str
        Name of the cluster for the merged hit.

    Returns
    -------
    dict
        Merged hit record inheriting properties from highest-scoring component.
    """
    # Find highest scoring hit
    best_hit = max(hits, key=lambda h: float(h.get('score_float', 0)))

    return {
        'model': cluster_name,
        'target': best_hit['target'],
        'hitStart': str(merged_start),
        'hitEnd': str(merged_end),
        'strand': best_hit['strand'],
        'evalue': best_hit['evalue'],
        'score': best_hit['score'],
        'bias': best_hit['bias'],
        'hmmStart': best_hit['hmmStart'],
        'hmmEnd': best_hit['hmmEnd'],
    }


def check_cross_cluster_overlaps(
    hit_table: pd.DataFrame,
    cluster_map: Dict[str, List[str]],
    min_overlap: int = 1,
) -> None:
    """
    Check for and warn about overlapping hits between different clusters.

    Parameters
    ----------
    hit_table : pandas.DataFrame
        Hit table with 'model', 'target', 'hitStart', 'hitEnd', 'strand' columns.
    cluster_map : dict
        Dictionary mapping cluster names to component lists.
    min_overlap : int, default 1
        Minimum overlap to report.
    """
    if hit_table.empty:
        return

    component_to_cluster = build_component_to_cluster_map(cluster_map)

    # Add cluster column
    hit_table = hit_table.copy()
    hit_table['cluster'] = hit_table['model'].map(component_to_cluster)
    hit_table['hitStart_int'] = hit_table['hitStart'].astype(int)
    hit_table['hitEnd_int'] = hit_table['hitEnd'].astype(int)

    # Check each target/strand combination
    warnings_reported = 0
    for (target, _strand), group in hit_table.groupby(['target', 'strand']):
        if len(group) < 2:
            continue

        # Sort by position
        group = group.sort_values('hitStart_int')
        hits_list = list(group.iterrows())

        for i, (_, hit1) in enumerate(hits_list):
            for _, hit2 in hits_list[i + 1 :]:
                # Skip if same cluster
                if hit1['cluster'] == hit2['cluster']:
                    continue

                # Check overlap
                start1 = int(hit1['hitStart_int'])
                end1 = int(hit1['hitEnd_int'])
                start2 = int(hit2['hitStart_int'])
                end2 = int(hit2['hitEnd_int'])

                overlap_start = max(start1, start2)
                overlap_end = min(end1, end2)
                overlap = overlap_end - overlap_start

                if overlap >= min_overlap:
                    if warnings_reported < 10:  # Limit warnings
                        logging.warning(
                            f"Cross-cluster overlap detected: {target}:{start1}-{end1} ({hit1['cluster']}) "
                            f"overlaps with {target}:{start2}-{end2} ({hit2['cluster']})"
                        )
                    warnings_reported += 1

    if warnings_reported > 10:
        logging.warning(f'... and {warnings_reported - 10} more cross-cluster overlaps')


# -----------------------------------------------------------------------------
# Nested Hit Removal
# -----------------------------------------------------------------------------


def remove_nested_paired_hits(
    hit_table: pd.DataFrame,
    pairing_map: Dict[str, str],
    min_score_ratio: float = 1.5,
) -> pd.DataFrame:
    """
    Remove weak hits nested within stronger paired terminal hits.

    Parameters
    ----------
    hit_table : pandas.DataFrame
        Hit table with 'model', 'target', 'hitStart', 'hitEnd', 'strand', 'score' columns.
    pairing_map : dict
        Dictionary mapping left feature names to right feature names.
    min_score_ratio : float, default 1.5
        Minimum score ratio for keeping nested hit (nested_score / enclosing_score).

    Returns
    -------
    pandas.DataFrame
        Hit table with nested weak hits removed.
    """
    if hit_table.empty or not pairing_map:
        return hit_table

    # Build reverse mapping (right -> left) - not used currently but may be useful
    # reverse_pairing = {v: k for k, v in pairing_map.items()}
    all_paired_features = set(pairing_map.keys()) | set(pairing_map.values())

    hit_table = hit_table.copy()
    hit_table['hitStart_int'] = hit_table['hitStart'].astype(int)
    hit_table['hitEnd_int'] = hit_table['hitEnd'].astype(int)
    hit_table['score_float'] = pd.to_numeric(hit_table['score'], errors='coerce').fillna(0)

    hits_to_remove: Set[int] = set()

    # Check each target/strand combination
    for (_target, _strand), group in hit_table.groupby(['target', 'strand']):
        if len(group) < 2:
            continue

        group_indices = list(group.index)

        for i, idx1 in enumerate(group_indices):
            hit1 = hit_table.loc[idx1]
            model1 = hit1['model']

            if model1 not in all_paired_features:
                continue

            for idx2 in group_indices[i + 1 :]:
                hit2 = hit_table.loc[idx2]
                model2 = hit2['model']

                if model2 not in all_paired_features:
                    continue

                # Check if they are paired (left-right or right-left)
                is_paired = (
                    (model1 in pairing_map and pairing_map[model1] == model2)
                    or (model2 in pairing_map and pairing_map[model2] == model1)
                )

                if not is_paired:
                    continue

                # Check if one is nested within the other
                start1 = int(hit1['hitStart_int'])
                end1 = int(hit1['hitEnd_int'])
                start2 = int(hit2['hitStart_int'])
                end2 = int(hit2['hitEnd_int'])
                score1 = float(hit1['score_float'])
                score2 = float(hit2['score_float'])

                # Check if hit2 is nested in hit1
                if start1 <= start2 and end2 <= end1:
                    # hit2 is nested in hit1
                    if score1 > 0 and (score2 / score1) < min_score_ratio:
                        hits_to_remove.add(idx2)
                        logging.debug(
                            f'Removing nested hit {model2} ({start2}-{end2}, score={score2:.1f}) '
                            f'within {model1} ({start1}-{end1}, score={score1:.1f})'
                        )

                # Check if hit1 is nested in hit2
                elif start2 <= start1 and end1 <= end2:
                    # hit1 is nested in hit2
                    if score2 > 0 and (score1 / score2) < min_score_ratio:
                        hits_to_remove.add(idx1)
                        logging.debug(
                            f'Removing nested hit {model1} ({start1}-{end1}, score={score1:.1f}) '
                            f'within {model2} ({start2}-{end2}, score={score2:.1f})'
                        )

    # Remove marked hits
    if hits_to_remove:
        logging.info(f'Removed {len(hits_to_remove)} nested weak hits between paired features')
        hit_table = hit_table.drop(index=list(hits_to_remove))

    # Clean up temporary columns
    hit_table = hit_table.drop(columns=['hitStart_int', 'hitEnd_int', 'score_float'])

    return hit_table.reset_index(drop=True)


# -----------------------------------------------------------------------------
# Output Functions
# -----------------------------------------------------------------------------


def write_hits_table(hit_table: pd.DataFrame, output_file: Path) -> None:
    """
    Write hit table to BLAST-compatible tabular format.

    Parameters
    ----------
    hit_table : pandas.DataFrame
        Hit table to write.
    output_file : Path
        Path to output file.
    """
    # Write header
    with open(output_file, 'w') as f:
        f.write('# TIRmite ensemble search output\n')
        f.write('# Format: model target hitStart hitEnd strand evalue score bias hmmStart hmmEnd\n')

    # Write data
    hit_table.to_csv(output_file, sep='\t', index=False, header=False, mode='a')

    logging.info(f'Wrote {len(hit_table)} hits to {output_file}')


# -----------------------------------------------------------------------------
# Search Execution Functions
# -----------------------------------------------------------------------------


def run_nhmmer_search(
    hmm_file: Path,
    target_file: Path,
    output_dir: Path,
    evalue: Optional[float] = None,
    threads: int = 1,
) -> Path:
    """
    Run nhmmer search with an HMM model against a target genome.

    Parameters
    ----------
    hmm_file : Path
        Path to HMM model file.
    target_file : Path
        Path to target genome FASTA file.
    output_dir : Path
        Directory to write output files.
    evalue : float, optional
        E-value threshold for reporting hits.
    threads : int, default 1
        Number of CPU threads to use.

    Returns
    -------
    Path
        Path to nhmmer output file.

    Raises
    ------
    EnsembleSearchError
        If nhmmer is not available or search fails.
    """
    if not shutil.which('nhmmer'):
        raise EnsembleSearchError('nhmmer not found in PATH. Please install HMMER.')

    command, results_dir = build_nhmmer_command(
        model_path=hmm_file,
        genome_path=target_file,
        output_dir=output_dir,
        evalue=evalue,
        cores=threads,
    )

    logging.info(f'Running nhmmer search: {hmm_file.name} vs {target_file.name}')

    try:
        result = run_command(command, verbose=True)
        if result.returncode != 0:
            raise EnsembleSearchError(f'nhmmer search failed: {result.stderr}')

        # Find output file
        output_file = results_dir / f'{hmm_file.stem}.out'
        if not output_file.exists():
            raise EnsembleSearchError(f'nhmmer output file not created: {output_file}')

        return output_file

    except Exception as e:
        raise EnsembleSearchError(f'nhmmer search failed: {e}') from e


def run_blastn_search(
    query_file: Path,
    target_file: Path,
    output_file: Path,
    evalue: float = 1e-3,
    identity: float = 60.0,
    threads: int = 1,
) -> Path:
    """
    Run BLAST search with query sequences against a target genome.

    Parameters
    ----------
    query_file : Path
        Path to query sequences FASTA file.
    target_file : Path
        Path to target genome FASTA file.
    output_file : Path
        Path to output file.
    evalue : float, default 1e-3
        E-value threshold for BLAST.
    identity : float, default 60.0
        Minimum percent identity threshold.
    threads : int, default 1
        Number of CPU threads to use.

    Returns
    -------
    Path
        Path to BLAST output file.

    Raises
    ------
    EnsembleSearchError
        If BLAST fails.
    """
    logging.info(f'Running BLAST search: {query_file.name} vs {target_file.name}')

    try:
        run_blastn(
            query=query_file,
            subject=target_file,
            output=output_file,
            word_size=4,
            perc_identity=identity,
            num_threads=threads,
            additional_args=['-evalue', str(evalue), '-max_target_seqs', '10000'],
        )

        return output_file

    except BlastError as e:
        raise EnsembleSearchError(f'BLAST search failed: {e}') from e


# -----------------------------------------------------------------------------
# Argument Validation
# -----------------------------------------------------------------------------


def validate_evalue(value: str) -> float:
    """Validate e-value argument."""
    try:
        fvalue = float(value)
        if fvalue <= 0:
            raise argparse.ArgumentTypeError(f'E-value must be positive: {value}')
        return fvalue
    except ValueError as err:
        raise argparse.ArgumentTypeError(f'Invalid e-value: {value}') from err


def validate_identity(value: str) -> float:
    """Validate identity percentage argument."""
    try:
        fvalue = float(value)
        if not 0 <= fvalue <= 100:
            raise argparse.ArgumentTypeError(
                f'Identity must be between 0 and 100: {value}'
            )
        return fvalue
    except ValueError as err:
        raise argparse.ArgumentTypeError(f'Invalid identity value: {value}') from err


def validate_threads(value: str) -> int:
    """Validate threads argument."""
    try:
        ivalue = int(value)
        if ivalue < 1:
            raise argparse.ArgumentTypeError(
                f'Threads must be at least 1: {value}'
            )
        return ivalue
    except ValueError as err:
        raise argparse.ArgumentTypeError(f'Invalid threads value: {value}') from err


# -----------------------------------------------------------------------------
# CLI Parser
# -----------------------------------------------------------------------------


def add_search_parser(subparsers: Any) -> argparse.ArgumentParser:
    """
    Add search subcommand parser.

    Parameters
    ----------
    subparsers : argparse._SubParsersAction
        Subparser object to add search command to.

    Returns
    -------
    argparse.ArgumentParser
        The configured search subcommand parser.
    """
    parser = cast(
        argparse.ArgumentParser,
        subparsers.add_parser(
            'search',
            help='Ensemble search: run searches and merge clustered hits',
            description='Run BLAST/nhmmer searches or load precomputed results, '
            'optionally merge overlapping hits from clustered features, '
            'and output merged hits for pairing.',
        ),
    )
    _configure_search_parser(parser)
    return parser


def _configure_search_parser(parser: argparse.ArgumentParser) -> None:
    """Configure parser with search command arguments."""
    # Input mode: either run searches or load precomputed results
    input_group = parser.add_argument_group('Input Options')

    input_group.add_argument(
        '--fasta',
        type=Path,
        nargs='+',
        help='FASTA file(s) containing query sequences for BLAST search',
    )
    input_group.add_argument(
        '--hmm',
        type=Path,
        nargs='+',
        help='HMM model file(s) for nhmmer search',
    )
    input_group.add_argument(
        '--genome',
        type=Path,
        help='Target genome FASTA file (required if running searches)',
    )
    input_group.add_argument(
        '--blast-results',
        type=Path,
        nargs='+',
        help='Precomputed BLAST tabular output file(s)',
    )
    input_group.add_argument(
        '--nhmmer-results',
        type=Path,
        nargs='+',
        help='Precomputed nhmmer tabular output file(s)',
    )

    # Clustering options
    cluster_group = parser.add_argument_group('Clustering Options')
    cluster_group.add_argument(
        '--cluster-map',
        type=Path,
        help='Tab-delimited file mapping cluster names to component features. '
        'Format: cluster_name<TAB>component1<TAB>component2...',
    )
    cluster_group.add_argument(
        '--pairing-map',
        type=Path,
        help='Tab-delimited file linking left and right terminal features for nested hit removal. '
        'Format: left_feature<TAB>right_feature',
    )

    # Filter options
    filter_group = parser.add_argument_group('Filter Options')
    filter_group.add_argument(
        '--max-evalue',
        type=validate_evalue,
        default=1e-3,
        help='Maximum e-value threshold for filtering hits (default: 1e-3)',
    )
    filter_group.add_argument(
        '--min-identity',
        type=validate_identity,
        default=60.0,
        help='Minimum percent identity for BLAST searches (default: 60.0)',
    )

    # Output options
    output_group = parser.add_argument_group('Output Options')
    output_group.add_argument(
        '--outdir',
        type=Path,
        default=Path.cwd(),
        help='Output directory (default: current directory)',
    )
    output_group.add_argument(
        '--prefix',
        type=str,
        default='ensemble_search',
        help='Prefix for output files (default: ensemble_search)',
    )

    # Processing options
    proc_group = parser.add_argument_group('Processing Options')
    proc_group.add_argument(
        '--threads',
        type=validate_threads,
        default=1,
        help='Number of CPU threads for searches (default: 1)',
    )
    proc_group.add_argument(
        '--loglevel',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='Set logging level (default: INFO)',
    )


def create_search_parser() -> argparse.ArgumentParser:
    """
    Create standalone search parser for direct invocation.

    Returns
    -------
    argparse.ArgumentParser
        Configured search command parser.
    """
    parser = argparse.ArgumentParser(
        prog='tirmite search',
        description='Ensemble search for transposon terminal detection',
    )
    _configure_search_parser(parser)
    return parser


# -----------------------------------------------------------------------------
# Main Function
# -----------------------------------------------------------------------------


def main(args: Optional[argparse.Namespace] = None) -> int:
    """
    Main function for ensemble search workflow.

    Parameters
    ----------
    args : argparse.Namespace, optional
        Parsed command-line arguments. If None, parses from sys.argv.

    Returns
    -------
    int
        Exit code (0 for success, 1 for error).
    """
    if args is None:
        parser = create_search_parser()
        args = parser.parse_args()

    assert args is not None

    # Set up logging
    init_logging(loglevel=args.loglevel)

    try:
        # Validate inputs
        has_search_inputs = args.fasta or args.hmm
        has_precomputed = args.blast_results or args.nhmmer_results

        if not has_search_inputs and not has_precomputed:
            raise EnsembleSearchError(
                'Must provide either search inputs (--fasta/--hmm + --genome) '
                'or precomputed results (--blast-results/--nhmmer-results)'
            )

        if has_search_inputs and not args.genome:
            raise EnsembleSearchError(
                '--genome is required when running searches with --fasta or --hmm'
            )

        # Create output directory
        args.outdir.mkdir(parents=True, exist_ok=True)

        # Collect result files
        blast_files: List[Path] = list(args.blast_results) if args.blast_results else []
        nhmmer_files: List[Path] = list(args.nhmmer_results) if args.nhmmer_results else []

        # Run searches if needed
        if has_search_inputs:
            with temporary_directory(prefix='tirmite_search_') as temp_dir:
                temp_path = Path(temp_dir)

                # Run BLAST searches
                if args.fasta:
                    for fasta_file in args.fasta:
                        if not fasta_file.exists():
                            logging.warning(f'FASTA file not found: {fasta_file}')
                            continue

                        output_file = temp_path / f'{fasta_file.stem}_blast.tab'
                        result_file = run_blastn_search(
                            query_file=fasta_file,
                            target_file=args.genome,
                            output_file=output_file,
                            evalue=args.max_evalue,
                            identity=args.min_identity,
                            threads=args.threads,
                        )
                        blast_files.append(result_file)

                # Run nhmmer searches
                if args.hmm:
                    for hmm_file in args.hmm:
                        if not hmm_file.exists():
                            logging.warning(f'HMM file not found: {hmm_file}')
                            continue

                        result_file = run_nhmmer_search(
                            hmm_file=hmm_file,
                            target_file=args.genome,
                            output_dir=temp_path,
                            evalue=args.max_evalue,
                            threads=args.threads,
                        )
                        nhmmer_files.append(result_file)

                # Load and process hits within temp context
                hit_table = _process_hits(args, blast_files, nhmmer_files)

        else:
            # Just load precomputed results
            hit_table = _process_hits(args, blast_files, nhmmer_files)

        # Write final output
        output_file = args.outdir / f'{args.prefix}_hits.tab'
        write_hits_table(hit_table, output_file)

        logging.info('Ensemble search completed successfully')
        return 0

    except EnsembleSearchError as e:
        logging.error(f'Ensemble search failed: {e}')
        return 1
    except Exception as e:
        logging.error(f'Unexpected error: {e}')
        return 1


def _process_hits(
    args: argparse.Namespace,
    blast_files: List[Path],
    nhmmer_files: List[Path],
) -> pd.DataFrame:
    """
    Process loaded hits: filter, merge, and clean.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments.
    blast_files : list of Path
        BLAST result files to load.
    nhmmer_files : list of Path
        nhmmer result files to load.

    Returns
    -------
    pandas.DataFrame
        Processed hit table.
    """
    # Load hits
    hit_table = load_hits_from_files(
        blast_files=blast_files if blast_files else None,
        nhmmer_files=nhmmer_files if nhmmer_files else None,
    )

    if hit_table.empty:
        logging.warning('No hits loaded')
        return hit_table

    # Report initial statistics
    report_hit_statistics(hit_table, stage='(before filtering)')

    # Filter by e-value
    hit_table = filter_hits_by_evalue(hit_table, args.max_evalue)

    # Report post-filter statistics
    report_hit_statistics(hit_table, stage='(after filtering)')

    # Load and apply cluster mapping if provided
    if args.cluster_map:
        if not args.cluster_map.exists():
            raise EnsembleSearchError(
                f'Cluster mapping file not found: {args.cluster_map}'
            )

        cluster_map = parse_cluster_mapping(args.cluster_map)

        if cluster_map:
            # Get available feature names
            available_features = set(hit_table['model'].unique())

            # Validate cluster mapping
            is_valid, warnings = validate_cluster_mapping(cluster_map, available_features)
            for warning in warnings:
                logging.warning(warning)

            if not is_valid:
                raise EnsembleSearchError(
                    'Cluster mapping validation failed. Check warnings above.'
                )

            # Check for cross-cluster overlaps before merging
            check_cross_cluster_overlaps(hit_table, cluster_map)

            # Merge overlapping hits within clusters
            hit_table = merge_overlapping_cluster_hits(hit_table, cluster_map)

            # Report post-merge statistics
            report_hit_statistics(hit_table, stage='(after merging)')

    # Remove nested weak hits if pairing map provided
    if args.pairing_map:
        if not args.pairing_map.exists():
            raise EnsembleSearchError(
                f'Pairing map file not found: {args.pairing_map}'
            )

        pairing_map = parse_pairing_map(args.pairing_map)

        if pairing_map:
            hit_table = remove_nested_paired_hits(hit_table, pairing_map)

            # Report final statistics
            report_hit_statistics(hit_table, stage='(after nested hit removal)')

    return hit_table


if __name__ == '__main__':
    import sys

    sys.exit(main())
