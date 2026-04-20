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
from dataclasses import dataclass, field
import logging
from pathlib import Path
import shutil
from typing import Any, Dict, List, Optional, Set, Tuple, cast

from Bio import SeqIO  # type: ignore[import-not-found]
import pandas as pd  # type: ignore[import-untyped]

from tirmite.runners.hmmer_wrappers import build_nhmmer_command
from tirmite.runners.runBlastn import BlastError, run_blastn
from tirmite.runners.wrapping import run_command
import tirmite.tirmitetools as tirmite
from tirmite.utils.logs import init_logging
from tirmite.utils.utils import prepare_genome_file, temporary_directory


class EnsembleSearchError(Exception):
    """Custom exception for ensemble search errors."""

    pass


# -----------------------------------------------------------------------------
# Filter Summary Dataclass
# -----------------------------------------------------------------------------


@dataclass
class SearchFilterSummary:
    """Accumulates hit-filtering statistics across all pairing map pipeline steps.

    Attributes
    ----------
    excluded_not_in_map : dict
        Hits excluded at Step 0 (model not in pairing map).
        Mapping of ``{model_name: hit_count}``.
    nested_removed : dict
        Hits removed at Step 1 (nested within a stronger partner hit).
        Mapping of ``{removed_model: {container_model: count}}``.
    cross_model_removed : dict
        Hits removed at Step 2 (weaker overlapping cross-model hit).
        Mapping of ``{(removed_model, better_model): count}``.
    """

    excluded_not_in_map: Dict[str, int] = field(default_factory=dict)
    nested_removed: Dict[str, Dict[str, int]] = field(default_factory=dict)
    cross_model_removed: Dict[Tuple[str, str], int] = field(default_factory=dict)


# -----------------------------------------------------------------------------
# Query Length Functions
# -----------------------------------------------------------------------------


def get_fasta_sequence_lengths(fasta_file: Path) -> Dict[str, int]:
    """
    Get sequence lengths from a FASTA file.

    Parameters
    ----------
    fasta_file : Path
        Path to FASTA file.

    Returns
    -------
    dict
        Dictionary mapping sequence IDs to their lengths.

    Raises
    ------
    EnsembleSearchError
        If file cannot be read.
    """
    lengths: Dict[str, int] = {}

    try:
        for record in SeqIO.parse(str(fasta_file), 'fasta'):
            lengths[record.id] = len(record.seq)
        logging.debug(f'Loaded {len(lengths)} sequence lengths from {fasta_file.name}')
    except Exception as e:
        raise EnsembleSearchError(f'Failed to read FASTA file {fasta_file}: {e}') from e

    return lengths


def get_hmm_model_lengths(hmm_file: Path) -> Dict[str, int]:
    """
    Extract HMM model lengths from HMM file by parsing LENG fields.

    Parameters
    ----------
    hmm_file : Path
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
    model_lengths: Dict[str, int] = {}

    try:
        with open(hmm_file, 'r') as f:
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
        logging.debug(f'Loaded {len(model_lengths)} model lengths from {hmm_file.name}')

    except Exception as e:
        logging.error(f'Error reading HMM file {hmm_file}: {e}')

    return model_lengths


def load_lengths_file(lengths_file: Path) -> Dict[str, int]:
    """
    Load model/query lengths from tab-delimited text file.

    Parameters
    ----------
    lengths_file : Path
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
    model_lengths: Dict[str, int] = {}

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
        raise EnsembleSearchError(
            f'Error reading model lengths file {lengths_file}: {e}'
        ) from e

    return model_lengths


def parse_genome_list(genome_list_file: Path) -> List[Path]:
    """
    Parse a file containing a list of genome paths.

    Parameters
    ----------
    genome_list_file : Path
        Path to file with one genome path per line.

    Returns
    -------
    list of Path
        List of genome file paths.

    Raises
    ------
    EnsembleSearchError
        If file cannot be read.
    """
    genomes: List[Path] = []

    try:
        with open(genome_list_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                genome_path = Path(line)
                if not genome_path.exists():
                    logging.warning(
                        f'Genome file not found (line {line_num}): {genome_path}'
                    )
                    continue

                genomes.append(genome_path)

        logging.info(f'Loaded {len(genomes)} genome paths from {genome_list_file.name}')

    except Exception as e:
        raise EnsembleSearchError(
            f'Failed to parse genome list file {genome_list_file}: {e}'
        ) from e

    return genomes


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
        warnings.append(f'Cluster components not found in hit data: {sorted(missing)}')

    return is_valid, warnings


def build_component_to_cluster_map(cluster_map: Dict[str, List[str]]) -> Dict[str, str]:
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
    failed_files: List[str] = []
    loaded_files: List[str] = []

    # Load BLAST files
    if blast_files:
        for blast_file in blast_files:
            if not blast_file.exists():
                logging.warning(f'BLAST file not found, skipping: {blast_file}')
                failed_files.append(f'{blast_file} (not found)')
                continue

            logging.info(f'Loading BLAST hits from: {blast_file}')
            try:
                hit_table = tirmite.import_blast(str(blast_file), hitTable=hit_table)
                loaded_files.append(str(blast_file))
            except Exception as e:
                logging.error(f'Failed to load BLAST file {blast_file}: {e}')
                failed_files.append(f'{blast_file} ({e})')
                continue

    # Load nhmmer files
    if nhmmer_files:
        for nhmmer_file in nhmmer_files:
            if not nhmmer_file.exists():
                logging.warning(f'nhmmer file not found, skipping: {nhmmer_file}')
                failed_files.append(f'{nhmmer_file} (not found)')
                continue

            logging.info(f'Loading nhmmer hits from: {nhmmer_file}')
            try:
                hit_table = tirmite.import_nhmmer(str(nhmmer_file), hitTable=hit_table)
                loaded_files.append(str(nhmmer_file))
            except Exception as e:
                logging.error(f'Failed to load nhmmer file {nhmmer_file}: {e}')
                failed_files.append(f'{nhmmer_file} ({e})')
                continue

    # Report summary of file loading
    if failed_files:
        logging.warning(
            f'Failed to load {len(failed_files)} file(s): {", ".join(failed_files[:5])}'
            + (f' and {len(failed_files) - 5} more' if len(failed_files) > 5 else '')
        )

    if loaded_files:
        logging.info(f'Successfully loaded {len(loaded_files)} file(s)')

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


def compute_outer_edge_offset(
    hmm_start: int,
    hmm_end: int,
    model_len: int,
    strand: str,
    terminus_type: str,
) -> int:
    """
    Compute the offset (in bases) between the hit alignment and the outer edge of the query model.

    The "outer edge" is the external end of the terminus model:
    - For a left terminus the outer edge faces the upstream/left genomic direction.
    - For a right terminus the outer edge faces the downstream/right genomic direction.

    The offset represents how many model positions are not covered between the
    alignment boundary and the outer edge.  A value of 0 means the alignment
    reaches the outer tip of the model; larger values mean the alignment stops
    progressively further from the outer edge.

    Parameters
    ----------
    hmm_start : int
        1-based start position of the alignment in the query/HMM model.
    hmm_end : int
        1-based end position of the alignment in the query/HMM model.
    model_len : int
        Total length of the query/HMM model.
    strand : str
        Strand of the hit on the target sequence: '+' or '-'.
    terminus_type : str
        Whether this is a 'left' or 'right' terminus.

    Returns
    -------
    int
        Number of unaligned model positions between the hit and the outer edge.
        Zero means the hit reaches the outer edge exactly.

    Notes
    -----
    Offset rules (matching ``compute_flank_coordinates`` in tirmitetools.py):

    Left terminus:
      - ``+`` strand: outer edge = model position 1  → offset = ``hmm_start - 1``
      - ``-`` strand: outer edge = model position ``model_len`` → offset = ``model_len - hmm_end``

    Right terminus:
      - ``+`` strand: outer edge = model position ``model_len`` → offset = ``model_len - hmm_end``
      - ``-`` strand: outer edge = model position 1 → offset = ``hmm_start - 1``
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
    query_lengths: Dict[str, int],
    max_offset: int,
    orientation: str = 'F,R',
    pairing_map: Optional[Dict[str, str]] = None,
) -> pd.DataFrame:
    """
    Filter hits to only those anchored within ``max_offset`` bases of the outer edge of the query model.

    For each hit the terminus type (left/right) is determined from the pairing map
    (asymmetric case) or from the hit strand combined with the orientation setting
    (symmetric F,R or R,F case).  The offset from the outer edge is then computed
    using :func:`compute_outer_edge_offset` and hits whose offset exceeds
    ``max_offset`` are removed.

    Hits are kept unchanged when:
    - The model length is not available in ``query_lengths``.
    - The terminus type cannot be determined (e.g. same-strand F,F or R,R symmetric
      pairing without a pairing map).

    Parameters
    ----------
    hit_table : pandas.DataFrame
        Hit table with columns: model, strand, hmmStart, hmmEnd.
    query_lengths : dict
        Mapping of model name to model length (number of positions).
    max_offset : int
        Maximum allowed offset (bases) from the outer edge of the model.
        Hits with offset > max_offset are removed.
    orientation : str, default 'F,R'
        Comma-separated orientation codes used to assign left/right terminus by
        strand for symmetric (single-model) cases.
        F = Forward (+), R = Reverse (-).  Examples: 'F,R', 'F,F', 'R,F', 'R,R'.
    pairing_map : dict, optional
        Mapping of left model name to right model name.  When provided, terminus
        type is determined by model name rather than strand.

    Returns
    -------
    pandas.DataFrame
        Filtered hit table containing only anchored hits.
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

    # Build reverse map: model_name -> 'left' or 'right'
    model_terminus: Dict[str, str] = {}
    if pairing_map:
        for left_model, right_model in pairing_map.items():
            model_terminus[left_model] = 'left'
            model_terminus[right_model] = 'right'

    kept: List[bool] = []
    skipped_no_terminus = 0
    removed = 0
    removed_per_model: Dict[str, int] = {}
    missing_lengths: Set[str] = set()

    for _, row in hit_table.iterrows():
        model = row['model']
        strand = row['strand']

        # Determine terminus type
        if model in model_terminus:
            terminus_type: Optional[str] = model_terminus[model]
        elif strands_differ:
            # Use strand to distinguish left from right
            if strand == left_strand:
                terminus_type = 'left'
            elif strand == right_strand:
                terminus_type = 'right'
            else:
                terminus_type = None
        else:
            # Same-strand orientation (F,F or R,R) without a pairing map –
            # cannot determine which end is "outer", check both ends
            terminus_type = None

        if terminus_type is None:
            # Same-strand symmetric (F,F or R,R) without a pairing map:
            # check both ends of the query model
            model_len = query_lengths.get(model)
            if model_len is None:
                missing_lengths.add(model)
                kept.append(True)
                continue

            try:
                hmm_start = int(row['hmmStart'])
                hmm_end = int(row['hmmEnd'])
            except (ValueError, TypeError):
                kept.append(True)
                continue

            offset_start = hmm_start - 1
            offset_end = model_len - hmm_end
            if offset_start <= max_offset and offset_end <= max_offset:
                kept.append(True)
            else:
                kept.append(False)
                removed += 1
                removed_per_model[model] = removed_per_model.get(model, 0) + 1
            continue

        model_len = query_lengths.get(model)
        if model_len is None:
            # Length required but not available – collect and report as error
            missing_lengths.add(model)
            kept.append(True)
            continue

        try:
            hmm_start = int(row['hmmStart'])
            hmm_end = int(row['hmmEnd'])
        except (ValueError, TypeError):
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

    # Raise immediately if any model lengths were needed but missing
    if missing_lengths:
        missing_list = ', '.join(sorted(missing_lengths))
        raise EnsembleSearchError(
            f'Anchor filter requires model lengths for {len(missing_lengths)} model(s) '
            f'that are not available: {missing_list}. '
            'Provide lengths via --fasta, --hmm, or --lengths-file.'
        )

    if skipped_no_terminus:
        logging.debug(
            f'Anchor filter: {skipped_no_terminus} hit(s) kept without checking – '
            'terminus type could not be determined (same-strand symmetric pairing requires '
            '--pairing-map to identify left/right models).'
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


def filter_hits_to_pairing_map_models(
    hit_table: pd.DataFrame,
    pairing_map: Dict[str, str],
    summary: Optional['SearchFilterSummary'] = None,
) -> pd.DataFrame:
    """
    Retain only hits from models listed in the pairing map.

    When a pairing map is provided, hits from models that are not part of any
    left/right pair are irrelevant for downstream pairing and are excluded from
    the output.

    Parameters
    ----------
    hit_table : pandas.DataFrame
        Hit table with a 'model' column.
    pairing_map : dict
        Dictionary mapping left feature names to right feature names.
    summary : SearchFilterSummary, optional
        If provided, per-model exclusion counts are stored in
        ``summary.excluded_not_in_map``.

    Returns
    -------
    pandas.DataFrame
        Hit table restricted to models in the pairing map.
    """
    if hit_table.empty or not pairing_map:
        return hit_table

    paired_models = set(pairing_map.keys()) | set(pairing_map.values())
    mask = hit_table['model'].isin(paired_models)
    n_before = len(hit_table)
    result = hit_table[mask].copy().reset_index(drop=True)
    n_removed = n_before - len(result)

    if n_removed:
        removed_models = sorted(set(hit_table.loc[~mask, 'model'].unique()))
        logging.info(
            f'Pairing map filter: removed {n_removed} hit(s) from '
            f'{len(removed_models)} model(s) not in the pairing map: '
            f'{", ".join(removed_models)}'
        )
        if summary is not None:
            for model in removed_models:
                count = int((~mask & (hit_table['model'] == model)).sum())
                summary.excluded_not_in_map[model] = count

    return result


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
        current_start: int = 0
        current_end: int = 0
        is_first_hit: bool = True

        for _, hit in group.iterrows():
            hit_start = int(hit['hitStart_int'])
            hit_end = int(hit['hitEnd_int'])

            if is_first_hit:
                # First hit
                current_start = hit_start
                current_end = hit_end
                current_hits = [hit]
                is_first_hit = False
            elif hit_start <= current_end + min_overlap:
                # Overlapping or adjacent - extend region
                current_end = max(current_end, hit_end)
                current_hits.append(hit)
            else:
                # Non-overlapping - save current merged hit and start new one
                merged_records.append(
                    _create_merged_hit(
                        current_hits, current_start, current_end, cluster
                    )
                )
                current_start = hit_start
                current_end = hit_end
                current_hits = [hit]

        # Don't forget the last merged region
        if current_hits:
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
                            f'Cross-cluster overlap detected: {target}:{start1}-{end1} ({hit1["cluster"]}) '
                            f'overlaps with {target}:{start2}-{end2} ({hit2["cluster"]})'
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
    summary: Optional['SearchFilterSummary'] = None,
) -> pd.DataFrame:
    """
    Remove weak hits nested within stronger hits from the same left/right pair.

    In asymmetric terminus models (i.e. separate left and right models), the two
    models sometimes share a small region of homology.  This means that the shorter
    (or weaker) model may produce hits that are entirely contained within hits of its
    paired model at the same genomic locus.  These nested "cross-hits" are artefacts
    of the shared homology rather than genuine detections by the nested model, and
    retaining them would generate spurious pair calls.

    This function compares all hits on the same target sequence and strand.  For each
    pair of hits where:

    - Both models appear in the ``pairing_map`` (either as left or right features), AND
    - The two models form a direct left/right pair in the ``pairing_map``, AND
    - One hit's genomic coordinates are completely contained within the other's (nesting),

    the nested hit is removed when its alignment score is less than
    ``min_score_ratio`` × the enclosing hit's score.  If the nested hit scores
    comparably (ratio ≥ ``min_score_ratio``), both hits are kept, because the nested
    model may still represent a genuine detection.

    **Handling models that appear in multiple pairs**

    Because ``pairing_map`` is a plain Python dictionary it maps each left-feature name
    to exactly one right-feature name.  If a single model name must be paired with
    multiple partners, add separate left–right rows to the pairing map file; only the
    last entry for a given left feature will be stored, so each left model effectively
    has one canonical right partner for this filter.  To support many-to-many pairing
    relationships, use :func:`filter_best_model_per_locus` after this step, which
    removes lower-quality overlapping cross-model hits without requiring a direct
    pairing relationship.

    Parameters
    ----------
    hit_table : pandas.DataFrame
        Hit table with 'model', 'target', 'hitStart', 'hitEnd', 'strand', 'score' columns.
    pairing_map : dict
        Dictionary mapping left feature names to right feature names.
        Only hits from models listed in this map (as keys or values) are evaluated.
    min_score_ratio : float, default 1.5
        A nested hit is removed when ``nested_score / enclosing_score < min_score_ratio``.
        Increase this value to be more aggressive (remove more nested hits); decrease
        it to be more conservative (keep more nested hits).
    summary : SearchFilterSummary, optional
        If provided, per-model nesting details are stored in
        ``summary.nested_removed`` as ``{removed_model: {container_model: count}}``.

    Returns
    -------
    pandas.DataFrame
        Hit table with nested weak hits removed.

    See Also
    --------
    filter_best_model_per_locus : Remove lower-quality overlapping hits across all
        models in the pairing map, regardless of direct pairing relationship.
    """
    if hit_table.empty or not pairing_map:
        return hit_table

    all_paired_features = set(pairing_map.keys()) | set(pairing_map.values())

    hit_table = hit_table.copy()
    hit_table['hitStart_int'] = hit_table['hitStart'].astype(int)
    hit_table['hitEnd_int'] = hit_table['hitEnd'].astype(int)
    hit_table['score_float'] = pd.to_numeric(
        hit_table['score'], errors='coerce'
    ).fillna(0)

    hits_to_remove: Set[int] = set()
    # Map removed hit index → (removed_model, container_model) for summary reporting
    removed_hit_containers: Dict[int, Tuple[str, str]] = {}

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
                    model1 in pairing_map and pairing_map[model1] == model2
                ) or (model2 in pairing_map and pairing_map[model2] == model1)

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
                        removed_hit_containers[idx2] = (model2, model1)
                        logging.debug(
                            f'Removing nested hit {model2} ({start2}-{end2}, score={score2:.1f}) '
                            f'within {model1} ({start1}-{end1}, score={score1:.1f})'
                        )

                # Check if hit1 is nested in hit2
                elif start2 <= start1 and end1 <= end2:
                    # hit1 is nested in hit2
                    if score2 > 0 and (score1 / score2) < min_score_ratio:
                        hits_to_remove.add(idx1)
                        removed_hit_containers[idx1] = (model1, model2)
                        logging.debug(
                            f'Removing nested hit {model1} ({start1}-{end1}, score={score1:.1f}) '
                            f'within {model2} ({start2}-{end2}, score={score2:.1f})'
                        )

    # Remove marked hits
    if hits_to_remove:
        logging.info(
            f'Removed {len(hits_to_remove)} nested weak hits between paired features'
        )
        hit_table = hit_table.drop(index=list(hits_to_remove))

        if summary is not None:
            for _idx, (
                removed_model,
                container_model,
            ) in removed_hit_containers.items():
                per_container = summary.nested_removed.setdefault(removed_model, {})
                per_container[container_model] = (
                    per_container.get(container_model, 0) + 1
                )

    # Clean up temporary columns
    hit_table = hit_table.drop(columns=['hitStart_int', 'hitEnd_int', 'score_float'])

    return hit_table.reset_index(drop=True)


def filter_best_model_per_locus(
    hit_table: pd.DataFrame,
    pairing_map: Dict[str, str],
    min_overlap: int = 1,
    min_score_ratio: float = 1.5,
    summary: Optional['SearchFilterSummary'] = None,
) -> pd.DataFrame:
    """
    Retain only the best-scoring model's hits at each overlapping genomic locus.

    Closely related transposon families often produce HMM or BLAST models that match
    their own target family well but also produce weaker hits against related families.
    When many models from different element families are searched simultaneously and
    their hits are assigned to pairs via ``pairing_map``, overlapping hits at the same
    locus from competing (non-paired) models must be resolved so that only the
    best-matching model is represented at each locus.

    **Algorithm**

    For each target sequence and strand, all pairwise combinations of hits from
    *different* models are examined.  When two hits overlap by at least ``min_overlap``
    bases and the better-scoring hit's score is at least ``min_score_ratio`` times the
    weaker hit's score, the weaker hit is marked for removal.  Hits from the same
    model are never compared against each other, so multiple non-overlapping hits from
    the best model at distinct loci are always retained.

    Only models listed in ``pairing_map`` (as either left or right features) are
    subject to this filter; hits from models not in the pairing map pass through
    unchanged.

    This step is complementary to :func:`remove_nested_paired_hits`:

    - :func:`remove_nested_paired_hits` handles hits where one is *completely contained*
      within another from its direct paired partner (nesting within a left/right pair).
    - :func:`filter_best_model_per_locus` handles hits that merely *overlap* (including
      but not limited to nesting) between *any* two models in the pairing map,
      regardless of whether those models form a direct left/right pair.

    Parameters
    ----------
    hit_table : pandas.DataFrame
        Hit table with 'model', 'target', 'hitStart', 'hitEnd', 'strand', 'score' columns.
    pairing_map : dict
        Dictionary mapping left feature names to right feature names.
        All model names that appear as keys or values are eligible for filtering.
    min_overlap : int, default 1
        Minimum overlap in base pairs required to trigger cross-model comparison.
        Hits that share fewer than this many bases are considered non-overlapping
        and are always retained.
    min_score_ratio : float, default 1.5
        A weaker overlapping hit is removed when
        ``better_score / weaker_score >= min_score_ratio``.
        When the ratio is below this threshold, both hits are kept because neither
        model dominates enough to confidently discard the other.
    summary : SearchFilterSummary, optional
        If provided, per-model-pair removal counts are stored in
        ``summary.cross_model_removed`` as ``{(removed_model, better_model): count}``.

    Returns
    -------
    pandas.DataFrame
        Hit table with lower-quality overlapping cross-model hits removed.

    See Also
    --------
    remove_nested_paired_hits : Remove hits nested within stronger hits from the
        same direct left/right pair.
    """
    if hit_table.empty or not pairing_map:
        return hit_table

    all_paired_features = set(pairing_map.keys()) | set(pairing_map.values())

    hit_table = hit_table.copy()
    hit_table['hitStart_int'] = hit_table['hitStart'].astype(int)
    hit_table['hitEnd_int'] = hit_table['hitEnd'].astype(int)
    hit_table['score_float'] = pd.to_numeric(
        hit_table['score'], errors='coerce'
    ).fillna(0)

    hits_to_remove: Set[int] = set()
    removed_per_model: Dict[str, int] = {}
    # Track (removed_model, better_model) → count for summary reporting
    cross_model_removal_counts: Dict[Tuple[str, str], int] = {}

    # Check each target/strand combination
    for (_target, _strand), group in hit_table.groupby(['target', 'strand']):
        if len(group) < 2:
            continue

        group_indices = list(group.index)

        for i, idx1 in enumerate(group_indices):
            if idx1 in hits_to_remove:
                continue

            hit1 = hit_table.loc[idx1]
            model1 = hit1['model']

            if model1 not in all_paired_features:
                continue

            start1 = int(hit1['hitStart_int'])
            end1 = int(hit1['hitEnd_int'])
            score1 = float(hit1['score_float'])

            for idx2 in group_indices[i + 1 :]:
                if idx2 in hits_to_remove:
                    continue

                hit2 = hit_table.loc[idx2]
                model2 = hit2['model']

                if model2 not in all_paired_features:
                    continue

                # Only compare hits from different models
                if model1 == model2:
                    continue

                start2 = int(hit2['hitStart_int'])
                end2 = int(hit2['hitEnd_int'])
                score2 = float(hit2['score_float'])

                # Compute overlap
                overlap = min(end1, end2) - max(start1, start2)
                if overlap < min_overlap:
                    continue

                # Remove the weaker hit when the score ratio is decisive
                if score1 > 0 and score2 > 0:
                    if score1 / score2 >= min_score_ratio:
                        hits_to_remove.add(idx2)
                        removed_per_model[model2] = removed_per_model.get(model2, 0) + 1
                        pair_key = (model2, model1)
                        cross_model_removal_counts[pair_key] = (
                            cross_model_removal_counts.get(pair_key, 0) + 1
                        )
                        logging.debug(
                            f'Removing overlapping cross-model hit {model2} '
                            f'({start2}-{end2}, score={score2:.1f}) in favour of '
                            f'{model1} ({start1}-{end1}, score={score1:.1f})'
                        )
                    elif score2 / score1 >= min_score_ratio:
                        hits_to_remove.add(idx1)
                        removed_per_model[model1] = removed_per_model.get(model1, 0) + 1
                        pair_key = (model1, model2)
                        cross_model_removal_counts[pair_key] = (
                            cross_model_removal_counts.get(pair_key, 0) + 1
                        )
                        logging.debug(
                            f'Removing overlapping cross-model hit {model1} '
                            f'({start1}-{end1}, score={score1:.1f}) in favour of '
                            f'{model2} ({start2}-{end2}, score={score2:.1f})'
                        )

    if hits_to_remove:
        logging.info(
            f'Removed {len(hits_to_remove)} overlapping lower-quality cross-model hits'
        )
        for model_name, count in sorted(removed_per_model.items()):
            logging.info(
                f'  {model_name}: {count} hit(s) removed by cross-model filter'
            )
        hit_table = hit_table.drop(index=list(hits_to_remove))

        if summary is not None:
            for pair_key, count in cross_model_removal_counts.items():
                summary.cross_model_removed[pair_key] = (
                    summary.cross_model_removed.get(pair_key, 0) + count
                )

    hit_table = hit_table.drop(columns=['hitStart_int', 'hitEnd_int', 'score_float'])

    return hit_table.reset_index(drop=True)


def log_filter_summary(summary: 'SearchFilterSummary') -> None:
    """
    Emit a structured INFO-level summary of all pairing map filtering steps.

    The report covers three filtering stages applied when a ``--pairing-map``
    is provided:

    * **Step 0** – Models excluded because they are not listed in the pairing map.
    * **Step 1** – Nested hit removal within direct left/right pairs.
    * **Step 2** – Cross-model overlap filtering (weaker hit at shared locus).

    Parameters
    ----------
    summary : SearchFilterSummary
        Accumulated statistics from the pairing map filtering pipeline.
    """
    lines = [
        '',
        '=' * 60,
        'Pairing Map Filter Summary',
        '=' * 60,
    ]

    # Step 0: Models not in pairing map
    if summary.excluded_not_in_map:
        total_excluded = sum(summary.excluded_not_in_map.values())
        lines.append(
            f'Step 0 — Excluded {total_excluded} hit(s) from models not in the pairing map:'
        )
        for model in sorted(summary.excluded_not_in_map):
            lines.append(
                f'  {model}: {summary.excluded_not_in_map[model]} hit(s) excluded'
            )
    else:
        lines.append('Step 0 — No hits excluded (all models present in pairing map)')

    # Step 1: Nested hit removal
    if summary.nested_removed:
        total_nested = sum(sum(c.values()) for c in summary.nested_removed.values())
        lines.append(
            f'Step 1 — Removed {total_nested} nested hit(s) within direct left/right pairs:'
        )
        for removed_model in sorted(summary.nested_removed):
            per_container = summary.nested_removed[removed_model]
            total_for_model = sum(per_container.values())
            container_str = ', '.join(
                f'{c} ({n})' for c, n in sorted(per_container.items())
            )
            lines.append(
                f'  {removed_model}: {total_for_model} hit(s) nested within [{container_str}]'
            )
    else:
        lines.append('Step 1 — No nested hits removed')

    # Step 2: Cross-model overlap filtering
    if summary.cross_model_removed:
        total_cross = sum(summary.cross_model_removed.values())
        lines.append(
            f'Step 2 — Removed {total_cross} cross-model hit(s) at shared loci:'
        )
        for removed_model, better_model in sorted(summary.cross_model_removed):
            count = summary.cross_model_removed[(removed_model, better_model)]
            lines.append(f'  {removed_model} → {better_model}: {count} hit(s) removed')
    else:
        lines.append('Step 2 — No cross-model overlapping hits removed')

    lines.append('=' * 60)

    for line in lines:
        logging.info(line)


# -----------------------------------------------------------------------------
# Output Functions
# -----------------------------------------------------------------------------


def write_hits_table(hit_table: pd.DataFrame, output_file: Path) -> None:
    """
    Write hit table to BLAST-compatible tabular format (outfmt 6).

    Parameters
    ----------
    hit_table : pandas.DataFrame
        Hit table to write with columns: model, target, hitStart, hitEnd, strand,
        evalue, score, bias, hmmStart, hmmEnd.
    output_file : Path
        Path to output file.

    Notes
    -----
    Output format is BLAST tabular format 6:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

    The internal hit table format is mapped to BLAST format as follows:
    - qseqid = model
    - sseqid = target
    - pident = 100.0 (not available, set to 100)
    - length = abs(hitEnd - hitStart) + 1
    - mismatch = 0 (not available)
    - gapopen = 0 (not available)
    - qstart = hmmStart
    - qend = hmmEnd
    - sstart = hitStart (for + strand) or hitEnd (for - strand)
    - send = hitEnd (for + strand) or hitStart (for - strand)
    - evalue = evalue
    - bitscore = score
    """
    # Write header
    with open(output_file, 'w') as f:
        f.write('# TIRmite ensemble search output - BLAST tabular format 6\n')
        f.write(
            '# Fields: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\n'
        )

    # Convert hit table to BLAST format 6
    blast_records = []
    for _, row in hit_table.iterrows():
        hit_start = int(row['hitStart'])
        hit_end = int(row['hitEnd'])
        length = abs(hit_end - hit_start) + 1

        # Handle strand - BLAST uses swapped coordinates for reverse strand
        if row['strand'] == '-':
            sstart = hit_end
            send = hit_start
        else:
            sstart = hit_start
            send = hit_end

        blast_records.append(
            {
                'qseqid': row['model'],
                'sseqid': row['target'],
                'pident': 100.0,  # Not available from merged hits
                'length': length,
                'mismatch': 0,
                'gapopen': 0,
                'qstart': row['hmmStart'],
                'qend': row['hmmEnd'],
                'sstart': sstart,
                'send': send,
                'evalue': row['evalue'],
                'bitscore': row['score'],
            }
        )

    # Write BLAST format output
    blast_df = pd.DataFrame(blast_records)
    if not blast_df.empty:
        blast_df.to_csv(output_file, sep='\t', index=False, header=False, mode='a')

    logging.info(f'Wrote {len(hit_table)} hits to {output_file}')


def validate_split_paired_output(pairing_map: Dict[str, str]) -> None:
    """
    Validate that no model name appears in both the left and right columns of the pairing map.

    Parameters
    ----------
    pairing_map : dict
        Dictionary mapping left feature names to right feature names.

    Raises
    ------
    EnsembleSearchError
        If any model name appears in both the left and right columns.
    """
    left_models = set(pairing_map.keys())
    right_models = set(pairing_map.values())
    overlap = left_models & right_models

    if overlap:
        overlap_list = ', '.join(sorted(overlap))
        raise EnsembleSearchError(
            f'Cannot use --split-paired-output: model(s) {overlap_list} appear in both '
            'the left and right columns of the pairing map. Each model must be '
            'exclusively left or right when splitting output.'
        )


def write_split_hits(
    hit_table: pd.DataFrame,
    pairing_map: Dict[str, str],
    outdir: Path,
    prefix: str,
) -> Tuple[Path, Path]:
    """
    Write left and right model hits to separate output files based on the pairing map.

    Parameters
    ----------
    hit_table : pandas.DataFrame
        Hit table to split and write.
    pairing_map : dict
        Dictionary mapping left feature names to right feature names.
    outdir : Path
        Output directory.
    prefix : str
        Prefix for output file names.

    Returns
    -------
    tuple of (Path, Path)
        Paths to the left and right output files.
    """
    left_models = set(pairing_map.keys())
    right_models = set(pairing_map.values())

    left_file = outdir / f'{prefix}_left_hits.tab'
    right_file = outdir / f'{prefix}_right_hits.tab'

    if hit_table.empty:
        write_hits_table(hit_table, left_file)
        write_hits_table(hit_table, right_file)
        logging.info(
            f'Split output: 0 left hits -> {left_file}, 0 right hits -> {right_file}'
        )
        return left_file, right_file

    left_hits = hit_table[hit_table['model'].isin(left_models)]
    right_hits = hit_table[hit_table['model'].isin(right_models)]

    # Hits from models not in the pairing map are written to neither file
    unassigned_hits = hit_table[~hit_table['model'].isin(left_models | right_models)]

    write_hits_table(left_hits, left_file)
    write_hits_table(right_hits, right_file)

    if not unassigned_hits.empty:
        logging.warning(
            f'{len(unassigned_hits)} hit(s) from model(s) not in the pairing map were not '
            'written to either the left or right output file: '
            f'{", ".join(sorted(unassigned_hits["model"].unique()))}'
        )

    logging.info(
        f'Split output: {len(left_hits)} left hits -> {left_file}, '
        f'{len(right_hits)} right hits -> {right_file}'
    )

    return left_file, right_file


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

    # Log the full command
    cmd_str = ' '.join(command)
    logging.info(f'Running nhmmer search: {hmm_file.name} vs {target_file.name}')
    logging.info(f'nhmmer command: {cmd_str}')

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
    word_size: int = 4,
    blast_timeout: Optional[int] = None,
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
    word_size : int, default 4
        Word size for initial matches.
    blast_timeout : int or None, default None
        Maximum number of seconds to wait for BLAST to complete.
        If None, BLAST is allowed to run indefinitely.

    Returns
    -------
    Path
        Path to BLAST output file.

    Raises
    ------
    EnsembleSearchError
        If BLAST fails.
    """
    # Standard BLAST format 6 columns expected by import_blast:
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    blast_outfmt = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'

    logging.info(f'Running BLAST search: {query_file.name} vs {target_file.name}')

    try:
        run_blastn(
            query=query_file,
            subject=target_file,
            output=output_file,
            word_size=word_size,
            perc_identity=identity,
            num_threads=threads,
            outfmt=blast_outfmt,
            additional_args=['-evalue', str(evalue), '-max_target_seqs', '10000'],
            timeout=blast_timeout,
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
            raise argparse.ArgumentTypeError(f'Threads must be at least 1: {value}')
        return ivalue
    except ValueError as err:
        raise argparse.ArgumentTypeError(f'Invalid threads value: {value}') from err


def validate_word_size(value: str) -> int:
    """Validate word_size argument."""
    try:
        ivalue = int(value)
        if ivalue < 4:
            raise argparse.ArgumentTypeError(f'Word size must be at least 4: {value}')
        return ivalue
    except ValueError as err:
        raise argparse.ArgumentTypeError(f'Invalid word_size value: {value}') from err


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
        help='Target genome FASTA file (required if running nhmmer searches with --hmm)',
    )
    input_group.add_argument(
        '--genome-list',
        type=Path,
        help='File containing list of genome FASTA paths (one per line, supports gzipped files)',
    )
    input_group.add_argument(
        '--blast-db',
        type=Path,
        help='Pre-built BLAST database prefix (alternative to --genome for BLAST searches)',
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
    input_group.add_argument(
        '--lengths-file',
        type=Path,
        help='Tab-delimited file with query lengths (model_name<TAB>length). '
        'Required when using --blast-results or --nhmmer-results without query files.',
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
        help='Maximum e-value threshold for filtering precomputed hits (default: 1e-3)',
    )
    filter_group.add_argument(
        '--blast-max-evalue',
        type=validate_evalue,
        default=1e-3,
        help='Maximum e-value threshold for BLAST searches (default: 1e-3)',
    )
    filter_group.add_argument(
        '--hmm-max-evalue',
        type=validate_evalue,
        default=1e-3,
        help='Maximum e-value threshold for nhmmer searches (default: 1e-3)',
    )
    filter_group.add_argument(
        '--min-identity',
        type=validate_identity,
        default=60.0,
        help='Minimum percent identity for BLAST searches (default: 60.0)',
    )
    filter_group.add_argument(
        '--word-size',
        type=validate_word_size,
        default=4,
        help='BLAST word size for initial matches (default: 4)',
    )
    filter_group.add_argument(
        '--max-offset',
        type=int,
        default=None,
        dest='max_offset',
        help=(
            'Maximum allowed offset (in bases) between the hit alignment boundary and '
            'the outer edge of the query model. '
            'Hits that do not reach within this many bases of the outer edge are removed. '
            'Requires query lengths to be available (from --fasta, --hmm, or --lengths-file). '
            'Default: no limit.'
        ),
    )
    filter_group.add_argument(
        '--orientation',
        type=str,
        default='F,R',
        dest='orientation',
        help=(
            'Orientation of left and right terminus models, used to determine the outer '
            'edge for --max-offset filtering. '
            'F = Forward (+), R = Reverse (-). '
            'Options: F,R (canonical TIR), F,F (LTR-like), R,R, R,F. '
            'Default: F,R'
        ),
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
    output_group.add_argument(
        '--split-paired-output',
        action='store_true',
        default=False,
        dest='split_paired_output',
        help=(
            'Write left and right model hits to separate output files based on '
            'the pairing map. Requires --pairing-map. Output files will be named '
            '<prefix>_left_hits.tab and <prefix>_right_hits.tab. '
            'Models appearing in both left and right columns of the pairing map '
            'are not allowed when this option is enabled.'
        ),
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
        '--blast-timeout',
        type=int,
        default=None,
        dest='blast_timeout',
        metavar='SECONDS',
        help='Maximum number of seconds to allow each BLAST search to run. '
        'Default: no limit (runs until complete).',
    )
    proc_group.add_argument(
        '--keep-temp',
        action='store_true',
        default=False,
        help='Keep temporary files after completion (default: False)',
    )
    proc_group.add_argument(
        '--loglevel',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='Set logging level (default: INFO)',
    )
    proc_group.add_argument(
        '--logfile',
        action='store_true',
        default=False,
        help='Write log messages to file in output directory',
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

    # Create output directory first (needed for logfile)
    args.outdir.mkdir(parents=True, exist_ok=True)

    # Set up logging with optional logfile
    logfile_path = None
    if args.logfile:
        logfile_name = (
            f'{args.prefix}_tirmite_search.log' if args.prefix else 'tirmite_search.log'
        )
        logfile_path = args.outdir / logfile_name

    init_logging(loglevel=args.loglevel, logfile=logfile_path)

    try:
        # Validate inputs
        has_search_inputs = args.fasta or args.hmm
        has_precomputed = args.blast_results or args.nhmmer_results

        if not has_search_inputs and not has_precomputed:
            raise EnsembleSearchError(
                'Must provide either search inputs (--fasta/--hmm + --genome/--blast-db) '
                'or precomputed results (--blast-results/--nhmmer-results)'
            )

        # Validate genome requirements for searches
        if args.fasta and not args.genome and not args.blast_db:
            raise EnsembleSearchError(
                '--genome or --blast-db is required when running BLAST searches with --fasta'
            )

        if args.hmm and not args.genome:
            raise EnsembleSearchError(
                '--genome is required when running nhmmer searches with --hmm'
            )

        # Validate --split-paired-output requirements
        if args.split_paired_output:
            if not args.pairing_map:
                raise EnsembleSearchError(
                    '--split-paired-output requires --pairing-map to identify '
                    'left and right models.'
                )

        # Validate lengths file requirement when using precomputed results without query files
        if has_precomputed and not has_search_inputs:
            if not args.lengths_file:
                msg = 'No --lengths-file provided. Query coverage filtering will not be available.'
                if getattr(args, 'max_offset', None) is not None:
                    msg += ' --max-offset filtering requires model lengths; please supply --lengths-file.'
                logging.warning(msg)

        # Collect genome files
        genomes: List[Path] = []
        if args.genome:
            genomes.append(args.genome)
        if args.genome_list:
            genomes.extend(parse_genome_list(args.genome_list))

        # Collect result files
        blast_files: List[Path] = list(args.blast_results) if args.blast_results else []
        nhmmer_files: List[Path] = (
            list(args.nhmmer_results) if args.nhmmer_results else []
        )

        # Collect query lengths from input files
        query_lengths: Dict[str, int] = {}
        if args.fasta:
            for fasta_file in args.fasta:
                if fasta_file.exists():
                    query_lengths.update(get_fasta_sequence_lengths(fasta_file))
        if args.hmm:
            for hmm_file in args.hmm:
                if hmm_file.exists():
                    query_lengths.update(get_hmm_model_lengths(hmm_file))
        if args.lengths_file:
            if args.lengths_file.exists():
                query_lengths.update(load_lengths_file(args.lengths_file))
            else:
                raise EnsembleSearchError(
                    f'Lengths file not found: {args.lengths_file}'
                )

        if query_lengths:
            logging.info(f'Loaded query lengths for {len(query_lengths)} models')

        # Run searches if needed
        if has_search_inputs:
            cleanup_temp = not args.keep_temp
            with temporary_directory(
                prefix='tirmite_search_', cleanup=cleanup_temp
            ) as temp_dir:
                temp_path = Path(temp_dir)

                if args.keep_temp:
                    logging.info(f'Temporary files will be kept in: {temp_path}')

                # Process each genome (or just run once if only using blast-db)
                genome_list: List[Optional[Path]] = list(genomes) if genomes else [None]
                for genome_file in genome_list:
                    # Prepare genome file (decompress if gzipped)
                    prepared_genome: Optional[Path] = None
                    if genome_file:
                        prepared_genome = prepare_genome_file(genome_file, temp_path)

                    # Run BLAST searches
                    if args.fasta:
                        for fasta_file in args.fasta:
                            if not fasta_file.exists():
                                logging.warning(f'FASTA file not found: {fasta_file}')
                                continue

                            # Use blast-db if provided, otherwise use prepared genome
                            if args.blast_db:
                                target = args.blast_db
                            elif prepared_genome:
                                target = prepared_genome
                            else:
                                logging.warning('No target for BLAST search')
                                continue

                            genome_suffix = (
                                f'_{genome_file.stem}' if genome_file else ''
                            )
                            output_file = (
                                temp_path
                                / f'{fasta_file.stem}{genome_suffix}_blast.tab'
                            )
                            result_file = run_blastn_search(
                                query_file=fasta_file,
                                target_file=target,
                                output_file=output_file,
                                evalue=args.blast_max_evalue,
                                identity=args.min_identity,
                                threads=args.threads,
                                word_size=args.word_size,
                                blast_timeout=args.blast_timeout,
                            )
                            blast_files.append(result_file)

                    # Run nhmmer searches
                    if args.hmm and prepared_genome:
                        for hmm_file in args.hmm:
                            if not hmm_file.exists():
                                logging.warning(f'HMM file not found: {hmm_file}')
                                continue

                            genome_suffix = (
                                f'_{genome_file.stem}' if genome_file else ''
                            )
                            result_file = run_nhmmer_search(
                                hmm_file=hmm_file,
                                target_file=prepared_genome,
                                output_dir=temp_path,
                                evalue=args.hmm_max_evalue,
                                threads=args.threads,
                            )
                            # Rename output file to include genome suffix
                            if genome_file:
                                new_result_file = (
                                    result_file.parent
                                    / f'{hmm_file.stem}{genome_suffix}.out'
                                )
                                if result_file != new_result_file:
                                    shutil.move(str(result_file), str(new_result_file))
                                    result_file = new_result_file
                            nhmmer_files.append(result_file)

                # Load and process hits within temp context
                hit_table = _process_hits(
                    args, blast_files, nhmmer_files, query_lengths
                )

        else:
            # Just load precomputed results
            hit_table = _process_hits(args, blast_files, nhmmer_files, query_lengths)

        # Write final output
        output_file = args.outdir / f'{args.prefix}_hits.tab'
        write_hits_table(hit_table, output_file)

        # Write split output if requested
        if args.split_paired_output:
            pairing_map = parse_pairing_map(args.pairing_map)
            validate_split_paired_output(pairing_map)
            write_split_hits(hit_table, pairing_map, args.outdir, args.prefix)

        # Log completion message with logfile location if enabled
        if logfile_path and args.logfile:
            logging.info(f'Log file saved to: {logfile_path}')

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
    query_lengths: Optional[Dict[str, int]] = None,
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
    query_lengths : dict, optional
        Mapping of model name to model length.  Required for anchor filtering.

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
    report_hit_statistics(hit_table, stage='(after e-value filtering)')

    # Apply anchor (outer-edge) filter if requested
    max_offset = getattr(args, 'max_offset', None)
    if max_offset is not None:
        orientation = getattr(args, 'orientation', 'F,R')

        # Resolve pairing map for left/right model identification
        anchor_pairing_map: Optional[Dict[str, str]] = None
        if args.pairing_map:
            if not args.pairing_map.exists():
                raise EnsembleSearchError(
                    f'Pairing map file not found: {args.pairing_map}'
                )
            anchor_pairing_map = parse_pairing_map(args.pairing_map)

        hit_table = filter_hits_by_anchor(
            hit_table,
            query_lengths=query_lengths if query_lengths else {},
            max_offset=max_offset,
            orientation=orientation,
            pairing_map=anchor_pairing_map,
        )

        # Report post-anchor statistics
        report_hit_statistics(hit_table, stage='(after anchor filtering)')

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
            is_valid, warnings = validate_cluster_mapping(
                cluster_map, available_features
            )
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
            raise EnsembleSearchError(f'Pairing map file not found: {args.pairing_map}')

        pairing_map = parse_pairing_map(args.pairing_map)

        if pairing_map:
            filter_summary = SearchFilterSummary()

            # Step 0: restrict output to models listed in the pairing map only.
            hit_table = filter_hits_to_pairing_map_models(
                hit_table, pairing_map, summary=filter_summary
            )

            # Report statistics after pairing map model filter
            report_hit_statistics(hit_table, stage='(after pairing map model filter)')

            # Step 1: remove hits from a paired model that are completely nested within
            # hits of its direct left/right partner and score significantly worse.
            hit_table = remove_nested_paired_hits(
                hit_table, pairing_map, summary=filter_summary
            )

            # Report statistics after nested hit removal
            report_hit_statistics(hit_table, stage='(after nested hit removal)')

            # Step 2: remove lower-quality overlapping hits from competing models across
            # all pairs in the pairing map.  This handles the case where models from
            # related element families hit the same locus: at each overlapping locus the
            # best-scoring model is retained and weaker cross-model hits are discarded.
            hit_table = filter_best_model_per_locus(
                hit_table, pairing_map, summary=filter_summary
            )

            # Report final statistics
            report_hit_statistics(
                hit_table, stage='(after cross-model overlap filtering)'
            )

            # Emit consolidated summary report for all pairing map filtering steps
            log_filter_summary(filter_summary)

    return hit_table


if __name__ == '__main__':
    import sys

    sys.exit(main())
