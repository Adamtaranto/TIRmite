#!/usr/bin/env python3
"""
HMM building workflow from seed sequences for TIRmite.

This module implements the seed-based HMM building workflow:
1. BLAST seed sequences against target genome(s)
2. Filter and process BLAST hits for quality
3. Extract and align genomic sequences
4. Build HMM models from alignments

Supports both symmetric (single seed) and asymmetric (left/right seeds)
transposon terminus models.
"""

import argparse
import hashlib
import logging
import os
from pathlib import Path
import shutil
from typing import Any, Dict, List, Optional, Set, Tuple, cast

from Bio import SeqIO  # type: ignore[import-not-found]
from Bio.Seq import Seq  # type: ignore[import-not-found]
from Bio.SeqRecord import SeqRecord  # type: ignore[import-not-found]
import pandas as pd  # type: ignore[import-untyped]
from pyhmmer.easel import (  # type: ignore[import-not-found]
    Alphabet,
    DigitalMSA,
    MSAFile,
    SequenceFile,
)
from pyhmmer.plan7 import (  # type: ignore[import-not-found]
    Background,
    Builder,
)

from tirmite.cli._argtypes import (
    validate_coverage,
    validate_evalue,
    validate_identity,
    validate_threads,
)
from tirmite.core.parsers import import_nhmmer
from tirmite.runners.hmmer_wrappers import (
    build_hmmalign_command,
    build_hmmpress_command,
    build_nhmmer_command,
)
from tirmite.runners.mafft import MafftError, align_to_file
from tirmite.runners.runBlastn import BlastError, run_blastn
from tirmite.runners.wrapping import run_command
from tirmite.utils.extract import (
    BlastDBSource,
    FastaSource,
    check_ids,
    clamp_region,
    fetch_sequence,
)
from tirmite.utils.logs import init_logging
from tirmite.utils.utils import (
    cleanID,
    create_output_directory,
    indexGenome,
    prepare_genome_file,
    temporary_directory,
)

logger = logging.getLogger(__name__)


class HMMBuildError(Exception):
    """
    Custom exception for HMM building errors.
    """

    pass


class BlastHit:
    """
    Container for BLAST hit information.

    Parameters
    ----------
    query_id : str
        Query sequence identifier.
    subject_id : str
        Subject/database sequence identifier.
    query_start : int
        Start position in query sequence.
    query_end : int
        End position in query sequence.
    subject_start : int
        Start position in subject sequence.
    subject_end : int
        End position in subject sequence.
    length : int
        Length of alignment.
    identity : float
        Percent identity of alignment.
    query_len : int
        Total length of query sequence.
    subject_len : int
        Total length of subject sequence.
    query_frame : int, default 1
        Reading frame of query alignment.
    subject_frame : int, default 1
        Reading frame of subject alignment.

    Attributes
    ----------
    query_coverage : float
        Query coverage calculated as length / query_len.
    subject_coverage : float
        Subject coverage calculated as length / subject_len.
    strand : str
        Strand orientation ('+' or '-') determined from coordinates and frames.
    """

    def __init__(
        self,
        query_id: str,
        subject_id: str,
        query_start: int,
        query_end: int,
        subject_start: int,
        subject_end: int,
        length: int,
        identity: float,
        query_len: int,
        subject_len: int,
        query_frame: int = 1,
        subject_frame: int = 1,
    ):
        self.query_id = query_id
        self.subject_id = subject_id
        self.query_start = query_start
        self.query_end = query_end
        self.subject_start = subject_start
        self.subject_end = subject_end
        self.length = length
        self.identity = identity
        self.query_len = query_len
        self.subject_len = subject_len
        self.query_frame = query_frame
        self.subject_frame = subject_frame
        self.query_coverage = length / query_len
        self.subject_coverage = length / subject_len

        # Determine strand from coordinates and frames
        # In BLAST tabular output, when subject_start > subject_end, it indicates reverse strand
        self.strand = '+' if subject_start <= subject_end else '-'

        # Additional validation using frames if available
        if subject_frame < 0:
            self.strand = '-'
        elif subject_frame > 0 and subject_start <= subject_end:
            self.strand = '+'

    def __repr__(self) -> str:
        return (
            f'BlastHit({self.query_id}->{self.subject_id}, '
            f'cov={self.query_coverage:.2f}, id={self.identity:.1f}%, strand={self.strand})'
        )


def check_dependencies() -> List[str]:
    """
    Check availability of required external command-line tools.

    Returns
    -------
    list of str
        Names of missing tools. Empty list if all tools are available.

    Notes
    -----
    Checks for: blastn, makeblastdb, mafft, hmmalign, nhmmer, hmmpress.
    """
    required_tools = [
        'blastn',
        'makeblastdb',
        'mafft',
        'hmmalign',
        'nhmmer',
        'hmmpress',
    ]
    missing = []

    for tool in required_tools:
        if not shutil.which(tool):
            missing.append(tool)

    return missing


def parse_blast_output(blast_file: Path) -> List[BlastHit]:
    """
    Parse BLAST tabular output file into BlastHit objects.

    Parameters
    ----------
    blast_file : Path
        Path to BLAST tabular output file.

    Returns
    -------
    list of BlastHit
        Parsed hit objects from all valid BLAST alignment records.

    Raises
    ------
    HMMBuildError
        If file cannot be read or parsing fails critically.

    Notes
    -----
    Expected BLAST output format:
    qstart qend sstart send length positive pident qlen slen qframe sframe qseqid sseqid

    Skips comment lines (starting with #) and malformed records.
    """
    hits = []

    try:
        with open(blast_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) >= 13:  # Need at least 13 fields for qseqid and sseqid
                    try:
                        hit = BlastHit(
                            query_id=fields[11],  # qseqid
                            subject_id=fields[12],  # sseqid
                            query_start=int(fields[0]),
                            query_end=int(fields[1]),
                            subject_start=int(fields[2]),
                            subject_end=int(fields[3]),
                            length=int(fields[4]),
                            identity=float(fields[6]),
                            query_len=int(fields[7]),
                            subject_len=int(fields[8]),
                            query_frame=int(fields[9]),
                            subject_frame=int(fields[10]),
                        )
                        hits.append(hit)
                    except (ValueError, IndexError) as e:
                        logger.warning(
                            f'Skipping malformed BLAST line {line_num}: {line.strip()} - {e}'
                        )
                        continue
                else:
                    logger.warning(
                        f'Skipping BLAST line {line_num} with insufficient fields: {line.strip()}'
                    )

    except Exception as e:
        raise HMMBuildError(f'Failed to parse BLAST output {blast_file}: {e}') from e

    return hits


def compare_seeds(
    left_seed: Path,
    right_seed: Path,
    temp_dir: Path,
    min_length: int = 10,
    min_identity: float = 50.0,
    num_threads: int = 1,  # Add threading parameter
) -> List[Tuple[BlastHit, object]]:
    """
    Compare left and right seed sequences using BLAST.

    Parameters
    ----------
    left_seed : Path
        Path to left seed sequence FASTA file.
    right_seed : Path
        Path to right seed sequence FASTA file.
    temp_dir : Path
        Temporary directory for BLAST output files.
    min_length : int, default 10
        Minimum alignment length threshold in base pairs.
    min_identity : float, default 50.0
        Minimum sequence identity threshold as percentage (0-100).
    num_threads : int, default 1
        Number of CPU threads for BLAST.

    Returns
    -------
    list of tuple
        List of (BlastHit, alignment) tuples for hits passing thresholds.

    Notes
    -----
    Identifies regions of similarity between left and right termini
    which may represent inverted repeats or conserved motifs.
    """
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner  # type: ignore[import-not-found]

    logger.info('Comparing left and right seed sequences...')

    blast_output = temp_dir / 'seed_comparison.tab'

    try:
        run_blastn(
            query=left_seed,
            subject=right_seed,
            output=blast_output,
            word_size=4,
            perc_identity=30.0,  # Lower threshold for initial search
            verbose=True,
            num_threads=num_threads,  # Pass threading parameter
        )

        all_hits = parse_blast_output(blast_output)
        logger.info(f'Found {len(all_hits)} initial seed comparison hits')

        if not all_hits:
            logger.info('No similarity found between left and right seeds')
            return []

        # Filter hits by length and identity thresholds
        filtered_hits = []
        for hit in all_hits:
            passed = True
            reasons = []

            if hit.length < min_length:
                passed = False
                reasons.append(f'length {hit.length}bp < {min_length}bp')

            if hit.identity < min_identity:
                passed = False
                reasons.append(f'identity {hit.identity:.1f}% < {min_identity:.1f}%')

            if passed:
                filtered_hits.append(hit)
                logger.debug(f'PASSED: {hit}')
            else:
                logger.debug(f'FILTERED: {hit} - Failed: {", ".join(reasons)}')

        logger.info(
            f'Threshold filtering: {len(all_hits)} -> {len(filtered_hits)} hits '
            f'(length>={min_length}bp, id>={min_identity:.1f}%)'
        )

        if not filtered_hits:
            logger.info(
                f'No seed comparison hits passed thresholds (length>={min_length}bp, id>={min_identity:.1f}%)'
            )
            return []

        # Load sequences for alignment
        left_records = list(SeqIO.parse(left_seed, 'fasta'))
        right_records = list(SeqIO.parse(right_seed, 'fasta'))

        # Create sequence lookup dictionaries
        left_seqs = {rec.id: rec.seq for rec in left_records}
        right_seqs = {rec.id: rec.seq for rec in right_records}

        # Generate alignments for each filtered hit
        hit_alignments = []
        aligner = PairwiseAligner()
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5

        for hit in filtered_hits:
            try:
                # Get query (left) and subject (right) sequences
                if hit.query_id not in left_seqs:
                    logger.warning(
                        f'Query sequence {hit.query_id} not found in left seed file'
                    )
                    continue

                if hit.subject_id not in right_seqs:
                    logger.warning(
                        f'Subject sequence {hit.subject_id} not found in right seed file'
                    )
                    continue

                # Extract hit regions with proper strand handling
                query_seq = left_seqs[hit.query_id][hit.query_start - 1 : hit.query_end]

                # Handle subject sequence based on coordinates
                if hit.subject_start <= hit.subject_end:
                    # Forward strand
                    subject_seq = right_seqs[hit.subject_id][
                        hit.subject_start - 1 : hit.subject_end
                    ]
                else:
                    # Reverse strand - coordinates are swapped in BLAST output
                    subject_seq = right_seqs[hit.subject_id][
                        hit.subject_end - 1 : hit.subject_start
                    ]
                    subject_seq = subject_seq.reverse_complement()

                # Perform pairwise alignment
                pairwise_alignments = aligner.align(query_seq, subject_seq)

                if pairwise_alignments:
                    best_alignment = pairwise_alignments[0]
                    hit_alignments.append((hit, best_alignment))

                    logger.debug(f'Generated alignment for hit: {hit}')
                    logger.debug(f'Alignment score: {best_alignment.score}')
                    logger.debug(
                        f'Query seq length: {len(query_seq)}, Subject seq length: {len(subject_seq)}'
                    )
                else:
                    logger.warning(f'Failed to generate alignment for hit: {hit}')

            except Exception as e:
                logger.warning(f'Error generating alignment for hit {hit}: {e}')
                continue

        logger.info(
            f'Generated {len(hit_alignments)} seed comparison alignments from {len(filtered_hits)} filtered hits'
        )

        # Sort by combined score (length * identity)
        hit_alignments.sort(key=lambda x: x[0].length * x[0].identity, reverse=True)

        return hit_alignments

    except BlastError as e:
        logger.warning(f'Seed comparison failed: {e}')
        return []
    except Exception as e:
        logger.error(f'Unexpected error in seed comparison: {e}')
        return []


def create_blast_database(genome_file: Path, db_dir: Path) -> Path:
    """
    Create BLAST nucleotide database from genome file using makeblastdb.

    Parameters
    ----------
    genome_file : Path
        Path to genome FASTA file to convert to BLAST database.
    db_dir : Path
        Directory where BLAST database files will be created.

    Returns
    -------
    Path
        Path to created BLAST database (without extension).

    Raises
    ------
    HMMBuildError
        If makeblastdb command fails.
    """
    db_name = db_dir / f'{genome_file.stem}_db'

    cmd = [
        'makeblastdb',
        '-in',
        str(genome_file),
        '-dbtype',
        'nucl',
        '-out',
        str(db_name),
        '-title',
        f'{genome_file.stem} database',
    ]

    logger.info(f'Creating BLAST database for {genome_file.name}')

    try:
        result = run_command(cmd, verbose=True)
        if result.returncode != 0:
            raise HMMBuildError(f'makeblastdb failed: {result.stderr}')

        logger.info(f'BLAST database created: {db_name}')
        return db_name

    except Exception as e:
        raise HMMBuildError(f'Failed to create BLAST database: {e}') from e


def blast_seed_against_genome(
    seed_file: Path,
    blast_db: Path,
    output_file: Path,
    identity_threshold: float = 60.0,
    num_threads: int = 1,
    evalue: float = 1e-3,
) -> List[BlastHit]:
    """
    BLAST seed sequence against genome database using blastn with threading support.

    Parameters
    ----------
    seed_file : Path
        Path to seed sequence file in FASTA format.
    blast_db : Path
        Path to BLAST database (without extension).
    output_file : Path
        Path to write BLAST tabular output.
    identity_threshold : float, default 60.0
        Minimum percent identity threshold for BLAST hits.
    num_threads : int, default 1
        Number of CPU threads for BLAST to use.
    evalue : float, default 1e-3
        E-value threshold for BLAST hits.

    Returns
    -------
    list of BlastHit
        List of parsed BLAST hit objects from output file.

    Raises
    ------
    HMMBuildError
        If BLAST command fails or output cannot be parsed.
    """
    cmd = [
        'blastn',
        '-query',
        str(seed_file),
        '-db',
        str(blast_db),
        '-out',
        str(output_file),
        '-outfmt',
        '6 qstart qend sstart send length positive pident qlen slen qframe sframe qseqid sseqid',
        '-word_size',
        '4',
        '-perc_identity',
        str(identity_threshold),
        '-max_target_seqs',
        '10000',
        '-evalue',
        str(evalue),
        '-num_threads',
        str(num_threads),
    ]

    # Log the full blastn command being executed
    logger.info('Running blastn with the following parameters:')
    logger.info(f'  Command: {" ".join(cmd)}')
    logger.info(f'  Query: {seed_file}')
    logger.info(f'  Database: {blast_db}')
    logger.info(f'  Output: {output_file}')
    logger.info(f'  Identity threshold: {identity_threshold}%')
    logger.info(f'  E-value: {evalue}')
    logger.info(f'  Threads: {num_threads}')

    try:
        result = run_command(cmd, verbose=True)
        if result.returncode != 0:
            raise BlastError(f'BLAST search failed: {result.stderr}')

        # Add header to BLAST output file
        add_blast_header(output_file)

        hits = parse_blast_output(output_file)
        logger.info(f'BLAST search found {len(hits)} hits')
        return hits

    except Exception as e:
        raise HMMBuildError(f'BLAST search failed: {e}') from e


def add_blast_header(blast_file: Path) -> None:
    """
    Add header comment lines to BLAST output file.

    Parameters
    ----------
    blast_file : Path
        Path to BLAST tabular output file to add header to.

    Returns
    -------
    None
        No return value. Modifies file in place.
    """
    # Read existing content
    with open(blast_file, 'r') as f:
        content = f.read()

    # Write header + content
    with open(blast_file, 'w') as f:
        f.write('# BLAST tabular output format 6\n')
        f.write(
            '# qstart\tqend\tsstart\tsend\tlength\tpositive\tpident\tqlen\tslen\tqframe\tsframe\tqseqid\tsseqid\n'
        )
        f.write(content)


def resolve_overlapping_hits(
    hits: List[BlastHit], max_overlap: int = 50
) -> List[BlastHit]:
    """
    Resolve overlapping hits using RepeatMasker-like logic.

    Parameters
    ----------
    hits : list of BlastHit
        List of BLAST hits to filter for overlaps.
    max_overlap : int, default 50
        Maximum allowed overlap in base pairs between hits.

    Returns
    -------
    list of BlastHit
        Filtered list with overlapping hits removed, prioritized by:
        1) alignment length, 2) percent identity, 3) query coverage.
    """
    if not hits:
        return hits

    # Sort hits by genomic position for overlap detection
    sorted_hits = sorted(
        hits, key=lambda h: (h.subject_id, h.subject_start, h.subject_end)
    )

    filtered_hits: List[Any] = []

    for current_hit in sorted_hits:
        keep_current = True
        remove_indices = []

        for i, existing_hit in enumerate(filtered_hits):
            # Only check hits on same chromosome/scaffold
            if current_hit.subject_id != existing_hit.subject_id:
                continue

            # Normalise before comparing. BLAST reports sstart > send for a
            # minus-strand HSP and BlastHit stores them unswapped (the strand
            # is derived from exactly that inversion), so comparing raw
            # coordinates made the overlap of any minus-strand hit provably
            # <= 0. No minus-strand hit was ever de-overlapped: every
            # reverse-oriented element contributed redundant, near-identical
            # copies to the HMM training set. hits_overlap() already
            # normalises this way.
            cur_lo, cur_hi = sorted(
                (current_hit.subject_start, current_hit.subject_end)
            )
            exist_lo, exist_hi = sorted(
                (existing_hit.subject_start, existing_hit.subject_end)
            )

            # Coordinates are 1-based inclusive, so two hits sharing a single
            # base overlap by 1.
            overlap_start = max(cur_lo, exist_lo)
            overlap_end = min(cur_hi, exist_hi)
            overlap_length = max(0, overlap_end - overlap_start + 1)

            if overlap_length > max_overlap:
                # Determine which hit to keep - prioritize length over identity
                current_score = (
                    current_hit.length,
                    current_hit.identity,
                    current_hit.query_coverage,
                )
                existing_score = (
                    existing_hit.length,
                    existing_hit.identity,
                    existing_hit.query_coverage,
                )

                if current_score > existing_score:
                    # Current hit is better, mark existing for removal
                    remove_indices.append(i)
                else:
                    # Existing hit is better, don't keep current
                    keep_current = False
                    break

        # Remove inferior hits
        for i in sorted(remove_indices, reverse=True):
            filtered_hits.pop(i)

        # Add current hit if it should be kept
        if keep_current:
            filtered_hits.append(current_hit)

    logger.info(f'Resolved overlaps: {len(hits)} -> {len(filtered_hits)} hits')
    return filtered_hits


def filter_hits_by_thresholds(
    hits: List[BlastHit], min_coverage: float, min_identity: float
) -> List[BlastHit]:
    """
    Filter hits by coverage and identity thresholds.

    Parameters
    ----------
    hits : list of BlastHit
        List of BLAST hits to filter.
    min_coverage : float
        Minimum query coverage required (0.0 to 1.0).
    min_identity : float
        Minimum percent identity required (0.0 to 100.0).

    Returns
    -------
    list of BlastHit
        Filtered list containing only hits meeting both thresholds.
    """
    filtered = []

    for hit in hits:
        passed = True
        reasons = []

        if hit.query_coverage < min_coverage:
            passed = False
            reasons.append(f'coverage {hit.query_coverage:.3f} < {min_coverage:.3f}')

        if hit.identity < min_identity:
            passed = False
            reasons.append(f'identity {hit.identity:.1f}% < {min_identity:.1f}%')

        if passed:
            filtered.append(hit)
            logger.debug(f'PASSED: {hit}')
        else:
            logger.debug(f'FILTERED: {hit} - Failed: {", ".join(reasons)}')

    logger.info(
        f'Threshold filtering: {len(hits)} -> {len(filtered)} hits '
        f'(cov>={min_coverage:.3f}, id>={min_identity:.1f})'
    )
    return filtered


def chain_fragmented_hits(
    hits: List[BlastHit], max_gap: int = 10
) -> List[List[BlastHit]]:
    """
    Chain hits that may represent fragments of the same element.

    Parameters
    ----------
    hits : list of BlastHit
        List of BLAST hits to chain together.
    max_gap : int, default 10
        Maximum gap in base pairs between hits to be considered part of same chain.

    Returns
    -------
    list of list of BlastHit
        List of hit chains where each chain is a list of consecutive hits.

    Notes
    -----
    Hits are chained only if they meet ALL of the following criteria:
    - Belong to the same query sequence
    - Are sequential and non-overlapping in query coordinates
    - Are sequential and non-overlapping in target coordinates
    - Are all in the same orientation (strand)
    - Are separated by less than max_gap bases in target sequence
    """
    if not hits:
        return []

    # Sort by query_id, subject_id, then query_start position
    # This ensures we process hits in the order they appear in the query
    sorted_hits = sorted(hits, key=lambda h: (h.query_id, h.subject_id, h.query_start))

    chains = []
    current_chain = [sorted_hits[0]]

    for hit in sorted_hits[1:]:
        last_hit = current_chain[-1]

        # Check all chaining criteria
        can_chain = True

        # Criterion 1: Must belong to the same query
        if hit.query_id != last_hit.query_id:
            can_chain = False

        # Criterion 2: Must be on the same subject/target sequence
        elif hit.subject_id != last_hit.subject_id:
            can_chain = False

        # Criterion 3: Must be in the same orientation/strand
        elif hit.strand != last_hit.strand:
            can_chain = False

        else:
            # Criterion 4: Hits must be sequential and non-overlapping in query
            # Query coordinates are always in forward orientation (start < end)
            if hit.query_start <= last_hit.query_end:
                # Hits overlap or are out of order in query
                can_chain = False

            # Criterion 5: Hits must be sequential and non-overlapping in target
            # Need to handle both forward and reverse strand cases
            elif hit.strand == '+':
                # Forward strand: subject_start < subject_end for both hits
                # Check if current hit starts after or at the end of last hit
                if hit.subject_start <= last_hit.subject_end:
                    # Hits overlap in target
                    can_chain = False
                # Check if gap is within max_gap
                elif hit.subject_start - last_hit.subject_end > max_gap:
                    # Gap too large
                    can_chain = False
            else:
                # Reverse strand: subject_start > subject_end for both hits
                # On reverse strand, we still need to ensure hits don't overlap
                # and are in the correct order
                min_last = min(last_hit.subject_start, last_hit.subject_end)
                max_curr = max(hit.subject_start, hit.subject_end)

                # Check for overlap
                if max_curr >= min_last:
                    # Hits overlap in target
                    can_chain = False
                # Check if gap is within max_gap
                elif min_last - max_curr > max_gap:
                    # Gap too large
                    can_chain = False

        if can_chain:
            current_chain.append(hit)
        else:
            # Start new chain
            chains.append(current_chain)
            current_chain = [hit]

    # Add the last chain
    chains.append(current_chain)

    # Log chaining results
    single_hits = sum(1 for chain in chains if len(chain) == 1)
    chained_hits = len(chains) - single_hits

    if chained_hits > 0:
        logger.warning(
            f'Hit chaining: {single_hits} single hits, {chained_hits} chains detected. '
            f'WARNING: Chained hits suggest fragmented BLAST matches. '
            f'It is recommended to manually inspect and trim the alignment before building the HMM with hmmbuild '
            f'to ensure high-quality terminal repeat models.'
        )
    else:
        logger.info(f'Hit chaining: {single_hits} single hits, {chained_hits} chains')

    return chains


def warn_multiple_queries(hits: List[BlastHit], context: str = '') -> None:
    """
    Warn if hits contain multiple distinct query sequence identifiers.

    Parameters
    ----------
    hits : list of BlastHit
        BLAST hits to inspect.
    context : str, optional
        Label used in log messages to identify which hit table is being checked.

    Returns
    -------
    None
        Emits a warning; returns nothing.
    """
    query_ids: Set[str] = {h.query_id for h in hits}
    if len(query_ids) > 1:
        label = f' ({context})' if context else ''
        logger.warning(
            f'Multiple query names detected in blast hit table{label}: '
            f'{", ".join(sorted(query_ids))}. '
            'Consider using a hit table derived from a single seed sequence.'
        )


class _LazySource:
    """Build a sequence source on first use and reuse it thereafter.

    Parameters
    ----------
    blast_db : Path or None
        BLAST database to read from, if any.
    genome_files : list of Path
        Genome FASTAs to index when no BLAST database is given.

    Notes
    -----
    A workflow needs the same source twice: once to check that every hit's
    target ID resolves, and again to extract the sequences. Building it twice
    meant indexing the genome twice, or -- for a BLAST database -- discarding
    the contig-length cache that the ID check had just populated, so extraction
    re-resolved every length at a ``blastdbcmd`` subprocess apiece.

    Construction stays lazy so that a run which fails before extraction (no
    hits passed thresholds, say) never pays to index a genome it will not read.
    """

    def __init__(self, blast_db: Optional[Path], genome_files: List[Path]):
        self._blast_db = blast_db
        self._genome_files = genome_files
        self._source: Optional[Any] = None

    @property
    def available(self) -> bool:
        """
        Report whether a source could be built at all.

        Returns
        -------
        bool
            True if a BLAST database or at least one genome was supplied.

        Notes
        -----
        Lets a caller skip an optional pre-flight check without forcing
        construction. Forcing it would raise "No genome or BLAST database
        available for extraction" from the ID check, pre-empting the more
        specific errors that come later -- a run with no hits above threshold
        should say so, not complain about extraction it never reached.
        """
        return self._blast_db is not None or bool(self._genome_files)

    def get(self) -> Any:
        """
        Return the shared source, building it on first call.

        Returns
        -------
        FastaSource or BlastDBSource
            The sequence source for this workflow.
        """
        if self._source is None:
            self._source = _build_source(self._blast_db, self._genome_files)
        return self._source


def check_targets_in_blastdb(hits: List[BlastHit], source: Any) -> List[str]:
    """
    Check that all target sequences in hits are present in a sequence source.

    Parameters
    ----------
    hits : list of BlastHit
        BLAST hits whose subject IDs should be validated.
    source : FastaSource or BlastDBSource
        Open sequence source to resolve the IDs against.

    Returns
    -------
    list of str
        Subject IDs that could NOT be found in the source.

    Notes
    -----
    Takes an open source rather than a database path so the caller can reuse
    the one it will extract with. This previously built its own
    ``BlastDBSource``, populated its contig-length cache while checking every
    ID, and then discarded it -- so extraction re-resolved every one of those
    lengths, at a ``blastdbcmd`` subprocess apiece.
    """
    return check_ids(source, (h.subject_id for h in hits))


def _build_source(blast_db: Optional[Path], genome_files: List[Path]):  # type: ignore[no-untyped-def]
    """
    Build the sequence source for extraction from CLI inputs.

    Parameters
    ----------
    blast_db : Path or None
        BLAST database path. Takes precedence over ``genome_files``.
    genome_files : list of Path
        Genome FASTA files. Only the first is indexed for extraction.

    Returns
    -------
    FastaSource or BlastDBSource
        A source wrapping whichever backend was supplied.

    Raises
    ------
    HMMBuildError
        If neither source is available, or the genome cannot be indexed.
    """
    if blast_db is not None:
        return BlastDBSource(blast_db)

    if not genome_files:
        raise HMMBuildError('No genome or BLAST database available for extraction')

    # Index every searched genome, not just the first.
    #
    # This used to index genome_files[0] alone and warn that hits on other
    # genomes "will be skipped". They were not skipped: extraction looks up
    # hit.subject_id in that single index, so when two assemblies share a
    # contig name -- 'chr1' is near-universal -- hits found in genome 2 were
    # extracted from genome 1 at genome 2's coordinates, yielding arbitrary
    # sequence with no warning at all.
    combined: Dict[str, Any] = {}
    combined_descriptions: Dict[str, str] = {}
    contig_origin: Dict[str, str] = {}
    collisions: Dict[str, List[str]] = {}

    for genome_file in genome_files:
        try:
            genome_index, descriptions = indexGenome(genome_file)
        except Exception as e:
            raise HMMBuildError(
                f'Failed to index genome {genome_file} for sequence extraction: {e}'
            ) from e

        for contig in genome_index.keys():
            if contig in combined:
                collisions.setdefault(contig, [contig_origin[contig]]).append(
                    genome_file.name
                )
                continue
            combined[contig] = genome_index[contig]
            contig_origin[contig] = genome_file.name

        combined_descriptions.update(descriptions)

    if collisions:
        shown = ', '.join(
            f'{contig} (in {" and ".join(sources)})'
            for contig, sources in list(collisions.items())[:5]
        )
        raise HMMBuildError(
            f'{len(collisions)} contig name(s) appear in more than one genome: '
            f'{shown}. A hit cannot be attributed to a genome by contig name '
            'alone, so extraction would silently take sequence from the wrong '
            'assembly. Rename the contigs to be unique across genomes, or '
            'search one genome at a time.'
        )

    logger.info(
        f'Indexed {len(combined)} contigs across {len(genome_files)} genome(s) '
        'for sequence extraction'
    )

    return FastaSource(combined, combined_descriptions)


def extract_sequences_from_chains(  # type: ignore[no-untyped-def]
    chains: List[List[BlastHit]], source, model_name: str
) -> List[SeqRecord]:
    """
    Extract sequences from hit chains, concatenating fragments where needed.

    Parameters
    ----------
    chains : list of list of BlastHit
        List of hit chains where each chain is a list of consecutive hits.
    source : FastaSource or BlastDBSource
        Sequence source, from :func:`tirmite.utils.extract.make_source`.
    model_name : str
        Name of model for sequence ID generation.

    Returns
    -------
    list of Bio.SeqRecord.SeqRecord
        List of extracted sequences, one per chain.

    Notes
    -----
    Coordinates are 1-based inclusive; see :mod:`tirmite.utils.extract`.
    """
    sequences = []

    for chain_idx, chain in enumerate(chains):
        try:
            if len(chain) == 1:
                # Single hit
                hit = chain[0]

                start = min(hit.subject_start, hit.subject_end)
                end = max(hit.subject_start, hit.subject_end)

                seq_str = fetch_sequence(source, hit.subject_id, start, end)
                if seq_str is None:
                    logger.warning(
                        f'Could not extract {hit.subject_id}:{start}-{end}. '
                        'Skipping hit.'
                    )
                    continue

                # Reverse complement if on negative strand
                if hit.strand == '-':
                    seq_str = str(Seq(seq_str).reverse_complement())

                # Include strand in sequence ID
                seq_id = f'{model_name}_{hit.subject_id}_{start}_{end}_{hit.strand}'

            else:
                # Chained hits - concatenate fragments
                seq_parts = []
                locations = []
                valid_chain = True

                for hit in chain:
                    start = min(hit.subject_start, hit.subject_end)
                    end = max(hit.subject_start, hit.subject_end)

                    fragment = fetch_sequence(source, hit.subject_id, start, end)
                    if fragment is None:
                        logger.warning(
                            f'Could not extract {hit.subject_id}:{start}-{end}. '
                            'Skipping chain.'
                        )
                        valid_chain = False
                        break

                    if hit.strand == '-':
                        fragment = str(Seq(fragment).reverse_complement())

                    seq_parts.append(fragment)
                    locations.append(f'{start}_{end}_{hit.strand}')

                if not valid_chain:
                    continue

                seq_str = ('N' * 10).join(seq_parts)  # Join fragments with Ns
                seq_id = f'{model_name}_{chain[0].subject_id}_chain_{chain_idx}_{"_".join(locations)}'

            # Create sequence record
            seq_record = SeqRecord(
                Seq(seq_str),
                id=seq_id,
                description=f'Extracted from {len(chain)} BLAST hit(s)',
            )
            sequences.append(seq_record)

        except KeyError as e:
            logger.error(f'Failed to extract sequence for chain {chain_idx}: {e}')
            continue
        except Exception as e:
            logger.error(
                f'Unexpected error extracting sequence for chain {chain_idx}: {e}'
            )
            continue

    logger.info(f'Extracted {len(sequences)} sequences from {len(chains)} hit chains')
    return sequences


def extract_flanked_sequences_from_chains(  # type: ignore[no-untyped-def]
    chains: List[List[BlastHit]], source, model_name: str, flank_size: int
) -> List[SeqRecord]:
    """
    Extract sequences from hit chains with flanking sequence.

    Parameters
    ----------
    chains : list of list of BlastHit
        List of hit chains where each chain is a list of consecutive hits.
    source : FastaSource or BlastDBSource
        Sequence source, from :func:`tirmite.utils.extract.make_source`.
    model_name : str
        Name of model for sequence ID generation.
    flank_size : int
        Number of base pairs of flanking sequence to add on each side.

    Returns
    -------
    list of Bio.SeqRecord.SeqRecord
        List of extracted sequences with flanks, one per chain.

    Notes
    -----
    For chained hits, only adds flanks to the beginning of the first segment
    and end of the last segment. Sequence IDs report the flank coordinates
    after clamping to the contig, so both backends label the same region
    identically.
    """
    sequences = []

    for chain_idx, chain in enumerate(chains):
        try:
            if len(chain) == 1:
                # Single hit with flanking
                hit = chain[0]

                hit_start = min(hit.subject_start, hit.subject_end)
                hit_end = max(hit.subject_start, hit.subject_end)

                # Add flanking sequence, clamped to the contig
                region = clamp_region(
                    source,
                    hit.subject_id,
                    hit_start - flank_size,
                    hit_end + flank_size,
                )
                if region is None:
                    logger.warning(f'Could not resolve {hit.subject_id}. Skipping hit.')
                    continue
                flanked_start, flanked_end = region

                seq_str = fetch_sequence(
                    source, hit.subject_id, flanked_start, flanked_end
                )
                if seq_str is None:
                    logger.warning(
                        f'Could not extract {hit.subject_id}:'
                        f'{flanked_start}-{flanked_end}. Skipping hit.'
                    )
                    continue

                # Reverse complement if on negative strand
                if hit.strand == '-':
                    seq_str = str(Seq(seq_str).reverse_complement())

                # Include flanking info in sequence ID
                seq_id = f'{model_name}_{hit.subject_id}_{flanked_start}_{flanked_end}_{hit.strand}_flank{flank_size}'

            else:
                # Chained hits - add flanks only to first and last segments
                first_hit = chain[0]

                seq_parts = []
                locations = []
                valid_chain = True

                for i, hit in enumerate(chain):
                    hit_start = min(hit.subject_start, hit.subject_end)
                    hit_end = max(hit.subject_start, hit.subject_end)

                    # Add flanking sequence only to first and last segments
                    if i == 0:  # First segment
                        req_start = hit_start - flank_size
                        req_end = hit_end
                        suffix = f'_5flank{flank_size}'
                    elif i == len(chain) - 1:  # Last segment
                        req_start = hit_start
                        req_end = hit_end + flank_size
                        suffix = f'_3flank{flank_size}'
                    else:  # Middle segments - no flanking
                        req_start = hit_start
                        req_end = hit_end
                        suffix = ''

                    region = clamp_region(source, hit.subject_id, req_start, req_end)
                    if region is None:
                        logger.warning(
                            f'Could not resolve {hit.subject_id}. Skipping chain.'
                        )
                        valid_chain = False
                        break
                    start, end = region

                    fragment = fetch_sequence(source, hit.subject_id, start, end)
                    if fragment is None:
                        logger.warning(
                            f'Could not extract {hit.subject_id}:{start}-{end}. '
                            'Skipping chain.'
                        )
                        valid_chain = False
                        break

                    if hit.strand == '-':
                        fragment = str(Seq(fragment).reverse_complement())

                    locations.append(f'{start}_{end}_{hit.strand}{suffix}')
                    seq_parts.append(fragment)

                if not valid_chain:
                    continue

                seq_str = ('N' * 10).join(seq_parts)  # Join fragments with Ns
                seq_id = f'{model_name}_{first_hit.subject_id}_chain_{chain_idx}_flank{flank_size}_{"_".join(locations)}'

            # Create sequence record
            seq_record = SeqRecord(
                Seq(seq_str),
                id=seq_id,
                description=f'Extracted from {len(chain)} BLAST hit(s) with {flank_size}bp flanking',
            )
            sequences.append(seq_record)

        except KeyError as e:
            logger.error(
                f'Failed to extract flanked sequence for chain {chain_idx}: {e}'
            )
            continue
        except Exception as e:
            logger.error(
                f'Unexpected error extracting flanked sequence for chain {chain_idx}: {e}'
            )
            continue

    logger.info(
        f'Extracted {len(sequences)} flanked sequences from {len(chains)} hit chains (flank size: {flank_size}bp)'
    )
    return sequences


def deduplicate_sequences(sequences: List[SeqRecord]) -> List[SeqRecord]:
    """
    Remove identical sequences, keeping first occurrence and prioritizing seeds.

    Parameters
    ----------
    sequences : list of Bio.SeqRecord.SeqRecord
        List of sequence records potentially containing duplicates.

    Returns
    -------
    list of Bio.SeqRecord.SeqRecord
        Deduplicated list with seed sequences prioritized over extracted sequences.
    """
    seen_sequences = set()
    unique_sequences = []

    # Sort sequences to prioritize seeds (seed sequences typically don't have genomic coordinates in their IDs)
    def is_seed_sequence(seq_record: Any) -> bool:
        """
        Check if sequence is likely a seed sequence (no genomic coordinates).

        Parameters
        ----------
        seq_record : Bio.SeqRecord.SeqRecord
            Sequence record to check.

        Returns
        -------
        bool
            True if likely a seed sequence, False otherwise.
        """
        seq_id = seq_record.id.lower()
        # Seed sequences typically don't have chromosome coordinates or chain info
        return not any(
            keyword in seq_id for keyword in ['_chain_', '_contig', '_scaffold', '_chr']
        )

    # Sort so seed sequences come first
    sorted_sequences = sorted(sequences, key=lambda x: (not is_seed_sequence(x), x.id))

    for seq in sorted_sequences:
        # Ensure all sequences are uppercase for consistent comparison
        seq_str = str(seq.seq).upper()
        if seq_str not in seen_sequences:
            seen_sequences.add(seq_str)
            # Create a new record with uppercase sequence
            uppercase_seq = SeqRecord(
                Seq(seq_str), id=seq.id, description=seq.description
            )
            unique_sequences.append(uppercase_seq)
        else:
            # Log when we skip a duplicate, noting if we kept the seed version
            if is_seed_sequence(seq):
                logger.info(
                    f'Seed sequence {seq.id} is identical to already retained sequence'
                )
            else:
                logger.debug(f'Skipping duplicate extracted sequence: {seq.id}')

    logger.info(
        f'Deduplication: {len(sequences)} -> {len(unique_sequences)} unique sequences'
    )
    return unique_sequences


def run_mafft_alignment(
    sequences: List[SeqRecord], output_file: Path, threads: int = 1
) -> Path:
    """
    Align sequences with MAFFT, reporting failures as :class:`HMMBuildError`.

    Parameters
    ----------
    sequences : list of Bio.SeqRecord.SeqRecord
        Sequences to align.
    output_file : Path
        Path to write the alignment to, in FASTA format.
    threads : int, default 1
        Number of CPU threads to pass to MAFFT.

    Returns
    -------
    Path
        ``output_file``.

    Raises
    ------
    HMMBuildError
        If the alignment could not be produced.

    Notes
    -----
    Thin adapter over :func:`tirmite.runners.mafft.align_to_file`. It exists so
    that every failure inside the seed workflow surfaces as a single exception
    type, which ``main`` catches and turns into a clean CLI error.
    """
    try:
        return align_to_file(sequences, output_file, threads=threads)
    except MafftError as e:
        raise HMMBuildError(f'MAFFT alignment failed: {e}') from e


def calculate_pairwise_identity(alignment_file: Path) -> pd.DataFrame:
    """
    Calculate pairwise identity matrix from alignment.

    Parameters
    ----------
    alignment_file : Path
        Path to alignment file in FASTA format.

    Returns
    -------
    pandas.DataFrame
        Pairwise identity matrix with sequence IDs as index and columns.
    """
    try:
        # Use SeqIO.parse instead of Align.read to get SeqRecord objects
        alignment_records = list(SeqIO.parse(alignment_file, 'fasta'))
        sequences = [record.seq for record in alignment_records]
        seq_ids = [record.id for record in alignment_records]

        n_seqs = len(sequences)
        identity_matrix = pd.DataFrame(index=seq_ids, columns=seq_ids, dtype=float)

        for i in range(n_seqs):
            for j in range(n_seqs):
                if i == j:
                    identity_matrix.iloc[i, j] = 100.0
                else:
                    seq1, seq2 = sequences[i], sequences[j]
                    matches = sum(
                        1
                        for a, b in zip(seq1, seq2)
                        if a == b and a != '-' and b != '-'
                    )
                    valid_positions = sum(
                        1 for a, b in zip(seq1, seq2) if a != '-' and b != '-'
                    )

                    if valid_positions > 0:
                        identity = (matches / valid_positions) * 100
                    else:
                        identity = 0.0

                    identity_matrix.iloc[i, j] = identity

        return identity_matrix

    except Exception as e:
        logger.warning(f'Failed to calculate pairwise identity: {e}')
        return pd.DataFrame()


def clean_hmm_name(name: str) -> str:
    """
    Clean model name for HMM compatibility with strict validation.

    Parameters
    ----------
    name : str
        Raw model name to clean.

    Returns
    -------
    str
        Cleaned name suitable for HMM models.
    """
    if not name:
        return 'unnamed_model'

    # Convert to string and strip whitespace
    name = str(name).strip()

    # Replace problematic characters with underscores
    # Keep only alphanumeric, underscore, and hyphen
    cleaned = ''.join(c if c.isalnum() or c in '_-' else '_' for c in name)

    # Remove leading/trailing underscores and collapse multiple underscores
    cleaned = '_'.join(part for part in cleaned.split('_') if part)

    # Ensure it starts with a letter (not a number or underscore)
    if cleaned and not cleaned[0].isalpha():
        cleaned = f'model_{cleaned}'

    # Ensure minimum length
    if not cleaned:
        cleaned = 'unnamed_model'

    # HMMER imposes no practical limit on the NAME field, so the cap here is
    # only to keep filenames manageable. It was 20, which silently collided:
    # '<name>_left' and '<name>_right' both truncate to the same string for any
    # name of 15 characters or more, and since the HMM path is
    # f'{clean_model_name}.hmm' the right model OVERWROTE the left one and both
    # returned paths pointed at the same file.
    max_length = 60
    if len(cleaned) > max_length:
        # Truncation can still collide, so disambiguate with a short digest of
        # the full name. Deterministic, so re-running gives the same filename.
        digest = hashlib.sha1(cleaned.encode('utf-8')).hexdigest()[:8]
        keep = max_length - len(digest) - 1
        cleaned = f'{cleaned[:keep].rstrip("_")}_{digest}'
        logger.warning(
            f'Model name truncated to {max_length} characters: "{cleaned}". '
            'A digest of the full name is appended so distinct models cannot '
            'collide.'
        )

    # Final validation - ensure it's not empty after processing
    if not cleaned:
        cleaned = 'model'

    return cleaned


def write_seed_comparison_report(
    seed_comparisons: List[Tuple[Any, Any]],
    output_file: Path,
    model_name: str,
    left_seed_name: str,
    right_seed_name: str,
) -> Path:
    """
    Write the detailed left-vs-right seed similarity report.

    Parameters
    ----------
    seed_comparisons : list of tuple
        ``(BlastHit, alignment)`` pairs as returned by :func:`compare_seeds`.
    output_file : Path
        Path to write the report to.
    model_name : str
        Model name, for the report title.
    left_seed_name : str
        Filename of the left seed.
    right_seed_name : str
        Filename of the right seed.

    Returns
    -------
    Path
        ``output_file``.

    Notes
    -----
    Extracted from ``main`` so that it can be tested. It previously lived
    inline inside a branch of ``main`` that also binds ``alignment`` to a
    ``Path`` in its sibling ``--update`` branch; a rename of the loop variable
    left three stale references behind and the block raised ``NameError`` for
    every asymmetric run that found any seed similarity. Nothing covered it,
    because ``main`` has no test that reaches this point.
    """
    with open(output_file, 'w') as handle:
        handle.write(f'Seed Comparison Results for {model_name}\n')
        handle.write('=' * 50 + '\n\n')
        handle.write(f'Left seed: {left_seed_name}\n')
        handle.write(f'Right seed: {right_seed_name}\n')
        handle.write(f'Total similarities found: {len(seed_comparisons)}\n\n')

        for i, (hit, seed_alignment) in enumerate(seed_comparisons, 1):
            handle.write(f'Similarity {i}:\n')
            handle.write(f'  Length: {hit.length}bp\n')
            handle.write(f'  Identity: {hit.identity:.1f}%\n')
            handle.write(f'  Query coverage: {hit.query_coverage:.3f}\n')
            handle.write(f'  Subject coverage: {hit.subject_coverage:.3f}\n')
            handle.write(
                f'  Query: {hit.query_id}[{hit.query_start}:{hit.query_end}]\n'
            )
            handle.write(
                f'  Subject: {hit.subject_id}[{hit.subject_start}:{hit.subject_end}]\n'
            )
            handle.write(f'  Alignment score: {seed_alignment.score}\n')
            handle.write('  Alignment:\n')
            for line in str(seed_alignment).split('\n'):
                handle.write(f'    {line}\n')
            handle.write('\n')

    return output_file


def read_alignment_records(alignment_file: Path) -> List[SeqRecord]:
    """
    Read an alignment without assuming its format.

    Parameters
    ----------
    alignment_file : Path
        Path to a multiple sequence alignment.

    Returns
    -------
    list of Bio.SeqRecord.SeqRecord
        The aligned records. Empty if the file holds no sequences in any of
        the supported formats.

    Raises
    ------
    HMMBuildError
        If the file cannot be read at all.

    Notes
    -----
    TIRmite produces alignments in two formats and feeds both to the same HMM
    builder: the seed workflow aligns with MAFFT and writes **FASTA**, while
    ``--update`` aligns with ``hmmalign`` and writes **Stockholm**.

    This previously hard-coded ``'stockholm'``, so every seed run failed with
    ``ValueError: Did not find STOCKHOLM header`` before it reached the
    builder. The format is now detected by trying each candidate in turn,
    matching ``pyhmmer.easel.MSAFile``, which auto-detects and does the actual
    work a few lines later.

    Stockholm is tried first because its header is unambiguous, so a Stockholm
    file can never be misread as something else.
    """
    last_error: Optional[Exception] = None

    for fmt in ('stockholm', 'fasta', 'clustal', 'phylip'):
        try:
            records = list(SeqIO.parse(alignment_file, fmt))
        except (ValueError, IndexError) as e:
            # Wrong format for this parser; try the next one.
            last_error = e
            continue

        if records:
            logger.debug(f'Read {len(records)} records from {alignment_file} as {fmt}')
            return records

    if last_error is not None:
        logger.debug(f'No parser accepted {alignment_file}: {last_error}')
    return []


def build_hmm_from_alignment_pyhmmer(
    alignment_file: Path, model_name: str, output_dir: Path
) -> Path:
    """
    Build HMM from multiple sequence alignment using pyhmmer.

    Parameters
    ----------
    alignment_file : Path
        Path to multiple sequence alignment file.
    model_name : str
        Name for the HMM model.
    output_dir : Path
        Directory to write HMM output file.

    Returns
    -------
    Path
        Path to created HMM file.
    """
    # Clean the model name with stricter rules
    clean_model_name = clean_hmm_name(model_name)
    output_hmm = output_dir / f'{clean_model_name}.hmm'

    # Add debugging
    logger.debug(f'Original model name: "{model_name}"')
    logger.debug(f'Cleaned model name: "{clean_model_name}"')
    logger.debug(f'Cleaned name length: {len(clean_model_name)}')
    logger.debug(f'Name characters: {[ord(c) for c in clean_model_name]}')

    try:
        alphabet = Alphabet.dna()

        # Read alignment file first to validate it
        alignment_records = read_alignment_records(alignment_file)

        if not alignment_records:
            raise HMMBuildError(
                f'No sequences found in alignment file: {alignment_file}'
            )

        if len(alignment_records) < 2:
            raise HMMBuildError(
                f'Need at least 2 sequences for HMM building, got {len(alignment_records)}'
            )

        logger.debug(f'Found {len(alignment_records)} sequences in alignment')

        # Try pyhmmer's native MSA file reading first
        try:
            with MSAFile(
                str(alignment_file), digital=True, alphabet=alphabet
            ) as msa_file:
                msa = msa_file.read()

            # len(msa) is the alignment WIDTH in columns; msa.sequences is the
            # sequence list. A zero-width alignment and an empty one are both
            # unusable, so check the sequence count.
            if msa is None or len(msa.sequences) == 0:
                raise HMMBuildError(
                    f'No sequences found in alignment file: {alignment_file}'
                )

            logger.debug(
                f'Successfully read MSA with {len(msa.sequences)} sequences '
                f'({len(msa)} columns) using MSAFile'
            )

        except Exception as msa_error:
            logger.debug(
                f'MSA file reading failed: {msa_error}. Trying sequence file approach.'
            )

            # Read as sequences and manually create MSA
            sequences = []
            try:
                with SequenceFile(
                    str(alignment_file), digital=True, alphabet=alphabet
                ) as seq_file:
                    for seq in seq_file:
                        sequences.append(seq)
            except Exception as seq_error:
                logger.debug(
                    f'SequenceFile reading failed: {seq_error}. Creating MSA manually.'
                )

                # Manual creation from BioPython records
                import pyhmmer.easel  # type: ignore[import-not-found]

                text_sequences = []
                for i, record in enumerate(alignment_records):
                    seq_str = str(record.seq).upper()

                    # Clean sequence ID
                    seq_id = clean_hmm_name(record.id) if record.id else f'seq_{i + 1}'

                    try:
                        text_seq = pyhmmer.easel.TextSequence(
                            name=seq_id.encode('ascii'), sequence=seq_str
                        )
                        digital_seq = text_seq.digitize(alphabet)
                        text_sequences.append(digital_seq)
                    except Exception as e:
                        logger.error(
                            f'Failed to create digital sequence for {seq_id}: {e}'
                        )
                        raise

                sequences = text_sequences

            if not sequences:
                raise HMMBuildError(
                    f'No sequences found in alignment file: {alignment_file}'
                ) from None  # Add 'from None' since this is a new logical error, not chained from seq_error

            # Every argument must be passed by keyword. `sequences` given
            # positionally lands in `name`, which expects bytes, so this
            # fallback raised "TypeError: Expected bytes, got list" every time
            # it was reached. pyhmmer also warns that `alphabet` will stop
            # accepting a positional argument after v0.12.0.
            msa = DigitalMSA(alphabet=alphabet, sequences=sequences)
            logger.debug(f'Created MSA manually with {len(msa.sequences)} sequences')

        # len(msa) is the alignment WIDTH in columns, not the number of
        # sequences; msa.sequences is the sequence list.
        logger.info(
            f'Building HMM for {model_name} from {len(msa.sequences)} sequences '
            f'({len(msa)} alignment columns)'
        )

        # Try to set MSA name before building (this might help)
        try:
            msa.name = clean_model_name.encode('ascii')
            logger.debug(f'Set MSA name to: {clean_model_name}')
        except Exception as name_error:
            logger.warning(f'Could not set MSA name: {name_error}')

        # Build HMM using pyhmmer Builder
        builder = Builder(alphabet)
        background = Background(alphabet)

        logger.debug('Starting HMM building...')

        try:
            # Build the HMM from the DigitalMSA
            hmm, _, _ = builder.build_msa(msa, background)
            logger.debug('HMM building successful')

        except Exception as build_error:
            logger.error(f'HMM building failed at builder.build_msa(): {build_error}')
            # Diagnostics must not themselves raise: an exception here would
            # mask build_error and skip the retry below. msa[:3] slices a
            # DigitalMSA, which is not guaranteed to support slicing.
            try:
                logger.error(
                    f'MSA details: {len(msa.sequences)} sequences, {len(msa)} columns'
                )
            except Exception:  # pragma: no cover - diagnostics only
                logger.error('MSA details unavailable')

            # Try building without setting MSA name
            logger.debug('Retrying HMM building without MSA name...')
            msa.name = None
            hmm, _, _ = builder.build_msa(msa, background)
            logger.debug('HMM building successful without MSA name')

        # Set HMM name after building
        try:
            hmm.name = clean_model_name.encode('ascii')
            logger.debug(f'Set HMM name to: {clean_model_name}')
        except Exception as name_error:
            logger.warning(f'Could not set HMM name: {name_error}')
            # Try with an even simpler name
            simple_name = 'model'
            hmm.name = simple_name.encode('ascii')
            logger.warning(f'Using simple name: {simple_name}')

        # Write HMM to file
        with open(output_hmm, 'wb') as f:
            hmm.write(f)

        logger.info(f'HMM written to {output_hmm}')
        return output_hmm

    except Exception as e:
        logger.error(f'pyhmmer HMM building failed: {e}')
        logger.error(f'Model name: "{model_name}" -> "{clean_model_name}"')
        raise HMMBuildError(f'pyhmmer HMM building failed: {e}') from e


def save_all_blast_hits(all_hits: List[BlastHit], output_file: Path) -> None:
    """
    Save all BLAST hits in tab-delimited format 6.

    Parameters
    ----------
    all_hits : list of BlastHit
        List of BLAST hit objects to save.
    output_file : Path
        Path to output file.

    Returns
    -------
    None
        No return value. Writes hits to file.
    """
    with open(output_file, 'w') as f:
        # Write header
        f.write('# BLAST tabular output format 6\n')
        f.write(
            '# qstart\tqend\tsstart\tsend\tlength\tpositive\tpident\tqlen\tslen\tqframe\tsframe\tqseqid\tsseqid\n'
        )

        # Write all hits
        for hit in all_hits:
            f.write(
                f'{hit.query_start}\t{hit.query_end}\t{hit.subject_start}\t{hit.subject_end}\t'
                f'{hit.length}\t{hit.length}\t{hit.identity:.2f}\t{hit.query_len}\t{hit.subject_len}\t'
                f'{hit.query_frame}\t{hit.subject_frame}\t{hit.query_id}\t{hit.subject_id}\n'
            )

    logger.info(f'All BLAST hits saved to {output_file}')


def process_seed_sequences(
    seed_file: Path,
    model_name: str,
    genome_files: List[Path],
    temp_dir: Path,
    output_dir: Path,
    min_coverage: float,
    min_identity: float,
    save_blast_hits: bool = False,
    flank_size: Optional[int] = None,
    threads: int = 1,
    evalue: float = 1e-3,
    blast_db: Optional[Path] = None,
    blast_hits_file: Optional[Path] = None,
) -> Tuple[Path, Path, Optional[Path], Optional[Path]]:
    """
    Process a seed file against genomes to build HMM.

    Parameters
    ----------
    seed_file : Path
        Path to seed sequence file in FASTA format.
    model_name : str
        Name for the HMM model.
    genome_files : list of Path
        List of genome FASTA files to search against.  May be empty when
        ``blast_db`` is provided.
    temp_dir : Path
        Directory for temporary files.
    output_dir : Path
        Directory for output files.
    min_coverage : float
        Minimum coverage threshold (0.0 to 1.0).
    min_identity : float
        Minimum percent identity threshold (0.0 to 100.0).
    save_blast_hits : bool, default False
        If True, save all BLAST hits to file.
    flank_size : int, optional
        Size of flanking sequence to extract.
    threads : int, default 1
        Number of CPU threads for BLAST.
    evalue : float, default 1e-3
        E-value threshold for BLAST hits.
    blast_db : Path, optional
        Path to a pre-built BLAST database (without extension).  When provided,
        blastn is run against this database and hit sequences are extracted via
        ``blastdbcmd`` rather than from a FASTA index.  Mutually exclusive with
        ``genome_files`` as the search target; the database must have been
        created with ``-parse_seqids``.
    blast_hits_file : Path, optional
        Path to a pre-computed BLAST hit table in the TIRmite tabular format
        (generated with ``--save-blast-hits``).  When provided the BLAST search
        step is skipped and these hits are used directly.  Sequence extraction
        still requires either ``genome_files`` or ``blast_db``.

    Returns
    -------
    tuple
        (hmm_file, alignment_file, blast_hits_file, flanked_alignment_file)
        where blast_hits_file and flanked_alignment_file may be None.
    """
    logger.info(f'Processing {model_name} seed: {seed_file.name}')

    # One source per workflow, built on first use. Shared between the
    # target-ID check and sequence extraction so the genome is indexed once
    # and blastdbcmd's contig-length cache survives between them.
    lazy_source = _LazySource(blast_db, genome_files)

    # ------------------------------------------------------------------ #
    # Step 1 – Obtain BLAST hits                                           #
    # ------------------------------------------------------------------ #
    if blast_hits_file is not None:
        # Use pre-computed hit table
        logger.info(f'Loading pre-computed BLAST hits from {blast_hits_file}')
        all_hits = parse_blast_output(blast_hits_file)
        logger.info(f'Loaded {len(all_hits)} hits from {blast_hits_file}')

        if not all_hits:
            raise HMMBuildError(
                f'No hits found in provided blast-hits file: {blast_hits_file}'
            )

        # Warn about multiple query names
        warn_multiple_queries(all_hits, context=model_name)

        # Validate that all target sequences are reachable, using the same
        # source extraction will use. Building it here rather than later means
        # the genome is indexed once instead of twice, and the contig lengths
        # resolved by this check are still cached when extraction runs.
        if lazy_source.available:
            missing = check_targets_in_blastdb(all_hits, lazy_source.get())
            if missing:
                logger.warning(
                    f'The following target sequences from the blast-hits table '
                    f'could not be resolved in the sequence source: '
                    f'{", ".join(sorted(set(missing)))}'
                )

    elif blast_db is not None:
        # Run BLAST against user-supplied pre-built database
        logger.info(f'Running BLAST against pre-built database: {blast_db}')
        blast_output = temp_dir / f'{model_name}_blast.tab'
        all_hits = blast_seed_against_genome(
            seed_file,
            blast_db,
            blast_output,
            min_identity,
            num_threads=threads,
            evalue=evalue,
        )

    else:
        # Build databases from genome FASTA files and run BLAST
        all_hits = []
        for genome_file in genome_files:
            logger.info(f'Searching {model_name} against {genome_file.name}')

            db_dir = temp_dir / f'blast_dbs_{genome_file.stem}'
            db_dir.mkdir(exist_ok=True)
            db_path = create_blast_database(genome_file, db_dir)

            blast_output = temp_dir / f'{model_name}_{genome_file.stem}_blast.tab'
            hits = blast_seed_against_genome(
                seed_file,
                db_path,
                blast_output,
                min_identity,
                num_threads=threads,
                evalue=evalue,
            )
            all_hits.extend(hits)

    logger.info(f'Total BLAST hits for {model_name}: {len(all_hits)}')

    if not all_hits:
        raise HMMBuildError(f'No BLAST hits found for {model_name}')

    # Save all BLAST hits in tabular format if requested
    saved_blast_hits_file = None
    if save_blast_hits:
        saved_blast_hits_file = output_dir / f'{cleanID(model_name)}_all_blast_hits.tab'
        save_all_blast_hits(all_hits, saved_blast_hits_file)

    # ------------------------------------------------------------------ #
    # Step 2 – Filter and resolve hits                                     #
    # ------------------------------------------------------------------ #
    filtered_hits = filter_hits_by_thresholds(all_hits, min_coverage, min_identity)

    if not filtered_hits:
        raise HMMBuildError(f'No hits passed thresholds for {model_name}')

    resolved_hits = resolve_overlapping_hits(filtered_hits)

    # Pre-compute the sequence source and single-element chains
    source = lazy_source.get()
    # Each resolved hit forms its own independent chain (no fragmented-hit chaining)
    hit_chains = [[h] for h in resolved_hits]

    # ------------------------------------------------------------------ #
    # Step 3 – Extract sequences                                           #
    # ------------------------------------------------------------------ #
    try:
        sequences = extract_sequences_from_chains(hit_chains, source, model_name)
    except HMMBuildError:
        raise
    except Exception as e:
        raise HMMBuildError(f'Failed to extract sequences: {e}') from e

    # extract_sequences_from_chains logs and skips each per-hit failure, so a
    # wholesale failure -- a BLAST database built without -parse_seqids, or a
    # genome whose contig names do not match the hits -- returned an empty list
    # here. The seed records were then appended and, for a multi-record seed
    # file, len(unique_sequences) >= 2 passed: the run exited 0 having built an
    # "HMM" from the seeds aligned to each other, with no genomic evidence.
    if resolved_hits and not sequences:
        raise HMMBuildError(
            f'None of the {len(resolved_hits)} hit(s) for {model_name} could be '
            'extracted. Check that the sequence source contains the hit target '
            'names; a BLAST database must be built with -parse_seqids for '
            'blastdbcmd to retrieve by sequence ID.'
        )

    # Add original seed sequence(s) - convert to uppercase and add BEFORE deduplication
    seed_records = list(SeqIO.parse(seed_file, 'fasta'))
    logger.info(f'Adding {len(seed_records)} seed sequence(s) to extracted sequences')

    for seed_record in seed_records:
        uppercase_seed = SeqRecord(
            Seq(str(seed_record.seq).upper()),
            id=seed_record.id,
            description=seed_record.description,
        )
        sequences.append(uppercase_seed)

    logger.info(
        f'Total sequences before deduplication: {len(sequences)} '
        f'(including {len(seed_records)} seed sequences)'
    )

    unique_sequences = deduplicate_sequences(sequences)

    if len(unique_sequences) < 2:
        raise HMMBuildError(f'Not enough unique sequences for {model_name} alignment')

    # ------------------------------------------------------------------ #
    # Step 4 – Alignment and HMM building                                  #
    # ------------------------------------------------------------------ #
    alignment_file = output_dir / f'{cleanID(model_name)}_alignment.fasta'
    run_mafft_alignment(unique_sequences, alignment_file, threads=threads)

    # Create flanked alignment if requested
    flanked_alignment_file = None
    if flank_size is not None and flank_size > 0:
        logger.info(f'Creating flanked alignment with {flank_size}bp flanking sequence')

        try:
            flanked_sequences = extract_flanked_sequences_from_chains(
                hit_chains,
                source,
                model_name,
                flank_size,
            )

            for seed_record in seed_records:
                uppercase_seed = SeqRecord(
                    Seq(str(seed_record.seq).upper()),
                    id=f'{seed_record.id}_seed',
                    description=f'{seed_record.description} (seed sequence)',
                )
                flanked_sequences.append(uppercase_seed)

            unique_flanked_sequences = deduplicate_sequences(flanked_sequences)

            if len(unique_flanked_sequences) >= 2:
                flanked_alignment_file = (
                    output_dir
                    / f'{cleanID(model_name)}_flanked_{flank_size}bp_alignment.fasta'
                )
                run_mafft_alignment(
                    unique_flanked_sequences, flanked_alignment_file, threads=threads
                )
                logger.info(f'Flanked alignment written to {flanked_alignment_file}')
            else:
                logger.warning(
                    f'Not enough unique flanked sequences for {model_name} flanked alignment'
                )

        except Exception as e:
            logger.warning(f'Failed to create flanked alignment: {e}')

    # Calculate pairwise identity for standard alignment
    identity_matrix = calculate_pairwise_identity(alignment_file)
    if not identity_matrix.empty:
        identity_file = output_dir / f'{cleanID(model_name)}_pairwise_identity.csv'
        identity_matrix.to_csv(identity_file)
        logger.info(f'Pairwise identity matrix written to {identity_file}')

    # Build HMM from standard alignment (not flanked) using pyhmmer
    hmm_file = build_hmm_from_alignment_pyhmmer(alignment_file, model_name, output_dir)

    return hmm_file, alignment_file, saved_blast_hits_file, flanked_alignment_file


def process_asymmetric_seeds(
    left_seed: Path,
    right_seed: Path,
    model_name: str,
    genome_files: List[Path],
    temp_dir: Path,
    output_dir: Path,
    min_coverage: float,
    min_identity: float,
    save_blast_hits: bool = False,
    flank_size: Optional[int] = None,
    threads: int = 1,
    evalue: float = 1e-3,
    blast_db: Optional[Path] = None,
    left_blast_hits_file: Optional[Path] = None,
    right_blast_hits_file: Optional[Path] = None,
) -> Tuple[Path, Path, Path, Path]:
    """
    Process asymmetric left and right seeds together to avoid filtering conflicts.

    Parameters
    ----------
    left_seed : Path
        Path to left seed sequence file.
    right_seed : Path
        Path to right seed sequence file.
    model_name : str
        Base name for the models.
    genome_files : list of Path
        List of genome FASTA files to search against.  May be empty when
        ``blast_db`` is provided.
    temp_dir : Path
        Directory for temporary files.
    output_dir : Path
        Directory for output files.
    min_coverage : float
        Minimum coverage threshold (0.0 to 1.0).
    min_identity : float
        Minimum percent identity threshold (0.0 to 100.0).
    save_blast_hits : bool, default False
        If True, save all BLAST hits to file.
    flank_size : int, optional
        Size of flanking sequence to extract.
    threads : int, default 1
        Number of CPU threads for BLAST.
    evalue : float, default 1e-3
        E-value threshold for BLAST hits.
    blast_db : Path, optional
        Path to a pre-built BLAST database (without extension).  When provided,
        blastn is run against this database and hit sequences are extracted via
        ``blastdbcmd``.  The database must have been created with
        ``-parse_seqids``.
    left_blast_hits_file : Path, optional
        Pre-computed BLAST hit table for the left seed.  When provided the
        BLAST search for the left seed is skipped.
    right_blast_hits_file : Path, optional
        Pre-computed BLAST hit table for the right seed.  When provided the
        BLAST search for the right seed is skipped.

    Returns
    -------
    tuple
        (left_hmm_file, right_hmm_file, left_alignment_file, right_alignment_file).
    """
    logger.info(f'Processing asymmetric seeds for {model_name}')

    # One source for BOTH seeds, built on first use. The left and right ID
    # checks and both extraction passes previously built their own, so an
    # asymmetric run indexed the genome three times over.
    lazy_source = _LazySource(blast_db, genome_files)

    # ------------------------------------------------------------------ #
    # Step 1 – Obtain BLAST hits for left seed                             #
    # ------------------------------------------------------------------ #
    if left_blast_hits_file is not None:
        logger.info(f'Loading pre-computed left BLAST hits from {left_blast_hits_file}')
        all_left_hits = parse_blast_output(left_blast_hits_file)
        logger.info(f'Loaded {len(all_left_hits)} left hits')
        if not all_left_hits:
            raise HMMBuildError(
                f'No hits in provided left-blast-hits file: {left_blast_hits_file}'
            )
        warn_multiple_queries(all_left_hits, context=f'{model_name}_left')
        if lazy_source.available:
            missing = check_targets_in_blastdb(all_left_hits, lazy_source.get())
            if missing:
                logger.warning(
                    f'Left hit targets could not be resolved in the sequence '
                    f'source: {", ".join(sorted(set(missing)))}'
                )
    elif blast_db is not None:
        logger.info(f'Running left BLAST against pre-built database: {blast_db}')
        left_output = temp_dir / f'{model_name}_left_blast.tab'
        all_left_hits = blast_seed_against_genome(
            left_seed,
            blast_db,
            left_output,
            min_identity,
            num_threads=threads,
            evalue=evalue,
        )
    else:
        all_left_hits = []
        for genome_file in genome_files:
            db_dir = temp_dir / f'blast_dbs_{genome_file.stem}'
            db_dir.mkdir(exist_ok=True)
            db_path = create_blast_database(genome_file, db_dir)
            left_output = temp_dir / f'{model_name}_left_{genome_file.stem}_blast.tab'
            hits = blast_seed_against_genome(
                left_seed,
                db_path,
                left_output,
                min_identity,
                num_threads=threads,
                evalue=evalue,
            )
            all_left_hits.extend(hits)

    # ------------------------------------------------------------------ #
    # Step 2 – Obtain BLAST hits for right seed                            #
    # ------------------------------------------------------------------ #
    if right_blast_hits_file is not None:
        logger.info(
            f'Loading pre-computed right BLAST hits from {right_blast_hits_file}'
        )
        all_right_hits = parse_blast_output(right_blast_hits_file)
        logger.info(f'Loaded {len(all_right_hits)} right hits')
        if not all_right_hits:
            raise HMMBuildError(
                f'No hits in provided right-blast-hits file: {right_blast_hits_file}'
            )
        warn_multiple_queries(all_right_hits, context=f'{model_name}_right')
        if lazy_source.available:
            missing = check_targets_in_blastdb(all_right_hits, lazy_source.get())
            if missing:
                logger.warning(
                    f'Right hit targets could not be resolved in the sequence '
                    f'source: {", ".join(sorted(set(missing)))}'
                )
    elif blast_db is not None:
        logger.info(f'Running right BLAST against pre-built database: {blast_db}')
        right_output = temp_dir / f'{model_name}_right_blast.tab'
        all_right_hits = blast_seed_against_genome(
            right_seed,
            blast_db,
            right_output,
            min_identity,
            num_threads=threads,
            evalue=evalue,
        )
    else:
        all_right_hits = []
        for genome_file in genome_files:
            db_dir = temp_dir / f'blast_dbs_{genome_file.stem}'
            db_dir.mkdir(exist_ok=True)
            db_path = create_blast_database(genome_file, db_dir)
            right_output = temp_dir / f'{model_name}_right_{genome_file.stem}_blast.tab'
            hits = blast_seed_against_genome(
                right_seed,
                db_path,
                right_output,
                min_identity,
                num_threads=threads,
                evalue=evalue,
            )
            all_right_hits.extend(hits)

    if not all_left_hits:
        raise HMMBuildError(f'No BLAST hits found for {model_name}_left')
    if not all_right_hits:
        raise HMMBuildError(f'No BLAST hits found for {model_name}_right')

    # Save all BLAST hits if requested
    if save_blast_hits:
        left_blast_file = output_dir / f'{cleanID(model_name)}_left_all_blast_hits.tab'
        right_blast_file = (
            output_dir / f'{cleanID(model_name)}_right_all_blast_hits.tab'
        )
        save_all_blast_hits(all_left_hits, left_blast_file)
        save_all_blast_hits(all_right_hits, right_blast_file)

    # ------------------------------------------------------------------ #
    # Step 3 – Filter and resolve hits                                     #
    # ------------------------------------------------------------------ #
    filtered_left = filter_hits_by_thresholds(all_left_hits, min_coverage, min_identity)
    filtered_right = filter_hits_by_thresholds(
        all_right_hits, min_coverage, min_identity
    )

    if not filtered_left:
        raise HMMBuildError(f'No hits passed thresholds for {model_name}_left')
    if not filtered_right:
        raise HMMBuildError(f'No hits passed thresholds for {model_name}_right')

    resolved_left = resolve_overlapping_hits(filtered_left)
    resolved_right = resolve_overlapping_hits(filtered_right)

    # Resolve conflicts between left and right hits
    resolved_left, resolved_right = resolve_asymmetric_conflicts(
        resolved_left, resolved_right
    )

    # ------------------------------------------------------------------ #
    # Step 4 – Extract sequences                                           #
    # ------------------------------------------------------------------ #
    # Determine the sequence source for extraction
    source = lazy_source.get()
    # Each resolved hit forms its own independent chain (no fragmented-hit chaining)
    left_chains = [[h] for h in resolved_left]
    right_chains = [[h] for h in resolved_right]

    # Process left seed
    logger.info(f'Processing left seed sequences for {model_name}_left')
    try:
        left_sequences = extract_sequences_from_chains(
            left_chains, source, f'{model_name}_left'
        )
    except HMMBuildError:
        raise
    except Exception as e:
        raise HMMBuildError(f'Failed to extract left sequences: {e}') from e

    # See the equivalent guard in process_seed_sequences: a wholesale
    # extraction failure must not be papered over by the seed records.
    if resolved_left and not left_sequences:
        raise HMMBuildError(
            f'None of the {len(resolved_left)} left hit(s) for {model_name} '
            'could be extracted. Check that the sequence source contains the hit '
            'target names; a BLAST database must be built with -parse_seqids.'
        )

    left_seed_records = list(SeqIO.parse(left_seed, 'fasta'))
    for seed_record in left_seed_records:
        left_sequences.append(
            SeqRecord(
                Seq(str(seed_record.seq).upper()),
                id=seed_record.id,
                description=seed_record.description,
            )
        )

    unique_left_sequences = deduplicate_sequences(left_sequences)

    if len(unique_left_sequences) < 2:
        raise HMMBuildError(
            f'Not enough unique sequences for {model_name}_left alignment'
        )

    left_alignment_file = output_dir / f'{cleanID(model_name)}_left_alignment.fasta'
    run_mafft_alignment(unique_left_sequences, left_alignment_file, threads=threads)

    left_hmm_file = build_hmm_from_alignment_pyhmmer(
        left_alignment_file, f'{model_name}_left', output_dir
    )

    # Process right seed
    logger.info(f'Processing right seed sequences for {model_name}_right')
    try:
        right_sequences = extract_sequences_from_chains(
            right_chains,
            source,
            f'{model_name}_right',
        )
    except HMMBuildError:
        raise
    except Exception as e:
        raise HMMBuildError(f'Failed to extract right sequences: {e}') from e

    # See the equivalent guard in process_seed_sequences: a wholesale
    # extraction failure must not be papered over by the seed records.
    if resolved_right and not right_sequences:
        raise HMMBuildError(
            f'None of the {len(resolved_right)} right hit(s) for {model_name} '
            'could be extracted. Check that the sequence source contains the hit '
            'target names; a BLAST database must be built with -parse_seqids.'
        )

    right_seed_records = list(SeqIO.parse(right_seed, 'fasta'))
    for seed_record in right_seed_records:
        right_sequences.append(
            SeqRecord(
                Seq(str(seed_record.seq).upper()),
                id=seed_record.id,
                description=seed_record.description,
            )
        )

    unique_right_sequences = deduplicate_sequences(right_sequences)

    if len(unique_right_sequences) < 2:
        raise HMMBuildError(
            f'Not enough unique sequences for {model_name}_right alignment'
        )

    right_alignment_file = output_dir / f'{cleanID(model_name)}_right_alignment.fasta'
    run_mafft_alignment(unique_right_sequences, right_alignment_file, threads=threads)

    right_hmm_file = build_hmm_from_alignment_pyhmmer(
        right_alignment_file, f'{model_name}_right', output_dir
    )

    # ------------------------------------------------------------------ #
    # Step 5 – Pairwise identity matrices                                  #
    # ------------------------------------------------------------------ #
    left_identity_matrix = calculate_pairwise_identity(left_alignment_file)
    if not left_identity_matrix.empty:
        left_identity_file = (
            output_dir / f'{cleanID(model_name)}_left_pairwise_identity.csv'
        )
        left_identity_matrix.to_csv(left_identity_file)
        logger.info(f'Left pairwise identity matrix written to {left_identity_file}')

    right_identity_matrix = calculate_pairwise_identity(right_alignment_file)
    if not right_identity_matrix.empty:
        right_identity_file = (
            output_dir / f'{cleanID(model_name)}_right_pairwise_identity.csv'
        )
        right_identity_matrix.to_csv(right_identity_file)
        logger.info(f'Right pairwise identity matrix written to {right_identity_file}')

    # ------------------------------------------------------------------ #
    # Step 6 – Optional flanked alignments                                 #
    # ------------------------------------------------------------------ #
    if flank_size is not None and flank_size > 0:
        logger.info(
            f'Creating flanked alignments with {flank_size}bp flanking sequence'
        )

        try:
            left_flanked_sequences = extract_flanked_sequences_from_chains(
                left_chains, source, f'{model_name}_left', flank_size
            )
            right_flanked_sequences = extract_flanked_sequences_from_chains(
                right_chains, source, f'{model_name}_right', flank_size
            )

            for seed_record in left_seed_records:
                left_flanked_sequences.append(
                    SeqRecord(
                        Seq(str(seed_record.seq).upper()),
                        id=f'{seed_record.id}_seed',
                        description=f'{seed_record.description} (seed sequence)',
                    )
                )

            for seed_record in right_seed_records:
                right_flanked_sequences.append(
                    SeqRecord(
                        Seq(str(seed_record.seq).upper()),
                        id=f'{seed_record.id}_seed',
                        description=f'{seed_record.description} (seed sequence)',
                    )
                )

            unique_left_flanked = deduplicate_sequences(left_flanked_sequences)
            unique_right_flanked = deduplicate_sequences(right_flanked_sequences)

            if len(unique_left_flanked) >= 2:
                left_flanked_file = (
                    output_dir
                    / f'{cleanID(model_name)}_left_flanked_{flank_size}bp_alignment.fasta'
                )
                run_mafft_alignment(
                    unique_left_flanked, left_flanked_file, threads=threads
                )
                logger.info(f'Left flanked alignment written to {left_flanked_file}')

            if len(unique_right_flanked) >= 2:
                right_flanked_file = (
                    output_dir
                    / f'{cleanID(model_name)}_right_flanked_{flank_size}bp_alignment.fasta'
                )
                run_mafft_alignment(
                    unique_right_flanked, right_flanked_file, threads=threads
                )
                logger.info(f'Right flanked alignment written to {right_flanked_file}')

        except Exception as e:
            logger.warning(f'Failed to create flanked alignments: {e}')

    logger.info(f'Asymmetric processing completed for {model_name}')
    logger.info(f'Left HMM: {left_hmm_file}')
    logger.info(f'Right HMM: {right_hmm_file}')

    return left_hmm_file, right_hmm_file, left_alignment_file, right_alignment_file


def resolve_asymmetric_conflicts(
    left_hits: List[BlastHit], right_hits: List[BlastHit]
) -> Tuple[List[BlastHit], List[BlastHit]]:
    """
    Resolve conflicts between left and right seed hits.

    Parameters
    ----------
    left_hits : list of BlastHit
        BLAST hits from left seed.
    right_hits : list of BlastHit
        BLAST hits from right seed.

    Returns
    -------
    tuple
        (filtered_left_hits, filtered_right_hits) with conflicts resolved,
        prioritizing longer and higher-identity hits.
    """
    # Group hits by chromosome
    left_by_chrom: Dict[str, List[Any]] = {}
    right_by_chrom: Dict[str, List[Any]] = {}

    for hit in left_hits:
        if hit.subject_id not in left_by_chrom:
            left_by_chrom[hit.subject_id] = []
        left_by_chrom[hit.subject_id].append(hit)

    for hit in right_hits:
        if hit.subject_id not in right_by_chrom:
            right_by_chrom[hit.subject_id] = []
        right_by_chrom[hit.subject_id].append(hit)

    # For each chromosome, resolve overlaps between left and right
    filtered_left = []
    filtered_right = []

    all_chroms = set(left_by_chrom.keys()) | set(right_by_chrom.keys())

    for chrom in all_chroms:
        chrom_left = left_by_chrom.get(chrom, [])
        chrom_right = right_by_chrom.get(chrom, [])

        # Left hits surviving this chromosome's conflict resolution.
        kept_left: List[BlastHit] = []

        # Check for overlaps between left and right hits
        for left_hit in chrom_left:
            keep_left = True
            for right_hit in chrom_right:
                if hits_overlap(left_hit, right_hit):
                    # Decide which to keep based on quality
                    left_score = (
                        left_hit.length,
                        left_hit.identity,
                        left_hit.query_coverage,
                    )
                    right_score = (
                        right_hit.length,
                        right_hit.identity,
                        right_hit.query_coverage,
                    )

                    if right_score > left_score:
                        keep_left = False
                        logger.debug(
                            f'Left hit {left_hit} removed due to better right hit {right_hit}'
                        )
                        break

            if keep_left:
                kept_left.append(left_hit)

        filtered_left.extend(kept_left)

        # Similar logic for right hits.
        #
        # This iterates the left hits that SURVIVED the loop above, not the
        # original chrom_left. Using the original allowed a left hit to be
        # discarded in the first loop and still evict a right hit in the
        # second, so a legitimate right-terminus hit could be deleted by a hit
        # that no longer existed. Verified against the previous behaviour: a
        # right hit overlapping only a discarded left hit was silently lost.
        for right_hit in chrom_right:
            keep_right = True
            for left_hit in kept_left:
                if hits_overlap(left_hit, right_hit):
                    left_score = (
                        left_hit.length,
                        left_hit.identity,
                        left_hit.query_coverage,
                    )
                    right_score = (
                        right_hit.length,
                        right_hit.identity,
                        right_hit.query_coverage,
                    )

                    if left_score >= right_score:  # Note: >= to avoid double removal
                        keep_right = False
                        logger.debug(
                            f'Right hit {right_hit} removed due to better left hit {left_hit}'
                        )
                        break

            if keep_right:
                filtered_right.append(right_hit)

    return filtered_left, filtered_right


def hits_overlap(hit1: BlastHit, hit2: BlastHit, min_overlap: int = 50) -> bool:
    """
    Check if two hits overlap significantly.

    Parameters
    ----------
    hit1 : BlastHit
        First BLAST hit.
    hit2 : BlastHit
        Second BLAST hit.
    min_overlap : int, default 50
        Minimum overlap in base pairs to consider significant.

    Returns
    -------
    bool
        True if hits overlap by at least min_overlap bases, False otherwise.
    """
    if hit1.subject_id != hit2.subject_id:
        return False

    start1, end1 = (
        min(hit1.subject_start, hit1.subject_end),
        max(hit1.subject_start, hit1.subject_end),
    )
    start2, end2 = (
        min(hit2.subject_start, hit2.subject_end),
        max(hit2.subject_start, hit2.subject_end),
    )

    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    overlap_length = max(0, overlap_end - overlap_start)

    return overlap_length >= min_overlap


def update_hmm_with_genome_hits(
    hmm_file: Path,
    genome_files: List[Path],
    model_name: str,
    temp_dir: Path,
    output_dir: Path,
    flank_size: Optional[int] = None,
    evalue: float = 1e-3,
    threads: int = 1,
) -> Tuple[Path, Path, Optional[Path]]:
    """
    Update an existing HMM by searching genomes, extracting hits, and rebuilding.

    Parameters
    ----------
    hmm_file : Path
        Path to existing HMM file to update.
    genome_files : list of Path
        List of genome files to search against.
    model_name : str
        Base name for the updated HMM model.
    temp_dir : Path
        Temporary directory for intermediate files.
    output_dir : Path
        Output directory for final results.
    flank_size : int, optional
        Size of flanking sequence to extract with hits.
    evalue : float, default 1e-3
        E-value threshold for nhmmer search.
    threads : int, default 1
        Number of CPU threads for alignment.

    Returns
    -------
    tuple[Path, Path, Path or None]
        - Path to updated HMM file
        - Path to alignment file
        - Path to flanked alignment file (if flank_size specified)

    Raises
    ------
    HMMBuildError
        If HMM update workflow fails.
    FileNotFoundError
        If HMM or genome files don't exist.

    Notes
    -----
    Workflow:
    1. Run nhmmer with existing HMM against genomes
    2. Parse hits and extract unique sequences
    3. Align sequences to HMM using hmmalign
    4. Build updated HMM from alignment
    5. Optionally create flanked alignment
    """
    if not hmm_file.exists():
        raise FileNotFoundError(f'HMM file not found: {hmm_file}')

    logger.info(f'Starting HMM update workflow with {hmm_file.name}')
    logger.info(f'Searching {len(genome_files)} genome(s)')

    # Collect all hits from all genomes
    all_hit_sequences: List[SeqRecord] = []
    total_hits = 0

    for genome_file in genome_files:
        logger.info(f'Processing genome: {genome_file.name}')

        # Prepare genome (decompress if needed)
        try:
            prepared_genome = prepare_genome_file(genome_file, temp_dir)
        except Exception as e:
            logger.error(f'Failed to prepare genome {genome_file}: {e}')
            continue

        # Run nhmmer search
        try:
            hit_sequences = search_and_extract_hits(
                hmm_file=hmm_file,
                genome_file=prepared_genome,
                temp_dir=temp_dir,
                evalue=evalue,
                flank_size=flank_size,
            )

            genome_hit_count = len(hit_sequences)
            total_hits += genome_hit_count
            all_hit_sequences.extend(hit_sequences)

            logger.info(f'  Found {genome_hit_count} hits in {genome_file.name}')

        except Exception as e:
            logger.error(f'Failed to search genome {genome_file}: {e}')
            continue

    if not all_hit_sequences:
        raise HMMBuildError('No hits found in any genome. Cannot update HMM.')

    # Deduplicate sequences
    unique_sequences = deduplicate_sequences(all_hit_sequences)
    unique_count = len(unique_sequences)

    logger.info(f'Total hits found: {total_hits}')
    logger.info(f'Unique sequences after deduplication: {unique_count}')

    # Save sequences to file for alignment
    sequences_file = temp_dir / f'{cleanID(model_name)}_hits.fasta'
    SeqIO.write(unique_sequences, sequences_file, 'fasta')
    logger.info(f'Saved {unique_count} unique sequences to {sequences_file.name}')

    # Align sequences to HMM using hmmalign
    alignment_file = align_sequences_to_hmm(
        hmm_file=hmm_file,
        sequences_file=sequences_file,
        temp_dir=temp_dir,
        model_name=model_name,
    )

    # Build updated HMM from alignment
    updated_hmm = build_hmm_from_alignment_pyhmmer(
        alignment_file=alignment_file,
        model_name=model_name,
        output_dir=output_dir,
    )

    logger.info(f'Updated HMM saved to: {updated_hmm}')

    # Copy alignment to output directory
    output_alignment = output_dir / alignment_file.name
    shutil.copy2(alignment_file, output_alignment)
    logger.info(f'Alignment saved to: {output_alignment}')

    # Create flanked alignment if requested
    flanked_alignment = None
    if flank_size:
        flanked_alignment = create_flanked_alignment_output(
            unique_sequences=unique_sequences,
            model_name=model_name,
            output_dir=output_dir,
            temp_dir=temp_dir,
            threads=threads,
        )
        if flanked_alignment:
            logger.info(f'Flanked alignment saved to: {flanked_alignment}')

    return updated_hmm, output_alignment, flanked_alignment


def search_and_extract_hits(
    hmm_file: Path,
    genome_file: Path,
    temp_dir: Path,
    evalue: float = 1e-3,
    flank_size: Optional[int] = None,
) -> List[SeqRecord]:
    """
    Search genome with HMM and extract hit sequences.

    Parameters
    ----------
    hmm_file : Path
        Path to HMM file.
    genome_file : Path
        Path to genome file.
    temp_dir : Path
        Temporary directory for intermediate files.
    evalue : float, default 1e-3
        E-value threshold for nhmmer.
    flank_size : int, optional
        Size of flanking sequence to extract.

    Returns
    -------
    list of SeqRecord
        List of extracted hit sequences.

    Raises
    ------
    HMMBuildError
        If search or extraction fails.
    """
    # Press HMM (index it)
    press_command = build_hmmpress_command(hmm_file)
    try:
        result = run_command(press_command, verbose=False)
        if result.returncode != 0:
            raise HMMBuildError(f'hmmpress failed: {result.stderr}')
    except Exception as e:
        raise HMMBuildError(f'Failed to press HMM: {e}') from e

    # Run nhmmer search
    nhmmer_command, results_dir = build_nhmmer_command(
        model_path=hmm_file,
        genome_path=genome_file,
        output_dir=temp_dir,
        evalue=evalue,
    )

    try:
        result = run_command(nhmmer_command, verbose=False)
        if result.returncode != 0:
            raise HMMBuildError(f'nhmmer failed: {result.stderr}')
    except Exception as e:
        raise HMMBuildError(f'Failed to run nhmmer: {e}') from e

    # Find output file
    output_file = results_dir / f'{hmm_file.stem}.out'
    if not output_file.exists():
        logger.warning(f'nhmmer output file not found: {output_file}')
        return []

    # Parse nhmmer output
    try:
        hits_df = import_nhmmer(str(output_file))
    except Exception as e:
        logger.warning(f'Failed to parse nhmmer output: {e}')
        return []

    if hits_df is None or hits_df.empty:
        logger.debug(f'No hits found in {genome_file.name}')
        return []

    logger.debug(f'Parsed {len(hits_df)} hits from nhmmer output')

    # Index genome for sequence extraction
    try:
        genome_index, descriptions = indexGenome(genome_file)
    except Exception as e:
        raise HMMBuildError(f'Failed to index genome: {e}') from e

    source = FastaSource(genome_index, descriptions)

    # Extract sequences from hits
    hit_sequences: List[SeqRecord] = []

    for idx, row in hits_df.iterrows():
        try:
            seq_id = row['target']
            # nhmmer coordinates are 1-based inclusive, as fetch_sequence expects
            start = int(row['hitStart'])
            end = int(row['hitEnd'])
            strand = row['strand']

            # Adjust coordinates for flanking sequence
            if flank_size:
                flank_start = start - flank_size
                flank_end = end + flank_size
            else:
                flank_start = start
                flank_end = end

            # Extract sequence (clamped to the contig by fetch_sequence)
            seq_str = fetch_sequence(source, seq_id, flank_start, flank_end)
            if seq_str is None:
                logger.warning(
                    f'Failed to extract sequence for hit {idx} '
                    f'({seq_id}:{flank_start}-{flank_end})'
                )
                continue

            # Reverse complement if on minus strand
            if strand == '-':
                seq_obj = Seq(seq_str)
                seq_str = str(seq_obj.reverse_complement())

            # Create SeqRecord
            hit_id = f'{seq_id}_{start}_{end}_{strand}'
            record = SeqRecord(
                Seq(seq_str),
                id=hit_id,
                description=f'HMM_hit {seq_id}:{start}-{end}({strand})',
            )

            hit_sequences.append(record)

        except Exception as e:
            logger.warning(f'Failed to extract sequence for hit {idx}: {e}')
            continue

    return hit_sequences


def align_sequences_to_hmm(
    hmm_file: Path,
    sequences_file: Path,
    temp_dir: Path,
    model_name: str,
) -> Path:
    """
    Align sequences to HMM using hmmalign.

    Parameters
    ----------
    hmm_file : Path
        Path to HMM file.
    sequences_file : Path
        Path to FASTA file with sequences to align.
    temp_dir : Path
        Temporary directory for output.
    model_name : str
        Model name for output files.

    Returns
    -------
    Path
        Path to Stockholm format alignment file.

    Raises
    ------
    HMMBuildError
        If hmmalign fails.
    """
    # Output alignment file (Stockholm format)
    alignment_file = temp_dir / f'{cleanID(model_name)}_updated_alignment.sto'

    # Build hmmalign command
    align_command = build_hmmalign_command(
        hmm_file=hmm_file,
        sequence_file=sequences_file,
        output_file=alignment_file,
        trim=False,
    )

    logger.info('Aligning sequences to HMM with hmmalign...')

    try:
        result = run_command(align_command, verbose=False)
        if result.returncode != 0:
            raise HMMBuildError(f'hmmalign failed: {result.stderr}')
    except Exception as e:
        raise HMMBuildError(f'Failed to align sequences to HMM: {e}') from e

    if not alignment_file.exists():
        raise HMMBuildError(f'Alignment file not created: {alignment_file}')

    logger.info(f'Alignment complete: {alignment_file.name}')

    return alignment_file


def create_flanked_alignment_output(
    unique_sequences: List[SeqRecord],
    model_name: str,
    output_dir: Path,
    temp_dir: Path,
    threads: int = 1,
) -> Optional[Path]:
    """
    Create alignment with flanked sequences using MAFFT.

    Parameters
    ----------
    unique_sequences : list of SeqRecord
        Sequences with flanks already included.
    model_name : str
        Model name for output files.
    output_dir : Path
        Output directory.
    temp_dir : Path
        Temporary directory.
    threads : int, default 1
        Number of threads for MAFFT.

    Returns
    -------
    Path or None
        Path to flanked alignment file, or None if creation failed.
    """
    try:
        # Prepare output file path
        output_file = temp_dir / f'{cleanID(model_name)}_flanked_alignment.fasta'

        # Run MAFFT alignment
        flanked_alignment = run_mafft_alignment(
            sequences=unique_sequences,
            output_file=output_file,
            threads=threads,
        )

        # Copy to output directory
        output_flanked = output_dir / flanked_alignment.name
        shutil.copy2(flanked_alignment, output_flanked)

        return output_flanked

    except Exception as e:
        logger.warning(f'Failed to create flanked alignment: {e}')
        return None


def create_seed_parser() -> argparse.ArgumentParser:
    """
    Create standalone argument parser for seed command.

    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser for seed workflow options.
    """
    parser = argparse.ArgumentParser(
        prog='tirmite-seed',
        description='Build HMM models from seed sequences for TIRmite',
    )
    _configure_seed_parser(parser)
    return parser


def add_seed_parser(subparsers: Any) -> argparse.ArgumentParser:
    """
    Add seed subcommand parser.

    Parameters
    ----------
    subparsers : argparse._SubParsersAction
        Subparser object to add seed command to.

    Returns
    -------
    argparse.ArgumentParser
        The configured seed subcommand parser.
    """
    parser = cast(
        argparse.ArgumentParser,
        subparsers.add_parser(
            'seed',
            help='Build HMM models from seed sequences',
            description='Build HMM models from seed sequences for TIRmite',
        ),
    )
    _configure_seed_parser(parser)
    return parser


def _configure_seed_parser(parser: argparse.ArgumentParser) -> None:
    """
    Configure parser with seed command arguments.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser to configure.

    Returns
    -------
    None
        Modifies parser in place.
    """
    # Update mode arguments
    parser.add_argument(
        '--update',
        action='store_true',
        help='Update existing HMM mode: search genomes with existing HMM, extract hits, and rebuild HMM',
    )
    parser.add_argument(
        '--hmm-file',
        type=Path,
        help='Path to existing HMM file to update (required when --update is specified)',
    )

    # Input arguments
    parser.add_argument(
        '--left-seed',
        type=Path,
        help='FASTA file containing left terminal seed sequence(s) (required unless --update is specified)',
    )
    parser.add_argument(
        '--right-seed',
        type=Path,
        help='FASTA file containing right terminal seed sequence(s) (optional for asymmetric elements)',
    )
    parser.add_argument(
        '--model-name', type=str, required=True, help='Base name for the HMM model(s)'
    )
    parser.add_argument(
        '--outdir', type=Path, default=Path.cwd(), help='Output directory for results'
    )

    # Genome input (mutually exclusive with blastdb)
    genome_group = parser.add_mutually_exclusive_group(required=True)
    genome_group.add_argument(
        '--genome', type=Path, help='Single genome FASTA file to search'
    )
    genome_group.add_argument(
        '--genome-list',
        type=Path,
        help='File containing list of genome paths (one per line)',
    )
    genome_group.add_argument(
        '--blastdb',
        type=Path,
        help=(
            'Path to a pre-built BLAST nucleotide database (without file extension). '
            'Provide this instead of --genome / --genome-list when a BLAST database '
            'has already been built. The database must have been created with '
            '``makeblastdb -parse_seqids`` so that sequences can be extracted with '
            '``blastdbcmd``. Mutually exclusive with --genome and --genome-list.'
        ),
    )

    # Optional pre-computed BLAST hit tables
    parser.add_argument(
        '--blast-hits',
        type=Path,
        help=(
            'Pre-computed BLAST hit table for the seed sequence (TIRmite tabular '
            'format, generated with --save-blast-hits). When provided the BLAST '
            'search step is skipped and these hits are used directly. Sequence '
            'extraction still requires --genome, --genome-list, or --blastdb.'
        ),
    )
    parser.add_argument(
        '--left-blast-hits',
        type=Path,
        help=(
            'Pre-computed BLAST hit table for the left seed sequence (asymmetric '
            'elements). Used together with --right-blast-hits to skip the BLAST '
            'search step for both seeds.'
        ),
    )
    parser.add_argument(
        '--right-blast-hits',
        type=Path,
        help=(
            'Pre-computed BLAST hit table for the right seed sequence (asymmetric '
            'elements). Used together with --left-blast-hits to skip the BLAST '
            'search step for both seeds.'
        ),
    )

    # Optional parameters
    parser.add_argument(
        '--tempdir',
        type=Path,
        help='Directory for temporary files (default: system temp)',
    )
    parser.add_argument(
        '--keep-temp', action='store_true', help='Keep temporary files after completion'
    )
    parser.add_argument(
        '--min-coverage',
        type=validate_coverage,
        default=0.7,
        help='Minimum query coverage threshold as fraction 0.0-1.0 (default: 0.7)',
    )
    parser.add_argument(
        '--min-identity',
        type=validate_identity,
        default=70.0,
        help='Minimum sequence identity threshold as percentage (default: 70.0)',
    )
    parser.add_argument(
        '--save-blast-hits',
        action='store_true',
        help='Save all BLAST hits in tabular format',
    )
    parser.add_argument(
        '--flank-size',
        type=int,
        help='Extract BLAST hits with flanking sequence of N bases and create additional alignment (optional)',
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
    parser.add_argument(
        '--threads',
        type=validate_threads,
        default=1,
        help='Number of CPU threads to use for MAFFT alignment (default: 1)',
    )
    parser.add_argument(
        '--evalue',
        type=validate_evalue,
        default=1e-3,
        help='E-value threshold for BLAST hits (default: 1e-3)',
    )


def main(args: Optional[argparse.Namespace] = None) -> int:
    """
    Main function for HMM building workflow.

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
        parser = create_seed_parser()
        args = parser.parse_args()

    # Mypy assertion: args is guaranteed non-None after parsing
    assert args is not None

    # Check available CPU threads
    max_threads = os.cpu_count() or 1
    if args.threads > max_threads:
        logger.warning(
            f'Requested threads ({args.threads}) exceeds available CPUs ({max_threads}). '
            f'Setting threads to {max_threads}.'
        )
        args.threads = max_threads

    # Set up logging. The log file, when requested, lands in the output
    # directory alongside the run's other products.
    logfile_path = None
    if getattr(args, 'logfile', False):
        log_dir = Path(args.outdir) if args.outdir else Path.cwd()
        log_dir.mkdir(parents=True, exist_ok=True)
        model_name = getattr(args, 'model_name', None)
        logfile_name = (
            f'{cleanID(model_name)}_tirmite_seed.log'
            if model_name
            else 'tirmite_seed.log'
        )
        logfile_path = log_dir / logfile_name

    init_logging(loglevel=args.loglevel, logfile=logfile_path)

    try:
        # Check dependencies
        missing_tools = check_dependencies()
        if missing_tools:
            raise HMMBuildError(f'Missing required tools: {", ".join(missing_tools)}')

        # Validate mode-specific requirements
        if args.update:
            # Update mode requires --hmm-file
            if not args.hmm_file:
                raise ValueError('--hmm-file is required when --update is specified')
            if not args.hmm_file.exists():
                raise FileNotFoundError(f'HMM file not found: {args.hmm_file}')
            # Update mode doesn't use --left-seed or --right-seed
            if args.left_seed or args.right_seed:
                logger.warning(
                    '--left-seed and --right-seed are ignored in --update mode'
                )
        else:
            # Normal mode requires --left-seed
            if not args.left_seed:
                raise ValueError('--left-seed is required unless --update is specified')
            if not args.left_seed.exists():
                raise FileNotFoundError(f'Left seed file not found: {args.left_seed}')
            if args.right_seed and not args.right_seed.exists():
                raise FileNotFoundError(f'Right seed file not found: {args.right_seed}')
            # Normal mode doesn't use --hmm-file
            if args.hmm_file:
                logger.warning('--hmm-file is ignored unless --update is specified')

            # Validate pre-computed hit table paths
            for attr, label in [
                ('blast_hits', '--blast-hits'),
                ('left_blast_hits', '--left-blast-hits'),
                ('right_blast_hits', '--right-blast-hits'),
            ]:
                path = getattr(args, attr, None)
                if path is not None and not path.exists():
                    raise FileNotFoundError(f'{label} file not found: {path}')

        # Get genome files / blastdb
        blast_db: Optional[Path] = None
        genome_files: List[Path] = []

        if args.blastdb is not None:
            blast_db = args.blastdb
            logger.info(f'Using pre-built BLAST database: {blast_db}')
            # Verify at least one database file exists (e.g. .nsq, .nhr, .nin)
            db_files = list(blast_db.parent.glob(blast_db.name + '.*'))
            if not db_files:
                raise FileNotFoundError(
                    f'No BLAST database files found for: {blast_db}'
                )
        elif args.genome:
            genome_files = [args.genome]
            if not args.genome.exists():
                raise FileNotFoundError(f'Genome file not found: {args.genome}')
        else:
            with open(args.genome_list, 'r') as f:
                genome_paths = [Path(line.strip()) for line in f if line.strip()]

            for path in genome_paths:
                if not path.exists():
                    logger.warning(f'Genome file not found, skipping: {path}')
                else:
                    genome_files.append(path)

            if not genome_files:
                raise FileNotFoundError('No valid genome files found')

        # Create output directory
        output_dir = create_output_directory(args.outdir)
        logger.info(f'Output directory: {output_dir}')

        # Set up temporary directory
        if args.tempdir:
            # Create the specified temp directory if it doesn't exist
            temp_base_dir = Path(args.tempdir)
            try:
                temp_base_dir.mkdir(parents=True, exist_ok=True)
                logger.info(f'Created/using temporary base directory: {temp_base_dir}')
            except Exception as e:
                raise HMMBuildError(
                    f'Failed to create temporary base directory {temp_base_dir}: {e}'
                ) from e
        else:
            temp_base_dir = None
            logger.info('Using system default temporary directory')

        with temporary_directory(
            prefix='tirmite_build_', base_dir=temp_base_dir, cleanup=not args.keep_temp
        ) as temp_dir:
            logger.info(f'Temporary directory: {temp_dir}')
            logger.info(f'Using {args.threads} CPU threads for alignment')

            # Route to appropriate workflow
            if args.update:
                # Update mode: search with existing HMM, extract hits, rebuild HMM
                if blast_db is not None:
                    raise HMMBuildError(
                        '--blastdb is not supported with --update mode. '
                        'Provide --genome or --genome-list instead.'
                    )
                if not genome_files:
                    raise HMMBuildError(
                        '--update mode requires --genome or --genome-list'
                    )
                logger.info('Running HMM update workflow')

                updated_hmm, alignment, flanked_alignment = update_hmm_with_genome_hits(
                    hmm_file=args.hmm_file,
                    genome_files=genome_files,
                    model_name=args.model_name,
                    temp_dir=temp_dir,
                    output_dir=output_dir,
                    flank_size=args.flank_size,
                    evalue=args.evalue,
                    threads=args.threads,
                )

                logger.info('HMM update workflow completed:')
                logger.info(f'  Updated HMM: {updated_hmm}')
                logger.info(f'  Alignment: {alignment}')
                if flanked_alignment:
                    logger.info(f'  Flanked alignment: {flanked_alignment}')

            else:
                # Normal seed-based workflow
                # Compare seeds if both provided
                if args.right_seed:
                    seed_comparisons = compare_seeds(
                        args.left_seed,
                        args.right_seed,
                        temp_dir,
                        min_length=10,  # Minimum 10bp hits
                        min_identity=50.0,  # Minimum 50% identity
                        num_threads=args.threads,
                    )

                    if seed_comparisons:
                        logger.info(
                            f'Found {len(seed_comparisons)} significant similarities between seeds:'
                        )

                        # Named seed_alignment, not alignment: the --update
                        # branch above uses `alignment` for an output Path.
                        for i, (hit, seed_alignment) in enumerate(
                            seed_comparisons[:3], 1
                        ):  # Show top 3
                            logger.info(
                                f'  Similarity {i}: {hit.length}bp, {hit.identity:.1f}% identity'
                            )
                            logger.info(
                                f'    Query: {hit.query_id}[{hit.query_start}:{hit.query_end}]'
                            )
                            logger.info(
                                f'    Subject: {hit.subject_id}[{hit.subject_start}:{hit.subject_end}] ({hit.strand} strand)'
                            )
                            logger.info(f'    Alignment score: {seed_alignment.score}')

                            # Optionally print the alignment for the best hit
                            if i == 1 and args.loglevel == 'DEBUG':
                                logger.debug('Best seed alignment:')
                                for line in str(seed_alignment).split('\n'):
                                    logger.debug(f'    {line}')

                        # Save detailed seed comparison results
                        seed_comparison_file = (
                            output_dir
                            / f'{cleanID(args.model_name)}_seed_comparison.txt'
                        )
                        write_seed_comparison_report(
                            seed_comparisons,
                            seed_comparison_file,
                            model_name=args.model_name,
                            left_seed_name=args.left_seed.name,
                            right_seed_name=args.right_seed.name,
                        )

                        logger.info(
                            f'Detailed seed comparison saved to: {seed_comparison_file}'
                        )

                    else:
                        logger.info(
                            'No significant similarity found between left and right seeds'
                        )

                        # Still proceed with asymmetric processing
                        logger.info('Proceeding with fully asymmetric seed processing')

                # Process asymmetric seeds if both provided
                if args.right_seed:
                    left_hmm, right_hmm, left_aln, right_aln = process_asymmetric_seeds(
                        args.left_seed,
                        args.right_seed,
                        args.model_name,
                        genome_files,
                        temp_dir,
                        output_dir,
                        args.min_coverage,
                        args.min_identity,
                        args.save_blast_hits,
                        args.flank_size,
                        threads=args.threads,
                        evalue=args.evalue,
                        blast_db=blast_db,
                        left_blast_hits_file=getattr(args, 'left_blast_hits', None),
                        right_blast_hits_file=getattr(args, 'right_blast_hits', None),
                    )

                    logger.info('Asymmetric processing completed:')
                    logger.info(f'  Left HMM: {left_hmm}')
                    logger.info(f'  Right HMM: {right_hmm}')
                    logger.info(f'  Left alignment: {left_aln}')
                    logger.info(f'  Right alignment: {right_aln}')

                else:
                    # Process single seed (existing logic)
                    left_hmm, left_aln, left_blast, left_flanked = (
                        process_seed_sequences(
                            args.left_seed,
                            args.model_name,
                            genome_files,
                            temp_dir,
                            output_dir,
                            args.min_coverage,
                            args.min_identity,
                            args.save_blast_hits,
                            args.flank_size,
                            threads=args.threads,
                            evalue=args.evalue,
                            blast_db=blast_db,
                            blast_hits_file=getattr(args, 'blast_hits', None),
                        )
                    )
                    logger.info(f'Single seed processing completed: {left_hmm}')

    except Exception as e:
        logger.error(f'HMM building failed: {e}')
        return 1

    return 0


if __name__ == '__main__':
    exit(main())
