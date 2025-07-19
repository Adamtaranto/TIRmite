#!/usr/bin/env python3
"""
HMM building from seed sequences for TIRmite.

This module builds HMM models from seed sequences by:
1. Comparing left/right seeds if both provided
2. BLASTing seeds against target genome(s)
3. Filtering and processing hits
4. Extracting and aligning sequences
5. Building HMM models
"""

import argparse
import io
import logging
from pathlib import Path
import shutil
from typing import List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

from tirmite.hmmer_wrappers import build_hmmbuild_command
from tirmite.logs import init_logging
from tirmite.runBlastn import BlastError, run_blastn
from tirmite.utils import (
    cleanID,
    create_output_directory,
    indexGenome,
    temporary_directory,
)
from tirmite.wrapping import run_command


class HMMBuildError(Exception):
    """Custom exception for HMM building errors."""

    pass


class BlastHit:
    """Container for BLAST hit information."""

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

        # Determine strand from coordinates and frame
        self.strand = '+' if subject_start <= subject_end else '-'

    def __repr__(self):
        return (
            f'BlastHit({self.query_id}->{self.subject_id}, '
            f'cov={self.query_coverage:.2f}, id={self.identity:.1f}%, strand={self.strand})'
        )


def check_dependencies() -> List[str]:
    """Check if required external tools are available."""
    required_tools = ['blastn', 'makeblastdb', 'mafft', 'hmmbuild']
    missing = []

    for tool in required_tools:
        if not shutil.which(tool):
            missing.append(tool)

    return missing


def parse_blast_output(blast_file: Path) -> List[BlastHit]:
    """
    Parse BLAST tabular output into BlastHit objects.

    Expected format: qstart qend sstart send length positive pident qlen slen qframe sframe qseqid sseqid
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
                        logging.warning(
                            f'Skipping malformed BLAST line {line_num}: {line.strip()} - {e}'
                        )
                        continue
                else:
                    logging.warning(
                        f'Skipping BLAST line {line_num} with insufficient fields: {line.strip()}'
                    )

    except Exception as e:
        raise HMMBuildError(f'Failed to parse BLAST output {blast_file}: {e}') from e

    return hits


def compare_seeds(
    left_seed: Path, right_seed: Path, temp_dir: Path
) -> Optional[BlastHit]:
    """
    Compare left and right seeds using BLAST to identify similarity regions.
    """
    logging.info('Comparing left and right seed sequences...')

    blast_output = temp_dir / 'seed_comparison.tab'

    try:
        run_blastn(
            query=left_seed,
            subject=right_seed,
            output=blast_output,
            word_size=4,
            perc_identity=50.0,  # Lower threshold for seed comparison
            verbose=True,
        )

        hits = parse_blast_output(blast_output)

        if hits:
            best_hit = max(hits, key=lambda h: h.length * h.identity)
            logging.info(f'Best seed similarity: {best_hit}')
            return best_hit
        else:
            logging.info('No significant similarity found between left and right seeds')
            return None

    except BlastError as e:
        logging.warning(f'Seed comparison failed: {e}')
        return None


def create_blast_database(genome_file: Path, db_dir: Path) -> Path:
    """Create BLAST nucleotide database from genome file using makeblastdb."""
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
        '-parse_seqids',  # Parse sequence identifiers for better handling
    ]

    logging.info(f'Creating BLAST database for {genome_file.name}')

    try:
        result = run_command(cmd, verbose=True)
        if result.returncode != 0:
            raise HMMBuildError(f'makeblastdb failed: {result.stderr}')

        logging.info(f'BLAST database created: {db_name}')
        return db_name

    except Exception as e:
        raise HMMBuildError(f'Failed to create BLAST database: {e}') from e


def blast_seed_against_genome(
    seed_file: Path, blast_db: Path, output_file: Path, identity_threshold: float = 60.0
) -> List[BlastHit]:
    """BLAST seed sequence against genome database using blastn."""

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
        '10000',  # Allow many hits for comprehensive search
        '-evalue',
        '1e-5',  # Reasonable e-value threshold
    ]

    logging.info(f'Running BLAST search with identity threshold {identity_threshold}%')

    try:
        result = run_command(cmd, verbose=True)
        if result.returncode != 0:
            raise BlastError(f'BLAST search failed: {result.stderr}')

        # Add header to BLAST output file
        add_blast_header(output_file)

        hits = parse_blast_output(output_file)
        logging.info(f'BLAST search found {len(hits)} hits')
        return hits

    except Exception as e:
        raise HMMBuildError(f'BLAST search failed: {e}') from e


def add_blast_header(blast_file: Path) -> None:
    """Add header to BLAST output file."""
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
    Prioritize by: 1) length, 2) identity, 3) query coverage
    """
    if not hits:
        return hits

    # Sort hits by genomic position for overlap detection
    sorted_hits = sorted(
        hits, key=lambda h: (h.subject_id, h.subject_start, h.subject_end)
    )

    filtered_hits = []

    for current_hit in sorted_hits:
        keep_current = True
        remove_indices = []

        for i, existing_hit in enumerate(filtered_hits):
            # Only check hits on same chromosome/scaffold
            if current_hit.subject_id != existing_hit.subject_id:
                continue

            # Calculate overlap
            overlap_start = max(current_hit.subject_start, existing_hit.subject_start)
            overlap_end = min(current_hit.subject_end, existing_hit.subject_end)
            overlap_length = max(0, overlap_end - overlap_start)

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

    logging.info(f'Resolved overlaps: {len(hits)} -> {len(filtered_hits)} hits')
    return filtered_hits


def filter_hits_by_thresholds(
    hits: List[BlastHit], min_coverage: float, min_identity: float
) -> List[BlastHit]:
    """Filter hits by coverage and identity thresholds."""
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
            logging.debug(f'PASSED: {hit}')
        else:
            logging.debug(f'FILTERED: {hit} - Failed: {", ".join(reasons)}')

    logging.info(
        f'Threshold filtering: {len(hits)} -> {len(filtered)} hits '
        f'(cov>={min_coverage:.3f}, id>={min_identity:.1f})'
    )
    return filtered


def chain_fragmented_hits(
    hits: List[BlastHit], max_gap: int = 500
) -> List[List[BlastHit]]:
    """
    Chain hits that may represent fragments of the same element.
    Returns list of hit chains (each chain is a list of hits).
    """
    if not hits:
        return []

    # Sort by genomic position
    sorted_hits = sorted(hits, key=lambda h: (h.subject_id, h.subject_start))

    chains = []
    current_chain = [sorted_hits[0]]

    for hit in sorted_hits[1:]:
        last_hit = current_chain[-1]

        # Check if hits are on same chromosome and within gap distance
        if (
            hit.subject_id == last_hit.subject_id
            and hit.subject_start - last_hit.subject_end <= max_gap
        ):
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
        logging.warning(
            f'Hit chaining: {single_hits} single hits, {chained_hits} chains detected. '
            f'WARNING: Chained hits suggest fragmented BLAST matches. '
            f'It is recommended to manually inspect and trim the alignment before building the HMM with hmmbuild '
            f'to ensure high-quality terminal repeat models.'
        )
    else:
        logging.info(f'Hit chaining: {single_hits} single hits, {chained_hits} chains')

    return chains


def extract_sequences_from_chains(
    chains: List[List[BlastHit]], genome_index, model_name: str
) -> List[SeqRecord]:
    """
    Extract sequences from hit chains, concatenating fragments where needed.
    """
    sequences = []

    for chain_idx, chain in enumerate(chains):
        try:
            if len(chain) == 1:
                # Single hit
                hit = chain[0]

                # Validate that subject_id exists in genome
                if hit.subject_id not in genome_index:
                    logging.warning(
                        f'Subject sequence {hit.subject_id} not found in genome. Skipping hit.'
                    )
                    continue

                chrom = genome_index[hit.subject_id]

                # Handle strand orientation
                start = (
                    min(hit.subject_start, hit.subject_end) - 1
                )  # Convert to 0-based
                end = max(hit.subject_start, hit.subject_end)

                seq_str = str(chrom[start:end]).upper()  # Convert to uppercase

                # Reverse complement if on negative strand
                if hit.strand == '-':
                    seq_str = str(Seq(seq_str).reverse_complement())

                # Include strand in sequence ID
                seq_id = f'{model_name}_{hit.subject_id}_{start + 1}_{end}_{hit.strand}'

            else:
                # Chained hits - concatenate fragments
                seq_parts = []
                locations = []
                valid_chain = True

                for hit in chain:
                    # Validate that subject_id exists in genome
                    if hit.subject_id not in genome_index:
                        logging.warning(
                            f'Subject sequence {hit.subject_id} not found in genome. Skipping chain.'
                        )
                        valid_chain = False
                        break

                    chrom = genome_index[hit.subject_id]
                    start = min(hit.subject_start, hit.subject_end) - 1
                    end = max(hit.subject_start, hit.subject_end)

                    fragment = str(chrom[start:end]).upper()  # Convert to uppercase
                    if hit.strand == '-':
                        fragment = str(Seq(fragment).reverse_complement())

                    seq_parts.append(fragment)
                    locations.append(f'{start + 1}_{end}_{hit.strand}')

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
            logging.error(f'Failed to extract sequence for chain {chain_idx}: {e}')
            continue
        except Exception as e:
            logging.error(
                f'Unexpected error extracting sequence for chain {chain_idx}: {e}'
            )
            continue

    logging.info(f'Extracted {len(sequences)} sequences from {len(chains)} hit chains')
    return sequences


def extract_flanked_sequences_from_chains(
    chains: List[List[BlastHit]], genome_index, model_name: str, flank_size: int
) -> List[SeqRecord]:
    """
    Extract sequences from hit chains with flanking sequence.
    For chained hits, only add flanks to the beginning of the first segment and end of the last segment.
    """
    sequences = []

    for chain_idx, chain in enumerate(chains):
        try:
            if len(chain) == 1:
                # Single hit with flanking
                hit = chain[0]

                # Validate that subject_id exists in genome
                if hit.subject_id not in genome_index:
                    logging.warning(
                        f'Subject sequence {hit.subject_id} not found in genome. Skipping hit.'
                    )
                    continue

                chrom = genome_index[hit.subject_id]
                chrom_len = len(chrom)

                # Handle strand orientation for coordinates
                hit_start = (
                    min(hit.subject_start, hit.subject_end) - 1
                )  # Convert to 0-based
                hit_end = max(hit.subject_start, hit.subject_end)

                # Add flanking sequence
                flanked_start = max(0, hit_start - flank_size)
                flanked_end = min(chrom_len, hit_end + flank_size)

                seq_str = str(chrom[flanked_start:flanked_end]).upper()

                # Reverse complement if on negative strand
                if hit.strand == '-':
                    seq_str = str(Seq(seq_str).reverse_complement())

                # Include flanking info in sequence ID
                seq_id = f'{model_name}_{hit.subject_id}_{flanked_start + 1}_{flanked_end}_{hit.strand}_flank{flank_size}'

            else:
                # Chained hits - add flanks only to first and last segments
                first_hit = chain[0]
                chain[-1]

                # Validate all subject_ids exist in genome
                valid_chain = True
                for hit in chain:
                    if hit.subject_id not in genome_index:
                        logging.warning(
                            f'Subject sequence {hit.subject_id} not found in genome. Skipping chain.'
                        )
                        valid_chain = False
                        break

                if not valid_chain:
                    continue

                # All hits should be on the same chromosome for a valid chain
                chrom = genome_index[first_hit.subject_id]
                chrom_len = len(chrom)

                seq_parts = []
                locations = []

                for i, hit in enumerate(chain):
                    hit_start = min(hit.subject_start, hit.subject_end) - 1
                    hit_end = max(hit.subject_start, hit.subject_end)

                    # Add flanking sequence only to first and last segments
                    if i == 0:  # First segment
                        start = max(0, hit_start - flank_size)
                        end = hit_end
                        locations.append(
                            f'{start + 1}_{end}_{hit.strand}_5flank{flank_size}'
                        )
                    elif i == len(chain) - 1:  # Last segment
                        start = hit_start
                        end = min(chrom_len, hit_end + flank_size)
                        locations.append(
                            f'{start + 1}_{end}_{hit.strand}_3flank{flank_size}'
                        )
                    else:  # Middle segments - no flanking
                        start = hit_start
                        end = hit_end
                        locations.append(f'{start + 1}_{end}_{hit.strand}')

                    fragment = str(chrom[start:end]).upper()
                    if hit.strand == '-':
                        fragment = str(Seq(fragment).reverse_complement())

                    seq_parts.append(fragment)

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
            logging.error(
                f'Failed to extract flanked sequence for chain {chain_idx}: {e}'
            )
            continue
        except Exception as e:
            logging.error(
                f'Unexpected error extracting flanked sequence for chain {chain_idx}: {e}'
            )
            continue

    logging.info(
        f'Extracted {len(sequences)} flanked sequences from {len(chains)} hit chains (flank size: {flank_size}bp)'
    )
    return sequences


def deduplicate_sequences(sequences: List[SeqRecord]) -> List[SeqRecord]:
    """Remove identical sequences, keeping the first occurrence. Prioritize seed sequences over extracted sequences."""
    seen_sequences = set()
    unique_sequences = []

    # Sort sequences to prioritize seeds (seed sequences typically don't have genomic coordinates in their IDs)
    def is_seed_sequence(seq_record):
        """Check if sequence is likely a seed sequence (no genomic coordinates)."""
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
                logging.info(
                    f'Seed sequence {seq.id} is identical to already retained sequence'
                )
            else:
                logging.debug(f'Skipping duplicate extracted sequence: {seq.id}')

    logging.info(
        f'Deduplication: {len(sequences)} -> {len(unique_sequences)} unique sequences'
    )
    return unique_sequences


def run_mafft_alignment(sequences: List[SeqRecord], output_file: Path) -> Path:
    """Run MAFFT multiple sequence alignment with default parameters."""
    if len(sequences) < 2:
        raise HMMBuildError('Need at least 2 sequences for alignment')

    # Write sequences to temporary file
    temp_fasta = output_file.parent / f'{output_file.stem}_input.fasta'

    with open(temp_fasta, 'w') as f:
        SeqIO.write(sequences, f, 'fasta')

    # Run MAFFT with default parameters
    cmd = [
        'mafft',
        str(temp_fasta),
    ]

    logging.info(f'Running MAFFT alignment on {len(sequences)} sequences')

    try:
        # Run MAFFT and capture output
        result = run_command(cmd, verbose=True)

        if result.returncode != 0:
            raise HMMBuildError(f'MAFFT failed: {result.stderr}')

        # Parse MAFFT output and convert to uppercase
        alignment_records = []
        for record in SeqIO.parse(io.StringIO(result.stdout), 'fasta'):
            uppercase_record = SeqRecord(
                Seq(str(record.seq).upper()),
                id=record.id,
                description=record.description,
            )
            alignment_records.append(uppercase_record)

        # Write the uppercase alignment to output file
        with open(output_file, 'w') as outfile:
            SeqIO.write(alignment_records, outfile, 'fasta')

        # Clean up temporary file
        temp_fasta.unlink()

        logging.info(f'Alignment written to {output_file}')
        return output_file

    except Exception as e:
        raise HMMBuildError(f'MAFFT alignment failed: {e}') from e


def calculate_pairwise_identity(alignment_file: Path) -> pd.DataFrame:
    """Calculate pairwise identity matrix from alignment."""
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
        logging.warning(f'Failed to calculate pairwise identity: {e}')
        return pd.DataFrame()


def build_hmm_from_alignment(
    alignment_file: Path, model_name: str, output_dir: Path
) -> Path:
    """Build HMM from multiple sequence alignment."""

    clean_model_name = cleanID(model_name)

    try:
        command, output_hmm = build_hmmbuild_command(
            model_name=clean_model_name,
            input_alignment=alignment_file,
            output_dir=output_dir,
        )

        logging.info(f'Building HMM for {model_name}')
        result = run_command(command, verbose=True)

        if result.returncode != 0:
            raise HMMBuildError(f'hmmbuild failed: {result.stderr}')

        logging.info(f'HMM written to {output_hmm}')
        return output_hmm

    except Exception as e:
        raise HMMBuildError(f'HMM building failed: {e}') from e


def save_all_blast_hits(all_hits: List[BlastHit], output_file: Path) -> None:
    """Save all BLAST hits in tab-delimited format 6."""
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

    logging.info(f'All BLAST hits saved to {output_file}')


def process_seed_sequences(
    seed_file: Path,
    model_name: str,
    genome_files: List[Path],
    temp_dir: Path,
    output_dir: Path,
    min_coverage: float,
    min_identity: float,
    max_gap: int = 500,
    save_blast_hits: bool = False,
    flank_size: Optional[int] = None,
) -> Tuple[Path, Path, Optional[Path], Optional[Path]]:
    """
    Process a seed file against genomes to build HMM.

    Returns:
        Tuple of (hmm_file, alignment_file, blast_hits_file, flanked_alignment_file)
    """
    logging.info(f'Processing {model_name} seed: {seed_file.name}')

    all_hits = []

    # BLAST against each genome
    for genome_file in genome_files:
        logging.info(f'Searching {model_name} against {genome_file.name}')

        # Create BLAST database
        db_dir = temp_dir / f'blast_dbs_{genome_file.stem}'
        db_dir.mkdir(exist_ok=True)
        blast_db = create_blast_database(genome_file, db_dir)

        # Run BLAST search
        blast_output = temp_dir / f'{model_name}_{genome_file.stem}_blast.tab'
        hits = blast_seed_against_genome(
            seed_file, blast_db, blast_output, min_identity
        )
        all_hits.extend(hits)

    logging.info(f'Total BLAST hits for {model_name}: {len(all_hits)}')

    if not all_hits:
        raise HMMBuildError(f'No BLAST hits found for {model_name}')

    # Save all BLAST hits in tabular format if requested
    blast_hits_file = None
    if save_blast_hits:
        blast_hits_file = output_dir / f'{cleanID(model_name)}_all_blast_hits.tab'
        save_all_blast_hits(all_hits, blast_hits_file)

    # Filter hits by thresholds
    filtered_hits = filter_hits_by_thresholds(all_hits, min_coverage, min_identity)

    if not filtered_hits:
        raise HMMBuildError(f'No hits passed thresholds for {model_name}')

    # Resolve overlapping hits
    resolved_hits = resolve_overlapping_hits(filtered_hits)

    # Chain fragmented hits (with warning message)
    hit_chains = chain_fragmented_hits(resolved_hits, max_gap)

    # Extract sequences with proper error handling
    try:
        # Use first genome for extraction (assuming all genomes are provided)
        genome_index = indexGenome(genome_files[0])
        sequences = extract_sequences_from_chains(hit_chains, genome_index, model_name)
    except Exception as e:
        raise HMMBuildError(f'Failed to extract sequences from genome: {e}') from e

    # Add original seed sequence(s) - convert to uppercase and add BEFORE deduplication
    seed_records = list(SeqIO.parse(seed_file, 'fasta'))
    logging.info(f'Adding {len(seed_records)} seed sequence(s) to extracted sequences')

    for seed_record in seed_records:
        # Convert seed sequence to uppercase
        uppercase_seed = SeqRecord(
            Seq(str(seed_record.seq).upper()),
            id=seed_record.id,
            description=seed_record.description,
        )
        sequences.append(uppercase_seed)

    logging.info(
        f'Total sequences before deduplication: {len(sequences)} (including {len(seed_records)} seed sequences)'
    )

    # Deduplicate - this will now prioritize seed sequences over extracted sequences
    unique_sequences = deduplicate_sequences(sequences)

    if len(unique_sequences) < 2:
        raise HMMBuildError(f'Not enough unique sequences for {model_name} alignment')

    # Create standard alignment
    alignment_file = output_dir / f'{cleanID(model_name)}_alignment.fasta'
    run_mafft_alignment(unique_sequences, alignment_file)

    # Create flanked alignment if requested
    flanked_alignment_file = None
    if flank_size is not None and flank_size > 0:
        logging.info(
            f'Creating flanked alignment with {flank_size}bp flanking sequence'
        )

        # Extract flanked sequences
        try:
            flanked_sequences = extract_flanked_sequences_from_chains(
                hit_chains, genome_index, model_name, flank_size
            )

            # Add seed sequences to flanked sequences as well
            for seed_record in seed_records:
                uppercase_seed = SeqRecord(
                    Seq(str(seed_record.seq).upper()),
                    id=f'{seed_record.id}_seed',
                    description=f'{seed_record.description} (seed sequence)',
                )
                flanked_sequences.append(uppercase_seed)

            # Deduplicate flanked sequences
            unique_flanked_sequences = deduplicate_sequences(flanked_sequences)

            if len(unique_flanked_sequences) >= 2:
                flanked_alignment_file = (
                    output_dir
                    / f'{cleanID(model_name)}_flanked_{flank_size}bp_alignment.fasta'
                )
                run_mafft_alignment(unique_flanked_sequences, flanked_alignment_file)
                logging.info(f'Flanked alignment written to {flanked_alignment_file}')
            else:
                logging.warning(
                    f'Not enough unique flanked sequences for {model_name} flanked alignment'
                )

        except Exception as e:
            logging.warning(f'Failed to create flanked alignment: {e}')

    # Calculate pairwise identity for standard alignment
    identity_matrix = calculate_pairwise_identity(alignment_file)
    if not identity_matrix.empty:
        identity_file = output_dir / f'{cleanID(model_name)}_pairwise_identity.csv'
        identity_matrix.to_csv(identity_file)
        logging.info(f'Pairwise identity matrix written to {identity_file}')

    # Build HMM from standard alignment (not flanked)
    hmm_file = build_hmm_from_alignment(alignment_file, model_name, output_dir)

    return hmm_file, alignment_file, blast_hits_file, flanked_alignment_file


def process_asymmetric_seeds(
    left_seed: Path,
    right_seed: Path,
    model_name: str,
    genome_files: List[Path],
    temp_dir: Path,
    output_dir: Path,
    min_coverage: float,
    min_identity: float,
    max_gap: int = 500,
    save_blast_hits: bool = False,
    flank_size: Optional[int] = None,
) -> Tuple[Path, Path, Path, Path]:
    """
    Process asymmetric left and right seeds together to avoid filtering conflicts.
    """
    logging.info(f'Processing asymmetric seeds for {model_name}')

    all_left_hits = []
    all_right_hits = []

    # BLAST both seeds against all genomes first
    for genome_file in genome_files:
        logging.info(f'Searching both seeds against {genome_file.name}')

        # Create BLAST database
        db_dir = temp_dir / f'blast_dbs_{genome_file.stem}'
        db_dir.mkdir(exist_ok=True)
        blast_db = create_blast_database(genome_file, db_dir)

        # BLAST left seed
        left_output = temp_dir / f'{model_name}_left_{genome_file.stem}_blast.tab'
        left_hits = blast_seed_against_genome(
            left_seed, blast_db, left_output, min_identity
        )
        all_left_hits.extend(left_hits)

        # BLAST right seed
        right_output = temp_dir / f'{model_name}_right_{genome_file.stem}_blast.tab'
        right_hits = blast_seed_against_genome(
            right_seed, blast_db, right_output, min_identity
        )
        all_right_hits.extend(right_hits)

    # Now apply intelligent filtering considering both seed results
    filtered_left = filter_hits_by_thresholds(all_left_hits, min_coverage, min_identity)
    filtered_right = filter_hits_by_thresholds(
        all_right_hits, min_coverage, min_identity
    )

    # Resolve overlaps within each seed type first
    resolved_left = resolve_overlapping_hits(filtered_left)
    resolved_right = resolve_overlapping_hits(filtered_right)

    # Then resolve conflicts between left and right hits
    resolved_left, resolved_right = resolve_asymmetric_conflicts(
        resolved_left, resolved_right, max_gap
    )

    # Process each seed type separately from here
    # ... continue with existing logic for each seed


def resolve_asymmetric_conflicts(
    left_hits: List[BlastHit], right_hits: List[BlastHit], max_gap: int
) -> Tuple[List[BlastHit], List[BlastHit]]:
    """
    Resolve conflicts between left and right seed hits.
    Prioritize longer, higher-identity hits.
    """
    # Group hits by chromosome
    left_by_chrom = {}
    right_by_chrom = {}

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
                        logging.debug(
                            f'Left hit {left_hit} removed due to better right hit {right_hit}'
                        )
                        break

            if keep_left:
                filtered_left.append(left_hit)

        # Similar logic for right hits
        for right_hit in chrom_right:
            keep_right = True
            for left_hit in chrom_left:
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
                        logging.debug(
                            f'Right hit {right_hit} removed due to better left hit {left_hit}'
                        )
                        break

            if keep_right:
                filtered_right.append(right_hit)

    return filtered_left, filtered_right


def hits_overlap(hit1: BlastHit, hit2: BlastHit, min_overlap: int = 50) -> bool:
    """Check if two hits overlap significantly."""
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


def main():
    """Main function for HMM building workflow."""
    parser = argparse.ArgumentParser(
        description='Build HMM models from seed sequences for TIRmite',
        prog='tirmite-build',
    )

    # Input arguments
    parser.add_argument(
        '--left-seed',
        type=Path,
        required=True,
        help='FASTA file containing left terminal seed sequence(s)',
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

    # Genome input (mutually exclusive)
    genome_group = parser.add_mutually_exclusive_group(required=True)
    genome_group.add_argument(
        '--genome', type=Path, help='Single genome FASTA file to search'
    )
    genome_group.add_argument(
        '--genome-list',
        type=Path,
        help='File containing list of genome paths (one per line)',
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
        type=float,
        default=0.7,
        help='Minimum query coverage threshold as fraction 0.0-1.0 (default: 0.7)',
    )
    parser.add_argument(
        '--min-identity',
        type=float,
        default=70.0,
        help='Minimum sequence identity threshold as percentage (default: 70.0)',
    )
    parser.add_argument(
        '--save-blast-hits',
        action='store_true',
        help='Save all BLAST hits in tabular format',
    )
    parser.add_argument(
        '--max-gap',
        type=int,
        default=500,
        help='Maximum gap size for chaining fragmented hits (default: 500)',
    )
    parser.add_argument(
        '--flank-size',
        type=int,
        help='Extract BLAST hits with flanking sequence of N bases and create additional alignment (optional)',
    )
    parser.add_argument(
        '--loglevel',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='Set logging level',
    )

    args = parser.parse_args()

    # Check thresholds are in range
    if not (0.0 <= args.min_coverage <= 1.0):
        parser.error('--min-coverage must be between 0.0 and 1.0')

    if not (0.0 <= args.min_identity <= 100.0):
        parser.error('--min-identity must be between 0.0 and 100.0')

    # Set up logging
    init_logging(loglevel=args.loglevel)

    try:
        # Check dependencies
        missing_tools = check_dependencies()
        if missing_tools:
            raise HMMBuildError(f'Missing required tools: {", ".join(missing_tools)}')

        # Validate input files
        if not args.left_seed.exists():
            raise FileNotFoundError(f'Left seed file not found: {args.left_seed}')

        if args.right_seed and not args.right_seed.exists():
            raise FileNotFoundError(f'Right seed file not found: {args.right_seed}')

        # Get genome files
        if args.genome:
            genome_files = [args.genome]
            if not args.genome.exists():
                raise FileNotFoundError(f'Genome file not found: {args.genome}')
        else:
            with open(args.genome_list, 'r') as f:
                genome_paths = [Path(line.strip()) for line in f if line.strip()]

            genome_files = []
            for path in genome_paths:
                if not path.exists():
                    logging.warning(f'Genome file not found, skipping: {path}')
                else:
                    genome_files.append(path)

            if not genome_files:
                raise FileNotFoundError('No valid genome files found')

        # Create output directory
        output_dir = create_output_directory(args.outdir)
        logging.info(f'Output directory: {output_dir}')

        # Set up temporary directory
        if args.tempdir:
            # Create the specified temp directory if it doesn't exist
            temp_base_dir = Path(args.tempdir)
            try:
                temp_base_dir.mkdir(parents=True, exist_ok=True)
                logging.info(f'Created/using temporary base directory: {temp_base_dir}')
            except Exception as e:
                raise HMMBuildError(
                    f'Failed to create temporary base directory {temp_base_dir}: {e}'
                ) from e
        else:
            temp_base_dir = None
            logging.info('Using system default temporary directory')

        with temporary_directory(
            prefix='tirmite_build_', base_dir=temp_base_dir, cleanup=not args.keep_temp
        ) as temp_dir:
            logging.info(f'Temporary directory: {temp_dir}')

            # Compare seeds if both provided
            if args.right_seed:
                seed_similarity = compare_seeds(
                    args.left_seed, args.right_seed, temp_dir
                )
                if seed_similarity:
                    logging.info(
                        f'Seed similarity detected: {seed_similarity.length}bp, '
                        f'{seed_similarity.identity:.1f}% identity'
                    )

            # Process left seed
            left_model_name = (
                f'{args.model_name}_left' if args.right_seed else args.model_name
            )
            left_hmm, left_aln, left_blast, left_flanked = process_seed_sequences(
                args.left_seed,
                left_model_name,
                genome_files,
                temp_dir,
                output_dir,
                args.min_coverage,
                args.min_identity,
                args.max_gap,
                args.save_blast_hits,
                args.flank_size,
            )

            logging.info(f'Left model completed: {left_hmm}')
            if left_flanked:
                logging.info(f'Left flanked alignment: {left_flanked}')

            # Process right seed if provided
            if args.right_seed:
                right_model_name = f'{args.model_name}_right'
                right_hmm, right_aln, right_blast, right_flanked = (
                    process_seed_sequences(
                        args.right_seed,
                        right_model_name,
                        genome_files,
                        temp_dir,
                        output_dir,
                        args.min_coverage,
                        args.min_identity,
                        args.max_gap,
                        args.save_blast_hits,
                        args.flank_size,
                    )
                )
                logging.info(f'Right model completed: {right_hmm}')
                if right_flanked:
                    logging.info(f'Right flanked alignment: {right_flanked}')

            logging.info('HMM building workflow completed successfully!')

            # Summary
            logging.info(f'Output files in {output_dir}:')
            for file in output_dir.glob(f'{cleanID(args.model_name)}*'):
                logging.info(f'  {file.name}')

    except Exception as e:
        logging.error(f'HMM building failed: {e}')
        return 1

    return 0


if __name__ == '__main__':
    exit(main())
