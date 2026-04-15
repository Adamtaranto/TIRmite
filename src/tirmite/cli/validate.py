#!/usr/bin/env python3
"""
TIRmite-validate: Validate reconstructed target sites.

This module validates reconstructed target sites by searching for
matching intact target sites in a genome database:
1. Run blastn of reconstructed target sites against a validation database
2. Filter for hits spanning the target site junction
3. Extract hit regions and re-align with mafft
4. Validate TSD/DR length from alignment gaps
"""

import argparse
import csv
import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Dict, List, Optional, Tuple, cast

from Bio import AlignIO, SeqIO  # type: ignore[import-not-found]
from Bio.Seq import Seq  # type: ignore[import-not-found]
from Bio.SeqRecord import SeqRecord  # type: ignore[import-not-found]

from tirmite._version import __version__  # type: ignore[import-not-found]
from tirmite.utils.logs import init_logging


def add_validate_parser(subparsers: Any) -> argparse.ArgumentParser:
    """
    Add validate subcommand parser.

    Parameters
    ----------
    subparsers : argparse._SubParsersAction
        Subparser object to add validate command to.

    Returns
    -------
    argparse.ArgumentParser
        The configured validate subcommand parser.
    """
    parser = cast(
        argparse.ArgumentParser,
        subparsers.add_parser(
            'validate',
            help='Validate reconstructed target sites',
            description=(
                'Validate reconstructed target sites by searching for matching '
                'intact (empty) sites in a genome database.'
            ),
        ),
    )
    _configure_validate_parser(parser)
    return parser


def create_validate_parser() -> argparse.ArgumentParser:
    """
    Create standalone argument parser for validate command.

    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser for validate workflow options.
    """
    parser = argparse.ArgumentParser(
        prog='tirmite-validate',
        description='Validate reconstructed target sites against a genome database',
    )
    _configure_validate_parser(parser)
    return parser


def _configure_validate_parser(parser: argparse.ArgumentParser) -> None:
    """
    Configure parser with validate command arguments.

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
        '--target-sites',
        type=str,
        required=True,
        dest='target_sites',
        help=(
            'Path to reconstructed target site FASTA file from tirmite pair '
            'with target site reconstruction enabled.'
        ),
    )

    parser.add_argument(
        '--blastdb',
        type=str,
        required=True,
        help='Path to BLAST database for validation searches.',
    )

    parser.add_argument(
        '--blast-results',
        type=str,
        default=None,
        dest='blast_results',
        help=(
            'Path to pre-computed BLAST results in format 6. If provided, '
            'skip running blastn and use these results instead. '
            'Expected columns: qseqid sseqid pident length mismatch gapopen '
            'qstart qend sstart send evalue bitscore qlen slen sstrand.'
        ),
    )

    parser.add_argument(
        '--min-coverage',
        type=float,
        default=0.95,
        dest='min_coverage',
        help=(
            'Minimum fraction of the query that must be covered by a hit '
            'to pass filtering. Default: 0.95'
        ),
    )

    parser.add_argument(
        '--evalue',
        type=float,
        default=1e-5,
        help='E-value threshold for blastn searches. Default: 1e-5',
    )

    parser.add_argument(
        '--tsd-length',
        type=int,
        default=0,
        dest='tsd_length',
        help=(
            'Default TSD/DR length to validate. '
            'Used when --tsd-length-map is not provided.'
        ),
    )

    parser.add_argument(
        '--tsd-length-map',
        type=str,
        default=None,
        dest='tsd_length_map',
        help=(
            'Path to tab-delimited file mapping model pairs to TSD lengths. '
            'Format: left_model<TAB>right_model<TAB>tsd_length.'
        ),
    )

    parser.add_argument(
        '--tsd-in-model',
        action='store_true',
        default=False,
        dest='tsd_in_model',
        help='TSD/DR feature is part of the termini model.',
    )

    parser.add_argument(
        '--outdir',
        type=str,
        default=None,
        help='Directory for output files. Default: current directory.',
    )

    parser.add_argument(
        '--prefix',
        type=str,
        default=None,
        help='Prefix for output file names.',
    )

    parser.add_argument(
        '--logfile',
        action='store_true',
        default=False,
        help='Write log messages to file in output directory.',
    )


def run_blastn(
    query_fasta: str,
    blastdb: str,
    output_file: str,
    evalue: float = 1e-5,
) -> None:
    """
    Run blastn search of query sequences against a BLAST database.

    Parameters
    ----------
    query_fasta : str
        Path to query FASTA file.
    blastdb : str
        Path to BLAST database.
    output_file : str
        Path to output file for BLAST results.
    evalue : float, default 1e-5
        E-value threshold.

    Raises
    ------
    RuntimeError
        If blastn is not found or fails.
    """
    if not shutil.which('blastn'):
        raise RuntimeError(
            'blastn not found in PATH. Please install NCBI BLAST+.'
        )

    cmd = [
        'blastn',
        '-query', query_fasta,
        '-db', blastdb,
        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand',
        '-evalue', str(evalue),
        '-out', output_file,
        '-dust', 'no',
    ]

    logging.info(f'Running blastn: {" ".join(cmd)}')
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logging.error(f'blastn failed: {result.stderr}')
        raise RuntimeError(f'blastn failed with return code {result.returncode}')


def parse_blast_results(
    blast_file: str,
) -> List[Dict[str, Any]]:
    """
    Parse BLAST format 6 results with extended columns.

    Parameters
    ----------
    blast_file : str
        Path to BLAST format 6 output file.

    Returns
    -------
    list of dict
        List of hit dictionaries with keys: qseqid, sseqid, pident, length,
        mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore,
        qlen, slen, sstrand.
    """
    columns = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
        'qlen', 'slen', 'sstrand',
    ]

    hits: List[Dict[str, Any]] = []

    if not os.path.exists(blast_file):
        logging.warning(f'BLAST results file not found: {blast_file}')
        return hits

    with open(blast_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < len(columns):
                continue
            hit: Dict[str, Any] = {}
            for i, col in enumerate(columns):
                val = row[i]
                if col in ('qstart', 'qend', 'sstart', 'send', 'length',
                           'mismatch', 'gapopen', 'qlen', 'slen'):
                    hit[col] = int(val)
                elif col in ('pident', 'evalue', 'bitscore'):
                    hit[col] = float(val)
                else:
                    hit[col] = val
            hits.append(hit)

    logging.info(f'Parsed {len(hits)} hits from {blast_file}')
    return hits


def filter_junction_spanning(
    hits: List[Dict[str, Any]],
    min_coverage: float = 0.95,
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Filter BLAST hits for junction-spanning hits with sufficient query coverage.

    Parameters
    ----------
    hits : list of dict
        Parsed BLAST hits.
    min_coverage : float, default 0.95
        Minimum fraction of query covered by the hit.

    Returns
    -------
    dict
        Dictionary mapping query ID to list of passing hits.
    """
    filtered: Dict[str, List[Dict[str, Any]]] = {}

    for hit in hits:
        qlen = hit['qlen']
        qstart = hit['qstart']
        qend = hit['qend']
        midpoint = qlen / 2.0

        # Check if hit spans the midpoint (junction)
        hit_qstart = min(qstart, qend)
        hit_qend = max(qstart, qend)
        if hit_qstart > midpoint or hit_qend < midpoint:
            continue

        # Check coverage
        alignment_len = hit_qend - hit_qstart + 1
        coverage = alignment_len / qlen
        if coverage < min_coverage:
            continue

        qid = hit['qseqid']
        if qid not in filtered:
            filtered[qid] = []
        filtered[qid].append(hit)

    total_passing = sum(len(v) for v in filtered.values())
    logging.info(
        f'Filtered to {total_passing} junction-spanning hits '
        f'across {len(filtered)} queries'
    )
    return filtered


def extract_hit_sequence(
    blastdb: str,
    subject_id: str,
    start: int,
    end: int,
    strand: str,
) -> Optional[str]:
    """
    Extract a sequence region from a BLAST database using blastdbcmd.

    Parameters
    ----------
    blastdb : str
        Path to BLAST database.
    subject_id : str
        Subject sequence identifier.
    start : int
        1-based start coordinate.
    end : int
        1-based end coordinate.
    strand : str
        Strand: 'plus' or 'minus'.

    Returns
    -------
    str or None
        Extracted sequence, or None if extraction failed.
    """
    if not shutil.which('blastdbcmd'):
        logging.error('blastdbcmd not found in PATH')
        return None

    cmd = [
        'blastdbcmd',
        '-db', blastdb,
        '-entry', subject_id,
        '-range', f'{start}-{end}',
        '-strand', strand,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.warning(
            f'blastdbcmd failed for {subject_id}:{start}-{end}: {result.stderr}'
        )
        return None

    # Parse FASTA output
    lines = result.stdout.strip().split('\n')
    seq = ''.join(line for line in lines if not line.startswith('>'))
    return seq if seq else None


def run_mafft_alignment(
    sequences: List[SeqRecord],
    tmpdir: str,
) -> Optional[List[SeqRecord]]:
    """
    Align sequences using MAFFT.

    Parameters
    ----------
    sequences : list of SeqRecord
        Sequences to align. First should be the query.
    tmpdir : str
        Temporary directory for intermediate files.

    Returns
    -------
    list of SeqRecord or None
        Aligned sequences, or None if alignment failed.
    """
    if not shutil.which('mafft'):
        logging.error('mafft not found in PATH. Please install MAFFT.')
        return None

    input_file = os.path.join(tmpdir, 'mafft_input.fasta')
    output_file = os.path.join(tmpdir, 'mafft_output.fasta')

    with open(input_file, 'w') as handle:
        SeqIO.write(sequences, handle, 'fasta')

    cmd = [
        'mafft',
        '--auto',
        '--quiet',
        input_file,
    ]

    with open(output_file, 'w') as out_handle:
        result = subprocess.run(cmd, stdout=out_handle, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        logging.warning(f'MAFFT alignment failed: {result.stderr}')
        return None

    try:
        alignment = list(AlignIO.read(output_file, 'fasta'))
        return alignment
    except Exception as e:
        logging.warning(f'Failed to parse MAFFT output: {e}')
        return None


def check_tsd_gaps(
    query_aligned: str,
    target_aligned: str,
    tsd_length: int,
    query_len: int,
) -> Tuple[int, str]:
    """
    Check for gaps around the alignment midpoint indicating TSD length errors.

    Parameters
    ----------
    query_aligned : str
        Aligned query sequence (with gaps).
    target_aligned : str
        Aligned target sequence (with gaps).
    tsd_length : int
        User-specified TSD length.
    query_len : int
        Original unaligned query length.

    Returns
    -------
    predicted_error : int
        Predicted error in TSD length. Positive means user specified too long,
        negative means too short, 0 means consistent.
    message : str
        Human-readable message about the validation result.
    """
    # Find the midpoint of the query in the aligned coordinates
    query_pos = 0
    midpoint_aligned = 0
    target_midpoint = query_len // 2

    for i, c in enumerate(query_aligned):
        if c != '-':
            query_pos += 1
        if query_pos >= target_midpoint:
            midpoint_aligned = i
            break

    # Check a window around the midpoint for gaps
    window = max(tsd_length + 5, 10)
    start = max(0, midpoint_aligned - window)
    end = min(len(query_aligned), midpoint_aligned + window + 1)

    query_gaps = query_aligned[start:end].count('-')
    target_gaps = target_aligned[start:end].count('-')

    if query_gaps > 0 and target_gaps == 0:
        # Gaps in query at midpoint → user TSD was too long
        return query_gaps, (
            f'Query has {query_gaps} gap(s) near midpoint: '
            f'TSD length may be {query_gaps}bp too long'
        )
    elif target_gaps > 0 and query_gaps == 0:
        # Gaps in target at midpoint → user TSD was too short
        return -target_gaps, (
            f'Target has {target_gaps} gap(s) near midpoint: '
            f'TSD length may be {target_gaps}bp too short'
        )
    elif query_gaps == 0 and target_gaps == 0:
        return 0, 'TSD length appears consistent with validation data'
    else:
        return 0, (
            f'Both query ({query_gaps}) and target ({target_gaps}) have gaps '
            f'near midpoint; validation inconclusive'
        )


def main(args: Optional[argparse.Namespace] = None) -> int:
    """
    Main entry point for tirmite-validate.

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
        parser = create_validate_parser()
        args = parser.parse_args()

    assert args is not None

    try:
        # Set up output directory
        if args.outdir:
            outdir = Path(args.outdir)
        else:
            outdir = Path.cwd()
        outdir.mkdir(parents=True, exist_ok=True)

        # Set up logging
        logfile_path = None
        if args.logfile:
            prefix_str = f'{args.prefix}_' if args.prefix else ''
            logfile_path = outdir / f'{prefix_str}tirmite_validate.log'

        init_logging(loglevel=args.loglevel, logfile=logfile_path)

        logging.info(f'TIRmite-validate version {__version__}')
        logging.info(f'Output directory: {outdir}')

        # Validate inputs
        if not Path(args.target_sites).exists():
            logging.error(f'Target sites file not found: {args.target_sites}')
            sys.exit(1)

        # Load query sequences
        logging.info(f'Loading target sites from {args.target_sites}')
        queries = list(SeqIO.parse(args.target_sites, 'fasta'))
        if not queries:
            logging.error('No sequences found in target sites file')
            sys.exit(1)
        logging.info(f'Loaded {len(queries)} target site queries')

        # Load TSD length map if provided
        tsd_length_map: Dict[str, int] = {}
        if args.tsd_length_map:
            from tirmite.tirmitetools import load_tsd_length_map
            tsd_length_map = load_tsd_length_map(args.tsd_length_map)

        # Run or load BLAST results
        with tempfile.TemporaryDirectory() as tmpdir:
            if args.blast_results:
                logging.info(
                    f'Loading pre-computed BLAST results from {args.blast_results}'
                )
                blast_file = args.blast_results
            else:
                blast_file = os.path.join(tmpdir, 'blast_results.txt')
                run_blastn(
                    query_fasta=args.target_sites,
                    blastdb=args.blastdb,
                    output_file=blast_file,
                    evalue=args.evalue,
                )

            # Parse and filter results
            all_hits = parse_blast_results(blast_file)
            filtered_hits = filter_junction_spanning(
                all_hits, min_coverage=args.min_coverage
            )

            # Process each query
            prefix_str = f'{args.prefix}_' if args.prefix else ''
            summary_rows: List[Dict[str, Any]] = []
            all_alignments: Dict[str, List[SeqRecord]] = {}

            for query in queries:
                qid = query.id
                query_len = len(query.seq)
                logging.info(f'Processing query: {qid} (len={query_len})')

                # Parse model pair info from description
                desc_parts = query.description.split()
                left_model = ''
                right_model = ''
                for part in desc_parts:
                    if part.startswith('left_model='):
                        left_model = part.split('=', 1)[1]
                    elif part.startswith('right_model='):
                        right_model = part.split('=', 1)[1]

                # Resolve TSD length for this query
                query_tsd_length = args.tsd_length
                if tsd_length_map and left_model and right_model:
                    key = f'{left_model}\t{right_model}'
                    if key in tsd_length_map:
                        query_tsd_length = tsd_length_map[key]

                query_hits = filtered_hits.get(qid, [])
                logging.info(f'  {len(query_hits)} hits passing filters')

                if not query_hits:
                    summary_rows.append({
                        'query_id': qid,
                        'left_model': left_model,
                        'right_model': right_model,
                        'query_length': query_len,
                        'tsd_length': query_tsd_length,
                        'num_valid_hits': 0,
                        'predicted_tsd_error': 'N/A',
                        'validation_message': 'No valid empty site hits found',
                    })
                    continue

                # Extract hit sequences and re-align
                alignment_seqs = [query]
                tsd_errors: List[int] = []
                tsd_messages: List[str] = []

                for hit in query_hits:
                    sstart = hit['sstart']
                    send = hit['send']

                    # Handle strand
                    sstrand = hit.get('sstrand', 'plus')
                    if sstart > send:
                        strand = 'minus'
                        sstart, send = send, sstart
                    else:
                        strand = 'plus' if sstrand == 'plus' else 'minus'

                    seq_str = extract_hit_sequence(
                        args.blastdb, hit['sseqid'], sstart, send, strand
                    )
                    if seq_str is None:
                        continue

                    hit_id = f"{hit['sseqid']}_{sstart}_{send}_{strand}"
                    hit_record = SeqRecord(
                        Seq(seq_str),
                        id=hit_id,
                        name=hit_id,
                        description=f"evalue={hit['evalue']} pident={hit['pident']}",
                    )
                    alignment_seqs.append(hit_record)

                if len(alignment_seqs) > 1:
                    # Run MAFFT alignment
                    aligned = run_mafft_alignment(alignment_seqs, tmpdir)

                    if aligned:
                        all_alignments[qid] = aligned
                        query_aligned = str(aligned[0].seq)

                        for aln_rec in aligned[1:]:
                            target_aligned = str(aln_rec.seq)
                            error, msg = check_tsd_gaps(
                                query_aligned, target_aligned,
                                query_tsd_length, query_len
                            )
                            tsd_errors.append(error)
                            tsd_messages.append(msg)

                # Determine consensus TSD error
                if tsd_errors:
                    avg_error = sum(tsd_errors) / len(tsd_errors)
                    if abs(avg_error) < 0.5:
                        consensus_msg = 'TSD length appears consistent'
                    elif avg_error > 0:
                        consensus_msg = (
                            f'TSD may be ~{abs(round(avg_error))}bp too long'
                        )
                    else:
                        consensus_msg = (
                            f'TSD may be ~{abs(round(avg_error))}bp too short'
                        )

                    if abs(avg_error) >= 1:
                        logging.warning(
                            f'TSD length validation warning for {qid}: {consensus_msg}'
                        )
                else:
                    avg_error = 0
                    consensus_msg = 'No alignments available for validation'

                summary_rows.append({
                    'query_id': qid,
                    'left_model': left_model,
                    'right_model': right_model,
                    'query_length': query_len,
                    'tsd_length': query_tsd_length,
                    'num_valid_hits': len(query_hits),
                    'predicted_tsd_error': round(avg_error, 1),
                    'validation_message': consensus_msg,
                })

            # Write summary report
            summary_file = outdir / f'{prefix_str}validation_summary.tsv'
            with open(summary_file, 'w', newline='') as f:
                if summary_rows:
                    writer = csv.DictWriter(
                        f, fieldnames=summary_rows[0].keys(), delimiter='\t'
                    )
                    writer.writeheader()
                    writer.writerows(summary_rows)
            logging.info(f'Validation summary written to {summary_file}')

            # Write aligned FASTA files
            for qid, aligned_records in all_alignments.items():
                safe_qid = qid.replace('/', '_').replace('\\', '_')
                aln_file = outdir / f'{prefix_str}{safe_qid}_alignment.fasta'
                with open(aln_file, 'w') as handle:
                    SeqIO.write(aligned_records, handle, 'fasta')
                logging.info(f'Alignment written to {aln_file}')

        logging.info('TIRmite-validate analysis completed successfully')

    except KeyboardInterrupt:
        logging.info('Analysis interrupted by user')
        sys.exit(130)
    except Exception as e:
        logging.error(f'Unexpected error: {e}')
        logging.exception('Full traceback:')
        sys.exit(1)

    return 0


if __name__ == '__main__':
    main()
