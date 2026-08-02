"""Synthetic fixture data for the TIRmite benchmarks.

The repository deliberately ships no test data, so every benchmark input is
generated here. Generation is deterministic (a fixed ``random.Random`` seed)
so that a run measures code changes rather than input changes -- CodSpeed
compares runs against each other, and a varying workload would make every
comparison noise.

Fixtures are session-scoped: building a 10k-hit table is itself measurable
work, and it is setup, not the thing under test.
"""

import random
from typing import Dict, List

import pandas as pd
import pytest

# Fixed seed: the workload must be identical across runs for CodSpeed's
# run-to-run comparison to mean anything.
SEED = 20260802

# Nucleotides, weighted evenly. Real genomes are not uniform, but the
# extraction paths do not branch on base composition.
BASES = 'ACGT'


def make_genome_sequence(length: int, seed: int = SEED) -> str:
    """
    Build a deterministic pseudo-random nucleotide sequence.

    Parameters
    ----------
    length : int
        Number of bases to generate.
    seed : int, default SEED
        Seed for the local random generator.

    Returns
    -------
    str
        Uppercase nucleotide sequence of the requested length.
    """
    rng = random.Random(seed)
    return ''.join(rng.choice(BASES) for _ in range(length))


def make_hit_table(
    n_hits: int,
    n_models: int = 4,
    n_targets: int = 2,
    contig_length: int = 5_000_000,
    seed: int = SEED,
) -> pd.DataFrame:
    """
    Build a synthetic hit table shaped like importer output.

    Parameters
    ----------
    n_hits : int
        Total number of hit rows.
    n_models : int, default 4
        Number of distinct model names, cycled across hits.
    n_targets : int, default 2
        Number of distinct target (contig) names.
    contig_length : int, default 5_000_000
        Upper bound for generated hit coordinates.
    seed : int, default SEED
        Seed for the local random generator.

    Returns
    -------
    pandas.DataFrame
        Hit table with the columns the importers produce, all values as
        strings, matching the package's dtype contract.

    Notes
    -----
    Hits are clustered into loci rather than scattered uniformly. The filters
    under benchmark are O(n^2) *within a target group* and only do real work
    where hits overlap, so a uniform scatter would understate their cost and
    make the benchmark insensitive to the very changes it exists to catch.
    """
    rng = random.Random(seed)
    rows: List[Dict[str, str]] = []

    # Roughly eight hits per locus, so groups genuinely overlap.
    n_loci = max(1, n_hits // 8)
    locus_starts = [rng.randrange(1, contig_length) for _ in range(n_loci)]

    for i in range(n_hits):
        locus = locus_starts[i % n_loci]
        # Jitter within the locus so hits overlap partially, not identically.
        start = locus + rng.randrange(0, 400)
        length = rng.randrange(80, 220)
        model = f'Model{i % n_models}'
        target = f'chr{i % n_targets}'
        rows.append(
            {
                'model': model,
                'target': target,
                'hitStart': str(start),
                'hitEnd': str(start + length),
                'strand': '+' if (i % 2 == 0) else '-',
                'evalue': '1e-20',
                'score': str(rng.randrange(50, 400)),
                'bias': '0.0',
                'hmmStart': '1',
                'hmmEnd': str(length),
            }
        )
    return pd.DataFrame(rows)


def make_pairing_hit_rows(
    n_elements: int,
    spacing: int = 5000,
    element_length: int = 2000,
    left_model: str = 'LeftTIR',
    right_model: str = 'RightTIR',
) -> List[Dict[str, object]]:
    """
    Build hit rows describing ``n_elements`` well-formed elements.

    Parameters
    ----------
    n_elements : int
        Number of elements; each contributes one left and one right hit.
    spacing : int, default 5000
        Distance between successive elements.
    element_length : int, default 2000
        Distance from the left terminus start to the right terminus start.
    left_model : str, default 'LeftTIR'
        Model name for the left terminus.
    right_model : str, default 'RightTIR'
        Model name for the right terminus.

    Returns
    -------
    list of dict
        Rows in the shape ``_make_hit_table`` in the test suite uses.

    Notes
    -----
    Every element is pairable, which is the expensive case: the engine has to
    evaluate reciprocity for each candidate rather than rejecting it early.
    """
    rows: List[Dict[str, object]] = []
    for i in range(n_elements):
        base = 1000 + (i * spacing)
        rows.append(
            {
                'model': left_model,
                'target': 'chr1',
                'hit_start': base,
                'hit_end': base + 100,
                'strand': '+',
            }
        )
        rows.append(
            {
                'model': right_model,
                'target': 'chr1',
                'hit_start': base + element_length,
                'hit_end': base + element_length + 100,
                'strand': '-',
            }
        )
    return rows


def rows_to_hit_table(rows: List[Dict[str, object]]) -> pd.DataFrame:
    """
    Convert pairing rows into the string-valued hit table ``table2dict`` expects.

    Parameters
    ----------
    rows : list of dict
        Rows with keys ``model``, ``target``, ``hit_start``, ``hit_end``,
        ``strand``.

    Returns
    -------
    pandas.DataFrame
        Hit table with the full column set, all values as strings.
    """
    return pd.DataFrame(
        [
            {
                'model': r['model'],
                'target': r['target'],
                'hitStart': str(r['hit_start']),
                'hitEnd': str(r['hit_end']),
                'strand': r['strand'],
                'evalue': '1e-10',
                'score': '100',
                'bias': 'NA',
                'hmmStart': '1',
                'hmmEnd': '100',
            }
            for r in rows
        ]
    )


def make_nhmmer_file(path, n_hits: int, seed: int = SEED) -> str:
    """
    Write a synthetic nhmmer ``--tblout`` file.

    Parameters
    ----------
    path : pathlib.Path
        Destination file path.
    n_hits : int
        Number of hit lines to write.
    seed : int, default SEED
        Seed for the local random generator.

    Returns
    -------
    str
        The path written, as a string.

    Notes
    -----
    Column order matches real nhmmer output, which ``import_nhmmer`` indexes
    positionally: target, accession, query, accession, hmmfrom, hmm to,
    alifrom, ali to, envfrom, env to, sq len, strand, E-value, score, bias,
    description.
    """
    rng = random.Random(seed)
    lines = [
        '# target name  accession  query name  accession  hmmfrom  hmm to  '
        'alifrom  ali to  envfrom  env to  sq len  strand  E-value  score  '
        'bias  description\n'
    ]
    for i in range(n_hits):
        start = rng.randrange(1, 5_000_000)
        end = start + rng.randrange(80, 220)
        strand = '+' if i % 2 == 0 else '-'
        lines.append(
            f'chr{i % 4} - Model{i % 8} - 1 60 {start} {end} {start} {end} '
            f'5000000 {strand} 1.2e-10 45.2 0.0 -\n'
        )
    path.write_text(''.join(lines))
    return str(path)


@pytest.fixture(scope='session')
def small_hit_table() -> pd.DataFrame:
    """Return a 500-row hit table."""
    return make_hit_table(500)


@pytest.fixture(scope='session')
def large_hit_table() -> pd.DataFrame:
    """Return a 5000-row hit table."""
    return make_hit_table(5000)


@pytest.fixture(scope='session')
def pairing_map() -> Dict[str, str]:
    """Return a pairing map covering the synthetic model names."""
    return {'Model0': 'Model1', 'Model2': 'Model3'}


@pytest.fixture(scope='session')
def genome_fasta(tmp_path_factory):
    """Write a 2 Mb single-contig FASTA and return an indexed FastaSource.

    ``FastaSource`` wraps an already-indexed ``pyfaidx.Fasta``, not a path, so
    the file is indexed here first. The fixture is session-scoped because
    building and indexing 2 Mb is setup, not the operation being measured.
    """
    from pyfaidx import Fasta

    from tirmite.utils.extract import FastaSource

    directory = tmp_path_factory.mktemp('genome')
    fasta = directory / 'genome.fa'
    sequence = make_genome_sequence(2_000_000)
    with open(fasta, 'w') as handle:
        handle.write('>chr1\n')
        for i in range(0, len(sequence), 60):
            handle.write(sequence[i : i + 60] + '\n')
    return FastaSource(Fasta(str(fasta)))


@pytest.fixture(scope='session')
def nhmmer_results(tmp_path_factory) -> str:
    """Write a 20,000-line nhmmer --tblout file and return its path.

    20k rather than 50k: parsing scales linearly, and 50k made this single
    benchmark dominate the suite's runtime without telling us anything the
    smaller size does not.
    """
    directory = tmp_path_factory.mktemp('nhmmer')
    return make_nhmmer_file(directory / 'hits.out', 20_000)
