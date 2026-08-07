"""Benchmarks for sequence extraction and hit parsing.

Extraction runs once per emitted feature, so a slow path here shows up
directly in the wall time of `tirmite pair` on a large result set. Parsing
runs once per input file but over every line in it.
"""

import pytest

from tirmite.core.extraction import fetch_padded_hit
from tirmite.core.parsers import import_nhmmer
from tirmite.utils.extract import fetch_region_padded, fetch_sequence

# Region widths spanning a short terminus hit up to a whole element.
REGION_WIDTHS = [100, 1000, 10_000]


@pytest.mark.parametrize('width', REGION_WIDTHS)
def test_fetch_sequence(benchmark, genome_fasta, width):
    """Plain in-bounds fetch; the baseline every other path builds on."""
    start = 500_000
    result = benchmark(fetch_sequence, genome_fasta, 'chr1', start, start + width - 1)

    assert result is not None
    assert len(result) == width


@pytest.mark.parametrize('width', REGION_WIDTHS)
def test_fetch_region_padded(benchmark, genome_fasta, width):
    """Padded fetch well inside the contig, so no padding is actually added."""
    start = 500_000
    result = benchmark(
        fetch_region_padded, genome_fasta, 'chr1', start, start + width - 1
    )

    assert result is not None


def test_fetch_region_padded_at_contig_start(benchmark, genome_fasta):
    """Padded fetch overrunning the contig start, exercising the pad path."""
    result = benchmark(fetch_region_padded, genome_fasta, 'chr1', -500, 500)

    assert result is not None


@pytest.mark.parametrize('strand', ['+', '-'])
def test_fetch_padded_hit(benchmark, genome_fasta, strand):
    """Hit fetch with lowercase-marked flanks, on both strands.

    The minus strand additionally reverse-complements, so the two are worth
    tracking separately.
    """
    result = benchmark(
        fetch_padded_hit,
        genome_fasta,
        'chr1',
        500_000,
        502_000,
        strand,
        500,
    )

    assert result is not None


def test_fetch_padded_hit_no_padding(benchmark, genome_fasta):
    """Hit fetch without flanks, isolating the padding cost by comparison."""
    result = benchmark(
        fetch_padded_hit, genome_fasta, 'chr1', 500_000, 502_000, '+', None
    )

    assert result is not None


def test_import_nhmmer(benchmark, nhmmer_results):
    """Parse a 20,000-hit nhmmer --tblout file.

    Includes the numeric-coordinate sort, which is the dominant cost at this
    size.
    """
    result = benchmark(import_nhmmer, nhmmer_results)

    assert len(result) == 20_000
