#!/usr/bin/env python3
"""
Differential tests proving the two sequence-extraction backends agree.

TIRmite retrieves hit sequences from either a pyfaidx-indexed FASTA genome
(``--genome``) or a BLAST database (``--blastdb``). These are user-selectable
alternatives, so the same coordinates in the same genome must yield the same
sequence from either backend.

Every test here builds a FASTA file *and* a BLAST database from that same file,
queries both, and asserts byte-identical results. Also covers the specific
``blastdbcmd`` behaviours that previously caused silent divergence:

* a range starting past the contig end returns the whole sequence, exit 0
* soft-masking is discarded (output is always uppercase)
* accessions differ from pyfaidx keys for headers such as ``>lcl|contig_1``
"""

import os
import shutil
import subprocess
import tempfile

import pandas as pd
from pyfaidx import Fasta
import pytest

from tirmite.cli.hmm_build import (
    BlastHit,
    extract_flanked_sequences_from_chains,
    extract_sequences_from_chains,
)
import tirmite.tirmitetools as tirmite
from tirmite.tirmitetools import extractTIRs, fetch_padded_hit
from tirmite.utils.extract import (
    BlastDBSource,
    FastaSource,
    annotate,
    check_ids,
    clamp_region,
    fetch_region_padded,
    fetch_sequence,
)

# Every test in this module shells out to blastdbcmd, which costs roughly a
# second per invocation. Marked slow so it can be deselected during quick runs;
# it should always run in CI, since backend equivalence is what it guards.
pytestmark = pytest.mark.slow

requires_blast = pytest.mark.skipif(
    shutil.which('makeblastdb') is None or shutil.which('blastdbcmd') is None,
    reason='BLAST+ (makeblastdb/blastdbcmd) not available',
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# 200bp with no repeats, so any misalignment shows up as a content difference
# rather than an accidental match.
GENOME_SEQ = (
    'CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGT'
    'CCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACCGCGGTGTTAAGTGTCGAGCT'
    'ACATCACTTCTCATGTAGCCAGAAGGCTGCAACTCATCGACTCTATGTAGTGACCGCGTC'
)
SHORT_SEQ = 'ACGTACGTAA'


_COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')


def _revcomp(seq):
    """Reverse complement, preserving case."""
    return seq.translate(_COMPLEMENT)[::-1]


def _write_fasta(path, records):
    """Write {id: (description, sequence)} to a wrapped FASTA file."""
    with open(path, 'w') as fh:
        for seq_id, (desc, seq) in records.items():
            header = f'>{seq_id} {desc}' if desc else f'>{seq_id}'
            fh.write(header + '\n')
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + '\n')


def _make_pair(tmpdir, records, name='test'):
    """
    Build a FASTA and a BLAST database from identical content.

    Returns
    -------
    tuple of (FastaSource, BlastDBSource)
    """
    fasta_path = os.path.join(tmpdir, f'{name}.fa')
    db_path = os.path.join(tmpdir, f'{name}_db')

    _write_fasta(fasta_path, records)

    subprocess.run(
        [
            'makeblastdb',
            '-in',
            fasta_path,
            '-dbtype',
            'nucl',
            '-out',
            db_path,
            '-parse_seqids',
        ],
        check=True,
        capture_output=True,
    )

    return FastaSource(Fasta(fasta_path)), BlastDBSource(db_path)


@pytest.fixture(scope='module')
def sources(tmp_path_factory):
    """
    A FastaSource and BlastDBSource over the same two-contig genome.

    Module-scoped: building a BLAST database costs about a second, and every
    test here reads the same immutable genome.
    """
    tmpdir = str(tmp_path_factory.mktemp('equiv'))
    return _make_pair(
        tmpdir,
        {
            'chr1': ('test chromosome one', GENOME_SEQ),
            'short1': ('a short contig', SHORT_SEQ),
        },
    )


# ---------------------------------------------------------------------------
# fetch_sequence: coordinate matrix
# ---------------------------------------------------------------------------


@requires_blast
@pytest.mark.parametrize(
    'start,end',
    [
        (1, 10),  # at contig start
        (50, 60),  # interior
        (1, 1),  # single base, first
        (200, 200),  # single base, last
        (191, 200),  # flush against contig end
        (1, 200),  # whole contig
        (100, 101),  # two-base window
    ],
)
def test_fetch_sequence_matches_across_backends(sources, start, end):
    """Both backends return identical sequence for in-bounds coordinates."""
    fasta_src, blast_src = sources

    from_fasta = fetch_sequence(fasta_src, 'chr1', start, end)
    from_blast = fetch_sequence(blast_src, 'chr1', start, end)

    assert from_fasta == from_blast
    # And the result is actually the region asked for.
    assert from_fasta == GENOME_SEQ[start - 1 : end]
    assert len(from_fasta) == end - start + 1


@requires_blast
def test_contig_lengths_match(sources):
    """Both backends agree on contig length."""
    fasta_src, blast_src = sources

    assert fasta_src.contig_length('chr1') == blast_src.contig_length('chr1')
    assert fasta_src.contig_length('chr1') == len(GENOME_SEQ)
    assert fasta_src.contig_length('short1') == blast_src.contig_length('short1')


# ---------------------------------------------------------------------------
# Out-of-bounds handling
# ---------------------------------------------------------------------------


@requires_blast
@pytest.mark.parametrize(
    'start,end',
    [
        (190, 250),  # overruns the end
        (0, 10),  # start below 1
        (-20, 5),  # start well below 1
        (150, 400),  # overruns by a lot
    ],
)
def test_out_of_range_is_clamped_identically(sources, start, end):
    """Partially out-of-range windows clamp to the same sequence on both."""
    fasta_src, blast_src = sources

    from_fasta = fetch_sequence(fasta_src, 'chr1', start, end)
    from_blast = fetch_sequence(blast_src, 'chr1', start, end)

    assert from_fasta == from_blast
    expected_start = max(1, start)
    expected_end = min(len(GENOME_SEQ), end)
    assert from_fasta == GENOME_SEQ[expected_start - 1 : expected_end]


@requires_blast
@pytest.mark.parametrize('start,end', [(250, 300), (201, 210), (500, 600)])
def test_range_entirely_past_contig_end_returns_none(sources, start, end):
    """
    A window wholly beyond the contig yields None from both backends.

    blastdbcmd answers such a request with the ENTIRE sequence and a zero exit
    status, which would silently substitute a whole contig for the intended
    region. Guard against that specifically.
    """
    fasta_src, blast_src = sources

    from_fasta = fetch_sequence(fasta_src, 'chr1', start, end)
    from_blast = fetch_sequence(blast_src, 'chr1', start, end)

    assert from_fasta is None
    assert from_blast is None
    # The failure mode being guarded against: never the whole contig.
    assert from_blast != GENOME_SEQ


@requires_blast
def test_raw_blastdbcmd_returns_whole_contig_past_end(sources):
    """
    Document the underlying blastdbcmd behaviour this module defends against.

    If this test ever fails, blastdbcmd's out-of-range semantics have changed
    and the length check in fetch_sequence can be revisited.
    """
    _, blast_src = sources

    raw = blast_src.fetch_raw('chr1', 250, 300)

    assert raw == GENOME_SEQ, 'blastdbcmd no longer returns the whole contig'


@requires_blast
def test_inverted_coordinates_return_none(sources):
    """start > end yields None from both backends rather than a stray base."""
    fasta_src, blast_src = sources

    assert fetch_sequence(fasta_src, 'chr1', 60, 50) is None
    assert fetch_sequence(blast_src, 'chr1', 60, 50) is None


@requires_blast
def test_clamp_region_agrees(sources):
    """Clamped coordinates match, so records get identical coordinate labels."""
    fasta_src, blast_src = sources

    assert clamp_region(fasta_src, 'chr1', -5, 250) == clamp_region(
        blast_src, 'chr1', -5, 250
    )
    assert clamp_region(fasta_src, 'chr1', -5, 250) == (1, 200)
    assert clamp_region(fasta_src, 'chr1', 300, 400) is None
    assert clamp_region(blast_src, 'chr1', 300, 400) is None


# ---------------------------------------------------------------------------
# Soft-masking
# ---------------------------------------------------------------------------


@requires_blast
def test_soft_masked_genome_yields_identical_uppercase():
    """
    Soft-masking must not make the backends disagree.

    pyfaidx preserves lowercase from the FASTA; blastdbcmd discards it. Both are
    normalised to uppercase.
    """
    masked = 'ACGTacgtACGTacgtACGT'

    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_src, blast_src = _make_pair(
            tmpdir, {'sm1': ('soft masked', masked)}, name='masked'
        )

        from_fasta = fetch_sequence(fasta_src, 'sm1', 1, 20)
        from_blast = fetch_sequence(blast_src, 'sm1', 1, 20)

        assert from_fasta == from_blast
        assert from_fasta == masked.upper()


# ---------------------------------------------------------------------------
# Sequence descriptions
# ---------------------------------------------------------------------------


@requires_blast
def test_sequence_descriptions_match_across_backends(sources):
    """
    Both backends report the same FASTA header description.

    The FASTA side reads it from the index, the BLAST side from the database
    title. If only one supplied descriptions, output headers would differ.
    """
    fasta_src, blast_src = sources

    assert fasta_src.sequence_description('chr1') == 'test chromosome one'
    assert blast_src.sequence_description('chr1') == 'test chromosome one'
    assert fasta_src.sequence_description('short1') == blast_src.sequence_description(
        'short1'
    )


@requires_blast
def test_annotate_matches_and_omits_empty_descriptions():
    """
    A sequence with no description gets no trailing space from either backend.

    Appending an empty description left a trailing space on the FASTA path only,
    which made otherwise identical headers differ.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_src, blast_src = _make_pair(
            tmpdir,
            {'plain': ('', SHORT_SEQ), 'annotated': ('has a description', SHORT_SEQ)},
            name='desc',
        )

        for src in (fasta_src, blast_src):
            assert annotate(src, 'plain', '[coords]') == '[coords]'
            assert annotate(src, 'annotated', '[coords]') == (
                '[coords] has a description'
            )

        assert annotate(fasta_src, 'plain', '[c]') == annotate(
            blast_src, 'plain', '[c]'
        )
        assert annotate(fasta_src, 'annotated', '[c]') == annotate(
            blast_src, 'annotated', '[c]'
        )


# ---------------------------------------------------------------------------
# Missing / unusual sequence IDs
# ---------------------------------------------------------------------------


@requires_blast
def test_missing_contig_returns_none_from_both(sources):
    """An unknown ID fails the same way on both backends, without raising."""
    fasta_src, blast_src = sources

    assert fetch_sequence(fasta_src, 'no_such_contig', 1, 10) is None
    assert fetch_sequence(blast_src, 'no_such_contig', 1, 10) is None


@requires_blast
def test_check_ids_reports_missing_from_both(sources):
    """check_ids agrees on which IDs are unresolvable."""
    fasta_src, blast_src = sources
    query = ['chr1', 'short1', 'missing_a', 'missing_b']

    assert check_ids(fasta_src, query) == ['missing_a', 'missing_b']
    assert check_ids(blast_src, query) == ['missing_a', 'missing_b']


@requires_blast
def test_pipe_prefixed_header_id_divergence_is_detectable():
    """
    Headers such as '>lcl|contig_1' index differently in each backend.

    pyfaidx keys on the whole first token ('lcl|contig_1') while the BLAST
    database stores the accession ('contig_1'). check_ids surfaces this rather
    than letting extraction silently return nothing.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_src, blast_src = _make_pair(
            tmpdir, {'lcl|contig_1': ('desc one', SHORT_SEQ)}, name='piped'
        )

        # The pyfaidx key resolves on both (blastdbcmd accepts the full seqid).
        assert fetch_sequence(fasta_src, 'lcl|contig_1', 1, 5) == SHORT_SEQ[:5]
        assert fetch_sequence(blast_src, 'lcl|contig_1', 1, 5) == SHORT_SEQ[:5]

        # The bare accession resolves only in the BLAST database. check_ids
        # reports the difference instead of failing silently at extraction.
        assert check_ids(fasta_src, ['contig_1']) == ['contig_1']
        assert check_ids(blast_src, ['contig_1']) == []


# ---------------------------------------------------------------------------
# fetch_padded_hit: strand and case marking
# ---------------------------------------------------------------------------


@requires_blast
@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('padlen', [None, 0, 5, 20])
def test_padded_hit_matches_across_backends(sources, strand, padlen):
    """Padded and unpadded hits agree on both backends, on both strands."""
    fasta_src, blast_src = sources

    from_fasta = fetch_padded_hit(fasta_src, 'chr1', 50, 79, strand, padlen)
    from_blast = fetch_padded_hit(blast_src, 'chr1', 50, 79, strand, padlen)

    assert from_fasta == from_blast
    assert from_fasta is not None


@requires_blast
def test_pad_case_marking_is_on_the_correct_flank_for_minus_strand(sources):
    """
    Lowercase pad markers must follow their own flank through the revcomp.

    Marking a sequence blastdbcmd had already reverse-complemented put the 5'
    and 3' markers on the wrong ends. Case marking now happens on the plus
    strand, before reverse complement.
    """
    fasta_src, blast_src = sources
    pad = 6
    start, end = 50, 79

    for src in (fasta_src, blast_src):
        seq = fetch_padded_hit(src, 'chr1', start, end, '-', pad)

        assert len(seq) == (end - start + 1) + 2 * pad
        # Both flanks lowercase, hit uppercase, regardless of strand.
        assert seq[:pad].islower()
        assert seq[-pad:].islower()
        assert seq[pad:-pad].isupper()

    # And the minus-strand result is exactly the revcomp of the plus-strand one.
    plus = fetch_padded_hit(fasta_src, 'chr1', start, end, '+', pad)
    minus = fetch_padded_hit(fasta_src, 'chr1', start, end, '-', pad)
    assert minus == _revcomp(plus)


@requires_blast
def test_asymmetric_pad_case_marking_survives_revcomp(sources):
    """
    Case marking must be correct when the two flanks differ in length.

    With equal flanks the marker pattern is symmetric, so marking before or
    after the reverse complement looks the same. When one flank is clipped at a
    contig boundary the two orders give different answers, and only marking on
    the plus strand first is correct.
    """
    fasta_src, blast_src = sources

    # short1 is 10bp. A hit at 2-5 with 5bp padding requests [-3, 10]: the left
    # pad runs 4bp off the contig start while the right pad fits exactly. The
    # two flanks therefore differ in composition even though both are 5bp.
    start, end, pad = 2, 5, 5

    plus = fetch_padded_hit(fasta_src, 'short1', start, end, '+', pad)
    # 4bp of pad, then 1 real base, then the hit, then the trailing flank.
    assert plus == (
        ('N' * 4 + SHORT_SEQ[0:1]).lower()
        + SHORT_SEQ[1:5].upper()
        + SHORT_SEQ[5:10].lower()
    )

    for src in (fasta_src, blast_src):
        minus = fetch_padded_hit(src, 'short1', start, end, '-', pad)

        # The marking follows the sequence through the reverse complement: the
        # trailing flank leads and the padded flank trails.
        assert minus == _revcomp(plus)
        assert minus[:5].islower()
        assert minus[5:9].isupper()
        assert minus[9:].islower()
        # The hit itself is real sequence; only the flanks may be padded.
        assert 'n' not in minus[5:9].lower()


@requires_blast
def test_padlen_beyond_contig_is_padded_to_full_width(sources):
    """A --padlen window overrunning both contig ends is padded, not clipped."""
    fasta_src, blast_src = sources

    start, end, pad = 5, 8, 20
    expected_width = (end - start + 1) + 2 * pad

    from_fasta = fetch_padded_hit(fasta_src, 'short1', start, end, '+', pad)
    from_blast = fetch_padded_hit(blast_src, 'short1', start, end, '+', pad)

    assert from_fasta == from_blast
    # Fixed width regardless of how little contig was available.
    assert len(from_fasta) == expected_width
    # The real bases are still present and the hit is still uppercase.
    assert SHORT_SEQ[4:8].upper() in from_fasta
    assert from_fasta.upper().count('N') == expected_width - len(SHORT_SEQ)


@requires_blast
def test_padlen_truncates_when_padding_disabled(sources):
    """With pad=False the --padlen window is clipped, as it was previously."""
    fasta_src, blast_src = sources

    from_fasta = fetch_padded_hit(fasta_src, 'short1', 5, 8, '+', 20, pad=False)
    from_blast = fetch_padded_hit(blast_src, 'short1', 5, 8, '+', 20, pad=False)

    assert from_fasta == from_blast
    assert len(from_fasta) == len(SHORT_SEQ)
    assert 'N' not in from_fasta.upper()


# ---------------------------------------------------------------------------
# extractTIRs: full record equivalence
# ---------------------------------------------------------------------------


def _hit_table(rows):
    """Build a minimal hitTable DataFrame."""
    return pd.DataFrame(
        [
            {
                'model': r['model'],
                'target': r['target'],
                'hitStart': str(r['start']),
                'hitEnd': str(r['end']),
                'strand': r['strand'],
                'evalue': '1e-10',
                'score': '100',
                'bias': 'NA',
                'hmmStart': '1',
                'hmmEnd': '100',
            }
            for r in rows
        ]
    ).reset_index(drop=True)


@requires_blast
@pytest.mark.parametrize('padlen', [None, 10])
def test_extractTIRs_records_match_across_backends(sources, padlen):
    """Full SeqRecords, not just sequences, match between backends."""
    fasta_src, blast_src = sources

    hitTable = _hit_table(
        [
            {'model': 'TIR', 'target': 'chr1', 'start': 20, 'end': 49, 'strand': '+'},
            {'model': 'TIR', 'target': 'chr1', 'start': 120, 'end': 149, 'strand': '-'},
        ]
    )

    from_fasta, count_fasta = extractTIRs(
        model='TIR', hitTable=hitTable, source=fasta_src, padlen=padlen
    )
    from_blast, count_blast = extractTIRs(
        model='TIR', hitTable=hitTable, source=blast_src, padlen=padlen
    )

    assert count_fasta == count_blast == 2
    assert len(from_fasta) == len(from_blast) == 2

    for rec_a, rec_b in zip(from_fasta, from_blast):
        assert str(rec_a.seq) == str(rec_b.seq)
        assert rec_a.id == rec_b.id
        assert rec_a.name == rec_b.name
        assert rec_a.description == rec_b.description


@requires_blast
def test_extractTIRs_minus_strand_id_suffix_on_both(sources):
    """Minus-strand records get the '_rc' suffix from either backend."""
    fasta_src, blast_src = sources

    hitTable = _hit_table(
        [{'model': 'TIR', 'target': 'chr1', 'start': 20, 'end': 49, 'strand': '-'}]
    )

    for src in (fasta_src, blast_src):
        records, _ = extractTIRs(model='TIR', hitTable=hitTable, source=src)
        assert records[0].id.endswith('_rc')
        assert records[0].name == records[0].id


@requires_blast
def test_extractTIRs_skips_unresolvable_target(sources):
    """A hit on a missing contig is skipped, not fatal, on both backends."""
    fasta_src, blast_src = sources

    hitTable = _hit_table(
        [
            {'model': 'TIR', 'target': 'chr1', 'start': 20, 'end': 49, 'strand': '+'},
            {'model': 'TIR', 'target': 'ghost', 'start': 20, 'end': 49, 'strand': '+'},
        ]
    )

    for src in (fasta_src, blast_src):
        records, count = extractTIRs(model='TIR', hitTable=hitTable, source=src)
        assert count == 2  # both counted as eligible
        assert len(records) == 1  # only the resolvable one extracted


# ---------------------------------------------------------------------------
# writeTIRs / fetchElements end-to-end equivalence
# ---------------------------------------------------------------------------


def _paired_inputs(hitTable):
    """Build hitIndex and paired dicts for the first two hits."""
    _, hitIndex = tirmite.table2dict(hitTable)
    idx = list(hitTable.index)
    model = hitTable.iloc[0]['model']
    paired = {model: [{idx[0], idx[1]}]}
    hitIndex[hitTable.iloc[0]['model']][idx[0]]['partner'] = idx[1]
    hitIndex[hitTable.iloc[1]['model']][idx[1]]['partner'] = idx[0]
    return hitIndex, paired


@requires_blast
def test_fetchElements_matches_across_backends():
    """Element sequences are identical, and always on the plus strand."""
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_src, blast_src = _make_pair(
            tmpdir, {'chr1': ('chromosome one', GENOME_SEQ)}, name='ele'
        )

        # Left hit on the minus strand: the case that previously diverged.
        hitTable = _hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'start': 20,
                    'end': 49,
                    'strand': '-',
                },
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'start': 120,
                    'end': 149,
                    'strand': '-',
                },
            ]
        )
        hitIndex, paired = _paired_inputs(hitTable)

        from_fasta = tirmite.fetchElements(
            paired=paired, hitIndex=hitIndex, genome=fasta_src.genome
        )
        from_blast = tirmite.fetchElements(
            paired=paired, hitIndex=hitIndex, blastdb=blast_src.path
        )

        seq_fasta = str(from_fasta['TIR'][0].seq.seq)
        seq_blast = str(from_blast['TIR'][0].seq.seq)

        assert seq_fasta == seq_blast
        # Plus strand, matching the genomic coordinates reported in the GFF.
        assert seq_fasta == GENOME_SEQ[19:149]


# ---------------------------------------------------------------------------
# hmm_build extraction: coordinate contract and backend equivalence
# ---------------------------------------------------------------------------


def _blast_hit(subject_id, sstart, send, subject_len):
    """Build a BlastHit with the given subject coordinates."""
    return BlastHit(
        query_id='seed',
        subject_id=subject_id,
        query_start=1,
        query_end=abs(send - sstart) + 1,
        subject_start=sstart,
        subject_end=send,
        length=abs(send - sstart) + 1,
        identity=100.0,
        query_len=abs(send - sstart) + 1,
        subject_len=subject_len,
    )


@requires_blast
def test_extract_sequences_from_chains_uses_1based_coordinates(sources):
    """
    Hit coordinates are 1-based inclusive and must not be shifted.

    An off-by-one here silently drops the first base of every extracted hit,
    which is invisible in the output because the length still looks plausible.
    """
    fasta_src, blast_src = sources
    hit = _blast_hit('chr1', 20, 49, len(GENOME_SEQ))

    for src in (fasta_src, blast_src):
        records = extract_sequences_from_chains([[hit]], src, 'TIR')

        assert len(records) == 1
        # Exactly the 1-based inclusive region 20..49, no shift in either
        # direction.
        assert str(records[0].seq) == GENOME_SEQ[19:49]
        assert len(records[0].seq) == 30


@requires_blast
@pytest.mark.parametrize('sstart,send', [(20, 49), (49, 20)])
def test_extract_sequences_from_chains_matches_across_backends(sources, sstart, send):
    """Chain extraction agrees between backends on both strands."""
    fasta_src, blast_src = sources
    hit = _blast_hit('chr1', sstart, send, len(GENOME_SEQ))

    from_fasta = extract_sequences_from_chains([[hit]], fasta_src, 'TIR')
    from_blast = extract_sequences_from_chains([[hit]], blast_src, 'TIR')

    assert len(from_fasta) == len(from_blast) == 1
    assert str(from_fasta[0].seq) == str(from_blast[0].seq)
    assert from_fasta[0].id == from_blast[0].id

    expected = GENOME_SEQ[19:49]
    if sstart > send:  # minus strand
        expected = _revcomp(expected)
    assert str(from_fasta[0].seq) == expected


@requires_blast
def test_extract_flanked_sequences_from_chains_matches_across_backends(sources):
    """Flanked chain extraction agrees, including the coordinates in the ID."""
    fasta_src, blast_src = sources
    hit = _blast_hit('chr1', 20, 49, len(GENOME_SEQ))

    from_fasta = extract_flanked_sequences_from_chains([[hit]], fasta_src, 'TIR', 10)
    from_blast = extract_flanked_sequences_from_chains([[hit]], blast_src, 'TIR', 10)

    assert str(from_fasta[0].seq) == str(from_blast[0].seq)
    assert from_fasta[0].id == from_blast[0].id
    assert str(from_fasta[0].seq) == GENOME_SEQ[9:59]


@requires_blast
def test_flanked_chain_ids_report_clamped_coordinates(sources):
    """
    IDs report post-clamping coordinates, so both backends label regions alike.

    Previously the BLAST path reported the unclamped request while the FASTA
    path reported the clipped region, giving different IDs for the same
    sequence.
    """
    fasta_src, blast_src = sources
    # A hit at the very end of short1, with flanks that overrun both boundaries.
    hit = _blast_hit('short1', 8, 10, len(SHORT_SEQ))

    from_fasta = extract_flanked_sequences_from_chains([[hit]], fasta_src, 'TIR', 20)
    from_blast = extract_flanked_sequences_from_chains([[hit]], blast_src, 'TIR', 20)

    assert from_fasta[0].id == from_blast[0].id
    # Clamped to the contig, not the requested -12..30.
    assert from_fasta[0].id == 'TIR_short1_1_10_+_flank20'
    assert str(from_fasta[0].seq) == SHORT_SEQ.upper()


@requires_blast
@pytest.mark.parametrize('padlen', [None, 8])
def test_writeTIRs_output_files_match_across_backends(padlen):
    """The written FASTA files are byte-identical between backends."""
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_src, blast_src = _make_pair(
            tmpdir, {'chr1': ('chromosome one', GENOME_SEQ)}, name='wt'
        )

        hitTable = _hit_table(
            [
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'start': 20,
                    'end': 49,
                    'strand': '+',
                },
                {
                    'model': 'TIR',
                    'target': 'chr1',
                    'start': 120,
                    'end': 149,
                    'strand': '-',
                },
            ]
        )

        out_fasta = os.path.join(tmpdir, 'out_genome')
        out_blast = os.path.join(tmpdir, 'out_blastdb')

        tirmite.writeTIRs(
            outDir=out_fasta,
            hitTable=hitTable,
            genome=fasta_src.genome,
            padlen=padlen,
        )
        tirmite.writeTIRs(
            outDir=out_blast,
            hitTable=hitTable,
            blastdb=blast_src.path,
            padlen=padlen,
        )

        files_a = sorted(os.listdir(out_fasta))
        files_b = sorted(os.listdir(out_blast))
        assert files_a == files_b
        assert files_a, 'expected at least one output file'

        for name in files_a:
            with open(os.path.join(out_fasta, name)) as fh:
                content_a = fh.read()
            with open(os.path.join(out_blast, name)) as fh:
                content_b = fh.read()
            assert content_a == content_b, f'{name} differs between backends'


# ---------------------------------------------------------------------------
# Out-of-bounds policy across every extraction method
# ---------------------------------------------------------------------------
#
# The policy, which these tests pin:
#   * partial overlap  -> pad (fetch_region_padded) or clamp (fetch_sequence)
#   * no overlap       -> None, from every method
#   * inverted/empty   -> None, from every method
# Both backends must agree in every case.

# (label, start, end) relative to a 200bp contig.
OUT_OF_BOUNDS_CASES = [
    ('past_start', -10, 5),
    ('past_end', 195, 220),
    ('spans_whole_contig', -50, 400),
    ('wholly_before', -50, -10),
    ('wholly_after', 300, 400),
    ('inverted', 60, 50),
    ('first_base', 1, 1),
    ('last_base', 200, 200),
]


def _case(label):
    """Look up an out-of-bounds case by label."""
    return next((s, e) for lbl, s, e in OUT_OF_BOUNDS_CASES if lbl == label)


@requires_blast
@pytest.mark.parametrize('label,start,end', OUT_OF_BOUNDS_CASES)
def test_fetch_sequence_out_of_bounds_matches(sources, label, start, end):
    """fetch_sequence clamps or skips identically on both backends."""
    fasta_src, blast_src = sources

    from_fasta = fetch_sequence(fasta_src, 'chr1', start, end)
    from_blast = fetch_sequence(blast_src, 'chr1', start, end)

    assert from_fasta == from_blast

    overlaps = start <= end and start <= len(GENOME_SEQ) and end >= 1
    if overlaps:
        expected = GENOME_SEQ[max(1, start) - 1 : min(len(GENOME_SEQ), end)]
        assert from_fasta == expected
    else:
        assert from_fasta is None


@requires_blast
@pytest.mark.parametrize('label,start,end', OUT_OF_BOUNDS_CASES)
def test_fetch_region_padded_out_of_bounds_matches(sources, label, start, end):
    """fetch_region_padded pads partial overlap and skips the rest, on both."""
    fasta_src, blast_src = sources

    from_fasta = fetch_region_padded(fasta_src, 'chr1', start, end)
    from_blast = fetch_region_padded(blast_src, 'chr1', start, end)

    assert from_fasta == from_blast

    overlaps = start <= end and start <= len(GENOME_SEQ) and end >= 1
    if not overlaps:
        # No overlap means no information: never an all-N record.
        assert from_fasta is None
        return

    # Always exactly the requested width, with the real bases in place.
    assert len(from_fasta.seq) == end - start + 1
    assert from_fasta.left_pad == max(0, 1 - start)
    assert from_fasta.right_pad == max(0, end - len(GENOME_SEQ))
    real = GENOME_SEQ[from_fasta.start - 1 : from_fasta.end]
    assert from_fasta.seq == (
        'N' * from_fasta.left_pad + real + 'N' * from_fasta.right_pad
    )


@requires_blast
def test_padded_and_clamped_agree_on_real_bases(sources):
    """Padding adds to, and never alters, what clamping returns."""
    fasta_src, _ = sources

    for label, start, end in OUT_OF_BOUNDS_CASES:
        padded = fetch_region_padded(fasta_src, 'chr1', start, end)
        clamped = fetch_sequence(fasta_src, 'chr1', start, end)

        if padded is None:
            assert clamped is None, label
            continue

        stripped = padded.seq.strip('N')
        # GENOME_SEQ has no Ns, so stripping recovers exactly the real region.
        assert stripped == clamped, label


@requires_blast
@pytest.mark.parametrize('label', ['wholly_before', 'wholly_after', 'inverted'])
def test_no_overlap_skips_in_every_method(sources, label):
    """Every extraction entry point declines a region it cannot observe."""
    fasta_src, blast_src = sources
    start, end = _case(label)

    for src in (fasta_src, blast_src):
        assert fetch_sequence(src, 'chr1', start, end) is None
        assert fetch_region_padded(src, 'chr1', start, end) is None
        assert clamp_region(src, 'chr1', start, end) is None
        # A hit window that does not overlap yields nothing, padded or not.
        assert fetch_padded_hit(src, 'chr1', start, end, '+', None) is None
        assert fetch_padded_hit(src, 'chr1', start, end, '+', 5, pad=True) is None


@requires_blast
def test_hmm_build_chain_extractors_clamp_not_pad(sources):
    """
    hmm_build extraction must not introduce N at contig edges.

    These sequences become MAFFT input and then HMM training data, so padding
    would turn a contig boundary into an N column in the model. Ragged ends are
    correct here.
    """
    fasta_src, blast_src = sources
    # Hit flush against the end of the 10bp contig, flanks overrunning both ends.
    hit = _blast_hit('short1', 3, 9, len(SHORT_SEQ))

    for src in (fasta_src, blast_src):
        plain = extract_sequences_from_chains([[hit]], src, 'TIR')
        flanked = extract_flanked_sequences_from_chains([[hit]], src, 'TIR', 20)

        assert 'N' not in str(plain[0].seq)
        assert 'N' not in str(flanked[0].seq)
        # Clamped to what the contig actually holds.
        assert str(flanked[0].seq) == SHORT_SEQ.upper()
