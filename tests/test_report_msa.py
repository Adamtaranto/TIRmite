"""Tests for tirmite.report.msa — stacked terminus alignment panels."""

import pytest

from tirmite.report import msa
from tirmite.report.model import HitRecord, ModelInfo, PairingGroup


class FakeSource:
    """Sequence source over an in-memory dict of contigs."""

    def __init__(self, contigs):
        self.contigs = contigs


def fake_fetch_region_padded(source, seqid, start, end, pad_char='N'):
    """Stand-in for utils.extract.fetch_region_padded over FakeSource."""
    from tirmite.utils.extract import PaddedRegion

    seq = source.contigs.get(seqid)
    if seq is None:
        return None
    length = len(seq)
    clamped_start = max(1, start)
    clamped_end = min(length, end)
    if clamped_start > clamped_end:
        return None
    left_pad = clamped_start - start
    right_pad = end - clamped_end
    body = seq[clamped_start - 1 : clamped_end]
    return PaddedRegion(
        seq=pad_char * left_pad + body + pad_char * right_pad,
        left_pad=left_pad,
        right_pad=right_pad,
        start=clamped_start,
        end=clamped_end,
    )


@pytest.fixture(autouse=True)
def patched_fetch(monkeypatch):
    monkeypatch.setattr(
        'tirmite.utils.extract.fetch_region_padded', fake_fetch_region_padded
    )


def make_hit(
    uid,
    contig,
    start,
    end,
    strand='+',
    hmm=(1, 20),
    model='TIR',
    trunc=(0, 0),
    evalue=1e-20,
):
    return HitRecord(
        uid=uid,
        model=model,
        contig=contig,
        start=start,
        end=end,
        strand=strand,
        evalue=evalue,
        hmm_start=hmm[0],
        hmm_end=hmm[1],
        model_length=20,
        trunc_left=trunc[0],
        trunc_right=trunc[1],
        group_ids=['TIR__TIR'],
        role='left',
    )


GROUPS = [
    PairingGroup(
        group_id='TIR__TIR',
        label='TIR',
        left_model='TIR',
        right_model='TIR',
        orientation='F,R',
        colour='#2a78d6',
    )
]
MODELS = [ModelInfo(name='TIR', length=20)]


def build(hits, source, monkeypatch, mode='anchor', **kwargs):
    return msa.build_msa_panels(
        hits, GROUPS, MODELS, source=source, tmpdir='/unused', mode=mode, **kwargs
    )


def build_padded(hits, source, monkeypatch, mode='anchor', **kwargs):
    """Build panels with model padding on, as --padlen requests."""
    return build(hits, source, monkeypatch, mode=mode, pad_model=True, **kwargs)


class TestPadRuns:
    def test_model_pads_become_m_runs(self):
        runs = msa.pad_runs('ACGTACGT', [(0, 3)])
        assert runs == [[0, 3, 'm']]

    def test_gaps_become_g_runs(self):
        runs = msa.pad_runs('--ACGT---', [])
        assert runs == [[0, 2, 'g'], [6, 3, 'g']]

    def test_trailing_gap_run_is_closed(self):
        assert msa.pad_runs('ACGT--', []) == [[4, 2, 'g']]

    def test_both_kinds_are_reported_and_sorted(self):
        # Model-pad ranges are half-open, so (2, 4) is two columns.
        runs = msa.pad_runs('--ACGT', [(2, 4)])
        assert runs == [[0, 2, 'g'], [2, 2, 'm']]

    def test_no_padding_gives_no_runs(self):
        assert msa.pad_runs('ACGT', []) == []


class TestModelPadVersusGap:
    def test_unmatched_model_within_the_contig_is_grey_not_a_gap(self, monkeypatch):
        # The hit matched model positions 6-20 only, and the genome extends
        # past the hit at both ends: the shortfall is a model pad.
        source = FakeSource({'chr1': 'ACGT' * 100})
        hit = make_hit(0, 'chr1', 101, 115, hmm=(6, 20), trunc=(5, 0))
        panel = build_padded([hit], source, monkeypatch)[0]
        row = panel.rows[0]
        assert [r for r in row.pad if r[2] == 'm'] == [[0, 5, 'm']]
        assert not [r for r in row.pad if r[2] == 'g']
        # The padded positions still carry real sequence.
        assert '-' not in row.seq

    def test_truncation_by_contig_start_is_a_gap(self, monkeypatch):
        # The hit starts at position 1, so the 5 missing model positions have
        # no sequence at all: they must be gaps, not grey padding.
        source = FakeSource({'chr1': 'ACGT' * 25})
        hit = make_hit(0, 'chr1', 1, 15, hmm=(6, 20), trunc=(5, 0))
        panel = build([hit], source, monkeypatch)[0]
        row = panel.rows[0]
        assert [r for r in row.pad if r[2] == 'g'] == [[0, 5, 'g']]
        assert not [r for r in row.pad if r[2] == 'm']
        assert row.seq.startswith('-----')

    def test_truncation_by_contig_end_is_a_gap(self, monkeypatch):
        source = FakeSource({'chr1': 'ACGT' * 5})  # 20 bp
        hit = make_hit(0, 'chr1', 6, 20, hmm=(1, 15), trunc=(0, 5))
        panel = build([hit], source, monkeypatch)[0]
        row = panel.rows[0]
        assert [r for r in row.pad if r[2] == 'g'] == [[15, 5, 'g']]
        assert row.seq.endswith('-----')

    def test_partly_clipped_truncation_splits_into_both_kinds(self, monkeypatch):
        # 5 model positions missing but only 2 bp of contig before the hit:
        # 3 of them do not exist at all, the other 2 are real sequence the
        # model did not claim. The grey run sits inside the gap run.
        source = FakeSource({'chr1': 'ACGT' * 25})
        hit = make_hit(0, 'chr1', 3, 17, hmm=(6, 20), trunc=(5, 0))
        panel = build_padded([hit], source, monkeypatch)[0]
        kinds = {r[2]: r for r in panel.rows[0].pad}
        assert kinds['g'] == [0, 3, 'g']
        assert kinds['m'] == [3, 2, 'm']
        assert panel.rows[0].seq.startswith('---')


class TestModelPaddingIsOptIn:
    """By default an alignment shows only what the hit actually matched."""

    def test_unmatched_model_positions_are_gaps_by_default(self, monkeypatch):
        source = FakeSource({'chr1': 'ACGT' * 100})
        hit = make_hit(0, 'chr1', 101, 115, hmm=(6, 20), trunc=(5, 0))
        row = build([hit], source, monkeypatch)[0].rows[0]
        # No grey at all: the five unclaimed model positions are gaps.
        assert not [r for r in row.pad if r[2] == 'm']
        assert [r for r in row.pad if r[2] == 'g'] == [[0, 5, 'g']]
        assert row.seq.startswith('-----')

    def test_padding_shows_the_flanking_sequence_instead(self, monkeypatch):
        source = FakeSource({'chr1': 'ACGT' * 100})
        hit = make_hit(0, 'chr1', 101, 115, hmm=(6, 20), trunc=(5, 0))
        row = build_padded([hit], source, monkeypatch)[0].rows[0]
        assert [r for r in row.pad if r[2] == 'm'] == [[0, 5, 'm']]
        assert '-' not in row.seq

    def test_both_modes_span_the_same_columns(self, monkeypatch):
        # Rows must still line up: only the content of the unclaimed columns
        # differs, never the width.
        source = FakeSource({'chr1': 'ACGT' * 100})
        hit = make_hit(0, 'chr1', 101, 115, hmm=(6, 20), trunc=(5, 0))
        plain = build([hit], source, monkeypatch)[0]
        padded = build_padded([hit], source, monkeypatch)[0]
        assert plain.n_cols == padded.n_cols
        assert len(plain.rows[0].seq) == len(padded.rows[0].seq)

    def test_matched_sequence_is_identical_in_both_modes(self, monkeypatch):
        source = FakeSource({'chr1': 'ACGT' * 100})
        hit = make_hit(0, 'chr1', 101, 115, hmm=(6, 20), trunc=(5, 0))
        plain = build([hit], source, monkeypatch)[0].rows[0].seq
        padded = build_padded([hit], source, monkeypatch)[0].rows[0].seq
        assert plain[5:] == padded[5:]

    def test_minus_strand_gaps_land_on_the_right(self, monkeypatch):
        # The deficit is at the lower genomic coordinate, which is the model's
        # far end on the minus strand.
        source = FakeSource({'chr1': 'ACGT' * 100})
        hit = make_hit(0, 'chr1', 101, 115, strand='-', hmm=(1, 15), trunc=(5, 0))
        row = build([hit], source, monkeypatch)[0].rows[0]
        assert row.seq.endswith('-----')
        assert not [r for r in row.pad if r[2] == 'm']

    def test_contig_truncation_still_gaps_without_padding(self, monkeypatch):
        source = FakeSource({'chr1': 'ACGT' * 25})
        hit = make_hit(0, 'chr1', 1, 15, hmm=(6, 20), trunc=(5, 0))
        row = build([hit], source, monkeypatch)[0].rows[0]
        assert row.seq.startswith('-----')


class TestOrientation:
    def test_minus_strand_row_is_reverse_complemented(self, monkeypatch):
        source = FakeSource({'chr1': 'AAAACCCCGGGGTTTT'})
        plus = make_hit(0, 'chr1', 1, 8)
        minus = make_hit(1, 'chr1', 1, 8, strand='-')
        forward = build([plus], source, monkeypatch)[0].rows[0].seq
        reverse = build([minus], source, monkeypatch)[0].rows[0].seq
        assert forward == 'AAAACCCC'
        assert reverse == 'GGGGTTTT'

    def test_mafft_receives_strand_corrected_sequence(self, monkeypatch, tmp_path):
        # The rows handed to the aligner must already be in model orientation.
        # Feeding raw genomic sequence would ask MAFFT to align a terminus
        # against its own reverse complement, which it cannot do.
        seen = []

        def capture(records, tmpdir):
            seen.extend(str(r.seq) for r in records)
            return None  # force the fallback; only the input matters here

        monkeypatch.setattr('tirmite.runners.mafft.mafft_available', lambda: True)
        monkeypatch.setattr('tirmite.runners.mafft.align_in_memory', capture)

        source = FakeSource({'chr1': 'AAAACCCCGGGGTTTT'})
        hits = [make_hit(0, 'chr1', 1, 8), make_hit(1, 'chr1', 1, 8, strand='-')]
        msa.build_msa_panels(
            hits, GROUPS, MODELS, source=source, tmpdir=str(tmp_path), mode='auto'
        )
        assert seen == ['AAAACCCC', 'GGGGTTTT']

    def test_every_row_is_in_model_orientation(self, monkeypatch):
        # A mixed-strand panel must not contain a row that is the reverse
        # complement of its neighbours.
        source = FakeSource({'chr1': 'ACGTTTGCAAACGTAC' * 4})
        hits = [
            make_hit(0, 'chr1', 1, 16),
            make_hit(1, 'chr1', 17, 32, strand='-'),
            make_hit(2, 'chr1', 33, 48),
        ]
        panel = build(hits, source, monkeypatch)[0]
        rows = {row.uid: row.seq for row in panel.rows}
        genome = source.contigs['chr1']
        complement = str.maketrans('ACGT', 'TGCA')
        assert rows[0] == genome[0:16]
        assert rows[1] == genome[16:32].translate(complement)[::-1]
        assert rows[2] == genome[32:48]

    def test_minus_strand_model_pad_moves_to_the_other_end(self, monkeypatch):
        # Model positions 16-20 are unmatched, and on the minus strand they lie
        # at the *lower* genomic coordinate. Rows are shown in model
        # orientation, so the pad belongs at the right of the panel.
        source = FakeSource({'chr1': 'ACGT' * 100})
        hit = make_hit(0, 'chr1', 101, 115, strand='-', hmm=(1, 15), trunc=(5, 0))
        panel = build_padded([hit], source, monkeypatch)[0]
        assert [r for r in panel.rows[0].pad if r[2] == 'm'] == [[15, 5, 'm']]


class TestAlignerSelection:
    def test_anchor_mode_never_calls_mafft(self, monkeypatch):
        called = []
        monkeypatch.setattr(
            'tirmite.runners.mafft.mafft_available', lambda: called.append(1) or True
        )
        source = FakeSource({'chr1': 'ACGT' * 100})
        hits = [make_hit(i, 'chr1', 1 + i * 40, 20 + i * 40) for i in range(3)]
        panel = build(hits, source, monkeypatch, mode='anchor')[0]
        assert panel.aligner == 'anchor'
        assert not called

    def test_falls_back_when_mafft_is_absent(self, monkeypatch):
        monkeypatch.setattr('tirmite.runners.mafft.mafft_available', lambda: False)
        source = FakeSource({'chr1': 'ACGT' * 100})
        hits = [make_hit(i, 'chr1', 1 + i * 40, 20 + i * 40) for i in range(3)]
        panel = build(hits, source, monkeypatch, mode='auto')[0]
        assert panel.aligner == 'anchor'
        assert 'MAFFT was not available' in panel.note

    def test_falls_back_when_mafft_fails(self, monkeypatch):
        monkeypatch.setattr('tirmite.runners.mafft.mafft_available', lambda: True)
        monkeypatch.setattr(
            'tirmite.runners.mafft.align_in_memory', lambda seqs, tmpdir: None
        )
        source = FakeSource({'chr1': 'ACGT' * 100})
        hits = [make_hit(i, 'chr1', 1 + i * 40, 20 + i * 40) for i in range(3)]
        panel = build(hits, source, monkeypatch, mode='auto')[0]
        assert panel.aligner == 'anchor'

    def test_uses_mafft_when_it_succeeds(self, monkeypatch, tmp_path):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        def fake_align(records, tmpdir):
            # Return a plausible alignment: every row gapped to a common width.
            width = max(len(r.seq) for r in records) + 2
            return [
                SeqRecord(Seq(str(r.seq) + '-' * (width - len(r.seq))), id=r.id)
                for r in records
            ]

        monkeypatch.setattr('tirmite.runners.mafft.mafft_available', lambda: True)
        monkeypatch.setattr('tirmite.runners.mafft.align_in_memory', fake_align)

        source = FakeSource({'chr1': 'ACGT' * 100})
        hits = [make_hit(i, 'chr1', 1 + i * 40, 20 + i * 40) for i in range(3)]
        panels = msa.build_msa_panels(
            hits, GROUPS, MODELS, source=source, tmpdir=str(tmp_path), mode='auto'
        )
        assert panels[0].aligner == 'mafft'
        assert panels[0].n_cols == 22
        assert all(len(row.seq) == 22 for row in panels[0].rows)

    def test_single_row_does_not_go_to_mafft(self, monkeypatch, tmp_path):
        def boom(records, tmpdir):
            raise AssertionError('MAFFT must not be called for a single row')

        monkeypatch.setattr('tirmite.runners.mafft.mafft_available', lambda: True)
        monkeypatch.setattr('tirmite.runners.mafft.align_in_memory', boom)
        source = FakeSource({'chr1': 'ACGT' * 100})
        panels = msa.build_msa_panels(
            [make_hit(0, 'chr1', 1, 20)],
            GROUPS,
            MODELS,
            source=source,
            tmpdir=str(tmp_path),
            mode='auto',
        )
        assert panels[0].aligner == 'anchor'


class TestMissingModelCoordinates:
    def test_rows_are_left_aligned_with_a_caveat(self, monkeypatch):
        source = FakeSource({'chr1': 'ACGT' * 100})
        hits = [
            HitRecord(
                uid=i,
                model='TIR',
                contig='chr1',
                start=1 + i * 40,
                end=20 + i * 40,
                strand='+',
                hmm_start=None,
                hmm_end=None,
                trunc_left=None,
                trunc_right=None,
            )
            for i in range(2)
        ]
        panel = build(hits, source, monkeypatch, mode='anchor')[0]
        assert panel.aligner == 'left'
        assert 'no model coordinates' in panel.note


class TestPanelAssembly:
    def test_off_mode_builds_nothing(self, monkeypatch):
        source = FakeSource({'chr1': 'ACGT' * 100})
        assert (
            build([make_hit(0, 'chr1', 1, 20)], source, monkeypatch, mode='off') == []
        )

    def test_no_source_builds_nothing(self, monkeypatch):
        assert (
            msa.build_msa_panels(
                [make_hit(0, 'chr1', 1, 20)],
                GROUPS,
                MODELS,
                source=None,
                tmpdir='/unused',
                mode='anchor',
            )
            == []
        )

    def test_one_panel_per_model(self, monkeypatch):
        source = FakeSource({'chr1': 'ACGT' * 100})
        hits = [
            make_hit(0, 'chr1', 1, 20, model='TIR_A'),
            make_hit(1, 'chr1', 101, 120, model='TIR_B'),
        ]
        models = [
            ModelInfo(name='TIR_A', length=20),
            ModelInfo(name='TIR_B', length=20),
        ]
        panels = msa.build_msa_panels(
            hits, GROUPS, models, source=source, tmpdir='/unused', mode='anchor'
        )
        assert [p.model for p in panels] == ['TIR_A', 'TIR_B']

    def test_row_cap_keeps_the_best_hits_and_records_the_total(self, monkeypatch):
        source = FakeSource({'chr1': 'ACGT' * 200})
        hits = [
            make_hit(i, 'chr1', 1 + i * 40, 20 + i * 40, evalue=10.0 ** -(i + 1))
            for i in range(6)
        ]
        panel = build(hits, source, monkeypatch, max_rows=2)[0]
        assert panel.n_rows_shown == 2
        assert panel.n_rows_total == 6
        # Best e-value first, so the last two hits survive the cap.
        assert {row.uid for row in panel.rows} == {4, 5}

    def test_cell_budget_drops_rows(self, monkeypatch):
        source = FakeSource({'chr1': 'ACGT' * 200})
        hits = [make_hit(i, 'chr1', 1 + i * 40, 20 + i * 40) for i in range(6)]
        panel = build(hits, source, monkeypatch, max_cells=40)[0]
        assert panel.n_rows_shown == 2
        assert panel.n_cols == 20

    def test_unreadable_contig_is_skipped(self, monkeypatch):
        source = FakeSource({'chr1': 'ACGT' * 100})
        hits = [make_hit(0, 'chr1', 1, 20), make_hit(1, 'missing', 1, 20)]
        panel = build(hits, source, monkeypatch)[0]
        assert [row.uid for row in panel.rows] == [0]

    def test_model_with_no_readable_rows_yields_no_panel(self, monkeypatch):
        source = FakeSource({'chr1': 'ACGT' * 100})
        assert build([make_hit(0, 'missing', 1, 20)], source, monkeypatch) == []
