#!/usr/bin/env python3
"""
Tests for external flank extraction feature.

Validates:
1. compute_flank_coordinates: correct genomic coordinates for all
   strand/terminus combinations and offset correction.
2. _determine_terminus_type: correct inference from PairingConfig.
3. writeFlanks: end-to-end integration using an in-memory pyfaidx genome.
"""

from collections import namedtuple
import os
import tempfile

import pandas as pd

import tirmite.tirmitetools as tirmite
from tirmite.tirmitetools import (
    PairingConfig,
    _determine_terminus_type,
    compute_flank_coordinates,
    writeFlanks,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

HitRec = namedtuple(
    'HitRec', ['model', 'target', 'hitStart', 'hitEnd', 'strand', 'idx', 'evalue']
)


def _make_hitTable(rows):
    """Build a hitTable DataFrame from a list of row dicts."""
    records = []
    for r in rows:
        records.append(
            {
                'model': r['model'],
                'target': r['target'],
                'hitStart': str(r['hit_start']),
                'hitEnd': str(r['hit_end']),
                'strand': r['strand'],
                'evalue': '1e-10',
                'score': '100',
                'bias': 'NA',
                'hmmStart': str(r.get('hmm_start', 1)),
                'hmmEnd': str(r.get('hmm_end', 100)),
            }
        )
    df = pd.DataFrame(records)
    df = df.reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# compute_flank_coordinates tests
# ---------------------------------------------------------------------------


class TestComputeFlankCoordinates:
    """Unit tests for the coordinate maths of compute_flank_coordinates."""

    # ---- LEFT terminus, + strand -----------------------------------------

    def test_left_plus_full_hit(self):
        """Left terminus, + strand, full hit: external_pos == hitStart."""
        fs, fe, offset = compute_flank_coordinates(
            hit_start=1000,
            hit_end=1099,
            strand='+',
            is_left_terminus=True,
            hmm_start=1,
            hmm_end=100,
            model_len=100,
            flank_len=10,
        )
        assert offset == 0
        assert fe == 999  # external_pos - 1 = 1000 - 1
        assert fs == 990  # external_pos - flank_len = 1000 - 10

    def test_left_plus_partial_hit_offset(self):
        """Left terminus, + strand, hmmStart=3 means 2bp gap → offset=2."""
        fs, fe, offset = compute_flank_coordinates(
            hit_start=1000,
            hit_end=1097,
            strand='+',
            is_left_terminus=True,
            hmm_start=3,
            hmm_end=100,
            model_len=100,
            flank_len=10,
        )
        # external_pos = 1000 - (3-1) = 998
        assert offset == 2
        assert fe == 997  # 998 - 1
        assert fs == 988  # 998 - 10

    # ---- LEFT terminus, - strand -----------------------------------------

    def test_left_minus_full_hit(self):
        """Left terminus, - strand, full hit: hmmEnd==model_len → offset=0."""
        # For - strand: hmmStart aligns to hitEnd, hmmEnd aligns to hitStart
        # offset = model_len - hmmEnd = 100 - 100 = 0
        # external_pos = hitStart - 0 = 1000
        fs, fe, offset = compute_flank_coordinates(
            hit_start=1000,
            hit_end=1099,
            strand='-',
            is_left_terminus=True,
            hmm_start=1,
            hmm_end=100,
            model_len=100,
            flank_len=10,
        )
        assert offset == 0
        assert fe == 999
        assert fs == 990

    def test_left_minus_partial_hit_offset(self):
        """Left terminus, - strand, hmmEnd=98 means 2bp gap → offset=2."""
        # offset = model_len - hmmEnd = 100 - 98 = 2
        # external_pos = hitStart - 2 = 998
        fs, fe, offset = compute_flank_coordinates(
            hit_start=1000,
            hit_end=1097,
            strand='-',
            is_left_terminus=True,
            hmm_start=1,
            hmm_end=98,
            model_len=100,
            flank_len=10,
        )
        assert offset == 2
        assert fe == 997
        assert fs == 988

    # ---- RIGHT terminus, + strand ----------------------------------------

    def test_right_plus_full_hit(self):
        """Right terminus, + strand, full hit: external_pos == hitEnd."""
        fs, fe, offset = compute_flank_coordinates(
            hit_start=2000,
            hit_end=2099,
            strand='+',
            is_left_terminus=False,
            hmm_start=1,
            hmm_end=100,
            model_len=100,
            flank_len=10,
        )
        # offset = model_len - hmmEnd = 0
        # external_pos = hitEnd + 0 = 2099
        assert offset == 0
        assert fs == 2100  # external_pos + 1
        assert fe == 2109  # external_pos + flank_len

    def test_right_plus_partial_hit_offset(self):
        """Right terminus, + strand, hmmEnd=98 means 2bp gap → offset=2."""
        # offset = 100 - 98 = 2 → external_pos = 2099 + 2 = 2101
        fs, fe, offset = compute_flank_coordinates(
            hit_start=2000,
            hit_end=2099,
            strand='+',
            is_left_terminus=False,
            hmm_start=1,
            hmm_end=98,
            model_len=100,
            flank_len=10,
        )
        assert offset == 2
        assert fs == 2102
        assert fe == 2111

    # ---- RIGHT terminus, - strand ----------------------------------------

    def test_right_minus_full_hit(self):
        """Right terminus, - strand, full hit: hmmStart=1 → offset=0."""
        # offset = hmmStart - 1 = 0
        # external_pos = hitEnd + 0 = 2099
        fs, fe, offset = compute_flank_coordinates(
            hit_start=2000,
            hit_end=2099,
            strand='-',
            is_left_terminus=False,
            hmm_start=1,
            hmm_end=100,
            model_len=100,
            flank_len=10,
        )
        assert offset == 0
        assert fs == 2100
        assert fe == 2109

    def test_right_minus_partial_hit_offset(self):
        """Right terminus, - strand, hmmStart=3 means 2bp gap → offset=2."""
        # offset = 3 - 1 = 2 → external_pos = 2099 + 2 = 2101
        fs, fe, offset = compute_flank_coordinates(
            hit_start=2000,
            hit_end=2099,
            strand='-',
            is_left_terminus=False,
            hmm_start=3,
            hmm_end=100,
            model_len=100,
            flank_len=10,
        )
        assert offset == 2
        assert fs == 2102
        assert fe == 2111

    def test_flank_len_zero(self):
        """flank_len=0 should yield zero-length flank (start > end)."""
        fs, fe, offset = compute_flank_coordinates(
            hit_start=1000,
            hit_end=1099,
            strand='+',
            is_left_terminus=True,
            hmm_start=1,
            hmm_end=100,
            model_len=100,
            flank_len=0,
        )
        # external_pos = 1000, flank = [1000, 999] → zero-length
        assert fe < fs


# ---------------------------------------------------------------------------
# _determine_terminus_type tests
# ---------------------------------------------------------------------------


class TestDetermineTerminusType:
    """Tests for _determine_terminus_type terminus inference."""

    def _hit(self, model, strand):
        return HitRec(
            model=model,
            target='chr1',
            hitStart=100,
            hitEnd=200,
            strand=strand,
            idx=0,
            evalue='1e-10',
        )

    def test_asymmetric_left_model(self):
        config = PairingConfig(orientation='F,R', left_model='L', right_model='R')
        assert _determine_terminus_type(self._hit('L', '+'), config) == 'left'

    def test_asymmetric_right_model(self):
        config = PairingConfig(orientation='F,R', left_model='L', right_model='R')
        assert _determine_terminus_type(self._hit('R', '-'), config) == 'right'

    def test_asymmetric_unknown_model(self):
        config = PairingConfig(orientation='F,R', left_model='L', right_model='R')
        assert _determine_terminus_type(self._hit('Other', '+'), config) is None

    def test_symmetric_FR_left_strand(self):
        """For F,R orientation: + strand → left, - strand → right."""
        config = PairingConfig(orientation='F,R', single_model='M')
        assert _determine_terminus_type(self._hit('M', '+'), config) == 'left'

    def test_symmetric_FR_right_strand(self):
        config = PairingConfig(orientation='F,R', single_model='M')
        assert _determine_terminus_type(self._hit('M', '-'), config) == 'right'

    def test_symmetric_RF_left_strand(self):
        """For R,F orientation: - strand → left, + strand → right."""
        config = PairingConfig(orientation='R,F', single_model='M')
        assert _determine_terminus_type(self._hit('M', '-'), config) == 'left'
        assert _determine_terminus_type(self._hit('M', '+'), config) == 'right'

    def test_symmetric_FF_returns_none_for_unpaired(self):
        """F,F same-strand: cannot determine terminus type → None."""
        config = PairingConfig(orientation='F,F', single_model='M')
        assert _determine_terminus_type(self._hit('M', '+'), config) is None

    def test_symmetric_RR_returns_none_for_unpaired(self):
        """R,R same-strand: cannot determine terminus type → None."""
        config = PairingConfig(orientation='R,R', single_model='M')
        assert _determine_terminus_type(self._hit('M', '-'), config) is None


# ---------------------------------------------------------------------------
# writeFlanks integration tests using a mock pyfaidx-like genome
# ---------------------------------------------------------------------------


class MockChrom:
    """Minimal pyfaidx Fasta chromosome substitute."""

    def __init__(self, seq):
        self._seq = seq

    def __len__(self):
        return len(self._seq)

    def __getitem__(self, key):
        if isinstance(key, slice):
            return self._seq[key]
        return self._seq[key]

    def __str__(self):
        return self._seq


class MockGenome(dict):
    """Dict-based mock for pyfaidx.Fasta."""

    def __getitem__(self, key):
        return MockChrom(super().__getitem__(key))


def _make_paired_data(hitTable):
    """Create hitsDict, hitIndex and a simple paired dict for one pair."""
    hitsDict, hitIndex = tirmite.table2dict(hitTable)
    # Build paired: assume first two hits (indices 0 and 1) form a pair
    idx_list = list(hitTable.index)
    assert len(idx_list) >= 2, 'Need at least 2 hits for pairing test'
    model = hitTable.iloc[0]['model']
    paired = {model: [{idx_list[0], idx_list[1]}]}
    # Mark them as partners in hitIndex
    m0 = hitTable.iloc[0]['model']
    m1 = hitTable.iloc[1]['model']
    hitIndex[m0][idx_list[0]]['partner'] = idx_list[1]
    hitIndex[m1][idx_list[1]]['partner'] = idx_list[0]
    return hitsDict, hitIndex, paired


class TestWriteFlanks:
    """Integration tests for writeFlanks output."""

    GENOME_SEQ = 'A' * 500 + 'C' * 500  # 1000bp: pos 1-500=A, 501-1000=C

    def _genome(self):
        return MockGenome({'chr1': self.GENOME_SEQ})

    def test_FR_orientation_paired_flanks(self):
        """F,R symmetric model: left flank upstream, right flank downstream."""
        # Left hit: + strand at 200-299 (hmmStart=1, hmmEnd=100)
        # Right hit: - strand at 700-799 (hmmStart=1, hmmEnd=100)
        rows = [
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 200,
                'hit_end': 299,
                'strand': '+',
                'hmm_start': 1,
                'hmm_end': 100,
            },
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 700,
                'hit_end': 799,
                'strand': '-',
                'hmm_start': 1,
                'hmm_end': 100,
            },
        ]
        hitTable = _make_hitTable(rows)
        _, hitIndex, paired = _make_paired_data(hitTable)

        config = PairingConfig(orientation='F,R', single_model='TIR')
        model_lengths = {'TIR': 100}

        with tempfile.TemporaryDirectory() as tmpdir:
            writeFlanks(
                outDir=tmpdir,
                hitTable=hitTable,
                model_lengths=model_lengths,
                paired=paired,
                hitIndex=hitIndex,
                config=config,
                genome=self._genome(),
                flank_len=10,
            )

            left_files = [f for f in os.listdir(tmpdir) if 'left_flank' in f]
            right_files = [f for f in os.listdir(tmpdir) if 'right_flank' in f]

            assert len(left_files) == 1, 'Expected one left flank file'
            assert len(right_files) == 1, 'Expected one right flank file'

            # Left flank: upstream of pos 200 → bases 190-199 (1-based)
            # All 'A' in our mock genome
            with open(os.path.join(tmpdir, left_files[0])) as fh:
                lines = fh.read().splitlines()
            seq_lines = [line for line in lines if not line.startswith('>')]
            left_seq = ''.join(seq_lines)
            assert len(left_seq) == 10
            assert left_seq == 'A' * 10, f'Expected all-A left flank, got {left_seq!r}'

            # Right flank: downstream of pos 799 → bases 800-809 (1-based)
            # Our mock has C's from pos 501 onwards
            with open(os.path.join(tmpdir, right_files[0])) as fh:
                lines = fh.read().splitlines()
            seq_lines = [line for line in lines if not line.startswith('>')]
            right_seq = ''.join(seq_lines)
            assert len(right_seq) == 10
            assert right_seq == 'C' * 10, (
                f'Expected all-C right flank, got {right_seq!r}'
            )

    def test_FF_orientation_paired_flanks(self):
        """F,F asymmetric: left model + strand, right model + strand."""
        # left_model hit at 200-299 (+), right_model hit at 700-799 (+)
        rows = [
            {
                'model': 'LEFT',
                'target': 'chr1',
                'hit_start': 200,
                'hit_end': 299,
                'strand': '+',
                'hmm_start': 1,
                'hmm_end': 100,
            },
            {
                'model': 'RIGHT',
                'target': 'chr1',
                'hit_start': 700,
                'hit_end': 799,
                'strand': '+',
                'hmm_start': 1,
                'hmm_end': 100,
            },
        ]
        hitTable = _make_hitTable(rows)
        _, hitIndex, paired = _make_paired_data(hitTable)
        # paired is keyed by first model (LEFT); manually correct:
        paired = {'LEFT': [{0, 1}]}
        hitIndex['LEFT'][0]['partner'] = 1
        hitIndex['RIGHT'][1]['partner'] = 0

        config = PairingConfig(
            orientation='F,F', left_model='LEFT', right_model='RIGHT'
        )
        model_lengths = {'LEFT': 100, 'RIGHT': 100}

        with tempfile.TemporaryDirectory() as tmpdir:
            writeFlanks(
                outDir=tmpdir,
                hitTable=hitTable,
                model_lengths=model_lengths,
                paired=paired,
                hitIndex=hitIndex,
                config=config,
                genome=self._genome(),
                flank_len=10,
            )

            left_files = [f for f in os.listdir(tmpdir) if 'left_flank' in f]
            right_files = [f for f in os.listdir(tmpdir) if 'right_flank' in f]

            assert len(left_files) == 1
            assert len(right_files) == 1

    def test_offset_correction_left_plus(self):
        """hmmStart=3 shifts left flank 2bp further upstream."""
        rows = [
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 200,
                'hit_end': 297,  # 98bp match, hmmStart=3
                'strand': '+',
                'hmm_start': 3,
                'hmm_end': 100,
            },
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 700,
                'hit_end': 799,
                'strand': '-',
                'hmm_start': 1,
                'hmm_end': 100,
            },
        ]
        hitTable = _make_hitTable(rows)
        _, hitIndex, paired = _make_paired_data(hitTable)

        config = PairingConfig(orientation='F,R', single_model='TIR')
        model_lengths = {'TIR': 100}

        with tempfile.TemporaryDirectory() as tmpdir:
            writeFlanks(
                outDir=tmpdir,
                hitTable=hitTable,
                model_lengths=model_lengths,
                paired=paired,
                hitIndex=hitIndex,
                config=config,
                genome=self._genome(),
                flank_len=5,
            )

            left_files = [f for f in os.listdir(tmpdir) if 'left_flank' in f]
            assert len(left_files) == 1

            with open(os.path.join(tmpdir, left_files[0])) as fh:
                content = fh.read()
            # Left flank: external_pos = 200 - 2 = 198; flank is [193, 197]
            # All 'A' in mock genome
            seq = ''.join(
                line for line in content.splitlines() if not line.startswith('>')
            )
            assert len(seq) == 5

    def test_flank_max_offset_skips_hit(self):
        """Hits whose offset exceeds flank_max_offset should be skipped."""
        rows = [
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 200,
                'hit_end': 295,  # hmmStart=6 → offset=5
                'strand': '+',
                'hmm_start': 6,
                'hmm_end': 100,
            },
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 700,
                'hit_end': 799,
                'strand': '-',
                'hmm_start': 1,
                'hmm_end': 100,
            },
        ]
        hitTable = _make_hitTable(rows)
        _, hitIndex, paired = _make_paired_data(hitTable)

        config = PairingConfig(orientation='F,R', single_model='TIR')
        model_lengths = {'TIR': 100}

        with tempfile.TemporaryDirectory() as tmpdir:
            # max_offset=3 < actual offset=5 for left hit → no left flank file
            writeFlanks(
                outDir=tmpdir,
                hitTable=hitTable,
                model_lengths=model_lengths,
                paired=paired,
                hitIndex=hitIndex,
                config=config,
                genome=self._genome(),
                flank_len=10,
                flank_max_offset=3,
            )

            files = os.listdir(tmpdir)
            left_files = [f for f in files if 'left_flank' in f]
            # Left hit has offset=5 > 3, so no left flank file
            assert len(left_files) == 0

    def test_paired_only_flag_skips_unpaired(self):
        """When paired_only=True, unpaired hits must not produce flank files."""
        rows = [
            # Paired hits
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 200,
                'hit_end': 299,
                'strand': '+',
                'hmm_start': 1,
                'hmm_end': 100,
            },
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 700,
                'hit_end': 799,
                'strand': '-',
                'hmm_start': 1,
                'hmm_end': 100,
            },
            # Unpaired hit
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 400,
                'hit_end': 499,
                'strand': '+',
                'hmm_start': 1,
                'hmm_end': 100,
            },
        ]
        hitTable = _make_hitTable(rows)
        hitsDict, hitIndex = tirmite.table2dict(hitTable)
        paired = {'TIR': [{0, 1}]}
        hitIndex['TIR'][0]['partner'] = 1
        hitIndex['TIR'][1]['partner'] = 0
        # index 2 stays unpaired

        config = PairingConfig(orientation='F,R', single_model='TIR')
        model_lengths = {'TIR': 100}

        with tempfile.TemporaryDirectory() as tmpdir:
            writeFlanks(
                outDir=tmpdir,
                hitTable=hitTable,
                model_lengths=model_lengths,
                paired=paired,
                hitIndex=hitIndex,
                config=config,
                genome=self._genome(),
                flank_len=10,
                paired_only=True,
            )

            left_files = [f for f in os.listdir(tmpdir) if 'left_flank' in f]
            right_files = [f for f in os.listdir(tmpdir) if 'right_flank' in f]
            # Only the paired hits produce flanks → 1 left + 1 right
            assert len(left_files) == 1
            assert len(right_files) == 1

            with open(os.path.join(tmpdir, left_files[0])) as fh:
                seqs = [line for line in fh.read().splitlines() if line.startswith('>')]
            # Only 1 sequence (from the single pair), not 2
            assert len(seqs) == 1

    def test_all_hits_includes_unpaired(self):
        """When paired_only=False, unpaired hits also produce flanks."""
        rows = [
            # Paired
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 200,
                'hit_end': 299,
                'strand': '+',
                'hmm_start': 1,
                'hmm_end': 100,
            },
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 700,
                'hit_end': 799,
                'strand': '-',
                'hmm_start': 1,
                'hmm_end': 100,
            },
            # Unpaired +strand hit (will be treated as left terminus)
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 50,
                'hit_end': 149,
                'strand': '+',
                'hmm_start': 1,
                'hmm_end': 100,
            },
        ]
        hitTable = _make_hitTable(rows)
        hitsDict, hitIndex = tirmite.table2dict(hitTable)
        paired = {'TIR': [{0, 1}]}
        hitIndex['TIR'][0]['partner'] = 1
        hitIndex['TIR'][1]['partner'] = 0
        # index 2 is unpaired

        # F,R: + strand = left terminus → unpaired + hit should produce a left flank
        config = PairingConfig(orientation='F,R', single_model='TIR')
        model_lengths = {'TIR': 100}

        with tempfile.TemporaryDirectory() as tmpdir:
            writeFlanks(
                outDir=tmpdir,
                hitTable=hitTable,
                model_lengths=model_lengths,
                paired=paired,
                hitIndex=hitIndex,
                config=config,
                genome=self._genome(),
                flank_len=10,
                paired_only=False,
            )

            left_files = [f for f in os.listdir(tmpdir) if 'left_flank' in f]
            assert len(left_files) == 1

            with open(os.path.join(tmpdir, left_files[0])) as fh:
                seqs = [line for line in fh.read().splitlines() if line.startswith('>')]
            # 1 from paired + 1 from unpaired = 2 left flank seqs
            assert len(seqs) == 2

    def test_rr_orientation_unpaired_returns_none(self):
        """R,R same-strand: unpaired hits can't be classified, no extra flanks."""
        rows = [
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 200,
                'hit_end': 299,
                'strand': '-',
                'hmm_start': 1,
                'hmm_end': 100,
            },
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 700,
                'hit_end': 799,
                'strand': '-',
                'hmm_start': 1,
                'hmm_end': 100,
            },
            # Extra unpaired - strand hit
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 50,
                'hit_end': 149,
                'strand': '-',
                'hmm_start': 1,
                'hmm_end': 100,
            },
        ]
        hitTable = _make_hitTable(rows)
        hitsDict, hitIndex = tirmite.table2dict(hitTable)
        paired = {'TIR': [{0, 1}]}
        hitIndex['TIR'][0]['partner'] = 1
        hitIndex['TIR'][1]['partner'] = 0

        config = PairingConfig(orientation='R,R', single_model='TIR')
        model_lengths = {'TIR': 100}

        with tempfile.TemporaryDirectory() as tmpdir:
            writeFlanks(
                outDir=tmpdir,
                hitTable=hitTable,
                model_lengths=model_lengths,
                paired=paired,
                hitIndex=hitIndex,
                config=config,
                genome=self._genome(),
                flank_len=10,
                paired_only=False,  # allow unpaired, but can't classify them
            )

            all_files = os.listdir(tmpdir)
            flank_files = [f for f in all_files if 'flank' in f]
            # Only the paired hits produce flanks
            assert len(flank_files) == 2  # one left, one right for the pair

    def test_no_flanks_when_flank_len_zero(self):
        """flank_len=0 should produce no output files."""
        rows = [
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 200,
                'hit_end': 299,
                'strand': '+',
                'hmm_start': 1,
                'hmm_end': 100,
            },
            {
                'model': 'TIR',
                'target': 'chr1',
                'hit_start': 700,
                'hit_end': 799,
                'strand': '-',
                'hmm_start': 1,
                'hmm_end': 100,
            },
        ]
        hitTable = _make_hitTable(rows)
        _, hitIndex, paired = _make_paired_data(hitTable)
        config = PairingConfig(orientation='F,R', single_model='TIR')
        model_lengths = {'TIR': 100}

        with tempfile.TemporaryDirectory() as tmpdir:
            writeFlanks(
                outDir=tmpdir,
                hitTable=hitTable,
                model_lengths=model_lengths,
                paired=paired,
                hitIndex=hitIndex,
                config=config,
                genome=self._genome(),
                flank_len=0,
            )
            flank_files = [f for f in os.listdir(tmpdir) if 'flank' in f]
            assert len(flank_files) == 0
