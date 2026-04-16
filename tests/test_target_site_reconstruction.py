#!/usr/bin/env python3
"""
Tests for target site reconstruction and validation features.

Validates:
1. hamming_distance: correct distance computation.
2. load_tsd_length_map: parsing of TSD length map files.
3. reconstruct_target_site: correct TSD deduplication for both in-model and
   outside-model modes.
4. format_interleaved_flanks: correct gap padding.
5. writeTargetSites: end-to-end integration with mock genome.
6. validate CLI: parser configuration and argument validation.
"""

import os
import tempfile

import pytest

from tirmite.tirmitetools import (
    PairingConfig,
    compute_flank_coordinates,
    compute_inner_tsd_coordinates,
    format_interleaved_flanks,
    hamming_distance,
    load_tsd_length_map,
    reconstruct_target_site,
    writeTargetSites,
)

# ---------------------------------------------------------------------------
# hamming_distance tests
# ---------------------------------------------------------------------------


class TestHammingDistance:
    """Tests for hamming_distance function."""

    def test_identical(self):
        assert hamming_distance('ATCG', 'ATCG') == 0

    def test_all_different(self):
        assert hamming_distance('AAAA', 'TTTT') == 4

    def test_one_mismatch(self):
        assert hamming_distance('ATCG', 'AACG') == 1

    def test_empty_strings(self):
        assert hamming_distance('', '') == 0

    def test_unequal_length_raises(self):
        with pytest.raises(ValueError, match='equal length'):
            hamming_distance('AT', 'ATC')

    def test_case_sensitive(self):
        # Upper vs lower should be different
        assert hamming_distance('atcg', 'ATCG') == 4


# ---------------------------------------------------------------------------
# compute_inner_tsd_coordinates tests
# ---------------------------------------------------------------------------


class TestComputeInnerTsdCoordinates:
    """Tests for compute_inner_tsd_coordinates function.

    Verifies that the inner-boundary TSD coordinates are correctly computed
    for all four strand × terminus-type combinations.

    Coordinate convention (same as compute_flank_coordinates):
      +strand: hmmStart aligns to hit_start; hmmEnd aligns to hit_end.
      -strand: hmmStart aligns to hit_end;  hmmEnd aligns to hit_start.
    """

    # ---- Left terminus, + strand ----------------------------------------
    def test_left_plus_full_model_coverage(self):
        """Left terminus + strand, full model coverage (hmm_start=1, hmm_end=model_len)."""
        # hit: 100-150, model_len=50, tsd_length=5
        # inner_pos = hit_end + (model_len - hmm_end) = 150 + (50-50) = 150
        # TSD: [146, 150]
        tsd_start, tsd_end = compute_inner_tsd_coordinates(
            hit_start=100,
            hit_end=150,
            strand='+',
            is_left_terminus=True,
            hmm_start=1,
            hmm_end=50,
            model_len=50,
            tsd_length=5,
        )
        assert tsd_start == 146
        assert tsd_end == 150

    def test_left_plus_inner_offset(self):
        """Left terminus + strand, inner offset of 3 (hmm_end = model_len - 3)."""
        # hit: 100-147, model_len=50, hmm_end=47, tsd_length=5
        # inner_pos = 147 + (50-47) = 150
        # TSD: [146, 150]
        tsd_start, tsd_end = compute_inner_tsd_coordinates(
            hit_start=100,
            hit_end=147,
            strand='+',
            is_left_terminus=True,
            hmm_start=1,
            hmm_end=47,
            model_len=50,
            tsd_length=5,
        )
        assert tsd_start == 146
        assert tsd_end == 150

    # ---- Left terminus, - strand ----------------------------------------
    def test_left_minus_full_model_coverage(self):
        """Left terminus - strand, full model coverage.

        For -strand hits: hmmStart(1) aligns to hit_end (higher coord).
        inner_pos = hit_end + (hmm_start - 1) = 150 + 0 = 150.
        TSD occupies the last tsd_length positions from inner_pos: [146, 150].
        """
        tsd_start, tsd_end = compute_inner_tsd_coordinates(
            hit_start=100,
            hit_end=150,
            strand='-',
            is_left_terminus=True,
            hmm_start=1,
            hmm_end=50,
            model_len=50,
            tsd_length=5,
        )
        assert tsd_start == 146
        assert tsd_end == 150

    def test_left_minus_outer_offset(self):
        """Left terminus - strand, outer offset of 2 (hmm_start=3 means 2 positions missing at outer edge)."""
        # For -strand: hmmStart(3) aligns to hit_end(150)
        # inner_pos = 150 + (3 - 1) = 152
        # TSD: [148, 152]
        tsd_start, tsd_end = compute_inner_tsd_coordinates(
            hit_start=100,
            hit_end=150,
            strand='-',
            is_left_terminus=True,
            hmm_start=3,
            hmm_end=50,
            model_len=50,
            tsd_length=5,
        )
        assert tsd_start == 148
        assert tsd_end == 152

    # ---- Right terminus, + strand ---------------------------------------
    def test_right_plus_full_model_coverage(self):
        """Right terminus + strand, full model coverage.

        inner_pos = hit_start - (hmm_start - 1) = 200 - 0 = 200.
        TSD: [200, 204].
        """
        tsd_start, tsd_end = compute_inner_tsd_coordinates(
            hit_start=200,
            hit_end=250,
            strand='+',
            is_left_terminus=False,
            hmm_start=1,
            hmm_end=50,
            model_len=50,
            tsd_length=5,
        )
        assert tsd_start == 200
        assert tsd_end == 204

    def test_right_plus_inner_offset(self):
        """Right terminus + strand, inner offset of 3 (hmm_start=4 means 3 positions before alignment)."""
        # inner_pos = hit_start - (hmm_start - 1) = 200 - 3 = 197
        # TSD: [197, 201]
        tsd_start, tsd_end = compute_inner_tsd_coordinates(
            hit_start=200,
            hit_end=250,
            strand='+',
            is_left_terminus=False,
            hmm_start=4,
            hmm_end=50,
            model_len=50,
            tsd_length=5,
        )
        assert tsd_start == 197
        assert tsd_end == 201

    # ---- Right terminus, - strand ---------------------------------------
    def test_right_minus_full_model_coverage(self):
        """Right terminus - strand, full model coverage.

        For -strand: hmmStart(1) aligns to hit_end(250), hmmEnd(50) aligns to hit_start(200).
        inner_pos = hit_start - (model_len - hmm_end) = 200 - 0 = 200.
        TSD: [200, 204].
        """
        tsd_start, tsd_end = compute_inner_tsd_coordinates(
            hit_start=200,
            hit_end=250,
            strand='-',
            is_left_terminus=False,
            hmm_start=1,
            hmm_end=50,
            model_len=50,
            tsd_length=5,
        )
        assert tsd_start == 200
        assert tsd_end == 204

    def test_right_minus_inner_offset(self):
        """Right terminus - strand with inner offset (hmm_end < model_len on - strand)."""
        # hmm_end=47 means 3 model positions not covered at inner end
        # inner_pos = hit_start - (model_len - hmm_end) = 200 - (50-47) = 197
        # TSD: [197, 201]
        tsd_start, tsd_end = compute_inner_tsd_coordinates(
            hit_start=200,
            hit_end=250,
            strand='-',
            is_left_terminus=False,
            hmm_start=1,
            hmm_end=47,
            model_len=50,
            tsd_length=5,
        )
        assert tsd_start == 197
        assert tsd_end == 201

    def test_tsd_length_one(self):
        """Single-base TSD."""
        tsd_start, tsd_end = compute_inner_tsd_coordinates(
            hit_start=100,
            hit_end=150,
            strand='+',
            is_left_terminus=True,
            hmm_start=1,
            hmm_end=50,
            model_len=50,
            tsd_length=1,
        )
        assert tsd_start == 150
        assert tsd_end == 150

    def test_inner_outer_complementary_with_flank(self):
        """Inner TSD start should be adjacent to (or touching) the inner end of the hit.

        For left terminus + strand, full coverage:
          compute_flank_coordinates gives flank just LEFT of the OUTER boundary.
          compute_inner_tsd_coordinates gives TSD at the INNER boundary (right end).
          These should be at opposite ends of the hit region.
        """
        hit_start, hit_end = 100, 150
        hmm_start, hmm_end, model_len = 1, 50, 50

        flank_start, flank_end, offset = compute_flank_coordinates(
            hit_start=hit_start,
            hit_end=hit_end,
            strand='+',
            is_left_terminus=True,
            hmm_start=hmm_start,
            hmm_end=hmm_end,
            model_len=model_len,
            flank_len=10,
        )
        tsd_start, tsd_end = compute_inner_tsd_coordinates(
            hit_start=hit_start,
            hit_end=hit_end,
            strand='+',
            is_left_terminus=True,
            hmm_start=hmm_start,
            hmm_end=hmm_end,
            model_len=model_len,
            tsd_length=5,
        )
        # Outer flank ends just left of outer boundary (hit_start=100)
        assert flank_end == hit_start - 1  # 99
        # Inner TSD ends at the inner boundary (hit_end=150)
        assert tsd_end == hit_end  # 150


# ---------------------------------------------------------------------------
# load_tsd_length_map tests
# ---------------------------------------------------------------------------


class TestLoadTsdLengthMap:
    """Tests for load_tsd_length_map function."""

    def test_basic_load(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write('modelA\tmodelB\t5\n')
            f.write('modelC\tmodelD\t8\n')
        try:
            result = load_tsd_length_map(f.name)
            assert result['modelA\tmodelB'] == 5
            assert result['modelC\tmodelD'] == 8
        finally:
            os.unlink(f.name)

    def test_comment_and_blank_lines(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write('# comment\n')
            f.write('\n')
            f.write('modelA\tmodelB\t3\n')
        try:
            result = load_tsd_length_map(f.name)
            assert len(result) == 1
            assert result['modelA\tmodelB'] == 3
        finally:
            os.unlink(f.name)

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            load_tsd_length_map('/nonexistent/file.tsv')

    def test_invalid_format(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write('modelA\tmodelB\n')  # Missing TSD length
        try:
            with pytest.raises(ValueError, match='3 tab-delimited columns'):
                load_tsd_length_map(f.name)
        finally:
            os.unlink(f.name)

    def test_invalid_tsd_length(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write('modelA\tmodelB\tabc\n')
        try:
            with pytest.raises(ValueError, match='Invalid TSD length'):
                load_tsd_length_map(f.name)
        finally:
            os.unlink(f.name)

    def test_negative_tsd_length(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write('modelA\tmodelB\t-5\n')
        try:
            with pytest.raises(ValueError, match='non-negative'):
                load_tsd_length_map(f.name)
        finally:
            os.unlink(f.name)

    def test_empty_file(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write('# only comments\n')
        try:
            with pytest.raises(ValueError, match='No valid TSD lengths'):
                load_tsd_length_map(f.name)
        finally:
            os.unlink(f.name)


# ---------------------------------------------------------------------------
# reconstruct_target_site tests
# ---------------------------------------------------------------------------


class TestReconstructTargetSite:
    """Tests for reconstruct_target_site function."""

    def test_no_tsd(self):
        """Without TSD, simply concatenate flanks."""
        ts, l_tsd, r_tsd, ham = reconstruct_target_site('AAAAA', 'CCCCC', tsd_length=0)
        assert ts == 'AAAAACCCCC'
        assert l_tsd == ''
        assert r_tsd == ''
        assert ham == 0

    def test_tsd_outside_model_exact_match(self):
        """TSD outside model: last 3 of left flank = first 3 of right flank."""
        ts, l_tsd, r_tsd, ham = reconstruct_target_site(
            'AAAAATGA', 'TGACCCCC', tsd_length=3, tsd_in_model=False
        )
        assert ts == 'AAAAATGACCCCC'
        assert l_tsd == 'TGA'
        assert r_tsd == 'TGA'
        assert ham == 0

    def test_tsd_outside_model_mismatch(self):
        """TSD outside model with mismatches should report hamming distance."""
        ts, l_tsd, r_tsd, ham = reconstruct_target_site(
            'AAAAATGA', 'TGCCCCCC', tsd_length=3, tsd_in_model=False
        )
        assert l_tsd == 'TGA'
        assert r_tsd == 'TGC'
        assert ham == 1
        # Right TSD is trimmed, so target site is left_flank + right_flank[3:]
        assert ts == 'AAAAATGA' + 'CCCCC'

    def test_tsd_in_model_exact_match(self):
        """TSD inside model: last 3 of left flank = first 3 of right flank."""
        ts, l_tsd, r_tsd, ham = reconstruct_target_site(
            'AAAAATGA', 'TGACCCCC', tsd_length=3, tsd_in_model=True
        )
        assert ts == 'AAAAATGACCCCC'
        assert l_tsd == 'TGA'
        assert r_tsd == 'TGA'
        assert ham == 0

    def test_tsd_in_model_mismatch(self):
        """TSD inside model with mismatches."""
        ts, l_tsd, r_tsd, ham = reconstruct_target_site(
            'AAAAATGA', 'TGCCCCCC', tsd_length=3, tsd_in_model=True
        )
        assert l_tsd == 'TGA'
        assert r_tsd == 'TGC'
        assert ham == 1

    def test_tsd_length_exceeds_flank(self):
        """When TSD length exceeds flank length, use available sequence."""
        ts, l_tsd, r_tsd, ham = reconstruct_target_site(
            'AT', 'CG', tsd_length=5, tsd_in_model=False
        )
        assert l_tsd == 'AT'
        assert r_tsd == 'CG'
        # Both TSDs are 2bp (clamped to available length), hamming distance computed
        assert ham == 2

    def test_single_base_tsd(self):
        """Single base TSD."""
        ts, l_tsd, r_tsd, ham = reconstruct_target_site(
            'AAAAAG', 'GCCCCC', tsd_length=1, tsd_in_model=False
        )
        assert l_tsd == 'G'
        assert r_tsd == 'G'
        assert ham == 0
        assert ts == 'AAAAAGCCCCC'


# ---------------------------------------------------------------------------
# format_interleaved_flanks tests
# ---------------------------------------------------------------------------


class TestFormatInterleavedFlanks:
    """Tests for format_interleaved_flanks function."""

    def test_no_tsd(self):
        left, right = format_interleaved_flanks('AAAAA', 'CCCCC', tsd_length=0)
        # With no overlap, left is padded with gaps for full right length
        assert left == 'AAAAA-----'
        assert right == '-----CCCCC'
        assert len(left) == len(right)

    def test_with_tsd(self):
        left, right = format_interleaved_flanks('AAAAATSD', 'TSDCCCCC', tsd_length=3)
        # left_row = 'AAAAATSD' + '-' * (8 - 3) = 'AAAAATSD-----'
        # right_row = '-' * (8 - 3) + 'TSDCCCCC' = '-----TSDCCCCC'
        assert left == 'AAAAATSD-----'
        assert right == '-----TSDCCCCC'
        assert len(left) == len(right)

    def test_equal_length_flanks(self):
        left, right = format_interleaved_flanks('AAAT', 'TAAA', tsd_length=1)
        assert left == 'AAAT---'
        assert right == '---TAAA'
        assert len(left) == len(right)

    def test_zero_length_flanks(self):
        left, right = format_interleaved_flanks('', '', tsd_length=0)
        assert left == ''
        assert right == ''


# ---------------------------------------------------------------------------
# writeTargetSites integration tests
# ---------------------------------------------------------------------------


class MockChrom:
    """Minimal pyfaidx Fasta chromosome substitute."""

    def __init__(self, seq):
        self._seq = seq

    def __len__(self):
        return len(self._seq)

    def __getitem__(self, key):
        return self._seq[key]


class MockGenome(dict):
    """Dict-based mock for pyfaidx.Fasta."""

    def __init__(self, chrom_seqs):
        for name, seq in chrom_seqs.items():
            self[name] = MockChrom(seq)


class TestWriteTargetSites:
    """Integration tests for writeTargetSites."""

    # Genome layout:
    #   positions 1-200: A's (upstream flank region)
    #   positions 201-203: T's (TSD-like region)
    #   positions 204-213: G's (element region)
    #   positions 214-216: T's (TSD-like region)
    #   positions 217-416: C's (downstream flank region)
    UPSTREAM_LEN = 200
    TSD_LEN = 3
    ELEMENT_LEN = 10
    DOWNSTREAM_LEN = 200
    GENOME_SEQ = (
        'A' * UPSTREAM_LEN
        + 'T' * TSD_LEN
        + 'G' * ELEMENT_LEN
        + 'T' * TSD_LEN
        + 'C' * DOWNSTREAM_LEN
    )

    def _genome(self):
        return MockGenome({'chr1': self.GENOME_SEQ})

    def _make_test_data(self):
        """Create minimal paired data for testing target site reconstruction."""
        from collections import namedtuple

        import pandas as pd

        HitRec = namedtuple(
            'HitRec',
            ['model', 'target', 'hitStart', 'hitEnd', 'strand', 'idx', 'evalue'],
        )

        # Simulate a pair: left hit on +, right hit on -
        # left terminus at position 201-210 on + strand
        # right terminus at position 214-223 on - strand
        hit1 = HitRec('myModel', 'chr1', 201, 210, '+', 0, 1e-10)
        hit2 = HitRec('myModel', 'chr1', 214, 223, '-', 1, 1e-10)

        hitTable = pd.DataFrame(
            [
                {
                    'model': 'myModel',
                    'target': 'chr1',
                    'hitStart': 201,
                    'hitEnd': 210,
                    'strand': '+',
                    'evalue': 1e-10,
                    'hmmStart': 1,
                    'hmmEnd': 10,
                },
                {
                    'model': 'myModel',
                    'target': 'chr1',
                    'hitStart': 214,
                    'hitEnd': 223,
                    'strand': '-',
                    'evalue': 1e-10,
                    'hmmStart': 1,
                    'hmmEnd': 10,
                },
            ]
        )

        model_lengths = {'myModel': 10}

        hitIndex = {
            'myModel': {
                0: {'rec': hit1, 'partner': 1, 'candidates': []},
                1: {'rec': hit2, 'partner': 0, 'candidates': []},
            }
        }

        paired = {'myModel': [{0, 1}]}

        config = PairingConfig(orientation='F,R', single_model='myModel')

        return hitTable, model_lengths, hitIndex, paired, config

    def test_target_site_no_tsd(self):
        """Target site reconstruction without TSD."""
        hitTable, model_lengths, hitIndex, paired, config = self._make_test_data()

        with tempfile.TemporaryDirectory() as tmpdir:
            writeTargetSites(
                outDir=tmpdir,
                hitTable=hitTable,
                model_lengths=model_lengths,
                paired=paired,
                hitIndex=hitIndex,
                config=config,
                genome=self._genome(),
                flank_len=5,
                tsd_length=0,
            )

            # Check that target sites file was created
            ts_files = [f for f in os.listdir(tmpdir) if 'target_sites' in f]
            assert len(ts_files) == 1

            # Check interleaved flanks file was created
            il_files = [f for f in os.listdir(tmpdir) if 'interleaved_flanks' in f]
            assert len(il_files) == 1

    def test_target_site_with_tsd_outside_model(self):
        """Target site reconstruction with TSD outside model."""
        hitTable, model_lengths, hitIndex, paired, config = self._make_test_data()

        with tempfile.TemporaryDirectory() as tmpdir:
            writeTargetSites(
                outDir=tmpdir,
                hitTable=hitTable,
                model_lengths=model_lengths,
                paired=paired,
                hitIndex=hitIndex,
                config=config,
                genome=self._genome(),
                flank_len=5,
                tsd_length=3,
                tsd_in_model=False,
            )

            ts_files = [f for f in os.listdir(tmpdir) if 'target_sites' in f]
            assert len(ts_files) == 1

            # Read and check FASTA header contains metadata
            with open(os.path.join(tmpdir, ts_files[0])) as fh:
                content = fh.read()
            assert 'tsd_len=3' in content
            assert 'tsd_in_model=False' in content
            assert 'tsd_hamming=' in content

    def test_target_site_with_tsd_in_model(self):
        """Target site reconstruction with TSD inside model."""
        hitTable, model_lengths, hitIndex, paired, config = self._make_test_data()

        with tempfile.TemporaryDirectory() as tmpdir:
            writeTargetSites(
                outDir=tmpdir,
                hitTable=hitTable,
                model_lengths=model_lengths,
                paired=paired,
                hitIndex=hitIndex,
                config=config,
                genome=self._genome(),
                flank_len=5,
                tsd_length=3,
                tsd_in_model=True,
            )

            ts_files = [f for f in os.listdir(tmpdir) if 'target_sites' in f]
            assert len(ts_files) == 1

            with open(os.path.join(tmpdir, ts_files[0])) as fh:
                content = fh.read()
            assert 'tsd_in_model=True' in content

    def test_target_site_with_prefix(self):
        """Prefix is included in output filenames."""
        hitTable, model_lengths, hitIndex, paired, config = self._make_test_data()

        with tempfile.TemporaryDirectory() as tmpdir:
            writeTargetSites(
                outDir=tmpdir,
                hitTable=hitTable,
                model_lengths=model_lengths,
                paired=paired,
                hitIndex=hitIndex,
                config=config,
                genome=self._genome(),
                prefix='test_prefix',
                flank_len=5,
                tsd_length=0,
            )

            ts_files = [f for f in os.listdir(tmpdir) if 'target_sites' in f]
            assert len(ts_files) == 1
            assert 'test_prefix' in ts_files[0]

    def test_target_site_with_tsd_length_map(self):
        """TSD length map overrides default tsd_length."""
        hitTable, model_lengths, hitIndex, paired, config = self._make_test_data()

        tsd_map = {'myModel\tmyModel': 2}

        with tempfile.TemporaryDirectory() as tmpdir:
            writeTargetSites(
                outDir=tmpdir,
                hitTable=hitTable,
                model_lengths=model_lengths,
                paired=paired,
                hitIndex=hitIndex,
                config=config,
                genome=self._genome(),
                flank_len=5,
                tsd_length=0,  # default
                tsd_length_map=tsd_map,
            )

            ts_files = [f for f in os.listdir(tmpdir) if 'target_sites' in f]
            assert len(ts_files) == 1

            with open(os.path.join(tmpdir, ts_files[0])) as fh:
                content = fh.read()
            assert 'tsd_len=2' in content


# ---------------------------------------------------------------------------
# Validate CLI parser tests
# ---------------------------------------------------------------------------


class TestValidateParser:
    """Tests for the validate CLI parser configuration."""

    def test_parser_creation(self):
        from tirmite.cli.validate import create_validate_parser

        parser = create_validate_parser()
        assert parser is not None

    def test_required_arguments(self):
        from tirmite.cli.validate import create_validate_parser

        parser = create_validate_parser()
        # Should fail without required args
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_valid_arguments(self):
        from tirmite.cli.validate import create_validate_parser

        parser = create_validate_parser()
        args = parser.parse_args(
            [
                '--target-sites',
                'targets.fa',
                '--blastdb',
                'mydb',
            ]
        )
        assert args.target_sites == 'targets.fa'
        assert args.blastdb == 'mydb'
        assert args.min_coverage == 0.95
        assert args.evalue == 1e-5

    def test_optional_arguments(self):
        from tirmite.cli.validate import create_validate_parser

        parser = create_validate_parser()
        args = parser.parse_args(
            [
                '--target-sites',
                'targets.fa',
                '--blastdb',
                'mydb',
                '--min-coverage',
                '0.8',
                '--tsd-length',
                '5',
                '--tsd-in-model',
                '--prefix',
                'test',
            ]
        )
        assert args.min_coverage == 0.8
        assert args.tsd_length == 5
        assert args.tsd_in_model is True
        assert args.prefix == 'test'


# ---------------------------------------------------------------------------
# Validate helper function tests
# ---------------------------------------------------------------------------


class TestValidateHelpers:
    """Tests for validate module helper functions."""

    def test_parse_blast_results_empty_file(self):
        from tirmite.cli.validate import parse_blast_results

        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            pass  # empty file
        try:
            result = parse_blast_results(f.name)
            assert result == []
        finally:
            os.unlink(f.name)

    def test_parse_blast_results_valid(self):
        from tirmite.cli.validate import parse_blast_results

        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(
                'query1\tsubject1\t95.0\t100\t5\t0\t1\t100\t500\t600\t1e-30\t200\t100\t10000\tplus\n'
            )
        try:
            result = parse_blast_results(f.name)
            assert len(result) == 1
            assert result[0]['qseqid'] == 'query1'
            assert result[0]['sseqid'] == 'subject1'
            assert result[0]['pident'] == 95.0
            assert result[0]['qstart'] == 1
            assert result[0]['qend'] == 100
            assert result[0]['qlen'] == 100
        finally:
            os.unlink(f.name)

    def test_filter_junction_spanning(self):
        from tirmite.cli.validate import filter_junction_spanning

        hits = [
            # Hit that spans midpoint (50) and covers 100% of query
            {
                'qseqid': 'q1',
                'qlen': 100,
                'qstart': 1,
                'qend': 100,
                'sseqid': 's1',
                'sstart': 1,
                'send': 100,
                'evalue': 1e-30,
                'pident': 95.0,
            },
            # Hit that does NOT span midpoint (starts after midpoint)
            {
                'qseqid': 'q1',
                'qlen': 100,
                'qstart': 60,
                'qend': 100,
                'sseqid': 's2',
                'sstart': 60,
                'send': 100,
                'evalue': 1e-10,
                'pident': 90.0,
            },
            # Hit that spans midpoint but low coverage
            {
                'qseqid': 'q1',
                'qlen': 100,
                'qstart': 45,
                'qend': 55,
                'sseqid': 's3',
                'sstart': 45,
                'send': 55,
                'evalue': 1e-5,
                'pident': 80.0,
            },
        ]
        result = filter_junction_spanning(hits, min_coverage=0.95)
        assert 'q1' in result
        assert len(result['q1']) == 1
        assert result['q1'][0]['sseqid'] == 's1'

    def test_check_tsd_gaps_consistent(self):
        from tirmite.cli.validate import check_tsd_gaps

        # No gaps anywhere = consistent
        error, msg = check_tsd_gaps(
            'AAAAAACCCCCCC',
            'AAAAAACCCCCCC',
            tsd_length=3,
            query_len=13,
        )
        assert error == 0
        assert 'consistent' in msg

    def test_check_tsd_gaps_query_has_gaps(self):
        from tirmite.cli.validate import check_tsd_gaps

        # Gaps in query at midpoint = TSD too long
        error, msg = check_tsd_gaps(
            'AAAAAA---CCCCCC',
            'AAAAAATTTCCCCCC',
            tsd_length=3,
            query_len=12,
        )
        assert error > 0
        assert 'too long' in msg

    def test_check_tsd_gaps_target_has_gaps(self):
        from tirmite.cli.validate import check_tsd_gaps

        # Gaps in target at midpoint = TSD too short
        error, msg = check_tsd_gaps(
            'AAAAAATTTCCCCCC',
            'AAAAAA---CCCCCC',
            tsd_length=3,
            query_len=15,
        )
        assert error < 0
        assert 'too short' in msg


# ---------------------------------------------------------------------------
# Pair CLI TSD argument tests
# ---------------------------------------------------------------------------


class TestPairCliTsdArgs:
    """Tests for TSD-related CLI arguments in pair command."""

    def test_parser_has_tsd_args(self):
        from tirmite.cli.hmm_pair import create_pair_parser

        parser = create_pair_parser()
        # Parse with minimum required args plus TSD options
        args = parser.parse_args(
            [
                '--genome',
                'genome.fa',
                '--nhmmer-file',
                'hits.out',
                '--hmm-file',
                'model.hmm',
                '--flank-len',
                '20',
                '--tsd-length',
                '5',
                '--tsd-in-model',
            ]
        )
        assert args.tsd_length == 5
        assert args.tsd_in_model is True
        assert args.tsd_length_map is None

    def test_tsd_length_map_arg(self):
        from tirmite.cli.hmm_pair import create_pair_parser

        parser = create_pair_parser()
        args = parser.parse_args(
            [
                '--genome',
                'genome.fa',
                '--nhmmer-file',
                'hits.out',
                '--hmm-file',
                'model.hmm',
                '--flank-len',
                '20',
                '--tsd-length-map',
                'tsd_map.tsv',
            ]
        )
        assert args.tsd_length_map == 'tsd_map.tsv'

    def test_tsd_requires_flank_len(self):
        import argparse

        from tirmite.cli.hmm_pair import validate_arguments

        # Create mock args with TSD but no flank_len
        args = argparse.Namespace(
            genome='genome.fa',
            blastdb=None,
            left_nhmmer=None,
            right_nhmmer=None,
            left_blast=None,
            right_blast=None,
            left_model=None,
            right_model=None,
            nhmmer_file='hits.out',
            hmm_file='model.hmm',
            blast_file=None,
            query_len=None,
            lengths_file=None,
            pairing_map=None,
            flank_len=None,
            tsd_length=5,
            tsd_length_map=None,
        )
        with pytest.raises(ValueError, match='--flank-len'):
            validate_arguments(args)
