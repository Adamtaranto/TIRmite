"""Tests for HMM construction from an alignment, and the seed comparison report.

Both areas covered here were completely untested, and both were broken:

* `build_hmm_from_alignment_pyhmmer` validated its input with
  `SeqIO.parse(file, 'stockholm')` while every seed path feeds it the FASTA
  that MAFFT writes, so `tirmite seed` failed on every run that reached HMM
  building.
* The seed comparison report referenced a loop variable that had been renamed,
  raising `NameError` for any asymmetric run that found a seed similarity.
"""

from collections import namedtuple

import pytest

from tirmite.cli.hmm_build import (
    HMMBuildError,
    build_hmm_from_alignment_pyhmmer,
    read_alignment_records,
    write_seed_comparison_report,
)

# MAFFT writes FASTA; hmmalign writes Stockholm. Both reach the same builder.
FASTA_ALIGNMENT = (
    '>s1\nACGTACGTACGTACGT\n>s2\nACGTACGTACGTACGA\n>s3\nACGTACGTACGTTCGT\n'
)
STOCKHOLM_ALIGNMENT = (
    '# STOCKHOLM 1.0\n'
    's1 ACGTACGTACGTACGT\n'
    's2 ACGTACGTACGTACGA\n'
    's3 ACGTACGTACGTTCGT\n'
    '//\n'
)


class TestReadAlignmentRecords:
    """Alignment reading must not assume a single format."""

    def test_reads_fasta(self, tmp_path):
        """FASTA is what the seed workflow produces, via MAFFT."""
        path = tmp_path / 'aln.fasta'
        path.write_text(FASTA_ALIGNMENT)

        records = read_alignment_records(path)

        assert len(records) == 3
        assert [r.id for r in records] == ['s1', 's2', 's3']

    def test_reads_stockholm(self, tmp_path):
        """Stockholm is what the --update workflow produces, via hmmalign."""
        path = tmp_path / 'aln.sto'
        path.write_text(STOCKHOLM_ALIGNMENT)

        records = read_alignment_records(path)

        assert len(records) == 3

    def test_stockholm_is_not_misread_as_fasta(self, tmp_path):
        """Stockholm is tried first, so its header can never be misparsed."""
        path = tmp_path / 'aln.sto'
        path.write_text(STOCKHOLM_ALIGNMENT)

        records = read_alignment_records(path)

        assert all(not r.id.startswith('#') for r in records)

    def test_unparseable_file_returns_empty(self, tmp_path):
        """A file in no known format yields no records rather than raising."""
        path = tmp_path / 'junk.txt'
        path.write_text('this is not an alignment in any format\n')

        assert read_alignment_records(path) == []


class TestBuildHmmFromAlignment:
    """The builder must accept every format TIRmite actually produces."""

    def test_builds_from_mafft_fasta(self, tmp_path):
        """The regression test for the broken seed workflow.

        Before the fix this raised
        `HMMBuildError: ... Did not find STOCKHOLM header`, so every
        `tirmite seed` run failed at the final step.
        """
        alignment = tmp_path / 'aln.fasta'
        alignment.write_text(FASTA_ALIGNMENT)

        hmm_file = build_hmm_from_alignment_pyhmmer(alignment, 'MyTIR', tmp_path)

        assert hmm_file.exists()
        text = hmm_file.read_text()
        assert text.startswith('HMMER3/f')
        assert 'NAME  MyTIR' in text

    def test_builds_from_hmmalign_stockholm(self, tmp_path):
        """The --update path, which was already working, must keep working."""
        alignment = tmp_path / 'aln.sto'
        alignment.write_text(STOCKHOLM_ALIGNMENT)

        hmm_file = build_hmm_from_alignment_pyhmmer(alignment, 'StoTIR', tmp_path)

        assert hmm_file.exists()
        assert 'NAME  StoTIR' in hmm_file.read_text()

    def test_model_length_matches_alignment_width(self, tmp_path):
        """A sanity check that the built model describes the input."""
        alignment = tmp_path / 'aln.fasta'
        alignment.write_text(FASTA_ALIGNMENT)

        hmm_file = build_hmm_from_alignment_pyhmmer(alignment, 'LenTIR', tmp_path)

        assert 'LENG  16' in hmm_file.read_text()

    def test_single_sequence_is_rejected(self, tmp_path):
        """One sequence is not an alignment."""
        alignment = tmp_path / 'aln.fasta'
        alignment.write_text('>only\nACGTACGT\n')

        with pytest.raises(HMMBuildError, match='at least 2 sequences'):
            build_hmm_from_alignment_pyhmmer(alignment, 'OneTIR', tmp_path)

    def test_empty_alignment_is_rejected(self, tmp_path):
        """An empty file is reported, not silently turned into an empty model."""
        alignment = tmp_path / 'aln.fasta'
        alignment.write_text('')

        with pytest.raises(HMMBuildError, match='No sequences found'):
            build_hmm_from_alignment_pyhmmer(alignment, 'EmptyTIR', tmp_path)


# Minimal stand-ins for the (BlastHit, alignment) pairs compare_seeds returns.
_Hit = namedtuple(
    'Hit',
    'length identity query_coverage subject_coverage query_id query_start '
    'query_end subject_id subject_start subject_end strand',
)


class _FakeAlignment:
    """Stand-in for a Bio.Align alignment: has .score and a str form."""

    score = 42.5

    def __str__(self):
        """Return a two-line pseudo-alignment."""
        return 'ACGT\nACGA'


class TestWriteSeedComparisonReport:
    """The report writer, extracted from main so it can be tested at all."""

    def _comparisons(self, n=2):
        """Build n (hit, alignment) pairs."""
        return [
            (
                _Hit(
                    length=100,
                    identity=95.5,
                    query_coverage=0.9,
                    subject_coverage=0.8,
                    query_id='LeftSeed',
                    query_start=1,
                    query_end=100,
                    subject_id='RightSeed',
                    subject_start=5,
                    subject_end=104,
                    strand='+',
                ),
                _FakeAlignment(),
            )
            for _ in range(n)
        ]

    def test_writes_without_error(self, tmp_path):
        """The regression test for the NameError.

        Three stale references to a renamed loop variable survived here, so
        this block raised for every asymmetric run that found a similarity.
        """
        out = tmp_path / 'report.txt'

        write_seed_comparison_report(
            self._comparisons(),
            out,
            model_name='MyTIR',
            left_seed_name='left.fa',
            right_seed_name='right.fa',
        )

        assert out.exists()

    def test_report_contains_alignment_score(self, tmp_path):
        """The alignment score is what the stale references were reading."""
        out = tmp_path / 'report.txt'

        write_seed_comparison_report(
            self._comparisons(1),
            out,
            model_name='MyTIR',
            left_seed_name='left.fa',
            right_seed_name='right.fa',
        )

        assert 'Alignment score: 42.5' in out.read_text()

    def test_report_lists_every_similarity(self, tmp_path):
        """All comparisons are written, not just the top few."""
        out = tmp_path / 'report.txt'

        write_seed_comparison_report(
            self._comparisons(5),
            out,
            model_name='MyTIR',
            left_seed_name='left.fa',
            right_seed_name='right.fa',
        )

        text = out.read_text()
        assert 'Total similarities found: 5' in text
        assert 'Similarity 5:' in text

    def test_report_includes_seed_filenames(self, tmp_path):
        """The header records which seeds were compared."""
        out = tmp_path / 'report.txt'

        write_seed_comparison_report(
            self._comparisons(1),
            out,
            model_name='MyTIR',
            left_seed_name='my_left.fa',
            right_seed_name='my_right.fa',
        )

        text = out.read_text()
        assert 'Left seed: my_left.fa' in text
        assert 'Right seed: my_right.fa' in text
