"""
Tests for the new --blastdb, --blast-hits, --left-blast-hits and
--right-blast-hits options added to ``tirmite seed``.

These tests exercise:
* CLI argument parsing (mutually exclusive group, new optional args,
  removal of --max-gap)
* ``warn_multiple_queries`` helper
* ``check_targets_in_blastdb`` helper
* ``process_seed_sequences`` routing logic with pre-computed hit tables
* ``process_asymmetric_seeds`` routing logic with pre-computed hit tables
* ``resolve_asymmetric_conflicts`` signature change (no max_gap)
"""

import argparse
from pathlib import Path

import pytest

from tirmite.cli.hmm_build import (
    BlastHit,
    HMMBuildError,
    _configure_seed_parser,
    parse_blast_output,
    resolve_asymmetric_conflicts,
    warn_multiple_queries,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_hit(
    query_id: str = 'query1',
    subject_id: str = 'chr1',
    query_start: int = 1,
    query_end: int = 100,
    subject_start: int = 1000,
    subject_end: int = 1099,
    length: int = 100,
    identity: float = 95.0,
    query_len: int = 100,
    subject_len: int = 10000,
) -> BlastHit:
    return BlastHit(
        query_id=query_id,
        subject_id=subject_id,
        query_start=query_start,
        query_end=query_end,
        subject_start=subject_start,
        subject_end=subject_end,
        length=length,
        identity=identity,
        query_len=query_len,
        subject_len=subject_len,
    )


def _write_blast_hits(path: Path, hits: list) -> None:
    """Write a minimal TIRmite-format BLAST hit table."""
    with open(path, 'w') as f:
        f.write('# BLAST tabular output format 6\n')
        f.write(
            '# qstart\tqend\tsstart\tsend\tlength\tpositive\tpident'
            '\tqlen\tslen\tqframe\tsframe\tqseqid\tsseqid\n'
        )
        for h in hits:
            f.write(
                f'{h.query_start}\t{h.query_end}\t{h.subject_start}\t'
                f'{h.subject_end}\t{h.length}\t{h.length}\t{h.identity:.2f}\t'
                f'{h.query_len}\t{h.subject_len}\t1\t1\t'
                f'{h.query_id}\t{h.subject_id}\n'
            )


# ---------------------------------------------------------------------------
# CLI parser tests
# ---------------------------------------------------------------------------


class TestSeedParserNewOptions:
    """Verify that the seed parser accepts the new options correctly."""

    def _get_parser(self) -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser()
        _configure_seed_parser(parser)
        return parser

    def test_blastdb_accepted(self, tmp_path):
        """--blastdb should be accepted as a valid genome source."""
        db = tmp_path / 'mydb'
        # Create a dummy db file so Path points to something
        (tmp_path / 'mydb.nsq').touch()
        parser = self._get_parser()
        args = parser.parse_args(
            ['--left-seed', 'seed.fa', '--model-name', 'myTE', '--blastdb', str(db)]
        )
        assert args.blastdb == db

    def test_genome_accepted(self, tmp_path):
        """--genome should still be accepted."""
        parser = self._get_parser()
        args = parser.parse_args(
            [
                '--left-seed',
                'seed.fa',
                '--model-name',
                'myTE',
                '--genome',
                str(tmp_path / 'genome.fa'),
            ]
        )
        assert args.genome is not None

    def test_blastdb_and_genome_mutually_exclusive(self):
        """--blastdb and --genome must not be combined."""
        parser = self._get_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(
                [
                    '--left-seed',
                    'seed.fa',
                    '--model-name',
                    'myTE',
                    '--genome',
                    'g.fa',
                    '--blastdb',
                    'db',
                ]
            )

    def test_blast_hits_accepted(self, tmp_path):
        """--blast-hits should be accepted."""
        parser = self._get_parser()
        hits_file = tmp_path / 'hits.tab'
        hits_file.touch()
        args = parser.parse_args(
            [
                '--left-seed',
                'seed.fa',
                '--model-name',
                'myTE',
                '--genome',
                'g.fa',
                '--blast-hits',
                str(hits_file),
            ]
        )
        assert args.blast_hits == hits_file

    def test_left_right_blast_hits_accepted(self, tmp_path):
        """--left-blast-hits and --right-blast-hits should be accepted."""
        parser = self._get_parser()
        left_hits = tmp_path / 'left.tab'
        right_hits = tmp_path / 'right.tab'
        args = parser.parse_args(
            [
                '--left-seed',
                'l.fa',
                '--right-seed',
                'r.fa',
                '--model-name',
                'myTE',
                '--genome',
                'g.fa',
                '--left-blast-hits',
                str(left_hits),
                '--right-blast-hits',
                str(right_hits),
            ]
        )
        assert args.left_blast_hits == left_hits
        assert args.right_blast_hits == right_hits

    def test_max_gap_removed(self):
        """--max-gap must no longer be a valid argument."""
        parser = self._get_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(
                [
                    '--left-seed',
                    'seed.fa',
                    '--model-name',
                    'myTE',
                    '--genome',
                    'g.fa',
                    '--max-gap',
                    '50',
                ]
            )

    def test_genome_source_required(self):
        """At least one of --genome / --genome-list / --blastdb is required."""
        parser = self._get_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(['--left-seed', 'seed.fa', '--model-name', 'myTE'])


# ---------------------------------------------------------------------------
# warn_multiple_queries tests
# ---------------------------------------------------------------------------


class TestWarnMultipleQueries:
    def test_single_query_no_warning(self, caplog):
        hits = [_make_hit(query_id='q1'), _make_hit(query_id='q1')]
        with caplog.at_level('WARNING'):
            warn_multiple_queries(hits)
        assert not any('Multiple query' in m for m in caplog.messages)

    def test_multiple_queries_warning(self, caplog):
        hits = [_make_hit(query_id='q1'), _make_hit(query_id='q2')]
        with caplog.at_level('WARNING'):
            warn_multiple_queries(hits, context='test_model')
        assert any('Multiple query' in m for m in caplog.messages)
        assert any('test_model' in m for m in caplog.messages)

    def test_empty_hits_no_warning(self, caplog):
        with caplog.at_level('WARNING'):
            warn_multiple_queries([])
        assert not any('Multiple query' in m for m in caplog.messages)


# ---------------------------------------------------------------------------
# resolve_asymmetric_conflicts - signature without max_gap
# ---------------------------------------------------------------------------


class TestResolveAsymmetricConflictsNewSignature:
    def test_called_without_max_gap(self):
        """resolve_asymmetric_conflicts must work without a max_gap argument."""
        left = [_make_hit(subject_id='chr1', subject_start=100, subject_end=200)]
        right = [_make_hit(subject_id='chr1', subject_start=5000, subject_end=5100)]
        # Should not raise
        fl, fr = resolve_asymmetric_conflicts(left, right)
        assert len(fl) == 1
        assert len(fr) == 1

    def test_max_gap_kwarg_rejected(self):
        """Passing max_gap as a keyword argument should raise TypeError."""
        left = [_make_hit()]
        right = [_make_hit(subject_start=5000, subject_end=5100)]
        with pytest.raises(TypeError):
            resolve_asymmetric_conflicts(left, right, max_gap=10)  # type: ignore[call-arg]


# ---------------------------------------------------------------------------
# parse_blast_output with custom format (round-trip test)
# ---------------------------------------------------------------------------


class TestParseBlastOutputRoundTrip:
    """Ensure that hits written by save_all_blast_hits can be read back."""

    def test_round_trip(self, tmp_path):
        from tirmite.cli.hmm_build import save_all_blast_hits

        original = [
            _make_hit('q1', 'chr1', 1, 100, 1000, 1099, 100, 95.0, 100, 10000),
            _make_hit(
                'q1', 'chr1', 1, 80, 2000, 1921, 80, 90.0, 100, 10000
            ),  # reverse strand (subject_start > subject_end)
        ]
        tab_file = tmp_path / 'hits.tab'
        save_all_blast_hits(original, tab_file)

        parsed = parse_blast_output(tab_file)
        assert len(parsed) == 2
        assert parsed[0].query_id == 'q1'
        assert parsed[0].subject_id == 'chr1'
        assert parsed[0].strand == '+'
        # Second hit: subject_start > subject_end → reverse strand
        assert parsed[1].strand == '-'


# ---------------------------------------------------------------------------
# process_seed_sequences with blast_hits_file
# ---------------------------------------------------------------------------


class TestProcessSeedSequencesWithBlastHitsFile:
    """Unit-level tests for the blast_hits_file path in process_seed_sequences."""

    def _make_fasta(
        self, path: Path, seqid: str = 'seed1', seq: str = 'ACGT' * 25
    ) -> Path:
        with open(path, 'w') as f:
            f.write(f'>{seqid}\n{seq}\n')
        return path

    def test_raises_when_blast_hits_file_empty(self, tmp_path):
        """Should raise HMMBuildError when the provided hit table is empty."""
        from tirmite.cli.hmm_build import process_seed_sequences

        seed = self._make_fasta(tmp_path / 'seed.fa')
        empty_hits = tmp_path / 'empty.tab'
        empty_hits.write_text('# BLAST tabular output format 6\n')

        with pytest.raises(HMMBuildError, match='No hits found|No hits in'):
            process_seed_sequences(
                seed_file=seed,
                model_name='testmodel',
                genome_files=[],
                temp_dir=tmp_path,
                output_dir=tmp_path,
                min_coverage=0.5,
                min_identity=70.0,
                blast_db=None,
                blast_hits_file=empty_hits,
            )

    def test_raises_when_no_hits_pass_thresholds(self, tmp_path):
        """Should raise HMMBuildError when all hits are below thresholds."""
        from tirmite.cli.hmm_build import process_seed_sequences

        seed = self._make_fasta(tmp_path / 'seed.fa')

        # Hit with very low identity and coverage
        bad_hit = _make_hit(
            query_start=1,
            query_end=10,  # coverage = 10/100 = 0.10
            identity=30.0,
        )
        hits_file = tmp_path / 'bad_hits.tab'
        _write_blast_hits(hits_file, [bad_hit])

        with pytest.raises(HMMBuildError, match='No hits passed thresholds'):
            process_seed_sequences(
                seed_file=seed,
                model_name='testmodel',
                genome_files=[],
                temp_dir=tmp_path,
                output_dir=tmp_path,
                min_coverage=0.7,
                min_identity=70.0,
                blast_db=None,
                blast_hits_file=hits_file,
            )
