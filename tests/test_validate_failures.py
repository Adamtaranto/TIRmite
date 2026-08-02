"""Tests that `tirmite validate` fails loudly instead of reporting success.

Every case here previously produced a complete `validation_summary.tsv` with
`predicted_tsd_error 0.0` on every row and exit code 0, while validating
nothing at all:

* MAFFT not installed
* blastdbcmd not installed, or a database built without ``-parse_seqids``
* a standard 12-column ``-outfmt 6`` file passed to ``--blast-results``
* a mistyped ``--blast-results`` path
* every comparison inconclusive because both sequences carry gaps

That last one is the same failure mode fixed in 1.5.0 for
``tirmite.core.tsd.compare_tsds``, which returns ``Optional[int]`` precisely so
that an unverifiable result cannot be averaged in as agreement.
"""

import argparse

import pytest

from tirmite.cli.validate import (
    ValidationError,
    _check_required_tools,
    check_tsd_gaps,
    parse_blast_results,
)

# The 15-column extended format TIRmite requires.
EXTENDED_ROW = (
    'site1\tchr1\t99.0\t100\t1\t0\t1\t100\t500\t599\t1e-40\t200\t100\t100000\tplus\n'
)
# The standard 12-column format, which is what users reach for by default.
STANDARD_ROW = 'site1\tchr1\t99.0\t100\t1\t0\t1\t100\t500\t599\t1e-40\t200\n'


def _args(**overrides):
    """Build a Namespace with the fields _check_required_tools reads."""
    defaults = {'blast_results': None}
    defaults.update(overrides)
    return argparse.Namespace(**defaults)


class TestParseBlastResults:
    """Unusable BLAST input must not look like an empty result."""

    def test_parses_extended_format(self, tmp_path):
        """The 15-column format is accepted."""
        path = tmp_path / 'hits.tsv'
        path.write_text(EXTENDED_ROW)

        hits = parse_blast_results(str(path))

        assert len(hits) == 1
        assert hits[0]['qlen'] == 100
        assert hits[0]['sstrand'] == 'plus'

    def test_standard_12_column_format_is_rejected(self, tmp_path):
        """A 12-column file previously produced zero hits and exit 0.

        This is the common case: the extended format is documented only in
        the --blast-results help text.
        """
        path = tmp_path / 'hits.tsv'
        path.write_text(STANDARD_ROW * 3)

        with pytest.raises(ValidationError, match='could be parsed'):
            parse_blast_results(str(path))

    def test_rejection_message_names_the_required_format(self, tmp_path):
        """The error tells the user exactly what to run."""
        path = tmp_path / 'hits.tsv'
        path.write_text(STANDARD_ROW)

        with pytest.raises(ValidationError, match='qlen slen sstrand'):
            parse_blast_results(str(path))

    def test_missing_file_raises(self, tmp_path):
        """A mistyped path was previously a warning, then exit 0."""
        with pytest.raises(ValidationError, match='not found'):
            parse_blast_results(str(tmp_path / 'does-not-exist.tsv'))

    def test_genuinely_empty_file_is_not_an_error(self, tmp_path):
        """Zero hits is a legitimate result and must stay distinguishable."""
        path = tmp_path / 'hits.tsv'
        path.write_text('')

        assert parse_blast_results(str(path)) == []

    def test_comment_lines_are_ignored(self, tmp_path):
        """BLAST writes '#' headers with some output options."""
        path = tmp_path / 'hits.tsv'
        path.write_text('# BLASTN 2.12.0+\n# Query: site1\n' + EXTENDED_ROW)

        assert len(parse_blast_results(str(path))) == 1

    def test_mixed_good_and_short_rows_warns_but_parses(self, tmp_path, caplog):
        """Some usable rows means the run can proceed, with a warning."""
        import logging

        path = tmp_path / 'hits.tsv'
        path.write_text(EXTENDED_ROW + STANDARD_ROW)

        with caplog.at_level(logging.WARNING):
            hits = parse_blast_results(str(path))

        assert len(hits) == 1
        assert 'Skipped 1 short' in caplog.text

    def test_malformed_numeric_value_is_skipped_not_crashed(self, tmp_path):
        """A bad number previously escaped as an unhandled ValueError."""
        path = tmp_path / 'hits.tsv'
        bad = EXTENDED_ROW.replace('\t100\t100000\tplus', '\tNOTANUMBER\t100000\tplus')
        path.write_text(EXTENDED_ROW + bad)

        hits = parse_blast_results(str(path))

        assert len(hits) == 1


class TestRequiredTools:
    """Missing external tools must abort, not silently produce empty results."""

    def test_missing_mafft_raises(self, monkeypatch):
        """align_in_memory returns None per alignment; nothing aggregated it."""
        import tirmite.cli.validate as mod

        monkeypatch.setattr(mod, 'mafft_available', lambda: False)
        monkeypatch.setattr(mod, 'blastdbcmd_available', lambda: True)
        monkeypatch.setattr(mod.shutil, 'which', lambda _n: '/usr/bin/blastn')

        with pytest.raises(ValidationError, match='mafft'):
            _check_required_tools(_args())

    def test_missing_blastdbcmd_raises(self, monkeypatch):
        """extract_hit_sequence returns None per hit; nothing aggregated it."""
        import tirmite.cli.validate as mod

        monkeypatch.setattr(mod, 'mafft_available', lambda: True)
        monkeypatch.setattr(mod, 'blastdbcmd_available', lambda: False)
        monkeypatch.setattr(mod.shutil, 'which', lambda _n: '/usr/bin/blastn')

        with pytest.raises(ValidationError, match='blastdbcmd'):
            _check_required_tools(_args())

    def test_blastn_not_required_with_precomputed_results(self, monkeypatch):
        """--blast-results means we never invoke blastn ourselves."""
        import tirmite.cli.validate as mod

        monkeypatch.setattr(mod, 'mafft_available', lambda: True)
        monkeypatch.setattr(mod, 'blastdbcmd_available', lambda: True)
        monkeypatch.setattr(mod.shutil, 'which', lambda _n: None)

        _check_required_tools(_args(blast_results='hits.tsv'))

    def test_blastn_required_without_precomputed_results(self, monkeypatch):
        """Without them, blastn must be present."""
        import tirmite.cli.validate as mod

        monkeypatch.setattr(mod, 'mafft_available', lambda: True)
        monkeypatch.setattr(mod, 'blastdbcmd_available', lambda: True)
        monkeypatch.setattr(mod.shutil, 'which', lambda _n: None)

        with pytest.raises(ValidationError, match='blastn'):
            _check_required_tools(_args())

    def test_all_missing_are_listed_together(self, monkeypatch):
        """One error naming every missing tool, not one failure at a time."""
        import tirmite.cli.validate as mod

        monkeypatch.setattr(mod, 'mafft_available', lambda: False)
        monkeypatch.setattr(mod, 'blastdbcmd_available', lambda: False)
        monkeypatch.setattr(mod.shutil, 'which', lambda _n: None)

        with pytest.raises(ValidationError) as excinfo:
            _check_required_tools(_args())

        message = str(excinfo.value)
        assert 'mafft' in message
        assert 'blastdbcmd' in message
        assert 'blastn' in message


class TestCheckTsdGapsInconclusive:
    """Inconclusive must be distinguishable from confirmed."""

    def test_confirmed_returns_zero(self):
        """No gaps on either sequence is genuine agreement."""
        error, message = check_tsd_gaps('ACGTACGTAC', 'ACGTACGTAC', 0, 10)

        assert error == 0
        assert 'consistent' in message

    def test_query_gaps_report_too_long(self):
        """Gaps in the query mean the declared TSD was too long."""
        error, _ = check_tsd_gaps('ACGT--GTAC', 'ACGTACGTAC', 0, 10)

        assert error is not None
        assert error > 0

    def test_target_gaps_report_too_short(self):
        """Gaps in the target mean the declared TSD was too short."""
        error, _ = check_tsd_gaps('ACGTACGTAC', 'ACGT--GTAC', 0, 10)

        assert error is not None
        assert error < 0

    def test_gaps_on_both_sides_return_none(self):
        """The regression test: this used to return 0, i.e. "consistent".

        Gaps on both sequences carry no directional information. Returning 0
        meant an unverifiable comparison was averaged in as agreement and
        reported as "TSD length appears consistent".
        """
        error, message = check_tsd_gaps('AC-TACGTAC', 'ACGTAC-TAC', 0, 10)

        assert error is None
        assert 'inconclusive' in message

    def test_none_is_not_zero(self):
        """Guard against the two collapsing again.

        `if error:` would treat both as falsey; the aggregation must test
        `is None` explicitly.
        """
        inconclusive, _ = check_tsd_gaps('AC-TACGTAC', 'ACGTAC-TAC', 0, 10)
        confirmed, _ = check_tsd_gaps('ACGTACGTAC', 'ACGTACGTAC', 0, 10)

        assert inconclusive is None
        assert confirmed == 0
        assert inconclusive is not confirmed
