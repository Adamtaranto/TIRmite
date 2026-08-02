"""Tests for the shared argparse type validators.

These were previously duplicated across hmm_build.py and ensemble_search.py
and had no direct test coverage at all; a bad bound would only have surfaced
as a confusing CLI error.
"""

import argparse

import pytest

from tirmite.cli._argtypes import (
    validate_coverage,
    validate_evalue,
    validate_identity,
    validate_score_ratio,
    validate_threads,
    validate_word_size,
)


class TestValidateEvalue:
    """E-values must be strictly positive."""

    @pytest.mark.parametrize(
        'raw,expected',
        [
            ('1e-3', 1e-3),
            ('0.001', 0.001),
            ('1', 1.0),
            ('1e-300', 1e-300),
            ('10', 10.0),  # Values above 1 are unusual but legal.
        ],
    )
    def test_accepts_positive(self, raw, expected):
        """Any strictly positive number is accepted."""
        assert validate_evalue(raw) == expected

    @pytest.mark.parametrize('raw', ['0', '0.0', '-1', '-1e-5'])
    def test_rejects_non_positive(self, raw):
        """Zero and negatives are rejected: they would reject every hit."""
        with pytest.raises(argparse.ArgumentTypeError, match='must be positive'):
            validate_evalue(raw)

    @pytest.mark.parametrize('raw', ['abc', '', '1e', 'None'])
    def test_rejects_non_numeric(self, raw):
        """Non-numeric input reports a type error, not a traceback."""
        with pytest.raises(argparse.ArgumentTypeError, match='must be a number'):
            validate_evalue(raw)


class TestValidateIdentity:
    """Identity is a percentage in [0, 100]."""

    @pytest.mark.parametrize('raw,expected', [('0', 0.0), ('60', 60.0), ('100', 100.0)])
    def test_accepts_in_range(self, raw, expected):
        """Both bounds are inclusive."""
        assert validate_identity(raw) == expected

    @pytest.mark.parametrize('raw', ['-0.1', '100.1', '1000'])
    def test_rejects_out_of_range(self, raw):
        """Anything outside 0-100 is rejected."""
        with pytest.raises(argparse.ArgumentTypeError, match='between 0 and 100'):
            validate_identity(raw)

    def test_rejects_non_numeric(self):
        """Non-numeric input reports a type error."""
        with pytest.raises(argparse.ArgumentTypeError, match='must be a number'):
            validate_identity('high')


class TestValidateCoverage:
    """Coverage is a fraction in [0.0, 1.0], not a percentage."""

    @pytest.mark.parametrize('raw,expected', [('0', 0.0), ('0.5', 0.5), ('1', 1.0)])
    def test_accepts_in_range(self, raw, expected):
        """Both bounds are inclusive."""
        assert validate_coverage(raw) == expected

    @pytest.mark.parametrize('raw', ['-0.1', '1.1', '80'])
    def test_rejects_out_of_range(self, raw):
        """80 is a percentage, and is correctly rejected as a fraction."""
        with pytest.raises(argparse.ArgumentTypeError, match='between 0.0 and 1.0'):
            validate_coverage(raw)


class TestValidateScoreRatio:
    """Score ratios are better/weaker, so they cannot be below 1.0."""

    @pytest.mark.parametrize('raw,expected', [('1', 1.0), ('1.5', 1.5), ('10', 10.0)])
    def test_accepts_at_least_one(self, raw, expected):
        """1.0 means 'any difference is decisive' and is the lower bound."""
        assert validate_score_ratio(raw) == expected

    @pytest.mark.parametrize('raw', ['0.9', '0', '-2'])
    def test_rejects_below_one(self, raw):
        """A ratio below 1.0 would discard the better of two hits."""
        with pytest.raises(argparse.ArgumentTypeError, match='at least 1.0'):
            validate_score_ratio(raw)

    def test_rejects_non_numeric(self):
        """Non-numeric input reports a type error."""
        with pytest.raises(argparse.ArgumentTypeError, match='must be a number'):
            validate_score_ratio('strict')


class TestValidateThreads:
    """Thread counts must be integers of at least 1."""

    @pytest.mark.parametrize('raw,expected', [('1', 1), ('8', 8), ('128', 128)])
    def test_accepts_positive_int(self, raw, expected):
        """Any integer >= 1 is accepted."""
        assert validate_threads(raw) == expected

    @pytest.mark.parametrize('raw', ['0', '-1'])
    def test_rejects_below_one(self, raw):
        """Zero threads would mean no work is done."""
        with pytest.raises(argparse.ArgumentTypeError, match='at least 1'):
            validate_threads(raw)

    @pytest.mark.parametrize('raw', ['1.5', 'many', ''])
    def test_rejects_non_integer(self, raw):
        """Floats are rejected rather than silently truncated."""
        with pytest.raises(argparse.ArgumentTypeError, match='must be an integer'):
            validate_threads(raw)


class TestValidateWordSize:
    """BLAST refuses nucleotide word sizes below 4."""

    @pytest.mark.parametrize('raw,expected', [('4', 4), ('11', 11), ('28', 28)])
    def test_accepts_at_least_four(self, raw, expected):
        """4 is the inclusive lower bound imposed by blastn."""
        assert validate_word_size(raw) == expected

    @pytest.mark.parametrize('raw', ['3', '0', '-4'])
    def test_rejects_below_four(self, raw):
        """A word size blastn would reject is caught at argument-parse time."""
        with pytest.raises(argparse.ArgumentTypeError, match='at least 4'):
            validate_word_size(raw)

    def test_rejects_non_integer(self):
        """Non-integer input reports a type error."""
        with pytest.raises(argparse.ArgumentTypeError, match='must be an integer'):
            validate_word_size('4.5')


class TestValidatorsAreSharedAcrossSubcommands:
    """The seed and search subcommands must use the same validator objects."""

    def test_search_and_seed_share_implementations(self):
        """A single definition backs both subcommands, so messages agree."""
        from tirmite.cli import _argtypes, ensemble_search, hmm_build

        assert ensemble_search.validate_evalue is _argtypes.validate_evalue
        assert hmm_build.validate_evalue is _argtypes.validate_evalue
        assert ensemble_search.validate_identity is hmm_build.validate_identity
        assert ensemble_search.validate_threads is hmm_build.validate_threads
