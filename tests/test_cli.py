"""Tests for the top-level `tirmite` CLI dispatcher.

tirmite/cli/cli.py previously had no direct test coverage at all, despite
being the entry point every invocation passes through.
"""

import argparse
import subprocess
import sys

import pytest

from tirmite.cli.cli import _subcommand_parsers, create_parser

SUBCOMMANDS = ['legacy', 'seed', 'pair', 'search', 'validate']


def run_cli(*args, cwd=None):
    """
    Invoke the tirmite CLI in a subprocess.

    Parameters
    ----------
    *args : str
        Arguments to pass after the program name.
    cwd : Path, optional
        Working directory for the invocation.

    Returns
    -------
    subprocess.CompletedProcess
        The completed process, with stdout and stderr captured as text.
    """
    return subprocess.run(
        [sys.executable, '-m', 'tirmite.cli.cli', *args],
        capture_output=True,
        text=True,
        cwd=cwd,
    )


class TestParserConstruction:
    """The parser exposes all five subcommands."""

    def test_all_subcommands_registered(self):
        """Every documented subcommand is reachable."""
        parsers = _subcommand_parsers(create_parser())
        assert sorted(parsers) == sorted(SUBCOMMANDS)

    def test_subcommand_parsers_returns_parsers(self):
        """Each entry is a usable ArgumentParser."""
        for name, parser in _subcommand_parsers(create_parser()).items():
            assert isinstance(parser, argparse.ArgumentParser), name

    def test_empty_for_parser_without_subcommands(self):
        """A parser with no subcommands yields an empty mapping."""
        assert _subcommand_parsers(argparse.ArgumentParser()) == {}


class TestUsageOnNoArguments:
    """Running a command with no arguments prints its own usage."""

    def test_bare_tirmite_prints_top_level_help(self):
        """`tirmite` lists the available subcommands."""
        result = run_cli()

        assert result.returncode == 2
        assert 'Available subcommands' in result.stdout

    @pytest.mark.parametrize('subcommand', SUBCOMMANDS)
    def test_bare_subcommand_exits_two(self, subcommand):
        """`tirmite <subcommand>` is a usage error, not a runtime failure."""
        result = run_cli(subcommand)

        assert result.returncode == 2

    @pytest.mark.parametrize(
        'subcommand,marker',
        [
            ('search', '--cluster-map'),
            ('seed', '--model-name'),
            ('pair', '--genome'),
            ('validate', '--target-sites'),
            ('legacy', '--hmm-file'),
        ],
    )
    def test_bare_subcommand_prints_its_own_help(self, subcommand, marker):
        """The help shown is the subcommand's, not the top-level help.

        Previously only a bare `tirmite` was special-cased, so `tirmite search`
        parsed successfully and failed later inside its main().
        """
        result = run_cli(subcommand)
        combined = result.stdout + result.stderr

        assert marker in combined, f'{subcommand} help did not mention {marker}'
        assert 'Available subcommands' not in result.stdout

    @pytest.mark.parametrize('subcommand', SUBCOMMANDS)
    def test_bare_subcommand_creates_nothing(self, subcommand, tmp_path):
        """A usage error must not touch the filesystem."""
        run_cli(subcommand, cwd=tmp_path)

        assert list(tmp_path.iterdir()) == []


class TestNoOutputDirectorySideEffect:
    """Invalid arguments must not leave an output directory behind."""

    def test_missing_inputs_creates_no_outdir(self, tmp_path):
        """`--outdir` is not created when validation fails.

        The output directory used to be created, and logging initialised,
        before any input validation ran -- so an invocation that could never
        succeed still left an empty directory behind.
        """
        result = run_cli(
            'search', '--outdir', 'RESULTS', '--max-evalue', '0.1', cwd=tmp_path
        )

        assert result.returncode == 2
        assert not (tmp_path / 'RESULTS').exists()

    def test_error_message_goes_to_stderr(self, tmp_path):
        """Usage errors are reported on stderr, since logging is not up yet."""
        result = run_cli('search', '--outdir', 'RESULTS', cwd=tmp_path)

        assert 'Must provide either search inputs' in result.stderr

    def test_split_paired_without_pairing_map_creates_no_outdir(self, tmp_path):
        """A dependent-option failure is caught before the mkdir too."""
        result = run_cli(
            'search',
            '--blast-results',
            'hits.tab',
            '--split-paired-output',
            '--outdir',
            'RESULTS',
            cwd=tmp_path,
        )

        assert result.returncode == 2
        assert not (tmp_path / 'RESULTS').exists()
        assert '--split-paired-output requires --pairing-map' in result.stderr


class TestVersionAndUnknownCommands:
    """Standard argparse behaviours still hold."""

    def test_version_flag(self):
        """`--version` reports the package version and exits cleanly."""
        result = run_cli('--version')

        assert result.returncode == 0
        assert 'tirmite' in result.stdout

    def test_unknown_subcommand_is_a_usage_error(self):
        """An unrecognised subcommand exits 2 with an argparse message."""
        result = run_cli('definitely-not-a-command')

        assert result.returncode == 2
        assert 'invalid choice' in result.stderr

    @pytest.mark.parametrize('subcommand', SUBCOMMANDS)
    def test_help_flag_exits_zero(self, subcommand):
        """Explicitly asking for help is success, unlike a bare invocation."""
        result = run_cli(subcommand, '--help')

        assert result.returncode == 0


class TestDispatch:
    """Each subcommand routes to its own main()."""

    @pytest.mark.parametrize(
        'subcommand,module_name',
        [
            ('legacy', 'tirmite.cli.legacy'),
            ('seed', 'tirmite.cli.hmm_build'),
            ('pair', 'tirmite.cli.hmm_pair'),
            ('search', 'tirmite.cli.ensemble_search'),
            ('validate', 'tirmite.cli.validate'),
        ],
    )
    def test_dispatches_to_expected_main(self, subcommand, module_name, monkeypatch):
        """The dispatcher calls the main() belonging to the named subcommand."""
        import importlib

        from tirmite.cli import cli as cli_module

        module = importlib.import_module(module_name)
        called = {}

        def fake_main(args):
            called['command'] = args.command
            return 0

        monkeypatch.setattr(module, 'main', fake_main)
        # Enough argv entries that the bare-subcommand help guard does not fire.
        monkeypatch.setattr(sys, 'argv', ['tirmite', subcommand, '--loglevel', 'ERROR'])
        monkeypatch.setattr(
            cli_module.argparse.ArgumentParser,
            'parse_args',
            lambda self, *a, **k: argparse.Namespace(command=subcommand),
        )

        assert cli_module.main() == 0
        assert called['command'] == subcommand
