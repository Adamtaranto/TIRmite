"""Tests for TIRmite's logging configuration.

TIRmite previously logged through the bare ``logging.info(...)`` module
functions -- i.e. the root logger -- from all ~650 call sites, and
``init_logging`` cleared and reconfigured root. That made per-module verbosity
control impossible and meant that importing TIRmite and calling any ``main()``
silently hijacked the host application's logging.

These tests pin the replacement: named per-module loggers under ``tirmite``,
a NullHandler on library modules, and configuration confined to the package
logger.
"""

import importlib
import logging
import pkgutil
import subprocess
import sys

import pytest

from tirmite.utils.logs import PACKAGE_LOGGER_NAME, init_logging

# Every module that emits log records.
LOGGING_MODULES = [
    'tirmite.cli.cli',
    'tirmite.cli.ensemble_search',
    'tirmite.cli.hmm_build',
    'tirmite.cli.hmm_pair',
    'tirmite.cli.legacy',
    'tirmite.cli.validate',
    'tirmite.core.extraction',
    'tirmite.core.flanks',
    'tirmite.core.hit_filters',
    'tirmite.core.pairing',
    'tirmite.core.parsers',
    'tirmite.core.termini',
    'tirmite.core.tsd',
    'tirmite.runners.hmmer_wrappers',
    'tirmite.runners.mafft',
    'tirmite.runners.runBlastn',
    'tirmite.runners.wrapping',
    'tirmite.utils.extract',
    'tirmite.utils.utils',
]

# Library modules must not emit anything until a host configures logging.
LIBRARY_MODULES = [m for m in LOGGING_MODULES if not m.startswith('tirmite.cli.')]


class TestNamedLoggers:
    """Each module logs through its own named logger."""

    @pytest.mark.parametrize('module_name', LOGGING_MODULES)
    def test_module_defines_named_logger(self, module_name):
        """The module-level `logger` is named after the module."""
        module = importlib.import_module(module_name)
        assert hasattr(module, 'logger'), f'{module_name} has no module logger'
        assert module.logger.name == module_name

    @pytest.mark.parametrize('module_name', LOGGING_MODULES)
    def test_logger_is_a_child_of_the_package_logger(self, module_name):
        """Configuring `tirmite` must reach every module logger."""
        module = importlib.import_module(module_name)
        assert module.logger.name.startswith(f'{PACKAGE_LOGGER_NAME}.')

    @pytest.mark.parametrize('module_name', LIBRARY_MODULES)
    def test_library_modules_have_a_null_handler(self, module_name):
        """Library modules attach a NullHandler so imports stay silent."""
        module = importlib.import_module(module_name)
        assert any(
            isinstance(h, logging.NullHandler) for h in module.logger.handlers
        ), f'{module_name} has no NullHandler'

    def test_no_module_uses_the_root_logger_functions(self):
        """No source file calls logging.info() and friends directly.

        Those go to the root logger, which defeats per-module control and is
        what this migration removed.
        """
        import re

        import tirmite

        package_root = tirmite.__path__[0]
        offenders = []
        pattern = re.compile(
            r'\blogging\.(debug|info|warning|error|critical|exception)\('
        )
        for info in pkgutil.walk_packages([package_root], prefix='tirmite.'):
            module = importlib.import_module(info.name)
            source_file = getattr(module, '__file__', None)
            if not source_file or not source_file.endswith('.py'):
                continue
            with open(source_file) as handle:
                if pattern.search(handle.read()):
                    offenders.append(info.name)

        assert offenders == [], f'modules still using root logger calls: {offenders}'


class TestInitLoggingDoesNotHijackRoot:
    """init_logging must not disturb a host application's logging."""

    def test_preserves_existing_root_handlers(self):
        """A handler the host installed on root survives init_logging.

        The previous implementation called root_logger.handlers.clear().
        """
        root = logging.getLogger()
        sentinel = logging.NullHandler()
        root.addHandler(sentinel)
        try:
            init_logging(loglevel='INFO')
            assert sentinel in root.handlers
        finally:
            root.removeHandler(sentinel)

    def test_configures_the_package_logger(self):
        """Handlers land on the tirmite logger, not on root."""
        init_logging(loglevel='WARNING')
        package_logger = logging.getLogger(PACKAGE_LOGGER_NAME)

        assert package_logger.level == logging.WARNING
        assert package_logger.handlers, 'expected a console handler'

    def test_repeated_calls_do_not_duplicate_handlers(self):
        """Calling init_logging twice must not double every message."""
        init_logging(loglevel='INFO')
        first = len(logging.getLogger(PACKAGE_LOGGER_NAME).handlers)
        init_logging(loglevel='INFO')
        second = len(logging.getLogger(PACKAGE_LOGGER_NAME).handlers)

        assert first == second

    def test_propagation_stays_enabled(self):
        """Records must still reach host handlers and pytest's caplog."""
        init_logging(loglevel='INFO')
        assert logging.getLogger(PACKAGE_LOGGER_NAME).propagate is True

    def test_rejects_invalid_level(self):
        """An unknown level is a clear error, not a silent default."""
        with pytest.raises(ValueError, match='Invalid log level'):
            init_logging(loglevel='CHATTY')


class TestLogFile:
    """The optional log file."""

    def test_writes_records_to_the_requested_file(self, tmp_path):
        """Messages logged after init_logging reach the file."""
        logfile = tmp_path / 'run.log'
        init_logging(loglevel='INFO', logfile=logfile)

        logging.getLogger('tirmite.core.pairing').info('canary message')
        for handler in logging.getLogger(PACKAGE_LOGGER_NAME).handlers:
            handler.flush()

        assert logfile.exists()
        assert 'canary message' in logfile.read_text()

    def test_creates_missing_parent_directories(self, tmp_path):
        """A log path inside a not-yet-created output directory works."""
        logfile = tmp_path / 'nested' / 'deeper' / 'run.log'
        init_logging(loglevel='INFO', logfile=logfile)

        assert logfile.exists()

    def test_record_carries_module_qualified_name(self, tmp_path):
        """The file format includes the logger name, so source is visible."""
        logfile = tmp_path / 'run.log'
        init_logging(loglevel='INFO', logfile=logfile)

        logging.getLogger('tirmite.core.parsers').info('from the parser')
        for handler in logging.getLogger(PACKAGE_LOGGER_NAME).handlers:
            handler.flush()

        assert 'tirmite.core.parsers' in logfile.read_text()


class TestImportingTirmiteIsSilent:
    """Importing the package must not emit anything to stderr."""

    def test_import_and_warn_produces_no_stderr(self):
        """A library warning with no logging configured is swallowed.

        Run in a subprocess so the check is not affected by handlers pytest
        has already installed in this process.
        """
        code = (
            'import tirmite.core.pairing as m\n'
            "m.logger.warning('this must not appear')\n"
        )
        result = subprocess.run(
            [sys.executable, '-c', code], capture_output=True, text=True
        )

        assert result.returncode == 0, result.stderr
        assert result.stderr == '', f'unexpected stderr: {result.stderr!r}'


class TestSubcommandLoggingOptions:
    """All five subcommands expose the same logging controls."""

    @pytest.mark.parametrize(
        'module_name,configure',
        [
            ('tirmite.cli.legacy', '_configure_legacy_parser'),
            ('tirmite.cli.hmm_build', '_configure_seed_parser'),
            ('tirmite.cli.hmm_pair', '_configure_pair_parser'),
            ('tirmite.cli.ensemble_search', '_configure_search_parser'),
            ('tirmite.cli.validate', '_configure_validate_parser'),
        ],
    )
    def test_subcommand_has_loglevel_and_logfile(self, module_name, configure):
        """legacy and seed previously lacked --logfile entirely."""
        import argparse

        module = importlib.import_module(module_name)
        parser = argparse.ArgumentParser()
        getattr(module, configure)(parser)

        options = {
            option for action in parser._actions for option in action.option_strings
        }
        assert '--loglevel' in options
        assert '--logfile' in options
