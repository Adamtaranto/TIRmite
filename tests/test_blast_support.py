#!/usr/bin/env python3
"""
Tests for BLAST input support in TIRmite.

Tests BLAST parsing, format detection, and basic functionality.
"""

import os
import tempfile

import pandas as pd
import pytest

from tirmite.runners.runBlastn import BlastError, blast_db_exists, run_blastn
import tirmite.tirmitetools as tirmite


@pytest.fixture
def blast_file():
    """Create a temporary BLAST tabular output file."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.blast', delete=False) as f:
        f.write('# BLASTN 2.12.0+\n')
        f.write('# Query: TIR_query1\n')
        f.write('# Database: test_genome\n')
        f.write(
            '# Fields: query acc.ver, subject acc.ver, % identity, alignment length, '
            'mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score\n'
        )
        f.write('TIR_query1\tchr1\t95.00\t100\t5\t0\t1\t100\t1000\t1099\t1e-40\t200\n')
        f.write('TIR_query1\tchr1\t93.00\t98\t7\t0\t3\t100\t5500\t5597\t1e-35\t180\n')
        f.write(
            'TIR_query1\tchr2\t90.00\t95\t10\t0\t5\t99\t2000\t1906\t1e-30\t160\n'
        )  # Reverse strand
        fname = f.name
    yield fname
    os.unlink(fname)


@pytest.fixture
def nhmmer_file():
    """Create a temporary nhmmer output file."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.nhmmer', delete=False) as f:
        f.write(
            '#target name         accession query name           accession mdl mdl from   '
            'mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc '
            'description of target\n'
        )
        f.write(
            '#------------------- --------- -------------------- --------- --- -------- '
            '-------- -------- -------- ------ ----- ---- ---- ----- ------ --------- '
            '--- ---------------------\n'
        )
        f.write(
            'chr1                 -         TIR_model1           -          cm        1       '
            '60     1000     1059      +    no    1 0.45   0.0   45.2   1.2e-10 yes -\n'
        )
        fname = f.name
    yield fname
    os.unlink(fname)


def test_detect_blast_format(blast_file):
    """Test detection of BLAST format."""
    fmt = tirmite.detect_input_format(blast_file)
    assert fmt == 'blast', f"Expected 'blast', got '{fmt}'"


def test_detect_nhmmer_format(nhmmer_file):
    """Test detection of nhmmer format."""
    fmt = tirmite.detect_input_format(nhmmer_file)
    assert fmt == 'nhmmer', f"Expected 'nhmmer', got '{fmt}'"


def test_import_blast(blast_file):
    """Test importing BLAST tabular output."""
    hitTable = tirmite.import_blast(blast_file)

    # Check basic structure
    assert isinstance(hitTable, pd.DataFrame)
    assert len(hitTable) == 3, f'Expected 3 hits, got {len(hitTable)}'

    # Check columns
    expected_cols = [
        'model',
        'target',
        'hitStart',
        'hitEnd',
        'strand',
        'evalue',
        'score',
        'bias',
        'hmmStart',
        'hmmEnd',
    ]
    assert list(hitTable.columns) == expected_cols

    # Check first hit (forward strand)
    first = hitTable.iloc[0]
    assert first['model'] == 'TIR_query1'
    assert first['target'] == 'chr1'
    assert first['hitStart'] == '1000'
    assert first['hitEnd'] == '1099'
    assert first['strand'] == '+'
    assert first['evalue'] == '1e-40'
    assert first['hmmStart'] == '1'
    assert first['hmmEnd'] == '100'

    # Check reverse strand hit (coordinates should be swapped)
    third = hitTable.iloc[2]
    assert third['strand'] == '-'
    assert third['hitStart'] == '1906'  # Should be smaller
    assert third['hitEnd'] == '2000'  # Should be larger
    assert int(third['hitStart']) < int(third['hitEnd'])


def test_import_blast_with_query_name(blast_file):
    """Test importing BLAST with explicit query name."""
    hitTable = tirmite.import_blast(blast_file, query_name='MyCustomQuery')

    # All hits should have the custom query name as model
    assert all(hitTable['model'] == 'MyCustomQuery')


def test_import_blast_empty_file():
    """Test importing empty BLAST file."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.blast', delete=False) as f:
        f.write('# BLASTN 2.12.0+\n')
        f.write('# No hits found\n')
        fname = f.name

    try:
        hitTable = tirmite.import_blast(fname)
        assert len(hitTable) == 0
        # Should still have correct columns
        expected_cols = [
            'model',
            'target',
            'hitStart',
            'hitEnd',
            'strand',
            'evalue',
            'score',
            'bias',
            'hmmStart',
            'hmmEnd',
        ]
        assert list(hitTable.columns) == expected_cols
    finally:
        os.unlink(fname)


def test_import_blast_concatenate(blast_file):
    """Test concatenating multiple BLAST imports."""
    # First import
    hitTable1 = tirmite.import_blast(blast_file)
    assert len(hitTable1) == 3

    # Second import with concatenation
    hitTable2 = tirmite.import_blast(blast_file, hitTable=hitTable1)
    assert len(hitTable2) == 6  # Should have double the hits


def test_import_blast_multiple_queries():
    """Test importing BLAST file with multiple query sequences."""
    # Create BLAST file with multiple queries
    with tempfile.NamedTemporaryFile(mode='w', suffix='.blast', delete=False) as f:
        f.write('query1\tseq1\t100.000\t50\t0\t0\t1\t50\t100\t149\t1e-20\t100\n')
        f.write('query1\tseq2\t95.000\t50\t2\t0\t1\t50\t200\t249\t1e-18\t95\n')
        f.write('query2\tseq1\t98.000\t45\t1\t0\t1\t45\t300\t344\t1e-19\t97\n')
        f.write('query3\tseq3\t100.000\t40\t0\t0\t1\t40\t500\t539\t1e-17\t90\n')
        fname = f.name

    try:
        hitTable = tirmite.import_blast(fname)

        # Should have 4 hits
        assert len(hitTable) == 4

        # Should have 3 unique models (query1, query2, query3)
        unique_models = hitTable['model'].unique()
        assert len(unique_models) == 3
        assert 'query1' in unique_models
        assert 'query2' in unique_models
        assert 'query3' in unique_models

        # Check hit counts per model
        assert len(hitTable[hitTable['model'] == 'query1']) == 2
        assert len(hitTable[hitTable['model'] == 'query2']) == 1
        assert len(hitTable[hitTable['model'] == 'query3']) == 1

    finally:
        os.unlink(fname)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])


# -----------------------------------------------------------------------------
# Tests for blast_db_exists
# -----------------------------------------------------------------------------


def test_blast_db_exists_returns_false_for_missing_prefix():
    """blast_db_exists returns False when no DB files exist for prefix."""
    assert blast_db_exists('/nonexistent/path/no_db') is False


def test_blast_db_exists_returns_true_for_nhr_file(tmp_path):
    """blast_db_exists returns True when a .nhr file exists for the prefix."""
    db_prefix = tmp_path / 'mydb'
    (tmp_path / 'mydb.nhr').touch()
    assert blast_db_exists(db_prefix) is True


def test_blast_db_exists_returns_true_for_nin_file(tmp_path):
    """blast_db_exists returns True when a .nin file exists for the prefix."""
    db_prefix = tmp_path / 'mydb'
    (tmp_path / 'mydb.nin').touch()
    assert blast_db_exists(db_prefix) is True


def test_blast_db_exists_returns_true_for_nsq_file(tmp_path):
    """blast_db_exists returns True when a .nsq file exists for the prefix."""
    db_prefix = tmp_path / 'mydb'
    (tmp_path / 'mydb.nsq').touch()
    assert blast_db_exists(db_prefix) is True


def test_blast_db_exists_string_prefix(tmp_path):
    """blast_db_exists accepts a string path prefix."""
    db_prefix = str(tmp_path / 'mydb')
    (tmp_path / 'mydb.nhr').touch()
    assert blast_db_exists(db_prefix) is True


# -----------------------------------------------------------------------------
# Tests for run_blastn subject validation
# -----------------------------------------------------------------------------


def test_run_blastn_raises_for_missing_subject(tmp_path):
    """run_blastn raises FileNotFoundError when subject is neither a file nor a DB."""
    query_file = tmp_path / 'query.fa'
    query_file.write_text('>seq1\nACGT\n')
    output_file = tmp_path / 'out.tab'

    with pytest.raises(
        FileNotFoundError, match='Subject file or BLAST database not found'
    ):
        run_blastn(
            query=query_file,
            subject=tmp_path / 'nonexistent_db',
            output=output_file,
        )


def test_run_blastn_raises_for_missing_query(tmp_path):
    """run_blastn raises FileNotFoundError when query file is missing."""
    output_file = tmp_path / 'out.tab'

    with pytest.raises(FileNotFoundError, match='Query file not found'):
        run_blastn(
            query=tmp_path / 'nonexistent_query.fa',
            subject=tmp_path / 'subject.fa',
            output=output_file,
        )


def test_run_blastn_accepts_blast_db_prefix(tmp_path, monkeypatch):
    """run_blastn does not raise FileNotFoundError when subject is a valid BLAST DB prefix."""
    import tirmite.runners.runBlastn as runBlastn_mod

    query_file = tmp_path / 'query.fa'
    query_file.write_text('>seq1\nACGT\n')
    output_file = tmp_path / 'out.tab'

    # Create fake BLAST DB files
    db_prefix = tmp_path / 'mydb'
    (tmp_path / 'mydb.nhr').touch()
    (tmp_path / 'mydb.nin').touch()
    (tmp_path / 'mydb.nsq').touch()

    # Mock check_blast_available to avoid requiring BLAST installed
    monkeypatch.setattr(runBlastn_mod, 'check_blast_available', lambda: False)

    # Should raise BlastError (not available) rather than FileNotFoundError
    with pytest.raises(BlastError, match='blastn not found in PATH'):
        run_blastn(
            query=query_file,
            subject=db_prefix,
            output=output_file,
        )


# -----------------------------------------------------------------------------
# Tests for timeout parameter
# -----------------------------------------------------------------------------


def test_run_blastn_timeout_raises_blast_error(tmp_path, monkeypatch):
    """run_blastn raises BlastError when the configured timeout expires."""
    import subprocess
    import tirmite.runners.runBlastn as runBlastn_mod

    query_file = tmp_path / 'query.fa'
    query_file.write_text('>seq1\nACGT\n')
    subject_file = tmp_path / 'subject.fa'
    subject_file.write_text('>seq2\nACGT\n')
    output_file = tmp_path / 'out.tab'

    # Simulate a long-running process that never finishes within the timeout
    class _SlowPopen:
        """Mock Popen that never finishes (poll always returns None)."""
        returncode = None

        def poll(self):
            return None  # Never done

        def kill(self):
            pass

        def wait(self):
            pass

        def communicate(self):
            return ('', '')

    monkeypatch.setattr(runBlastn_mod, 'check_blast_available', lambda: True)
    monkeypatch.setattr(subprocess, 'Popen', lambda *args, **kwargs: _SlowPopen())

    with pytest.raises(BlastError, match='timed out after 1 second'):
        run_blastn(
            query=query_file,
            subject=subject_file,
            output=output_file,
            timeout=1,
        )


def test_run_blastn_no_timeout_by_default(tmp_path, monkeypatch):
    """run_blastn accepts timeout=None (the default) without raising."""
    import subprocess
    import tirmite.runners.runBlastn as runBlastn_mod

    query_file = tmp_path / 'query.fa'
    query_file.write_text('>seq1\nACGT\n')
    subject_file = tmp_path / 'subject.fa'
    subject_file.write_text('>seq2\nACGT\n')
    output_file = tmp_path / 'out.tab'

    # Create an output file so the existence check passes
    output_file.touch()

    class _FastPopen:
        """Mock Popen that finishes immediately."""
        returncode = 0
        _polled = False

        def poll(self):
            # Return 0 (finished) on first call
            return 0

        def communicate(self):
            return ('', '')

    monkeypatch.setattr(runBlastn_mod, 'check_blast_available', lambda: True)
    monkeypatch.setattr(subprocess, 'Popen', lambda *args, **kwargs: _FastPopen())

    # Should not raise – no timeout means no limit
    result = run_blastn(
        query=query_file,
        subject=subject_file,
        output=output_file,
        timeout=None,
    )
    assert result.returncode == 0


# -----------------------------------------------------------------------------
# Tests for --blast-timeout CLI option
# -----------------------------------------------------------------------------


def test_search_parser_accepts_blast_timeout():
    """The search subcommand parser accepts --blast-timeout."""
    from tirmite.cli.ensemble_search import create_search_parser

    parser = create_search_parser()
    args = parser.parse_args([
        '--blast-results', 'dummy.tab',
        '--blast-timeout', '3600',
    ])
    assert args.blast_timeout == 3600


def test_search_parser_blast_timeout_default_is_none():
    """The --blast-timeout option defaults to None (no limit)."""
    from tirmite.cli.ensemble_search import create_search_parser

    parser = create_search_parser()
    args = parser.parse_args(['--blast-results', 'dummy.tab'])
    assert args.blast_timeout is None
