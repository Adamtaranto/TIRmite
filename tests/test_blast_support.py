#!/usr/bin/env python3
"""
Tests for BLAST input support in TIRmite.

Tests BLAST parsing, format detection, and basic functionality.
"""

import os
import tempfile

import pandas as pd
import pytest

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


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
