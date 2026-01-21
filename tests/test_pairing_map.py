#!/usr/bin/env python3
"""
Tests for pairing map functionality.

Tests pairing map file parsing, validation, and multi-model pairing.
"""

import os
import tempfile
from pathlib import Path

import pytest

from tirmite.cli.hmm_pair import load_pairing_map, check_multiple_models


@pytest.fixture
def symmetric_pairing_map():
    """Create a temporary symmetric pairing map file."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("# Symmetric pairing map\n")
        f.write("model1\tmodel1\n")
        f.write("model2\tmodel2\n")
        fname = f.name
    yield fname
    os.unlink(fname)


@pytest.fixture
def asymmetric_pairing_map():
    """Create a temporary asymmetric pairing map file."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("# Asymmetric pairing map\n")
        f.write("left1\tright1\n")
        f.write("left2\tright2\n")
        fname = f.name
    yield fname
    os.unlink(fname)


@pytest.fixture
def mixed_pairing_map():
    """Create a pairing map with both symmetric and asymmetric pairings."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("# Mixed pairing map\n")
        f.write("model_A\tmodel_A\n")  # Symmetric
        f.write("left_B\tright_B\n")  # Asymmetric
        f.write("model_C\tmodel_C\n")  # Symmetric
        fname = f.name
    yield fname
    os.unlink(fname)


@pytest.fixture
def duplicate_feature_map():
    """Create a pairing map with features appearing multiple times."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("model1\tmodel2\n")
        f.write("model1\tmodel3\n")  # model1 appears again
        fname = f.name
    yield fname
    os.unlink(fname)


def test_load_symmetric_pairing_map(symmetric_pairing_map):
    """Test loading symmetric pairing map."""
    pairings = load_pairing_map(symmetric_pairing_map)
    
    assert len(pairings) == 2
    assert pairings[0] == ('model1', 'model1')
    assert pairings[1] == ('model2', 'model2')


def test_load_asymmetric_pairing_map(asymmetric_pairing_map):
    """Test loading asymmetric pairing map."""
    pairings = load_pairing_map(asymmetric_pairing_map)
    
    assert len(pairings) == 2
    assert pairings[0] == ('left1', 'right1')
    assert pairings[1] == ('left2', 'right2')


def test_load_mixed_pairing_map(mixed_pairing_map):
    """Test loading mixed symmetric and asymmetric pairings."""
    pairings = load_pairing_map(mixed_pairing_map)
    
    assert len(pairings) == 3
    assert pairings[0] == ('model_A', 'model_A')
    assert pairings[1] == ('left_B', 'right_B')
    assert pairings[2] == ('model_C', 'model_C')


def test_duplicate_features_warning(duplicate_feature_map, caplog):
    """Test that duplicate features generate warnings."""
    pairings = load_pairing_map(duplicate_feature_map)
    
    assert len(pairings) == 2
    # Check that warning was logged (feature appears in multiple lines)
    assert any('appears in multiple pairing combinations' in record.message for record in caplog.records)


def test_empty_pairing_map():
    """Test that empty pairing map raises error."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("# Only comments\n")
        f.write("\n")
        fname = f.name
    
    try:
        with pytest.raises(ValueError, match='No valid pairings found'):
            load_pairing_map(fname)
    finally:
        os.unlink(fname)


def test_malformed_pairing_map():
    """Test that malformed pairing map raises error."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("only_one_column\n")
        fname = f.name
    
    try:
        with pytest.raises(ValueError, match='expected 2 tab-delimited columns'):
            load_pairing_map(fname)
    finally:
        os.unlink(fname)


def test_missing_pairing_map():
    """Test that missing file raises appropriate error."""
    with pytest.raises(FileNotFoundError, match='Pairing map file not found'):
        load_pairing_map('/nonexistent/pairing_map.txt')


def test_check_multiple_models():
    """Test check_multiple_models function."""
    import pandas as pd
    
    # Test single model
    single_model_table = pd.DataFrame({
        'model': ['model1', 'model1', 'model1'],
        'target': ['chr1', 'chr1', 'chr2'],
    })
    models = check_multiple_models(single_model_table)
    assert len(models) == 1
    assert models[0] == 'model1'
    
    # Test multiple models
    multi_model_table = pd.DataFrame({
        'model': ['model1', 'model2', 'model3'],
        'target': ['chr1', 'chr1', 'chr2'],
    })
    models = check_multiple_models(multi_model_table)
    assert len(models) == 3
    assert set(models) == {'model1', 'model2', 'model3'}


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
