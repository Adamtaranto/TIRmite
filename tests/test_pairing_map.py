#!/usr/bin/env python3
"""
Tests for pairing map functionality.

Tests pairing map file parsing, validation, and multi-model pairing.
"""

from pathlib import Path

import pytest

from tirmite.cli.hmm_pair import check_multiple_models, load_pairing_map


@pytest.fixture
def symmetric_pairing_map(tmp_path: Path) -> str:
    """Create a temporary symmetric pairing map file."""
    path = tmp_path / 'symmetric.txt'
    path.write_text('# Symmetric pairing map\nmodel1\tmodel1\nmodel2\tmodel2\n')
    return str(path)


@pytest.fixture
def asymmetric_pairing_map(tmp_path: Path) -> str:
    """Create a temporary asymmetric pairing map file."""
    path = tmp_path / 'asymmetric.txt'
    path.write_text('# Asymmetric pairing map\nleft1\tright1\nleft2\tright2\n')
    return str(path)


@pytest.fixture
def mixed_pairing_map(tmp_path: Path) -> str:
    """Create a pairing map with both symmetric and asymmetric pairings."""
    path = tmp_path / 'mixed.txt'
    path.write_text(
        '# Mixed pairing map\n'
        'model_A\tmodel_A\n'  # Symmetric
        'left_B\tright_B\n'  # Asymmetric
        'model_C\tmodel_C\n'  # Symmetric
    )
    return str(path)


@pytest.fixture
def duplicate_feature_map(tmp_path: Path) -> str:
    """Create a pairing map with features appearing multiple times."""
    path = tmp_path / 'duplicate.txt'
    # model1 appears as the left feature of two different pairs.
    path.write_text('model1\tmodel2\nmodel1\tmodel3\n')
    return str(path)


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
    assert any(
        'appears in multiple pairing combinations' in record.message
        for record in caplog.records
    )


def test_empty_pairing_map(tmp_path: Path) -> None:
    """Test that empty pairing map raises error."""
    path = tmp_path / 'empty.txt'
    path.write_text('# Only comments\n\n')

    with pytest.raises(ValueError, match='No valid pairings found'):
        load_pairing_map(str(path))


def test_malformed_pairing_map(tmp_path: Path) -> None:
    """Test that malformed pairing map raises error."""
    path = tmp_path / 'malformed.txt'
    path.write_text('only_one_column\n')

    with pytest.raises(ValueError, match='expected 2 tab-delimited columns'):
        load_pairing_map(str(path))


def test_missing_pairing_map():
    """Test that missing file raises appropriate error."""
    with pytest.raises(FileNotFoundError, match='Pairing map file not found'):
        load_pairing_map('/nonexistent/pairing_map.txt')


def test_check_multiple_models():
    """Test check_multiple_models function."""
    import pandas as pd

    # Test single model
    single_model_table = pd.DataFrame(
        {
            'model': ['model1', 'model1', 'model1'],
            'target': ['chr1', 'chr1', 'chr2'],
        }
    )
    models = check_multiple_models(single_model_table)
    assert len(models) == 1
    assert models[0] == 'model1'

    # Test multiple models
    multi_model_table = pd.DataFrame(
        {
            'model': ['model1', 'model2', 'model3'],
            'target': ['chr1', 'chr1', 'chr2'],
        }
    )
    models = check_multiple_models(multi_model_table)
    assert len(models) == 3
    assert set(models) == {'model1', 'model2', 'model3'}


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
