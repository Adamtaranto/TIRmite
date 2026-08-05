"""Tests for gzip handling utilities in utils.py"""

import gzip
from pathlib import Path

import pytest

from tirmite.utils.utils import decompress_genome, is_gzipped_file, prepare_genome_file


def test_is_gzipped_file_with_gz_extension(tmp_path: Path) -> None:
    """Test detection of gzipped files by extension."""
    gz_file = tmp_path / 'sample.gz'
    with gzip.open(gz_file, 'wb') as gz:
        gz.write(b'test data')

    assert is_gzipped_file(gz_file) is True


def test_is_gzipped_file_with_magic_bytes(tmp_path: Path) -> None:
    """Test detection of gzipped files by magic bytes.

    The file deliberately carries a ``.fa`` suffix so that detection cannot
    succeed on the extension alone.
    """
    disguised = tmp_path / 'sample.fa'
    with gzip.open(disguised, 'wb') as gz:
        gz.write(b'test data')

    assert is_gzipped_file(disguised) is True


def test_is_gzipped_file_regular_file(tmp_path: Path) -> None:
    """Test detection returns False for regular files."""
    plain = tmp_path / 'sample.fa'
    plain.write_text('test data')

    assert is_gzipped_file(plain) is False


def test_decompress_genome(tmp_path: Path) -> None:
    """Test decompression of gzipped genome file."""
    input_gz = tmp_path / 'test.fa.gz'
    test_content = b'>seq1\nACGT\n>seq2\nTGCA\n'
    with gzip.open(input_gz, 'wb') as gz:
        gz.write(test_content)

    output_fa = decompress_genome(input_gz, tmp_path)

    # The ``.gz`` suffix is stripped and the payload round-trips byte for byte.
    assert output_fa.exists()
    assert output_fa.name == 'test.fa'
    assert output_fa.read_bytes() == test_content


def test_prepare_genome_file_gzipped(tmp_path: Path) -> None:
    """Test prepare_genome_file with gzipped input."""
    input_gz = tmp_path / 'test.fa.gz'
    with gzip.open(input_gz, 'wb') as gz:
        gz.write(b'>seq1\nACGT\n')

    prepared = prepare_genome_file(input_gz, tmp_path)

    # A gzipped input must be replaced by a decompressed copy.
    assert prepared != input_gz
    assert prepared.exists()
    assert not prepared.name.endswith('.gz')


def test_prepare_genome_file_regular(tmp_path: Path) -> None:
    """Test prepare_genome_file with regular file."""
    input_fa = tmp_path / 'test.fa'
    input_fa.write_text('>seq1\nACGT\n')

    prepared = prepare_genome_file(input_fa, tmp_path)

    # An uncompressed input is passed straight through, not copied.
    assert prepared == input_fa


def test_decompress_genome_missing_file(tmp_path: Path) -> None:
    """Test decompress_genome raises error for missing file."""
    missing_file = tmp_path / 'nonexistent.fa.gz'

    with pytest.raises(FileNotFoundError):
        decompress_genome(missing_file, tmp_path)
