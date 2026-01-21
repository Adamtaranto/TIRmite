"""Tests for gzip handling utilities in utils.py"""

import gzip
import tempfile
from pathlib import Path

import pytest

from tirmite.utils.utils import is_gzipped_file, decompress_genome, prepare_genome_file


def test_is_gzipped_file_with_gz_extension():
    """Test detection of gzipped files by extension."""
    # Create a temporary gzipped file
    with tempfile.NamedTemporaryFile(suffix='.gz', delete=False) as tmp:
        tmp_path = Path(tmp.name)
        with gzip.open(tmp_path, 'wb') as gz:
            gz.write(b'test data')
    
    try:
        assert is_gzipped_file(tmp_path) is True
    finally:
        tmp_path.unlink()


def test_is_gzipped_file_with_magic_bytes():
    """Test detection of gzipped files by magic bytes."""
    # Create a temporary gzipped file without .gz extension
    with tempfile.NamedTemporaryFile(suffix='.fa', delete=False) as tmp:
        tmp_path = Path(tmp.name)
        with gzip.open(tmp_path, 'wb') as gz:
            gz.write(b'test data')
    
    try:
        assert is_gzipped_file(tmp_path) is True
    finally:
        tmp_path.unlink()


def test_is_gzipped_file_regular_file():
    """Test detection returns False for regular files."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp:
        tmp_path = Path(tmp.name)
        tmp.write('test data')
    
    try:
        assert is_gzipped_file(tmp_path) is False
    finally:
        tmp_path.unlink()


def test_decompress_genome():
    """Test decompression of gzipped genome file."""
    # Create a temporary gzipped file with FASTA content
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        
        # Create gzipped input
        input_gz = tmpdir_path / 'test.fa.gz'
        test_content = b'>seq1\nACGT\n>seq2\nTGCA\n'
        with gzip.open(input_gz, 'wb') as gz:
            gz.write(test_content)
        
        # Decompress
        output_fa = decompress_genome(input_gz, tmpdir_path)
        
        # Verify output exists and has correct content
        assert output_fa.exists()
        assert output_fa.name == 'test.fa'
        
        with open(output_fa, 'rb') as f:
            assert f.read() == test_content


def test_prepare_genome_file_gzipped():
    """Test prepare_genome_file with gzipped input."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        
        # Create gzipped input
        input_gz = tmpdir_path / 'test.fa.gz'
        test_content = b'>seq1\nACGT\n'
        with gzip.open(input_gz, 'wb') as gz:
            gz.write(test_content)
        
        # Prepare genome
        prepared = prepare_genome_file(input_gz, tmpdir_path)
        
        # Should return decompressed file
        assert prepared != input_gz
        assert prepared.exists()
        assert not prepared.name.endswith('.gz')


def test_prepare_genome_file_regular():
    """Test prepare_genome_file with regular file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        
        # Create regular input
        input_fa = tmpdir_path / 'test.fa'
        input_fa.write_text('>seq1\nACGT\n')
        
        # Prepare genome
        prepared = prepare_genome_file(input_fa, tmpdir_path)
        
        # Should return original file unchanged
        assert prepared == input_fa


def test_decompress_genome_missing_file():
    """Test decompress_genome raises error for missing file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        missing_file = tmpdir_path / 'nonexistent.fa.gz'
        
        with pytest.raises(FileNotFoundError):
            decompress_genome(missing_file, tmpdir_path)
