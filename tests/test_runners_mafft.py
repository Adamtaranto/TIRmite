"""Tests for the shared MAFFT wrappers.

MAFFT is not installed in CI, so every test here fakes the subprocess call.
The point is to pin the two distinct result contracts -- align_to_file raises
and persists, align_in_memory returns None and does not -- which previously
lived as two separate functions that happened to share a name.
"""

from pathlib import Path
import subprocess

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest

from tirmite.runners import mafft
from tirmite.runners.mafft import MafftError, align_in_memory, align_to_file

# A two-sequence alignment as MAFFT would emit it on stdout. Lowercase is
# deliberate: align_to_file must uppercase it, align_in_memory must not.
MAFFT_STDOUT = '>seq1\nacgt-acgt\n>seq2\nacgtaacgt\n'


def _records(n=2):
    """Build ``n`` trivial SeqRecords for use as MAFFT input."""
    return [
        SeqRecord(Seq('ACGTACGT'), id=f'seq{i + 1}', description='') for i in range(n)
    ]


@pytest.fixture
def fake_mafft(monkeypatch):
    """Pretend MAFFT is installed and return a canned alignment.

    Yields the list that captures each invoked command line, so tests can
    assert on the flags passed.
    """
    calls = []

    def _fake_run(cmd, capture_output=False, text=False, check=False):
        calls.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout=MAFFT_STDOUT, stderr='')

    monkeypatch.setattr(mafft.shutil, 'which', lambda _name: '/usr/bin/mafft')
    monkeypatch.setattr(mafft.subprocess, 'run', _fake_run)
    return calls


@pytest.fixture
def failing_mafft(monkeypatch):
    """Pretend MAFFT is installed but exits non-zero."""

    def _fake_run(cmd, capture_output=False, text=False, check=False):
        return subprocess.CompletedProcess(cmd, 1, stdout='', stderr='mafft: boom')

    monkeypatch.setattr(mafft.shutil, 'which', lambda _name: '/usr/bin/mafft')
    monkeypatch.setattr(mafft.subprocess, 'run', _fake_run)


@pytest.fixture
def missing_mafft(monkeypatch):
    """Pretend MAFFT is not installed."""
    monkeypatch.setattr(mafft.shutil, 'which', lambda _name: None)


class TestAlignToFile:
    """align_to_file persists an uppercased alignment or raises."""

    def test_writes_uppercased_alignment(self, fake_mafft, tmp_path: Path):
        """Lowercase MAFFT output is uppercased before being written.

        hmmbuild treats lowercase residues as masked, so a soft-masked region
        would silently drop out of the model if the case were preserved.
        """
        out = tmp_path / 'aln.fasta'
        result = align_to_file(_records(), out, threads=1)

        assert result == out
        text = out.read_text()
        assert 'ACGT-ACGT' in text
        assert 'acgt' not in text

    def test_passes_thread_count(self, fake_mafft, tmp_path: Path):
        """The requested thread count reaches the MAFFT command line."""
        align_to_file(_records(), tmp_path / 'aln.fasta', threads=8)

        assert fake_mafft[0][:3] == ['mafft', '--thread', '8']

    def test_removes_temporary_input(self, fake_mafft, tmp_path: Path):
        """The MAFFT input file does not survive a successful run."""
        out = tmp_path / 'aln.fasta'
        align_to_file(_records(), out, threads=1)

        assert not (tmp_path / 'aln_input.fasta').exists()

    def test_removes_temporary_input_on_failure(self, failing_mafft, tmp_path: Path):
        """The MAFFT input file does not survive a failed run either."""
        out = tmp_path / 'aln.fasta'
        with pytest.raises(MafftError):
            align_to_file(_records(), out, threads=1)

        assert not (tmp_path / 'aln_input.fasta').exists()

    def test_raises_when_mafft_fails(self, failing_mafft, tmp_path: Path):
        """A non-zero exit is reported as MafftError with the stderr text."""
        with pytest.raises(MafftError, match='boom'):
            align_to_file(_records(), tmp_path / 'aln.fasta', threads=1)

    def test_raises_when_mafft_missing(self, missing_mafft, tmp_path: Path):
        """An absent binary is reported clearly rather than as FileNotFoundError."""
        with pytest.raises(MafftError, match='not found in PATH'):
            align_to_file(_records(), tmp_path / 'aln.fasta', threads=1)

    def test_raises_for_single_sequence(self, fake_mafft, tmp_path: Path):
        """One sequence cannot be aligned; MAFFT is never invoked."""
        with pytest.raises(MafftError, match='at least 2 sequences'):
            align_to_file(_records(1), tmp_path / 'aln.fasta', threads=1)

        assert fake_mafft == []


class TestAlignInMemory:
    """align_in_memory returns records, or None on any failure."""

    def test_returns_aligned_records(self, fake_mafft, tmp_path: Path):
        """The parsed alignment is returned to the caller."""
        aligned = align_in_memory(_records(), str(tmp_path))

        assert aligned is not None
        assert len(aligned) == 2
        assert str(aligned[0].seq) == 'acgt-acgt'

    def test_uses_auto_and_quiet(self, fake_mafft, tmp_path: Path):
        """Validation alignments run with --auto --quiet, not --thread."""
        align_in_memory(_records(), str(tmp_path))

        assert fake_mafft[0][:3] == ['mafft', '--auto', '--quiet']

    def test_returns_none_when_mafft_fails(self, failing_mafft, tmp_path: Path):
        """A failed alignment is inconclusive, not fatal, so None is returned.

        The validation workflow carries on with the remaining junctions rather
        than aborting the whole run.
        """
        assert align_in_memory(_records(), str(tmp_path)) is None

    def test_returns_none_when_mafft_missing(self, missing_mafft, tmp_path: Path):
        """An absent binary also yields None rather than raising."""
        assert align_in_memory(_records(), str(tmp_path)) is None

    def test_returns_none_for_single_sequence(self, fake_mafft, tmp_path: Path):
        """Too few sequences is handled the same way as any other failure."""
        assert align_in_memory(_records(1), str(tmp_path)) is None


class TestMafftAvailable:
    """mafft_available reflects PATH lookup."""

    def test_true_when_present(self, fake_mafft):
        """Reports True when shutil.which resolves mafft."""
        assert mafft.mafft_available() is True

    def test_false_when_absent(self, missing_mafft):
        """Reports False when shutil.which returns None."""
        assert mafft.mafft_available() is False
