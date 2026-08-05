"""Tests that BLAST database sources are built once and reused.

`BlastDBSource` caches contig lengths and descriptions because each cache miss
costs a `blastdbcmd` subprocess. Constructing one per hit discarded that cache
every time, so every hit paid two spawns -- one to resolve the contig length
for coordinate clamping, one to fetch the sequence -- instead of amortising
the first across all hits on the same contig.

These tests count subprocess invocations rather than timing anything, so they
are deterministic and need no BLAST installation.
"""

import subprocess

import pytest

from tirmite.cli.hmm_build import _LazySource
from tirmite.cli.validate import extract_hit_sequence
import tirmite.utils.extract as extract_module


@pytest.fixture
def blastdbcmd_counter(monkeypatch):
    """Count blastdbcmd invocations and serve canned responses.

    Returns a dict whose ``count`` key tracks how many times ``blastdbcmd``
    was spawned. Avoids needing BLAST+ installed.
    """
    state = {'count': 0}
    sequence = 'ACGT' * 500

    def fake_run(cmd, *args, **kwargs):
        if not cmd or 'blastdbcmd' not in str(cmd[0]):
            raise AssertionError(f'unexpected subprocess: {cmd}')
        state['count'] += 1

        parts = [str(c) for c in cmd]

        # A metadata lookup asks for length and title; anything else is a
        # range fetch returning FASTA.
        if '%l' in ' '.join(parts):
            stdout = f'{len(sequence)}\tchr1 test contig\n'
        else:
            # Honour -range, because fetch_raw verifies the returned length
            # against the request -- blastdbcmd returns the WHOLE contig with
            # a zero exit status when a range starts past the end, and that
            # check exists to reject exactly that.
            start, end = 1, len(sequence)
            if '-range' in parts:
                start_s, _, end_s = parts[parts.index('-range') + 1].partition('-')
                start, end = int(start_s), int(end_s)
            stdout = f'>chr1\n{sequence[start - 1 : end]}\n'

        return subprocess.CompletedProcess(cmd, 0, stdout=stdout, stderr='')

    monkeypatch.setattr(extract_module.subprocess, 'run', fake_run)
    return state


class TestValidateReusesOneSource:
    """`tirmite validate` must not rebuild the source per hit."""

    def test_shared_source_amortises_the_length_lookup(self, blastdbcmd_counter):
        """One shared source costs one spawn per hit, plus one for the length.

        A fresh source per hit costs two per hit, because the contig-length
        cache is discarded before it can ever be hit.
        """
        source = extract_module.BlastDBSource('dummy_db')

        n_hits = 20
        for i in range(n_hits):
            result = extract_hit_sequence(
                source, 'chr1', 100 + i * 10, 200 + i * 10, 'plus'
            )
            assert result is not None

        # n fetches + a single metadata lookup for chr1.
        assert blastdbcmd_counter['count'] == n_hits + 1

    def test_fresh_source_per_hit_would_double_the_spawns(self, blastdbcmd_counter):
        """Pins the old cost, so a regression back to it is visible."""
        n_hits = 20
        for i in range(n_hits):
            extract_hit_sequence(
                extract_module.BlastDBSource('dummy_db'),
                'chr1',
                100 + i * 10,
                200 + i * 10,
                'plus',
            )

        assert blastdbcmd_counter['count'] == n_hits * 2

    def test_extract_hit_sequence_takes_a_source_not_a_path(self):
        """The signature change is what makes reuse possible.

        Passing a path forced construction inside the function, which is why
        the cache could never survive between calls.
        """
        import inspect

        first_param = list(inspect.signature(extract_hit_sequence).parameters)[0]
        assert first_param == 'source'


class TestLazySource:
    """`tirmite seed` builds one source per workflow, on first use."""

    def test_builds_only_once(self, monkeypatch):
        """Repeated get() calls return the same object.

        The ID pre-flight and the extraction pass both need the source; each
        used to build its own, so the genome was indexed twice (three times
        for an asymmetric run) and any populated cache was thrown away.
        """
        builds = {'count': 0}

        def fake_build(blast_db, genome_files):
            builds['count'] += 1
            return object()

        monkeypatch.setattr('tirmite.cli.hmm_build._build_source', fake_build)

        lazy = _LazySource(None, ['genome.fa'])
        first = lazy.get()
        second = lazy.get()

        assert builds['count'] == 1
        assert first is second

    def test_is_lazy(self, monkeypatch):
        """Constructing the holder must not index anything.

        A run that fails before extraction -- no hits above threshold, say --
        should never pay to index a genome it will not read.
        """
        builds = {'count': 0}

        def fake_build(blast_db, genome_files):
            builds['count'] += 1
            return object()

        monkeypatch.setattr('tirmite.cli.hmm_build._build_source', fake_build)

        _LazySource(None, ['genome.fa'])

        assert builds['count'] == 0

    @pytest.mark.parametrize(
        'blast_db,genome_files,expected',
        [
            (None, [], False),
            (None, ['genome.fa'], True),
            ('db', [], True),
            ('db', ['genome.fa'], True),
        ],
    )
    def test_available_reports_whether_a_source_could_be_built(
        self, blast_db, genome_files, expected
    ):
        """`available` lets callers skip an optional check without building.

        Forcing construction in the ID pre-flight would raise "No genome or
        BLAST database available for extraction" there, pre-empting the more
        specific error a run is actually entitled to -- a run with no hits
        above threshold should say so, not complain about extraction it never
        reached.
        """
        assert _LazySource(blast_db, genome_files).available is expected
