"""Correctness tests for `tirmite seed` hit resolution and model naming.

Each class here covers a bug that silently degraded the HMM training set
rather than producing an error.
"""

import pytest

from tirmite.cli.hmm_build import (
    BlastHit,
    clean_hmm_name,
    hits_overlap,
    resolve_asymmetric_conflicts,
    resolve_overlapping_hits,
)


def make_hit(subject_start, subject_end, identity=90.0, subject_id='chr1'):
    """
    Build a BlastHit with the given subject coordinates.

    BLAST reports ``sstart > send`` for a minus-strand HSP, and BlastHit stores
    them unswapped, so passing ``start > end`` here produces a minus-strand hit
    exactly as the parser would.
    """
    length = abs(subject_end - subject_start) + 1
    return BlastHit(
        query_id='seed',
        subject_id=subject_id,
        query_start=1,
        query_end=length,
        subject_start=subject_start,
        subject_end=subject_end,
        length=length,
        identity=identity,
        query_len=length,
        subject_len=100000,
    )


class TestResolveOverlappingHitsStrand:
    """Overlap resolution must normalise coordinates before comparing.

    BLAST inverts sstart/send for minus-strand HSPs. Comparing them raw made
    the computed overlap of any minus-strand hit provably <= 0, so no
    minus-strand hit was ever de-overlapped and every reverse-oriented element
    contributed redundant near-identical copies to the HMM training set.
    """

    def test_minus_strand_overlaps_are_resolved(self):
        """Two overlapping minus-strand hits collapse to one."""
        hits = [make_hit(1000, 801), make_hit(990, 791)]

        assert len(resolve_overlapping_hits(hits)) == 1

    def test_plus_strand_control(self):
        """The equivalent plus-strand pair behaves the same way."""
        hits = [make_hit(801, 1000), make_hit(791, 990)]

        assert len(resolve_overlapping_hits(hits)) == 1

    def test_mixed_strand_overlap_is_resolved(self):
        """A plus and a minus hit at one locus are still competitors."""
        hits = [make_hit(801, 1000), make_hit(990, 791)]

        assert len(resolve_overlapping_hits(hits)) == 1

    def test_distant_minus_hits_both_survive(self):
        """Normalising must not make the filter over-eager."""
        hits = [make_hit(1000, 801), make_hit(9000, 8801)]

        assert len(resolve_overlapping_hits(hits)) == 2

    def test_hits_on_different_contigs_never_conflict(self):
        """Only hits on the same subject are compared."""
        hits = [make_hit(1000, 801), make_hit(1000, 801, subject_id='chr2')]

        assert len(resolve_overlapping_hits(hits)) == 2

    def test_longer_hit_wins(self):
        """The retained hit is the longer of an overlapping pair."""
        hits = [make_hit(1000, 801), make_hit(1200, 791)]

        kept = resolve_overlapping_hits(hits)

        assert len(kept) == 1
        assert kept[0].length == 410


class TestResolveAsymmetricConflicts:
    """Left/right conflict resolution must not cascade."""

    def test_discarded_left_hit_cannot_evict_a_right_hit(self):
        """The regression test for cascading eviction.

        L1 loses to R1, so L1 is discarded. R2 overlaps only L1. The second
        loop iterated the ORIGINAL left list, so the already-discarded L1 still
        evicted R2 -- a legitimate right-terminus hit lost to a hit that no
        longer existed.
        """
        left_hit = make_hit(1000, 1600)  # 601 bp
        better_right = make_hit(1500, 2200)  # 701 bp, beats left_hit
        innocent_right = make_hit(900, 1100)  # overlaps left_hit only

        kept_left, kept_right = resolve_asymmetric_conflicts(
            [left_hit], [better_right, innocent_right]
        )

        assert kept_left == []
        assert len(kept_right) == 2
        assert any(h.subject_start == 900 for h in kept_right)

    def test_left_wins_ties(self):
        """The >= tie-break deliberately favours the left model.

        With > on both sides an exact tie would survive in both sets and
        duplicate the sequence; with >= on both it would vanish from both.
        """
        left_hit = make_hit(1000, 1200)
        right_hit = make_hit(1000, 1200)

        kept_left, kept_right = resolve_asymmetric_conflicts([left_hit], [right_hit])

        assert len(kept_left) == 1
        assert kept_right == []

    def test_non_overlapping_hits_both_survive(self):
        """Hits that do not conflict are untouched."""
        kept_left, kept_right = resolve_asymmetric_conflicts(
            [make_hit(1000, 1200)], [make_hit(5000, 5200)]
        )

        assert len(kept_left) == 1
        assert len(kept_right) == 1

    def test_hits_overlap_is_symmetric(self):
        """The overlap predicate must not depend on argument order."""
        a, b = make_hit(1000, 1200), make_hit(1100, 1300)

        assert hits_overlap(a, b) == hits_overlap(b, a)

    def test_hits_overlap_normalises_strand(self):
        """A minus-strand hit overlaps its plus-strand twin."""
        assert hits_overlap(make_hit(1000, 1200), make_hit(1300, 1100))


class TestCleanHmmName:
    """Model names must stay distinct after cleaning."""

    @pytest.mark.parametrize(
        'base',
        [
            'hAT_Tag1_Arabidopsis',
            'MuDR_Mutator_Zea_mays',
            'a_very_long_transposon_family_name_indeed',
        ],
    )
    def test_left_and_right_never_collide(self, base):
        """The regression test.

        The cap was 20 characters, so '<name>_left' and '<name>_right'
        truncated to the same string for any name of 15+ characters. Since the
        HMM path is f'{clean_model_name}.hmm', the right model OVERWROTE the
        left one and both returned paths pointed at the same file.
        """
        assert clean_hmm_name(f'{base}_left') != clean_hmm_name(f'{base}_right')

    def test_typical_names_are_not_truncated(self):
        """Ordinary names pass through recognisably."""
        assert clean_hmm_name('hAT_Tag1_Arabidopsis_left') == (
            'hAT_Tag1_Arabidopsis_left'
        )

    def test_very_long_names_stay_distinct(self):
        """Truncation appends a digest, so collisions remain impossible."""
        long_left = 'X' * 70 + '_left'
        long_right = 'X' * 70 + '_right'

        assert clean_hmm_name(long_left) != clean_hmm_name(long_right)

    def test_truncation_is_deterministic(self):
        """Re-running must produce the same filename."""
        name = 'Y' * 80

        assert clean_hmm_name(name) == clean_hmm_name(name)

    def test_truncated_name_respects_the_cap(self):
        """The result stays within the length limit."""
        assert len(clean_hmm_name('Z' * 200)) <= 60

    def test_empty_name_gets_a_placeholder(self):
        """A name that cleans away entirely still yields something usable."""
        assert clean_hmm_name('---') != ''


class TestAddBlastHeader:
    """The header is prepended by streaming, not by rewriting in memory."""

    def _row(self, i):
        """One BLAST tabular row in the format the seed search emits."""
        return (
            f'1\t100\t{i * 100}\t{i * 100 + 99}\t100\t95\t95.0\t100\t100000'
            f'\t1\t1\tquery\tchr1\n'
        )

    def test_prepends_header_and_preserves_body(self, tmp_path):
        """Every original byte survives, after the header."""
        from tirmite.cli.hmm_build import add_blast_header

        path = tmp_path / 'hits.tab'
        body = ''.join(self._row(i) for i in range(5))
        path.write_text(body)

        add_blast_header(path)
        text = path.read_text()

        assert text.startswith('# BLAST tabular output format 6\n')
        assert text.endswith(body)

    def test_result_still_parses(self, tmp_path):
        """parse_blast_output skips the comment lines and finds every hit."""
        from tirmite.cli.hmm_build import add_blast_header, parse_blast_output

        path = tmp_path / 'hits.tab'
        path.write_text(''.join(self._row(i) for i in range(5)))

        add_blast_header(path)

        assert len(parse_blast_output(path)) == 5

    def test_empty_file_gets_only_a_header(self, tmp_path):
        """A zero-hit search is legitimate and must not error."""
        from tirmite.cli.hmm_build import add_blast_header, parse_blast_output

        path = tmp_path / 'hits.tab'
        path.write_text('')

        add_blast_header(path)

        assert path.read_text().startswith('# BLAST tabular output format 6')
        assert parse_blast_output(path) == []

    def test_peak_memory_does_not_scale_with_file_size(self, tmp_path):
        """The regression test for the 2x-file-size peak.

        The previous implementation read the whole table into a string and
        wrote header + string, so peak memory was ~2x the file. That matters
        because the genome search runs with -word_size 4 and
        -max_target_seqs 10000, which can make the table very large.
        """
        import tracemalloc

        from tirmite.cli.hmm_build import add_blast_header

        def peak_for(n_rows):
            path = tmp_path / f'hits{n_rows}.tab'
            with open(path, 'w') as handle:
                for i in range(n_rows):
                    handle.write(self._row(i))
            tracemalloc.start()
            add_blast_header(path)
            _, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            return path.stat().st_size, peak

        small_size, small_peak = peak_for(5_000)
        large_size, large_peak = peak_for(50_000)

        # The file grew ~10x; peak allocation must not follow it.
        assert large_size > small_size * 5
        assert large_peak < small_peak * 2
        # And peak must stay well under the file size itself.
        assert large_peak < large_size

    def test_original_survives_a_failure_midway(self, tmp_path, monkeypatch):
        """An interrupted rewrite must not destroy the search results.

        The old version opened the destination with 'w', truncating the
        original *before* writing, so a failure anywhere in the write lost the
        table outright. The replacement writes a sibling temp file and swaps it
        in with os.replace.
        """
        from tirmite.cli import hmm_build

        path = tmp_path / 'hits.tab'
        body = ''.join(self._row(i) for i in range(20))
        path.write_text(body)

        def explode(*args, **kwargs):
            raise OSError('disk full')

        monkeypatch.setattr(hmm_build.shutil, 'copyfileobj', explode)

        with pytest.raises(OSError, match='disk full'):
            hmm_build.add_blast_header(path)

        assert path.read_text() == body

    def test_no_temporary_files_are_left_behind(self, tmp_path, monkeypatch):
        """A failed rewrite cleans up its scratch file."""
        from tirmite.cli import hmm_build

        path = tmp_path / 'hits.tab'
        path.write_text(''.join(self._row(i) for i in range(5)))

        def explode(*args, **kwargs):
            raise OSError('disk full')

        monkeypatch.setattr(hmm_build.shutil, 'copyfileobj', explode)

        with pytest.raises(OSError):
            hmm_build.add_blast_header(path)

        assert [p.name for p in tmp_path.iterdir()] == ['hits.tab']
