"""Integration tests for the `tirmite search` hit-processing pipeline.

These exercise `_process_hits` end to end rather than a single filter, because
several of the bugs fixed here only appear from the *interaction* between
steps -- notably the anchor filter running before cluster merging renames hit
models.
"""

import argparse

import pytest

from tirmite.cli.ensemble_search import _process_hits

# BLAST outfmt 6:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
BLAST_COLUMNS = 12


def blast_row(model, target, sstart, send, qstart, qend, bitscore, evalue='1e-40'):
    """
    Build one BLAST tabular row.

    Parameters
    ----------
    model : str
        Query id, which becomes the hit's model name.
    target : str
        Subject id.
    sstart, send : int
        Subject (genomic) coordinates. ``sstart > send`` marks a minus-strand
        hit, matching how BLAST reports one.
    qstart, qend : int
        Query (model) coordinates, which become hmmStart/hmmEnd.
    bitscore : float
        Alignment score.
    evalue : str, default '1e-40'
        E-value field.

    Returns
    -------
    str
        A single tab-delimited BLAST line, newline-terminated.
    """
    fields = [
        model,
        target,
        '95.0',
        '100',
        '5',
        '0',
        str(qstart),
        str(qend),
        str(sstart),
        str(send),
        evalue,
        str(bitscore),
    ]
    assert len(fields) == BLAST_COLUMNS
    return '\t'.join(fields) + '\n'


def make_args(**overrides):
    """
    Build an argparse.Namespace with the fields `_process_hits` reads.

    Parameters
    ----------
    **overrides
        Values to override on the default namespace.

    Returns
    -------
    argparse.Namespace
        Arguments suitable for passing to `_process_hits`.
    """
    defaults = {
        'max_evalue': 1.0,
        'max_offset': None,
        'orientation': 'F,R',
        'cluster_map': None,
        'pairing_map': None,
    }
    defaults.update(overrides)
    return argparse.Namespace(**defaults)


class TestAnchorFilterWithClusterAndPairingMaps:
    """--max-offset, --cluster-map and --pairing-map used together.

    The anchor filter runs *before* cluster merging renames hit models to
    cluster names, because it measures each hit against its own model's
    length and those lengths only exist per component. A cluster-level pairing
    map therefore matched no hit at that point, and terminus assignment fell
    back to strand.

    That fallback only *differs* from the model-name assignment for the
    same-strand orientations F,F and R,R. There, an unmatched model has no
    identifiable outer edge, so the filter falls through to requiring BOTH
    model ends to be anchored -- much stricter than intended, silently
    discarding valid termini. Under F,R the strand fallback happens to agree
    with the role-based assignment, so those runs were unaffected; the tests
    below cover both to pin that down.
    """

    @pytest.fixture
    def cluster_file(self, tmp_path):
        """Cluster-first format: cluster<TAB>component1<TAB>component2..."""
        path = tmp_path / 'clusters.tsv'
        path.write_text('LEFT_TIR\tL1\nRIGHT_TIR\tR1\n')
        return path

    @pytest.fixture
    def pairing_file(self, tmp_path):
        """The pairing map names CLUSTERS, not component models."""
        path = tmp_path / 'pairing.tsv'
        path.write_text('LEFT_TIR\tRIGHT_TIR\n')
        return path

    @pytest.fixture
    def hits_file(self, tmp_path):
        """One left-terminus hit anchored at the model START but not the END.

        A left terminus only needs its outer edge (model position 1) covered;
        the inner edge trailing off is normal and expected.
        """
        path = tmp_path / 'hits.blast'
        path.write_text(blast_row('L1', 'chr1', 1000, 1100, 1, 60, 200))
        return path

    def test_same_strand_orientation_keeps_outer_anchored_hit(
        self, hits_file, cluster_file, pairing_file
    ):
        """F,F: the cluster-level map must still resolve terminus roles.

        Without expansion the model matched nothing, so the filter demanded
        both model ends be anchored and dropped this hit even though its outer
        edge sits exactly at model position 1.
        """
        args = make_args(
            max_offset=5,
            orientation='F,F',
            cluster_map=cluster_file,
            pairing_map=pairing_file,
        )

        result = _process_hits(
            args,
            blast_files=[hits_file],
            nhmmer_files=[],
            query_lengths={'L1': 100, 'R1': 100},
        )

        assert len(result) == 1
        assert result.iloc[0]['model'] == 'LEFT_TIR'

    def test_same_strand_orientation_still_drops_unanchored_hit(
        self, tmp_path, cluster_file, pairing_file
    ):
        """The expansion must not simply disable the filter."""
        hits = tmp_path / 'hits.blast'
        # Aligned from model position 40: 39 short of the outer edge.
        hits.write_text(blast_row('L1', 'chr1', 1000, 1100, 40, 100, 200))

        args = make_args(
            max_offset=5,
            orientation='F,F',
            cluster_map=cluster_file,
            pairing_map=pairing_file,
        )

        result = _process_hits(
            args,
            blast_files=[hits],
            nhmmer_files=[],
            query_lengths={'L1': 100, 'R1': 100},
        )

        assert len(result) == 0

    def test_canonical_fr_orientation_is_unchanged(
        self, hits_file, cluster_file, pairing_file
    ):
        """F,R was already correct, because strand alone identifies the terminus."""
        args = make_args(
            max_offset=5,
            orientation='F,R',
            cluster_map=cluster_file,
            pairing_map=pairing_file,
        )

        result = _process_hits(
            args,
            blast_files=[hits_file],
            nhmmer_files=[],
            query_lengths={'L1': 100, 'R1': 100},
        )

        assert len(result) == 1

    def test_without_max_offset_the_anchor_filter_does_not_run(
        self, tmp_path, cluster_file, pairing_file
    ):
        """A hit anchored at neither end survives when --max-offset is unset."""
        hits = tmp_path / 'hits.blast'
        hits.write_text(blast_row('L1', 'chr1', 1000, 1100, 40, 60, 200))

        args = make_args(
            orientation='F,F', cluster_map=cluster_file, pairing_map=pairing_file
        )

        result = _process_hits(
            args,
            blast_files=[hits],
            nhmmer_files=[],
            query_lengths={'L1': 100, 'R1': 100},
        )

        assert len(result) == 1


class TestPreParsedMapsAreHonoured:
    """_process_hits accepts already-parsed maps so files are read once."""

    def test_supplied_maps_are_used_without_reparsing(self, tmp_path):
        """A caller-supplied map is used even when the path does not exist.

        Proves the pre-parsed value is used rather than silently re-read.
        """
        hits = tmp_path / 'hits.blast'
        hits.write_text(blast_row('L1', 'chr1', 1000, 1100, 1, 60, 200))

        args = make_args(
            cluster_map=tmp_path / 'does-not-exist.tsv',
            pairing_map=tmp_path / 'also-missing.tsv',
        )

        result = _process_hits(
            args,
            blast_files=[hits],
            nhmmer_files=[],
            query_lengths={'L1': 100, 'R1': 100},
            cluster_map={'LEFT_TIR': ['L1']},
            pairing_map={'LEFT_TIR': 'RIGHT_TIR'},
        )

        assert list(result['model'].unique()) == ['LEFT_TIR']

    def test_missing_map_file_is_reported(self, tmp_path):
        """Without a pre-parsed map, a missing file is a clear error."""
        from tirmite.cli.ensemble_search import EnsembleSearchError

        hits = tmp_path / 'hits.blast'
        hits.write_text(blast_row('L1', 'chr1', 1000, 1100, 1, 60, 200))

        args = make_args(cluster_map=tmp_path / 'nope.tsv')

        with pytest.raises(EnsembleSearchError, match='Cluster mapping file not found'):
            _process_hits(
                args,
                blast_files=[hits],
                nhmmer_files=[],
                query_lengths={'L1': 100},
            )
