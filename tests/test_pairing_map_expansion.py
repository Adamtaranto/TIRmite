"""Tests for expanding a cluster-level pairing map to component model names.

``tirmite search`` applies the anchor filter before cluster merging renames
hit models to cluster names, because anchor offsets are measured against
per-model query lengths that only exist at component level. A cluster-level
pairing map therefore matches nothing at that point unless it is expanded.
"""

from tirmite.core.hit_filters import expand_pairing_map_to_components


class TestExpandPairingMap:
    """Cross-product expansion of cluster pairs into component pairs."""

    def test_expands_cluster_pair_to_cross_product(self):
        """Each left component is paired with each right component."""
        result = expand_pairing_map_to_components(
            {'LeftCluster': 'RightCluster'},
            {'LeftCluster': ['L1', 'L2'], 'RightCluster': ['R1', 'R2']},
        )
        assert result == [('L1', 'R1'), ('L1', 'R2'), ('L2', 'R1'), ('L2', 'R2')]

    def test_accepts_tuple_list_shape(self):
        """A list of tuples expands identically to the dict form."""
        result = expand_pairing_map_to_components(
            [('LeftCluster', 'RightCluster')],
            {'LeftCluster': ['L1'], 'RightCluster': ['R1']},
        )
        assert result == [('L1', 'R1')]

    def test_unclustered_feature_passes_through(self):
        """A feature absent from the cluster map is used as-is.

        This is what lets a pairing map mix cluster names with bare model
        names without the user having to declare single-member clusters.
        """
        result = expand_pairing_map_to_components(
            {'LeftCluster': 'PlainRightModel'},
            {'LeftCluster': ['L1', 'L2']},
        )
        assert result == [('L1', 'PlainRightModel'), ('L2', 'PlainRightModel')]

    def test_self_paired_cluster_yields_symmetric_rows(self):
        """A self-paired cluster produces left == right rows.

        Those rows are exactly what filter_hits_by_anchor treats as symmetric,
        so a symmetric element expressed at cluster level stays symmetric at
        component level.
        """
        result = expand_pairing_map_to_components(
            {'SymCluster': 'SymCluster'},
            {'SymCluster': ['S1', 'S2']},
        )
        assert ('S1', 'S1') in result
        assert ('S2', 'S2') in result
        assert ('S1', 'S2') in result

    def test_no_cluster_map_returns_pairs_unchanged(self):
        """Without a cluster map the pairing map is already component-level."""
        result = expand_pairing_map_to_components({'L': 'R'}, None)
        assert result == [('L', 'R')]

    def test_empty_pairing_map_returns_empty(self):
        """No pairs in, no pairs out."""
        assert expand_pairing_map_to_components(None, {'C': ['a']}) == []
        assert expand_pairing_map_to_components({}, {'C': ['a']}) == []

    def test_deduplicates_repeated_pairs(self):
        """Overlapping clusters must not yield the same pair twice.

        A duplicated pair would be harmless but would inflate the map handed
        to the anchor filter, so it is collapsed here.
        """
        result = expand_pairing_map_to_components(
            [('CA', 'CB'), ('CA', 'CB')],
            {'CA': ['x'], 'CB': ['y']},
        )
        assert result == [('x', 'y')]

    def test_preserves_order(self):
        """Expansion order follows the pairing map, then cluster membership."""
        result = expand_pairing_map_to_components(
            [('C2', 'C1')],
            {'C1': ['b1', 'b2'], 'C2': ['a1']},
        )
        assert result == [('a1', 'b1'), ('a1', 'b2')]


class TestExpandedMapDrivesAnchorFilter:
    """The expanded map must actually resolve terminus roles for components."""

    def test_component_hits_get_terminus_roles(self):
        """Component-named hits are filtered using cluster-level pairing.

        Before expansion the cluster-level map named LeftCluster/RightCluster,
        which matched no hit, so terminus assignment silently fell back to
        strand.
        """
        import pandas as pd

        from tirmite.core.hit_filters import filter_hits_by_anchor

        cluster_map = {'LeftCluster': ['L1'], 'RightCluster': ['R1']}
        expanded = expand_pairing_map_to_components(
            {'LeftCluster': 'RightCluster'}, cluster_map
        )

        hits = pd.DataFrame(
            [
                # L1 on '+' anchored at model start: a valid left terminus.
                {
                    'model': 'L1',
                    'target': 'chr1',
                    'hitStart': '100',
                    'hitEnd': '200',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '1',
                    'hmmEnd': '60',
                },
                # L1 on '+' starting 39 positions into the model: interior.
                {
                    'model': 'L1',
                    'target': 'chr1',
                    'hitStart': '500',
                    'hitEnd': '600',
                    'strand': '+',
                    'evalue': '1e-10',
                    'hmmStart': '40',
                    'hmmEnd': '100',
                },
            ]
        )

        result = filter_hits_by_anchor(
            hits,
            {'L1': 100, 'R1': 100},
            max_offset=5,
            orientation='F,R',
            pairing_map=expanded,
        )

        assert len(result) == 1
        assert result.iloc[0]['hitStart'] == '100'
