"""Tests for hit-table sorting in the importers.

Hit coordinates are stored as strings, so sorting on them directly compares
lexicographically -- '1000' sorts before '999'. Row order is not cosmetic:
`filter_best_model_per_locus` walks each group in order and removes greedily,
so ordering decides which of two mutually-eliminating hits survives.
"""

import pandas as pd
import pytest

from tirmite.core.parsers import import_blast, import_nhmmer, sort_hit_table

NHMMER_HEADER = (
    '#target name         accession query name           accession mdl mdl from   '
    'mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc '
    'description of target\n'
)


def nhmmer_row(target, model, seq_from, seq_to, strand='+'):
    """Build one nhmmer ``--tblout`` row.

    Column order, which `import_nhmmer` indexes positionally:
    0 target name, 1 accession, 2 query name, 3 accession, 4 hmmfrom,
    5 hmm to, 6 alifrom, 7 ali to, 8 envfrom, 9 env to, 10 sq len,
    11 strand, 12 E-value, 13 score, 14 bias, 15 description.
    """
    fields = [
        target,
        '-',
        model,
        '-',  # target/query names and accessions
        '1',
        '60',  # hmm from / to
        str(seq_from),
        str(seq_to),  # ali from / to
        str(seq_from),
        str(seq_to),  # env from / to
        '100000',  # sq len
        strand,
        '1.2e-10',
        '45.2',
        '0.0',  # E-value, score, bias
        '-',  # description
    ]
    return ' '.join(fields) + '\n'


def blast_row(model, target, sstart, send):
    """Build one BLAST outfmt 6 row."""
    return f'{model}\t{target}\t95.0\t100\t5\t0\t1\t100\t{sstart}\t{send}\t1e-40\t200\n'


class TestSortHitTable:
    """The shared sort helper."""

    def _table(self, starts):
        """Build a minimal hit table with the given string hitStart values."""
        return pd.DataFrame(
            {
                'model': ['M'] * len(starts),
                'target': ['chr1'] * len(starts),
                'hitStart': [str(s) for s in starts],
                'hitEnd': [str(s + 50) for s in starts],
                'strand': ['+'] * len(starts),
            }
        )

    def test_sorts_coordinates_numerically(self):
        """999 precedes 1000, which lexicographic ordering gets backwards."""
        result = sort_hit_table(self._table([1000, 999, 10000, 99]))

        assert list(result['hitStart']) == ['99', '999', '1000', '10000']

    def test_preserves_string_dtype(self):
        """The all-string dtype contract must survive sorting.

        table2dict, the pairing engine and the GFF writer all rely on it.
        """
        result = sort_hit_table(self._table([1000, 999]))

        assert result['hitStart'].map(type).eq(str).all()
        assert result['hitEnd'].map(type).eq(str).all()

    def test_drops_temporary_sort_columns(self):
        """No helper columns leak into the returned table."""
        result = sort_hit_table(self._table([1000, 999]))

        assert list(result.columns) == [
            'model',
            'target',
            'hitStart',
            'hitEnd',
            'strand',
        ]

    def test_resets_the_index(self):
        """Downstream code indexes positionally after sorting."""
        result = sort_hit_table(self._table([1000, 999, 500]))

        assert list(result.index) == [0, 1, 2]

    def test_empty_table_is_returned_unchanged(self):
        """An empty result is common and must not raise."""
        empty = self._table([]).iloc[0:0]

        assert sort_hit_table(empty).empty

    def test_malformed_coordinate_does_not_raise(self):
        """A non-numeric coordinate sorts last rather than aborting the run."""
        table = self._table([1000, 999])
        table.loc[0, 'hitStart'] = 'not-a-number'

        result = sort_hit_table(table)

        assert len(result) == 2
        assert result.iloc[-1]['hitStart'] == 'not-a-number'

    def test_sorts_by_model_then_target_then_position(self):
        """The full key order is model, target, start, end, strand."""
        table = pd.DataFrame(
            {
                'model': ['B', 'A', 'A'],
                'target': ['chr1', 'chr2', 'chr1'],
                'hitStart': ['100', '100', '2000'],
                'hitEnd': ['200', '200', '2100'],
                'strand': ['+', '+', '+'],
            }
        )

        result = sort_hit_table(table)

        assert list(result['model']) == ['A', 'A', 'B']
        assert list(result['target']) == ['chr1', 'chr2', 'chr1']


class TestImportersSortNumerically:
    """Both importers apply the numeric sort."""

    def test_import_nhmmer_sorts_numerically(self, tmp_path):
        """999 precedes 1000 on the same target."""
        path = tmp_path / 'hits.nhmmer'
        path.write_text(
            NHMMER_HEADER
            + nhmmer_row('chr1', 'ModelA', 1000, 1060)
            + nhmmer_row('chr1', 'ModelA', 999, 1059)
        )

        result = import_nhmmer(str(path))

        assert list(result['hitStart']) == ['999', '1000']

    def test_import_blast_sorts_numerically(self, tmp_path):
        """Same for the BLAST importer."""
        path = tmp_path / 'hits.blast'
        path.write_text(
            blast_row('ModelA', 'chr1', 1000, 1100)
            + blast_row('ModelA', 'chr1', 999, 1099)
        )

        result = import_blast(str(path))

        assert list(result['hitStart']) == ['999', '1000']

    def test_import_blast_orders_a_wide_coordinate_range(self, tmp_path):
        """Ordering holds across several magnitudes."""
        path = tmp_path / 'hits.blast'
        path.write_text(
            ''.join(
                blast_row('ModelA', 'chr1', start, start + 100)
                for start in (10000, 99, 1000, 999, 2)
            )
        )

        result = import_blast(str(path))

        assert list(result['hitStart']) == ['2', '99', '999', '1000', '10000']

    def test_concatenated_import_is_sorted_as_a_whole(self, tmp_path):
        """Appending to an existing table re-sorts the combined result."""
        first = tmp_path / 'a.blast'
        first.write_text(blast_row('ModelA', 'chr1', 1000, 1100))
        second = tmp_path / 'b.blast'
        second.write_text(blast_row('ModelA', 'chr1', 999, 1099))

        table = import_blast(str(first))
        combined = import_blast(str(second), hitTable=table)

        assert list(combined['hitStart']) == ['999', '1000']


class TestSortAffectsGreedyFiltering:
    """Row order determines which hit a greedy filter keeps."""

    def test_ordering_is_deterministic_and_numeric(self):
        """The filter sees hits in genomic order, not lexicographic order.

        With lexicographic ordering a hit at 10000 preceded one at 999, so a
        greedy pass could eliminate the wrong member of an overlapping pair.
        """
        table = pd.DataFrame(
            {
                'model': ['M'] * 3,
                'target': ['chr1'] * 3,
                'hitStart': ['10000', '999', '2000'],
                'hitEnd': ['10100', '1099', '2100'],
                'strand': ['+'] * 3,
            }
        )

        result = sort_hit_table(table)

        starts = [int(s) for s in result['hitStart']]
        assert starts == sorted(starts)


@pytest.mark.parametrize('importer', [import_nhmmer, import_blast])
def test_importers_expose_the_documented_columns(importer, tmp_path):
    """Both importers return the same normalised column set."""
    if importer is import_nhmmer:
        path = tmp_path / 'hits.nhmmer'
        path.write_text(NHMMER_HEADER + nhmmer_row('chr1', 'ModelA', 100, 200))
    else:
        path = tmp_path / 'hits.blast'
        path.write_text(blast_row('ModelA', 'chr1', 100, 200))

    result = importer(str(path))

    assert list(result.columns) == [
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
