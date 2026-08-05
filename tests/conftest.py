"""Shared pytest fixtures and test doubles for the TIRmite test suite.

This module holds helpers that were previously duplicated across individual
test modules.  Keeping a single definition means a change to the fake genome
interface (for example, when ``pyfaidx`` changes its slicing semantics) has to
be made in exactly one place.

Notes
-----
Hit tables throughout TIRmite are ``pandas`` DataFrames whose values are
**strings**, not integers.  This mirrors what the nhmmer/BLAST importers in
:mod:`tirmite.core.parsers` produce, and downstream consumers (the pairing
engine, the GFF writer) rely on that dtype contract.  The ``hit_table_factory``
fixture below preserves it so tests exercise the same types as production.
"""

from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional

import pandas as pd
import pytest

# Column order produced by ``import_nhmmer`` / ``import_blast``.  Kept here so
# tests that build synthetic hit tables stay in sync with the importers.
HIT_TABLE_COLUMNS = [
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

# Values used for any hit-table column a test does not care about.  Chosen to
# be individually harmless: a strong e-value, a mid-range score, no bias.
_HIT_DEFAULTS: Dict[str, str] = {
    'model': 'TIR',
    'target': 'chr1',
    'hitStart': '100',
    'hitEnd': '200',
    'strand': '+',
    'evalue': '1e-10',
    'score': '100.0',
    'bias': '0.0',
    'hmmStart': '1',
    'hmmEnd': '100',
}


class MockChrom:
    """Minimal substitute for a ``pyfaidx`` chromosome record.

    Supports the three operations TIRmite performs on a contig: taking its
    length, slicing out a region, and coercing the whole thing to ``str``.

    Parameters
    ----------
    seq : str
        Nucleotide sequence for this contig.
    """

    def __init__(self, seq: str) -> None:
        self._seq = seq

    def __len__(self) -> int:
        """Return the contig length in bases."""
        return len(self._seq)

    def __getitem__(self, key: Any) -> str:
        """Return a base or a slice of the contig as a plain string.

        ``pyfaidx`` returns a ``Sequence`` object here whose ``str()`` is the
        bases; returning ``str`` directly is close enough for every call site
        in TIRmite, all of which immediately stringify the result.
        """
        return self._seq[key]

    def __str__(self) -> str:
        """Return the full contig sequence."""
        return self._seq


class MockGenome(dict):
    """Dict-backed stand-in for :class:`pyfaidx.Fasta`.

    Parameters
    ----------
    chrom_seqs : mapping of str to str
        Mapping of contig name to nucleotide sequence.  Each value is wrapped
        in a :class:`MockChrom` at construction so that ``genome[name]``
        behaves like a ``pyfaidx`` record.
    """

    def __init__(self, chrom_seqs: Mapping[str, str]) -> None:
        super().__init__()
        for name, seq in chrom_seqs.items():
            self[name] = MockChrom(seq)


@pytest.fixture
def mock_genome() -> Callable[[Mapping[str, str]], MockGenome]:
    """Return a factory building a :class:`MockGenome` from name/sequence pairs.

    Returns
    -------
    callable
        Function taking a ``{contig_name: sequence}`` mapping and returning a
        :class:`MockGenome`.
    """
    return MockGenome


@pytest.fixture
def hit_table_factory() -> Callable[
    [Optional[Iterable[Mapping[str, Any]]]], pd.DataFrame
]:
    """Return a factory building hit-table DataFrames with sensible defaults.

    Each row supplied is merged over :data:`_HIT_DEFAULTS`, so a test need only
    state the fields relevant to what it is checking.  All values are coerced
    to ``str`` to match the importers' dtype contract.

    Returns
    -------
    callable
        Function taking an iterable of partial row mappings and returning a
        ``pandas.DataFrame`` with :data:`HIT_TABLE_COLUMNS`.

    Examples
    --------
    >>> table = hit_table_factory([{'model': 'LeftA', 'hitStart': '10'}])
    """

    def _factory(rows: Optional[Iterable[Mapping[str, Any]]] = None) -> pd.DataFrame:
        if rows is None:
            rows = []
        built: List[Dict[str, str]] = []
        for row in rows:
            merged = dict(_HIT_DEFAULTS)
            merged.update({key: str(value) for key, value in row.items()})
            built.append(merged)
        return pd.DataFrame(built, columns=HIT_TABLE_COLUMNS)

    return _factory


@pytest.fixture(autouse=True)
def reset_logging():
    """Restore the logging configuration after every test.

    ``tirmite.utils.logs.init_logging`` installs handlers and sets levels on
    the ``tirmite`` logger.  Without this fixture, a test that calls it leaks
    that configuration into every subsequent test in the session, which both
    duplicates log output and can break ``caplog`` assertions elsewhere.
    """
    import logging

    root = logging.getLogger()
    tirmite_logger = logging.getLogger('tirmite')

    saved = {
        'root_handlers': list(root.handlers),
        'root_level': root.level,
        'tirmite_handlers': list(tirmite_logger.handlers),
        'tirmite_level': tirmite_logger.level,
        'tirmite_propagate': tirmite_logger.propagate,
    }

    yield

    root.handlers[:] = saved['root_handlers']
    root.setLevel(saved['root_level'])
    tirmite_logger.handlers[:] = saved['tirmite_handlers']
    tirmite_logger.setLevel(saved['tirmite_level'])
    tirmite_logger.propagate = saved['tirmite_propagate']
