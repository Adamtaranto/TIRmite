"""
Build report data from the internal results of a ``tirmite pair`` run.

This is the only report module that knows about TIRmite's internal structures
(``hitIndex``, ``paired``, the element namedtuples and the hit table). Keeping
that knowledge in one place means the rest of the report -- the layout, the
renderer, the figures -- works against the stable dataclasses in
:mod:`tirmite.report.model` and does not have to change when the pairing
internals do.

Usage follows the shape of ``hmm_pair.main``: construct the accumulator once,
call :meth:`PairReportAccumulator.add_group` after each pairing procedure, then
:meth:`PairReportAccumulator.finalise` at the end.
"""

from dataclasses import replace
from datetime import datetime, timezone
import logging
import math
from statistics import median
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from tirmite.core.termini import _pair_roles, flipTIRs
from tirmite.report.layout import stack_contig
from tirmite.report.model import (
    ContigInfo,
    ElementRecord,
    HitColumns,
    HitRecord,
    ModelInfo,
    PairingGroup,
    ReportData,
    SequenceStore,
)
from tirmite.report.palette import group_colours
from tirmite.report.stats import PairSummary

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = ['PairReportAccumulator', 'model_truncation']

# Per-element ceiling on embedded sequence. A single multi-megabase "element"
# is almost always a mis-pairing, and embedding it would blow the whole budget
# on one bogus record.
MAX_EMBEDDED_ELEMENT_BP = 50_000


def _to_float(value: Any) -> Optional[float]:
    """
    Coerce a hit-table value to float, returning None when it is not numeric.

    Parameters
    ----------
    value : any
        Value from the hit table. These columns are typed as strings and may
        hold 'NA' for formats that do not supply the statistic.

    Returns
    -------
    float or None
        The numeric value, or None if it is missing, non-numeric, NaN or
        infinite. Non-finite values are rejected because they cannot be
        represented in strict JSON.
    """
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(result) or math.isinf(result):
        return None
    return result


def _to_int(value: Any) -> Optional[int]:
    """
    Coerce a hit-table value to int, returning None when it is not numeric.

    Parameters
    ----------
    value : any
        Value from the hit table.

    Returns
    -------
    int or None
        The integer value, or None if missing or non-numeric.
    """
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def model_truncation(
    hmm_start: Optional[int],
    hmm_end: Optional[int],
    model_length: Optional[int],
    strand: str,
) -> Tuple[Optional[int], Optional[int]]:
    """
    Compute how many model positions are uncovered at each genomic end.

    Parameters
    ----------
    hmm_start, hmm_end : int or None
        One-based inclusive alignment coordinates on the model. Ascending even
        for a minus-strand hit.
    model_length : int or None
        Declared model length.
    strand : str
        '+' or '-'.

    Returns
    -------
    tuple of (int or None, int or None)
        Base pairs of model missing at the lower and higher genomic coordinate
        ends. ``(None, None)`` when the inputs are insufficient.

    Notes
    -----
    On the plus strand the deficit before `hmm_start` sits at the lower
    genomic coordinate; on the minus strand the mapping is mirrored. This
    matches :func:`tirmite.core.extraction._model_extension`, which uses the
    same quantities to extend hits to full model length, so the jagged edges
    drawn in the report agree with the sequence that ``--extend-hits-to-model``
    would extract.

    Negative deficits are clamped to zero. They mean the alignment ran past the
    declared end of the model, which indicates a mismatched ``--lengths-file``
    or ``--query-len`` applied across several queries rather than a real
    truncation.
    """
    if hmm_start is None or hmm_end is None or model_length is None:
        return None, None

    start_deficit = max(0, hmm_start - 1)
    end_deficit = max(0, model_length - hmm_end)

    if strand == '-':
        return end_deficit, start_deficit
    return start_deficit, end_deficit


class PairReportAccumulator:
    """
    Collect the results of a ``tirmite pair`` run into report data.

    Parameters
    ----------
    tirmite_version : str
        Version string to record in the report.
    command : str
        Command line that produced the run.
    title : str
        Report heading.
    params : dict
        Run parameters worth surfacing to the reader.
    filter_stats : dict
        Hit-filtering statistics as accumulated by the CLI.
    hit_table : pandas.DataFrame
        The filtered hit table. Used for the per-hit statistics that the
        pairing structures do not carry, notably the model coordinates.
    model_lengths : dict, optional
        Mapping of model name to declared length.
    contig_length : callable, optional
        Called with a contig name, returning its length or None. Normally
        ``SequenceSource.contig_length``. When absent or returning None,
        lengths are inferred from the hits and flagged as such.
    max_hits : int, default 200000
        Cap on the number of hits included in the report.
    max_rows : int, default 30
        Cap on stacked annotation rows per contig.

    Notes
    -----
    ``add_group`` copies the pair sets it is given. The pairing-map loop in
    ``hmm_pair.main`` reuses one hit index across every row, so the pairing
    state it holds is overwritten by later rows; snapshotting at call time is
    what keeps each group's membership correct.
    """

    def __init__(
        self,
        *,
        tirmite_version: str = '',
        command: str = '',
        title: str = 'TIRmite pair report',
        params: Optional[Dict[str, Any]] = None,
        filter_stats: Optional[Dict[str, Any]] = None,
        hit_table: Any = None,
        model_lengths: Optional[Dict[str, int]] = None,
        contig_length: Optional[Callable[[str], Optional[int]]] = None,
        max_hits: int = 200_000,
        max_rows: int = 30,
    ) -> None:
        self.tirmite_version = tirmite_version
        self.command = command
        self.title = title
        self.params = dict(params) if params else {}
        self.filter_stats = dict(filter_stats) if filter_stats else {}
        self.model_lengths = dict(model_lengths) if model_lengths else {}
        self.contig_length = contig_length
        self.max_hits = max_hits
        self.max_rows = max_rows
        self.warnings: List[str] = []

        # Per-hit facts that only the hit table carries, keyed by row index.
        self._hit_meta = self._index_hit_table(hit_table)

        self._groups: List[PairingGroup] = []
        self._group_ids: List[str] = []
        # uid -> ordered list of group ids. A model listed in several
        # pairing-map rows produces hits belonging to more than one group.
        self._hit_groups: Dict[int, List[str]] = {}
        self._hit_records: Dict[int, HitRecord] = {}
        # pair_id -> (group_id, left_hitTup, right_hitTup, element or None)
        self._pairs: List[Tuple[str, str, Any, Any, Any]] = []
        self._element_seqs: Dict[str, str] = {}
        # Pairing configuration per group, needed to assign terminus roles.
        self._config_by_group: Dict[str, Any] = {}
        # Contig each pair sits on, filled while building the elements.
        self._pair_contig: Dict[str, str] = {}

    # -- input indexing ---------------------------------------------------

    @staticmethod
    def _index_hit_table(hit_table: Any) -> Dict[int, Dict[str, Any]]:
        """
        Extract the per-hit fields the pairing structures do not carry.

        Parameters
        ----------
        hit_table : pandas.DataFrame or None
            The filtered hit table.

        Returns
        -------
        dict
            Mapping of row index to a dict with 'hmm_start', 'hmm_end',
            'score' and 'span_coverage'.
        """
        meta: Dict[int, Dict[str, Any]] = {}
        if hit_table is None or not len(hit_table):
            return meta

        has_coverage = 'coverage' in hit_table.columns
        for idx, row in hit_table.iterrows():
            meta[idx] = {
                'hmm_start': _to_int(row.get('hmmStart')),
                'hmm_end': _to_int(row.get('hmmEnd')),
                'score': _to_float(row.get('score')),
                'span_coverage': (
                    _to_float(row.get('coverage')) if has_coverage else None
                ),
            }
        return meta

    # -- accumulation -----------------------------------------------------

    def add_group(
        self,
        *,
        left_feature: str,
        right_feature: str,
        config: Any,
        paired: Dict[str, List[Set[int]]],
        hit_index: Dict[str, Dict[int, Dict[str, Any]]],
        elements: Optional[Dict[str, List[Any]]] = None,
        summary: Optional[PairSummary] = None,
    ) -> None:
        """
        Record the outcome of one pairing procedure.

        Parameters
        ----------
        left_feature, right_feature : str
            Models forming this pairing.
        config : PairingConfig
            Configuration the pairing used. Supplies the orientation and, for
            asymmetric pairings, the terminus roles.
        paired : dict
            ``paired[model] = [{idxA, idxB}, ...]`` for this pairing.
        hit_index : dict
            Hit index used by this pairing, holding the hit records.
        elements : dict, optional
            Extracted elements keyed by model, as returned by ``fetchElements``.
            Absent when element extraction was disabled.
        summary : PairSummary, optional
            Counts for this pairing, so the report reproduces exactly the
            numbers written to the text summary.

        Returns
        -------
        None
            Updates the accumulator in place.
        """
        group_id = f'{left_feature}__{right_feature}'
        label = (
            f'{left_feature} <-> {right_feature}'
            if left_feature != right_feature
            else left_feature
        )
        orientation = ','.join(getattr(config, 'orientation', []) or [])

        # Element records are keyed by their terminus hit indices so each pair
        # can find the element it produced. fetchElements iterates a set, so
        # the two hits are re-identified by coordinate, never by tuple order.
        element_by_pair = self._index_elements(elements)

        n_pairs = 0
        for model_pairs in paired.values():
            for pair_set in model_pairs:
                members = self._resolve_pair(pair_set, hit_index)
                if members is None:
                    continue
                left_rec, right_rec = members
                pair_id = f'{group_id}:{n_pairs}'
                key = frozenset((left_rec.idx, right_rec.idx))
                self._pairs.append(
                    (
                        pair_id,
                        group_id,
                        left_rec,
                        right_rec,
                        element_by_pair.get(key),
                    )
                )
                n_pairs += 1

        # Every hit in this pairing's index belongs to the group, paired or
        # not: an unpaired terminus is still evidence about this model pair.
        for model in (left_feature, right_feature):
            for uid, entry in hit_index.get(model, {}).items():
                groups = self._hit_groups.setdefault(uid, [])
                if group_id not in groups:
                    groups.append(group_id)
                if uid not in self._hit_records:
                    self._hit_records[uid] = self._make_hit(entry['rec'])

        n_elements = sum(len(v) for v in elements.values()) if elements else 0
        n_unpaired = summary.total_unpaired if summary else 0
        hits_per_model = dict(summary.hits_per_model) if summary else {}

        self._groups.append(
            PairingGroup(
                group_id=group_id,
                label=label,
                left_model=left_feature,
                right_model=right_feature,
                orientation=orientation,
                colour='',  # assigned in finalise, once all groups are known
                n_pairs=n_pairs,
                n_elements=n_elements,
                n_unpaired=n_unpaired,
                hits_per_model=hits_per_model,
            )
        )
        self._group_ids.append(group_id)
        self._config_by_group[group_id] = config

    def add_unpaired(self, hit_index: Dict[str, Dict[int, Dict[str, Any]]]) -> None:
        """
        Record hits that no pairing procedure claimed.

        Parameters
        ----------
        hit_index : dict
            The run's full hit index.

        Returns
        -------
        None
            Updates the accumulator in place.

        Notes
        -----
        Hits for a model absent from the pairing map never enter any group, but
        they were still found and belong on the tracks.
        """
        for entries in hit_index.values():
            for uid, entry in entries.items():
                if uid not in self._hit_records:
                    self._hit_records[uid] = self._make_hit(entry['rec'])
                    self._hit_groups.setdefault(uid, [])

    @staticmethod
    def _index_elements(
        elements: Optional[Dict[str, List[Any]]],
    ) -> Dict[frozenset, Any]:
        """
        Key extracted elements by the pair of hit indices that produced them.

        Parameters
        ----------
        elements : dict or None
            Elements keyed by model, as returned by ``fetchElements``.

        Returns
        -------
        dict
            Mapping of ``frozenset({left_idx, right_idx})`` to the element.
        """
        indexed: Dict[frozenset, Any] = {}
        if not elements:
            return indexed
        for records in elements.values():
            for element in records:
                left = getattr(element, 'leftHit', None)
                right = getattr(element, 'rightHit', None)
                if left is None or right is None:
                    continue
                indexed[frozenset((left.idx, right.idx))] = element
        return indexed

    @staticmethod
    def _resolve_pair(
        pair_set: Iterable[int],
        hit_index: Dict[str, Dict[int, Dict[str, Any]]],
    ) -> Optional[Tuple[Any, Any]]:
        """
        Look up a pair's two hit records and order them genomically.

        Parameters
        ----------
        pair_set : iterable of int
            The pair's two hit indices. This is a set, so its iteration order
            carries no meaning.
        hit_index : dict
            Hit index holding the records.

        Returns
        -------
        tuple or None
            ``(left_record, right_record)`` in genomic order, or None when
            either member could not be found.
        """
        found = []
        for uid in pair_set:
            for entries in hit_index.values():
                if uid in entries:
                    found.append(entries[uid]['rec'])
                    break
        if len(found) != 2:
            logger.debug(
                f'Skipping a pair whose members could not both be resolved: '
                f'{sorted(pair_set)}'
            )
            return None
        return flipTIRs(found[0], found[1])

    def _make_hit(self, rec: Any) -> HitRecord:
        """
        Convert a pairing hit record into a :class:`HitRecord`.

        Parameters
        ----------
        rec : namedtuple
            Hit record from ``table2dict``.

        Returns
        -------
        HitRecord
            The report's representation, without pairing or layout fields,
            which are filled in during :meth:`finalise`.
        """
        meta = self._hit_meta.get(rec.idx, {})
        hmm_start = meta.get('hmm_start')
        hmm_end = meta.get('hmm_end')
        model_length = self.model_lengths.get(rec.model)

        model_coverage = None
        if hmm_start is not None and hmm_end is not None and model_length:
            covered = max(0, hmm_end - hmm_start + 1)
            model_coverage = round(min(1.0, covered / model_length), 4)

        trunc_left, trunc_right = model_truncation(
            hmm_start, hmm_end, model_length, rec.strand
        )

        return HitRecord(
            uid=rec.idx,
            model=rec.model,
            contig=rec.target,
            start=int(rec.hitStart),
            end=int(rec.hitEnd),
            strand=rec.strand,
            evalue=_to_float(rec.evalue),
            score=meta.get('score'),
            hmm_start=hmm_start,
            hmm_end=hmm_end,
            model_length=model_length,
            model_coverage=model_coverage,
            span_coverage=meta.get('span_coverage'),
            trunc_left=trunc_left,
            trunc_right=trunc_right,
        )

    # -- finalisation -----------------------------------------------------

    def finalise(
        self,
        *,
        embed_sequences: bool = True,
        max_seq_bytes: int = 20 * 1024 * 1024,
        elements_fasta_path: str = '',
        generated: Optional[str] = None,
    ) -> ReportData:
        """
        Produce the finished report data.

        Parameters
        ----------
        embed_sequences : bool, default True
            Whether to embed element sequences for the copy-to-clipboard
            popup.
        max_seq_bytes : int, default 20971520
            Byte budget for embedded sequences.
        elements_fasta_path : str, optional
            Path shown to the reader for elements whose sequence was not
            embedded.
        generated : str, optional
            ISO-8601 timestamp. Defaults to now, in UTC.

        Returns
        -------
        ReportData
            Ready to render.
        """
        colours = group_colours(self._group_ids)
        groups = [
            PairingGroup(
                group_id=g.group_id,
                label=g.label,
                left_model=g.left_model,
                right_model=g.right_model,
                orientation=g.orientation,
                colour=colours[g.group_id],
                n_pairs=g.n_pairs,
                n_elements=g.n_elements,
                n_unpaired=g.n_unpaired,
                hits_per_model=g.hits_per_model,
            )
            for g in self._groups
        ]
        group_index = {g.group_id: i for i, g in enumerate(groups)}

        hits = self._select_hits()
        contig_lengths = self._resolve_contig_lengths(hits)
        self._apply_clipping(hits, contig_lengths)

        pairs_by_contig = self._pairs_on_contigs(hits)
        elements, element_index = self._build_elements(
            hits, group_index, contig_lengths
        )
        self._annotate_pairing(hits, elements, element_index)

        hits = sorted(hits, key=lambda h: (h.contig, h.start, h.end, h.uid))
        hits = self._stack(hits, pairs_by_contig, contig_lengths)

        contigs = self._build_contigs(hits, elements, contig_lengths)
        contig_index = {c.name: i for i, c in enumerate(contigs)}
        models = self._build_models(hits)
        model_index = {m.name: i for i, m in enumerate(models)}

        # Element contig indices were placeholders until the contig list
        # existed; rebuild them now that the ordering is fixed.
        elements = [
            ElementRecord(
                pair_id=e.pair_id,
                group_i=e.group_i,
                element_id=e.element_id,
                contig_i=contig_index[self._pair_contig[e.pair_id]],
                start=e.start,
                end=e.end,
                strand=e.strand,
                length=e.length,
                left_uid=e.left_uid,
                right_uid=e.right_uid,
                inner_distance=e.inner_distance,
            )
            for e in elements
        ]

        sequences = self._build_sequences(
            elements, embed_sequences, max_seq_bytes, elements_fasta_path
        )

        return ReportData(
            kind='pair',
            tirmite_version=self.tirmite_version,
            generated=generated or datetime.now(timezone.utc).isoformat(),
            command=self.command,
            title=self.title,
            params=self.params,
            filter_stats=self.filter_stats,
            models=models,
            groups=groups,
            contigs=contigs,
            hits=HitColumns.from_records(
                hits, model_index, contig_index, group_index, element_index
            ),
            elements=elements,
            sequences=sequences,
            warnings=list(self.warnings),
        )

    def _select_hits(self) -> List[HitRecord]:
        """
        Return the hits to include, applying the hit cap.

        Returns
        -------
        list of HitRecord
            At most `max_hits` records. When the cap bites, paired hits are
            kept in preference to unpaired ones, and the omission is recorded
            as a warning rather than passing silently.
        """
        records = list(self._hit_records.values())
        if len(records) <= self.max_hits:
            return records

        paired_uids = set()
        for _, _, left, right, _ in self._pairs:
            paired_uids.add(left.idx)
            paired_uids.add(right.idx)

        ranked = sorted(
            records,
            key=lambda h: (
                0 if h.uid in paired_uids else 1,
                h.evalue if h.evalue is not None else float('inf'),
                h.uid,
            ),
        )
        dropped = len(records) - self.max_hits
        self.warnings.append(
            f'{dropped:,} of {len(records):,} hits were omitted to keep the '
            f'report responsive (--report-max-hits {self.max_hits:,}). Paired '
            'hits and the most significant unpaired hits were kept.'
        )
        return ranked[: self.max_hits]

    def _resolve_contig_lengths(
        self, hits: Sequence[HitRecord]
    ) -> Dict[str, Tuple[int, str]]:
        """
        Determine each contig's length and where that length came from.

        Parameters
        ----------
        hits : sequence of HitRecord
            The hits to be reported.

        Returns
        -------
        dict
            Mapping of contig name to ``(length, source)`` where source is
            'source' or 'inferred'.
        """
        furthest: Dict[str, int] = {}
        for hit in hits:
            if hit.end > furthest.get(hit.contig, 0):
                furthest[hit.contig] = hit.end

        lengths: Dict[str, Tuple[int, str]] = {}
        inferred = 0
        for name, max_end in furthest.items():
            length = None
            if self.contig_length is not None:
                try:
                    length = self.contig_length(name)
                except Exception as exc:  # noqa: BLE001 - a length is optional
                    logger.debug(f'Could not read length of contig {name}: {exc}')
            if length:
                lengths[name] = (int(length), 'source')
            else:
                # Without a sequence source the axis can only be approximate.
                # Pad past the last hit so the final annotation is not flush
                # against the edge of the track.
                lengths[name] = (int(max_end * 1.05) + 1, 'inferred')
                inferred += 1

        if inferred:
            self.warnings.append(
                f'Contig lengths were unavailable for {inferred:,} sequence(s) '
                'and have been estimated from the hits. Track axes for those '
                'contigs are approximate.'
            )
        return lengths

    def _apply_clipping(
        self,
        hits: List[HitRecord],
        contig_lengths: Dict[str, Tuple[int, str]],
    ) -> None:
        """
        Mark truncations that are caused by a contig end.

        Parameters
        ----------
        hits : list of HitRecord
            Hits to update. Replaced in place.
        contig_lengths : dict
            Contig name to ``(length, source)``.

        Returns
        -------
        None
            Rewrites the list in place.

        Notes
        -----
        A hit is short for one of two reasons: it only partially matched its
        model, or the sequence ran out. The report draws those differently -- a
        jagged edge versus a contig-end cap -- so they must be distinguished
        here. Only a truncation that would need bases beyond the contig counts
        as clipped.
        """
        for i, hit in enumerate(hits):
            length, source = contig_lengths.get(hit.contig, (0, 'inferred'))
            clip_left = bool(hit.trunc_left and hit.start - hit.trunc_left < 1)
            clip_right = bool(
                source == 'source'
                and hit.trunc_right
                and hit.end + hit.trunc_right > length
            )
            if clip_left or clip_right:
                hits[i] = _replace_hit(hit, clip_left=clip_left, clip_right=clip_right)

    def _pairs_on_contigs(
        self, hits: Sequence[HitRecord]
    ) -> Dict[str, List[Tuple[int, int]]]:
        """
        Group pair memberships by contig for stacking.

        Parameters
        ----------
        hits : sequence of HitRecord
            The hits to be reported.

        Returns
        -------
        dict
            Contig name to a list of ``(left_uid, right_uid)``.
        """
        included = {h.uid for h in hits}
        by_contig: Dict[str, List[Tuple[int, int]]] = {}
        for _, _, left, right, _ in self._pairs:
            if left.idx not in included or right.idx not in included:
                continue
            by_contig.setdefault(left.target, []).append((left.idx, right.idx))
        return by_contig

    def _build_elements(
        self,
        hits: Sequence[HitRecord],
        group_index: Dict[str, int],
        contig_lengths: Dict[str, Tuple[int, str]],
    ) -> Tuple[List[ElementRecord], Dict[str, int]]:
        """
        Build element records for every resolvable pair.

        Parameters
        ----------
        hits : sequence of HitRecord
            The hits to be reported.
        group_index : dict
            Group id to index.
        contig_lengths : dict
            Contig name to ``(length, source)``. Unused for coordinates, but
            confirms the contig is present.

        Returns
        -------
        elements : list of ElementRecord
            One per pair, in the order pairs were accumulated.
        element_index : dict
            Pair id to index in `elements`.
        """
        included = {h.uid for h in hits}
        elements: List[ElementRecord] = []
        element_index: Dict[str, int] = {}
        self._pair_contig.clear()

        for pair_id, group_id, left, right, element in self._pairs:
            if left.idx not in included or right.idx not in included:
                continue
            if left.target not in contig_lengths:
                continue

            start = int(left.hitStart)
            end = int(right.hitEnd)
            element_id = getattr(element, 'id', None)
            if element_id is None:
                element_id = pair_id.rsplit(':', 1)[-1]

            element_index[pair_id] = len(elements)
            self._pair_contig[pair_id] = left.target
            elements.append(
                ElementRecord(
                    pair_id=pair_id,
                    group_i=group_index[group_id],
                    element_id=str(element_id),
                    contig_i=0,  # rewritten once the contig ordering is fixed
                    start=start,
                    end=end,
                    strand=getattr(element, 'strand', left.strand),
                    length=end - start + 1,
                    left_uid=left.idx,
                    right_uid=right.idx,
                    inner_distance=int(right.hitStart) - int(left.hitEnd) - 1,
                )
            )

            seq = getattr(element, 'seq', None)
            if seq is not None:
                self._element_seqs[pair_id] = str(getattr(seq, 'seq', seq))

        return elements, element_index

    def _annotate_pairing(
        self,
        hits: List[HitRecord],
        elements: Sequence[ElementRecord],
        element_index: Dict[str, int],
    ) -> None:
        """
        Attach pairing group, role and element membership to each hit.

        Parameters
        ----------
        hits : list of HitRecord
            Hits to update. Replaced in place.
        elements : sequence of ElementRecord
            Element records.
        element_index : dict
            Pair id to element index.

        Returns
        -------
        None
            Rewrites the list in place.
        """
        roles: Dict[int, str] = {}
        pair_of: Dict[int, str] = {}
        configs = getattr(self, '_config_by_group', {})

        for pair_id, group_id, left, right, _ in self._pairs:
            if pair_id not in element_index:
                continue
            left_role, right_role = _pair_roles_safe(left, right, configs.get(group_id))
            roles[left.idx] = left_role
            roles[right.idx] = right_role
            pair_of[left.idx] = pair_id
            pair_of[right.idx] = pair_id

        for i, hit in enumerate(hits):
            hits[i] = _replace_hit(
                hit,
                group_ids=list(self._hit_groups.get(hit.uid, [])),
                role=roles.get(hit.uid),
                pair_id=pair_of.get(hit.uid),
            )

    def _stack(
        self,
        hits: List[HitRecord],
        pairs_by_contig: Dict[str, List[Tuple[int, int]]],
        contig_lengths: Dict[str, Tuple[int, str]],
    ) -> List[HitRecord]:
        """
        Assign a track row to every hit.

        Parameters
        ----------
        hits : list of HitRecord
            Hits sorted by contig then start.
        pairs_by_contig : dict
            Pair memberships per contig.
        contig_lengths : dict
            Contig name to ``(length, source)``.

        Returns
        -------
        list of HitRecord
            The same hits with `row` and `row_overflow` set.
        """
        by_contig: Dict[str, List[HitRecord]] = {}
        for hit in hits:
            by_contig.setdefault(hit.contig, []).append(hit)

        rows: Dict[int, int] = {}
        overflow: Set[int] = set()
        for contig, contig_hits in by_contig.items():
            length, _ = contig_lengths.get(contig, (1, 'inferred'))
            contig_rows, contig_overflow = stack_contig(
                contig_hits,
                pairs_by_contig.get(contig, []),
                length,
                max_rows=self.max_rows,
            )
            rows.update(contig_rows)
            overflow.update(contig_overflow)

        return [
            _replace_hit(
                hit,
                row=rows.get(hit.uid, 0),
                row_overflow=hit.uid in overflow,
            )
            for hit in hits
        ]

    def _build_contigs(
        self,
        hits: Sequence[HitRecord],
        elements: Sequence[ElementRecord],
        contig_lengths: Dict[str, Tuple[int, str]],
    ) -> List[ContigInfo]:
        """
        Summarise each contig that carries at least one hit.

        Parameters
        ----------
        hits : sequence of HitRecord
            Hits sorted by contig then start.
        elements : sequence of ElementRecord
            Element records.
        contig_lengths : dict
            Contig name to ``(length, source)``.

        Returns
        -------
        list of ContigInfo
            Ordered by descending pair count then descending hit count, so the
            contigs most worth looking at appear first.
        """
        slices: Dict[str, Tuple[int, int]] = {}
        n_hits: Dict[str, int] = {}
        n_rows: Dict[str, int] = {}
        for i, hit in enumerate(hits):
            lo, _ = slices.get(hit.contig, (i, i))
            slices[hit.contig] = (lo, i + 1)
            n_hits[hit.contig] = n_hits.get(hit.contig, 0) + 1
            n_rows[hit.contig] = max(n_rows.get(hit.contig, 1), hit.row + 1)

        n_elements: Dict[str, int] = {}
        element_bp: Dict[str, int] = {}
        for element in elements:
            contig = self._pair_contig[element.pair_id]
            n_elements[contig] = n_elements.get(contig, 0) + 1
            element_bp[contig] = element_bp.get(contig, 0) + element.length

        infos = []
        for name, (lo, hi) in slices.items():
            length, source = contig_lengths.get(name, (1, 'inferred'))
            infos.append(
                ContigInfo(
                    name=name,
                    length=length,
                    length_source=source,
                    n_hits=n_hits.get(name, 0),
                    n_pairs=n_elements.get(name, 0),
                    n_elements=n_elements.get(name, 0),
                    n_rows=n_rows.get(name, 1),
                    hit_lo=lo,
                    hit_hi=hi,
                    element_bp=element_bp.get(name, 0),
                )
            )

        # Hits keep their (contig, start) ordering, so the slices stay valid
        # only if the contig list is not reordered. Sort a copy of the display
        # order instead and let the browser use it for presentation.
        return sorted(infos, key=lambda c: c.hit_lo)

    def _build_models(self, hits: Sequence[HitRecord]) -> List[ModelInfo]:
        """
        Summarise each model that produced at least one hit.

        Parameters
        ----------
        hits : sequence of HitRecord
            The hits to be reported.

        Returns
        -------
        list of ModelInfo
            Ordered by model name.
        """
        by_model: Dict[str, List[HitRecord]] = {}
        for hit in hits:
            by_model.setdefault(hit.model, []).append(hit)

        infos = []
        for name in sorted(by_model):
            model_hits = by_model[name]
            coverages = [
                h.model_coverage for h in model_hits if h.model_coverage is not None
            ]
            n_paired = sum(1 for h in model_hits if h.pair_id is not None)
            clipped = sum(1 for h in model_hits if h.clip_left or h.clip_right)
            infos.append(
                ModelInfo(
                    name=name,
                    length=self.model_lengths.get(name),
                    n_hits=len(model_hits),
                    n_contigs=len({h.contig for h in model_hits}),
                    n_paired=n_paired,
                    n_unpaired=len(model_hits) - n_paired,
                    median_model_coverage=(
                        round(median(coverages), 4) if coverages else None
                    ),
                    frac_full_length=(
                        round(sum(1 for c in coverages if c >= 1.0) / len(coverages), 4)
                        if coverages
                        else None
                    ),
                    frac_clipped=round(clipped / len(model_hits), 4),
                )
            )
        return infos

    def _build_sequences(
        self,
        elements: Sequence[ElementRecord],
        embed: bool,
        max_bytes: int,
        fasta_path: str,
    ) -> SequenceStore:
        """
        Embed element sequences within the byte budget.

        Parameters
        ----------
        elements : sequence of ElementRecord
            Element records.
        embed : bool
            Whether embedding is enabled at all.
        max_bytes : int
            Byte budget.
        fasta_path : str
            Path shown for elements whose sequence was not embedded.

        Returns
        -------
        SequenceStore
            The embedded sequences and a count of what did not fit.

        Notes
        -----
        Shorter elements are embedded first. A handful of very long elements
        would otherwise consume the whole budget and leave the great majority
        of the copy buttons non-functional.
        """
        store = SequenceStore(embedded=embed, limit_bytes=max_bytes)
        if not embed or not self._element_seqs:
            store.omitted = len(elements) if not embed else 0
            return store

        candidates = sorted(
            (
                (len(self._element_seqs[e.pair_id]), e.pair_id)
                for e in elements
                if e.pair_id in self._element_seqs
            ),
        )
        total = 0
        for size, pair_id in candidates:
            if size > MAX_EMBEDDED_ELEMENT_BP or total + size > max_bytes:
                store.omitted += 1
                continue
            store.seq[pair_id] = self._element_seqs[pair_id]
            total += size
        store.total_bytes = total

        if store.omitted:
            where = f' They can be read from {fasta_path}.' if fasta_path else ''
            self.warnings.append(
                f'{store.omitted:,} element sequence(s) were too large to embed '
                f'and are not available from the copy button.{where}'
            )
        return store


def _pair_roles_safe(left: Any, right: Any, config: Any) -> Tuple[str, str]:
    """
    Assign terminus roles to a pair, tolerating a missing configuration.

    Parameters
    ----------
    left, right : namedtuple
        The pair's hits in genomic order.
    config : PairingConfig or None
        Pairing configuration.

    Returns
    -------
    tuple of (str, str)
        Roles of the left and right hits.
    """
    try:
        return _pair_roles(left, right, config)
    except Exception as exc:  # noqa: BLE001 - labelling must never abort a run
        logger.debug(f'Falling back to genomic terminus roles: {exc}')
        return 'left', 'right'


def _replace_hit(hit: HitRecord, **changes: Any) -> HitRecord:
    """
    Return a copy of a hit with some fields replaced.

    Parameters
    ----------
    hit : HitRecord
        The record to copy.
    **changes : any
        Fields to override.

    Returns
    -------
    HitRecord
        The updated copy. :class:`HitRecord` is frozen, so updates are copies.
    """
    return replace(hit, **changes)
