"""
Build report data from a ``tirmite search`` run.

A search run answers a different question from a pairing run: not "which
termini pair into elements" but "which model matched where, and where did the
models disagree". So the report carries no elements and no terminus pairs, and
its groups are clusters rather than pairing procedures.

Everything that does not depend on pairing -- contig lengths, contig-end
clipping, row stacking, per-contig and per-model summaries -- comes from
:class:`tirmite.report.collect.ReportAccumulatorBase`, so the annotation
tracks, the alignment panels and the statistics tables behave identically in
both reports.

The three things a search report exists to show, and where each comes from:

*Per-query statistics*
    The hit table itself, summarised per model.
*Cross-matches between terminus models*
    ``check_cross_cluster_overlaps``, which returns every contested locus.
*Clustering decisions*
    The cluster map, the per-component hit counts recorded before merging, and
    the pairing-map filter statistics that record which model won each
    contested locus.
"""

from datetime import datetime, timezone
import logging
from typing import Any, Callable, Dict, List, Optional, Sequence

from tirmite.report.collect import (
    ReportAccumulatorBase,
    _to_float,
    _to_int,
    model_truncation,
)
from tirmite.report.model import (
    HitColumns,
    HitRecord,
    PairingGroup,
    ReportData,
)
from tirmite.report.palette import group_colours

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = ['SearchReportAccumulator']


class SearchReportAccumulator(ReportAccumulatorBase):
    """
    Turn the results of a ``tirmite search`` run into report data.

    Parameters
    ----------
    tirmite_version : str
        Version string to record in the report.
    command : str
        Command line that produced the run.
    title : str
        Report heading.
    params : dict, optional
        Run parameters worth surfacing to the reader.
    model_lengths : dict, optional
        Model or query name to declared length.
    contig_length : callable, optional
        Called with a contig name, returning its length or None.
    cluster_map : dict, optional
        Cluster name to its component model names.
    pairing_map : dict, optional
        Left model name to right model name.
    max_hits : int, default 200000
        Cap on the number of hits included in the report.
    max_rows : int, default 30
        Cap on stacked annotation rows per contig.

    Notes
    -----
    Hits are grouped by cluster when a cluster map was used, and by model
    otherwise. Either way a group is "the thing whose colour you are looking
    at", which keeps the legend meaningful without the report having to know
    whether clustering ran.
    """

    def __init__(
        self,
        *,
        tirmite_version: str = '',
        command: str = '',
        title: str = 'TIRmite search report',
        params: Optional[Dict[str, Any]] = None,
        model_lengths: Optional[Dict[str, int]] = None,
        contig_length: Optional[Callable[[str], Optional[int]]] = None,
        cluster_map: Optional[Dict[str, List[str]]] = None,
        pairing_map: Optional[Dict[str, str]] = None,
        max_hits: int = 200_000,
        max_rows: int = 30,
    ) -> None:
        self.tirmite_version = tirmite_version
        self.command = command
        self.title = title
        self.params = dict(params) if params else {}
        self.model_lengths = dict(model_lengths) if model_lengths else {}
        self.contig_length = contig_length
        self.cluster_map = dict(cluster_map) if cluster_map else {}
        self.pairing_map = dict(pairing_map) if pairing_map else {}
        self.max_hits = max_hits
        self.max_rows = max_rows
        self.warnings: List[str] = []
        self._check_cluster_membership()
        self._resolve_cluster_lengths()

    def _check_cluster_membership(self) -> None:
        """
        Warn when a query model belongs to more than one cluster.

        Returns
        -------
        None
            Appends to `warnings`.

        Notes
        -----
        A model in two clusters makes the run ambiguous rather than merely
        untidy: merging assigns each of its hits to whichever cluster claims
        them, and the overlap heatmap has to place the model in one block or
        the other. The report names the offenders and states which cluster it
        used, so the resulting layout is explicable rather than arbitrary.
        """
        owners: Dict[str, List[str]] = {}
        for cluster, components in self.cluster_map.items():
            for component in components:
                owners.setdefault(component, []).append(cluster)

        shared = {
            model: sorted(clusters)
            for model, clusters in owners.items()
            if len(clusters) > 1
        }
        if not shared:
            return

        detail = '; '.join(
            f'{model} in {", ".join(clusters)}'
            for model, clusters in sorted(shared.items())
        )
        self.warnings.append(
            f'{len(shared)} query model(s) belong to more than one cluster '
            f'({detail}). Hits from a shared model are assigned to whichever '
            'cluster claimed them during merging, and the overlap heatmap '
            'places each such model under the first cluster alphabetically.'
        )

    def _resolve_cluster_lengths(self) -> None:
        """
        Give each cluster a model length where its components agree on one.

        Returns
        -------
        None
            Extends `model_lengths` in place.

        Notes
        -----
        Cluster merging renames every hit to its cluster, but the lengths file
        is keyed by component, so model coverage and the jagged edges that
        depend on it would be unavailable for every clustered run.

        A length is adopted only when every component declares the same one.
        When they differ there is no single denominator: a merged hit's
        alignment coordinates come from whichever component scored best, and
        that component is not recoverable from the merged row. Inventing a
        length -- the longest, say -- would silently understate coverage for
        hits from the shorter components, so those clusters are left without
        one and the report says so.
        """
        if not self.cluster_map or not self.model_lengths:
            return

        ambiguous = []
        for cluster, components in self.cluster_map.items():
            if cluster in self.model_lengths:
                continue
            lengths = {
                self.model_lengths[c] for c in components if c in self.model_lengths
            }
            if len(lengths) == 1:
                self.model_lengths[cluster] = lengths.pop()
            elif len(lengths) > 1:
                ambiguous.append(cluster)

        if ambiguous:
            self.warnings.append(
                f'Model coverage is unavailable for {len(ambiguous)} cluster(s) '
                f'({", ".join(sorted(ambiguous))}): their component queries '
                'declare different lengths, so a merged hit has no single '
                'model length to be measured against.'
            )

    # -- hit records ------------------------------------------------------

    def _hits_from_table(self, hit_table: Any) -> List[HitRecord]:
        """
        Convert the final hit table into report hit records.

        Parameters
        ----------
        hit_table : pandas.DataFrame
            The processed hit table, as written to ``*_hits.tab``.

        Returns
        -------
        list of HitRecord
            One per row, without layout or clipping, which are filled in
            during :meth:`finalise`.

        Notes
        -----
        Unlike the pairing report there is no hit index to read from: a search
        run never builds one. Everything comes from the table, which is also
        what the user gets on disk, so the two cannot disagree.
        """
        if hit_table is None or not len(hit_table):
            return []

        records: List[HitRecord] = []
        for idx, row in hit_table.iterrows():
            model = str(row['model'])
            start = _to_int(row['hitStart'])
            end = _to_int(row['hitEnd'])
            if start is None or end is None:
                logger.debug(f'Skipping a {model} hit with unreadable coordinates')
                continue
            # The importers guarantee start <= end, but a merged hit is built
            # by hand, so the order is re-established rather than assumed.
            if start > end:
                start, end = end, start

            hmm_start = _to_int(row.get('hmmStart'))
            hmm_end = _to_int(row.get('hmmEnd'))
            model_length = self.model_lengths.get(model)

            model_coverage = None
            if hmm_start is not None and hmm_end is not None and model_length:
                covered = max(0, hmm_end - hmm_start + 1)
                model_coverage = round(min(1.0, covered / model_length), 4)

            strand = str(row.get('strand', '+'))
            trunc_left, trunc_right = model_truncation(
                hmm_start, hmm_end, model_length, strand
            )

            records.append(
                HitRecord(
                    uid=int(idx) if isinstance(idx, int) else len(records),
                    model=model,
                    contig=str(row['target']),
                    start=start,
                    end=end,
                    strand=strand,
                    evalue=_to_float(row.get('evalue')),
                    score=_to_float(row.get('score')),
                    hmm_start=hmm_start,
                    hmm_end=hmm_end,
                    model_length=model_length,
                    model_coverage=model_coverage,
                    trunc_left=trunc_left,
                    trunc_right=trunc_right,
                    group_ids=[self._group_for(model)],
                )
            )
        return records

    def _group_for(self, model: str) -> str:
        """
        Return the group a model belongs to.

        Parameters
        ----------
        model : str
            Model name as it appears in the final hit table. After cluster
            merging this is already a cluster name.

        Returns
        -------
        str
            The cluster name if the model is one, or the model name itself.
        """
        if model in self.cluster_map:
            return model
        for cluster, components in self.cluster_map.items():
            if model in components:
                return cluster
        return model

    # -- finalisation -----------------------------------------------------

    def finalise(
        self,
        hit_table: Any,
        *,
        summary: Any = None,
        generated: Optional[str] = None,
        source: Any = None,
        tempdir: Optional[str] = None,
        msa_mode: str = 'off',
        msa_max_rows: int = 500,
        msa_max_cells: int = 2_000_000,
    ) -> ReportData:
        """
        Produce the finished report data.

        Parameters
        ----------
        hit_table : pandas.DataFrame
            The processed hit table.
        summary : SearchFilterSummary, optional
            Per-stage counts, cross-cluster overlaps and pairing-map filter
            statistics collected during the run.
        generated : str, optional
            ISO-8601 timestamp. Defaults to now, in UTC.
        source : SequenceSource, optional
            Sequence source for the terminus alignment panels.
        tempdir : str, optional
            Directory for MAFFT's intermediate files.
        msa_mode : {'off', 'auto', 'mafft', 'anchor'}, default 'off'
            Alignment strategy for the panels.
        msa_max_rows, msa_max_cells : int
            Panel caps.

        Returns
        -------
        ReportData
            Ready to render.
        """
        hits = self._hits_from_table(hit_table)
        hits = self._apply_cap(hits)

        contig_lengths = self._resolve_contig_lengths(hits)
        self._apply_clipping(hits, contig_lengths)

        hits = sorted(hits, key=lambda h: (h.contig, h.start, h.end, h.uid))
        # No pairs, so every hit stacks as its own interval.
        hits = self._stack(hits, {}, contig_lengths)

        contigs = self._build_contigs(hits, [], contig_lengths)
        contig_index = {c.name: i for i, c in enumerate(contigs)}
        models = self._build_models(hits)
        model_index = {m.name: i for i, m in enumerate(models)}

        groups = self._build_groups(hits)
        group_index = {g.group_id: i for i, g in enumerate(groups)}

        msa = self._build_msa(
            hits,
            groups,
            models,
            source=source,
            tempdir=tempdir,
            mode=msa_mode,
            max_rows=msa_max_rows,
            max_cells=msa_max_cells,
        )

        return ReportData(
            kind='search',
            tirmite_version=self.tirmite_version,
            generated=generated or datetime.now(timezone.utc).isoformat(),
            command=self.command,
            title=self.title,
            params=self.params,
            models=models,
            groups=groups,
            contigs=contigs,
            hits=HitColumns.from_records(
                hits, model_index, contig_index, group_index, {}
            ),
            elements=[],
            msa=msa,
            stats=self._build_stats(summary),
            warnings=list(self.warnings),
        )

    def _apply_cap(self, hits: List[HitRecord]) -> List[HitRecord]:
        """
        Apply the hit cap, keeping the most significant hits.

        Parameters
        ----------
        hits : list of HitRecord
            All hits.

        Returns
        -------
        list of HitRecord
            At most `max_hits` records.
        """
        if len(hits) <= self.max_hits:
            return hits

        ranked = sorted(
            hits,
            key=lambda h: (h.evalue if h.evalue is not None else float('inf'), h.uid),
        )
        dropped = len(hits) - self.max_hits
        self.warnings.append(
            f'{dropped:,} of {len(hits):,} hits were omitted to keep the report '
            f'responsive (--report-max-hits {self.max_hits:,}). The most '
            'significant hits were kept.'
        )
        return ranked[: self.max_hits]

    def _build_groups(self, hits: Sequence[HitRecord]) -> List[PairingGroup]:
        """
        Build the colour groups: clusters when clustering ran, else models.

        Parameters
        ----------
        hits : sequence of HitRecord
            The hits to be reported.

        Returns
        -------
        list of PairingGroup
            One per group, ordered by descending hit count so the busiest
            group takes the first palette slot.
        """
        counts: Dict[str, int] = {}
        members: Dict[str, set] = {}
        for hit in hits:
            for group_id in hit.group_ids:
                counts[group_id] = counts.get(group_id, 0) + 1
                members.setdefault(group_id, set()).add(hit.model)

        ordered = sorted(counts, key=lambda g: (-counts[g], g))
        colours = group_colours(ordered)

        groups = []
        for group_id in ordered:
            components = self.cluster_map.get(group_id)
            label = group_id
            if components:
                plural = '' if len(components) == 1 else 's'
                label = f'{group_id} ({len(components)} component{plural})'
            groups.append(
                PairingGroup(
                    group_id=group_id,
                    label=label,
                    left_model=group_id,
                    right_model=self.pairing_map.get(group_id, ''),
                    orientation='',
                    colour=colours[group_id][0],
                    colour_dark=colours[group_id][1],
                    n_pairs=0,
                    n_elements=0,
                    n_unpaired=counts[group_id],
                    hits_per_model={
                        model: sum(1 for h in hits if h.model == model)
                        for model in sorted(members[group_id])
                    },
                )
            )
        return groups

    def _build_stats(self, summary: Any) -> Dict[str, Any]:
        """
        Collect the search-specific statistics into the report payload.

        Parameters
        ----------
        summary : SearchFilterSummary or None
            Statistics collected during the run.

        Returns
        -------
        dict
            With keys 'stages', 'cross_matches', 'cluster_map',
            'hits_per_model_before_merge', 'nested_removed',
            'cross_model_removed' and 'excluded_not_in_map'.
        """
        stats: Dict[str, Any] = {
            'cluster_map': {k: list(v) for k, v in self.cluster_map.items()},
            'pairing_map': dict(self.pairing_map),
        }
        if summary is None:
            return stats

        stats['stages'] = [
            {'label': label, 'hits': count}
            for label, count in getattr(summary, 'stages', [])
        ]
        stats['cross_matches'] = list(getattr(summary, 'cross_cluster_overlaps', []))
        stats['model_overlaps'] = list(getattr(summary, 'model_overlaps', []))
        stats['hits_per_model_before_merge'] = dict(
            getattr(summary, 'hits_per_model_before_merge', {})
        )
        stats['excluded_not_in_map'] = dict(getattr(summary, 'excluded_not_in_map', {}))
        stats['nested_removed'] = {
            removed: dict(containers)
            for removed, containers in getattr(summary, 'nested_removed', {}).items()
        }
        # Tuple keys cannot be JSON object keys, so the pair becomes a record.
        stats['cross_model_removed'] = [
            {'removed': removed, 'kept': kept, 'hits': count}
            for (removed, kept), count in getattr(
                summary, 'cross_model_removed', {}
            ).items()
        ]
        return stats

    def _build_msa(
        self,
        hits: Sequence[HitRecord],
        groups: Sequence[PairingGroup],
        models: Sequence[Any],
        *,
        source: Any,
        tempdir: Optional[str],
        mode: str,
        max_rows: int,
        max_cells: int,
    ) -> List[Any]:
        """
        Build the terminus alignment panels, tolerating failure.

        Parameters
        ----------
        hits : sequence of HitRecord
            The hits in the report.
        groups : sequence of PairingGroup
            Colour groups.
        models : sequence of ModelInfo
            Models in the report.
        source : SequenceSource or None
            Sequence source.
        tempdir : str or None
            Directory for MAFFT's intermediate files.
        mode : str
            Alignment mode.
        max_rows, max_cells : int
            Panel caps.

        Returns
        -------
        list of MsaPanel
            The panels, or an empty list when they could not be built.
        """
        if mode == 'off' or source is None or tempdir is None:
            return []

        from tirmite.report.msa import build_msa_panels

        try:
            return list(
                build_msa_panels(
                    hits,
                    groups,
                    models,
                    source=source,
                    tmpdir=tempdir,
                    mode=mode,
                    max_rows=max_rows,
                    max_cells=max_cells,
                )
            )
        except Exception as exc:  # noqa: BLE001 - panels are not worth a run
            logger.warning(f'Could not build terminus alignments: {exc}')
            self.warnings.append(
                'Terminus alignment panels could not be built for this run.'
            )
            return []
