"""
Serialisable data model for TIRmite HTML reports.

Everything the browser needs is described here as plain frozen dataclasses, so
``dataclasses.asdict`` is enough to serialise a report and no custom JSON
encoder is involved. The same structures serve both the ``tirmite pair`` and
``tirmite search`` reports; ``ReportData.kind`` says which.

Hits are stored columnar rather than as a list of objects. At a few hundred
thousand hits the columnar form is several times smaller and much faster for
the browser to parse, and the JavaScript rehydrates a per-hit view on demand.
"""

from dataclasses import asdict, dataclass, field
import json
from typing import Any, Dict, List, Optional

__all__ = [
    'SCHEMA_VERSION',
    'ContigInfo',
    'ElementRecord',
    'FigureSpec',
    'HitColumns',
    'HitRecord',
    'MsaPanel',
    'MsaRow',
    'ModelInfo',
    'PairingGroup',
    'ReportData',
    'SequenceStore',
    'to_json',
]

# Bumped whenever the embedded JSON changes shape. The JavaScript checks it and
# refuses to render a document it does not understand, which is a clearer
# failure than a blank page.
SCHEMA_VERSION = 1

# Role codes as embedded in the JSON. Kept numeric to keep the payload small.
ROLE_NONE = 0
ROLE_LEFT = 1
ROLE_RIGHT = 2

ROLE_CODES = {None: ROLE_NONE, 'left': ROLE_LEFT, 'right': ROLE_RIGHT}


@dataclass(frozen=True)
class HitRecord:
    """
    One terminus model hit, before serialisation into columns.

    Attributes
    ----------
    uid : int
        Stable identifier for this hit, taken from the hit table's row index.
    model : str
        Name of the terminus model that produced the hit.
    contig : str
        Name of the sequence the hit lies on.
    start, end : int
        Inclusive one-based genomic coordinates, always with ``start <= end``
        regardless of strand.
    strand : str
        '+' or '-'.
    evalue, score : float, optional
        Hit statistics. None when the source format did not supply them or
        they were not numeric.
    hmm_start, hmm_end : int, optional
        Coordinates of the alignment on the model (or on the BLAST query).
        None for BED input, which carries no model coordinates.
    model_length : int, optional
        Declared length of the model, when known.
    model_coverage : float, optional
        Fraction of the model covered by the alignment, derived from
        `hmm_start`/`hmm_end`. This is the quantity the jagged-edge rendering
        uses.
    span_coverage : float, optional
        Genomic span of the hit divided by the model length, as computed by
        the ``--mincov`` filter. Can exceed 1.0 when the hit contains
        insertions, so it is reported separately rather than conflated with
        `model_coverage`.
    trunc_left, trunc_right : int, optional
        Base pairs of the model not covered at the lower and higher genomic
        coordinate ends respectively.
    clip_left, clip_right : bool
        True when the corresponding truncation is because the contig ends
        there. A hit that is short because the chromosome stopped is a
        different fact from one that partially matched its model, and the two
        are drawn differently.
    group_ids : list of str
        Pairing groups this hit participates in. Usually one, but a model that
        appears in several pairing-map rows can belong to more than one.
    role : str, optional
        'left' or 'right' terminus of its element, or None when unpaired or
        undeterminable.
    pair_id : str, optional
        Identifier of the element this hit belongs to, or None when unpaired.
    row : int
        Track row assigned by :mod:`tirmite.report.layout`.
    row_overflow : bool
        True when the track row cap was reached and this hit was forced onto
        the last row.
    """

    uid: int
    model: str
    contig: str
    start: int
    end: int
    strand: str
    evalue: Optional[float] = None
    score: Optional[float] = None
    hmm_start: Optional[int] = None
    hmm_end: Optional[int] = None
    model_length: Optional[int] = None
    model_coverage: Optional[float] = None
    span_coverage: Optional[float] = None
    trunc_left: Optional[int] = None
    trunc_right: Optional[int] = None
    clip_left: bool = False
    clip_right: bool = False
    group_ids: List[str] = field(default_factory=list)
    role: Optional[str] = None
    pair_id: Optional[str] = None
    row: int = 0
    row_overflow: bool = False

    @property
    def length(self) -> int:
        """
        Return the genomic span of the hit in base pairs.

        Returns
        -------
        int
            Inclusive length, so a single-base hit has length 1.
        """
        return self.end - self.start + 1


@dataclass
class HitColumns:
    """
    Column-oriented serialisation of a list of :class:`HitRecord`.

    Attributes
    ----------
    n : int
        Number of hits.
    uid, model_i, contig_i, start, end : list
        Per-hit values. `model_i` and `contig_i` index into the report's
        `models` and `contigs` lists.
    strand : str
        One character per hit, packed into a single string.
    evalue, score, hmm_start, hmm_end, model_cov, span_cov : list
        Per-hit values, with None where unavailable.
    trunc_l, trunc_r : list
        Model truncation in base pairs at each genomic end.
    clip_l, clip_r : list of int
        1 when the corresponding truncation is a contig end.
    row : list of int
        Track row.
    overflow : list of int
        1 when the hit was forced onto the last track row.
    group_ix : list of list of int
        Indices into the report's `groups` list.
    role : list of int
        0 none, 1 left, 2 right.
    pair_ix : list
        Index into the report's `elements` list, or None when unpaired.
    """

    n: int = 0
    uid: List[int] = field(default_factory=list)
    model_i: List[int] = field(default_factory=list)
    contig_i: List[int] = field(default_factory=list)
    start: List[int] = field(default_factory=list)
    end: List[int] = field(default_factory=list)
    strand: str = ''
    evalue: List[Optional[float]] = field(default_factory=list)
    score: List[Optional[float]] = field(default_factory=list)
    hmm_start: List[Optional[int]] = field(default_factory=list)
    hmm_end: List[Optional[int]] = field(default_factory=list)
    model_cov: List[Optional[float]] = field(default_factory=list)
    span_cov: List[Optional[float]] = field(default_factory=list)
    trunc_l: List[Optional[int]] = field(default_factory=list)
    trunc_r: List[Optional[int]] = field(default_factory=list)
    clip_l: List[int] = field(default_factory=list)
    clip_r: List[int] = field(default_factory=list)
    row: List[int] = field(default_factory=list)
    overflow: List[int] = field(default_factory=list)
    group_ix: List[List[int]] = field(default_factory=list)
    role: List[int] = field(default_factory=list)
    pair_ix: List[Optional[int]] = field(default_factory=list)

    @classmethod
    def from_records(
        cls,
        records: List[HitRecord],
        model_index: Dict[str, int],
        contig_index: Dict[str, int],
        group_index: Dict[str, int],
        element_index: Dict[str, int],
    ) -> 'HitColumns':
        """
        Build columns from hit records.

        Parameters
        ----------
        records : list of HitRecord
            Hits, already sorted by ``(contig, start)`` so that the contig
            slices recorded in :class:`ContigInfo` are contiguous.
        model_index : dict
            Model name to index in the report's `models` list.
        contig_index : dict
            Contig name to index in the report's `contigs` list.
        group_index : dict
            Group id to index in the report's `groups` list.
        element_index : dict
            Pair id to index in the report's `elements` list.

        Returns
        -------
        HitColumns
            The columnar form.
        """
        cols = cls(n=len(records))
        strands: List[str] = []
        for rec in records:
            cols.uid.append(rec.uid)
            cols.model_i.append(model_index[rec.model])
            cols.contig_i.append(contig_index[rec.contig])
            cols.start.append(rec.start)
            cols.end.append(rec.end)
            strands.append(rec.strand if rec.strand in '+-' else '.')
            cols.evalue.append(rec.evalue)
            cols.score.append(rec.score)
            cols.hmm_start.append(rec.hmm_start)
            cols.hmm_end.append(rec.hmm_end)
            cols.model_cov.append(rec.model_coverage)
            cols.span_cov.append(rec.span_coverage)
            cols.trunc_l.append(rec.trunc_left)
            cols.trunc_r.append(rec.trunc_right)
            cols.clip_l.append(1 if rec.clip_left else 0)
            cols.clip_r.append(1 if rec.clip_right else 0)
            cols.row.append(rec.row)
            cols.overflow.append(1 if rec.row_overflow else 0)
            cols.group_ix.append(
                [group_index[g] for g in rec.group_ids if g in group_index]
            )
            cols.role.append(ROLE_CODES.get(rec.role, ROLE_NONE))
            cols.pair_ix.append(element_index.get(rec.pair_id) if rec.pair_id else None)
        cols.strand = ''.join(strands)
        return cols


@dataclass(frozen=True)
class ElementRecord:
    """
    A predicted element spanning one pair of termini.

    Attributes
    ----------
    pair_id : str
        Unique identifier for the pair, of the form ``'<group_id>:<n>'``.
    group_i : int
        Index into the report's `groups` list.
    element_id : str
        Element name as written to the FASTA and GFF3 outputs.
    contig_i : int
        Index into the report's `contigs` list.
    start, end : int
        Inclusive one-based genomic coordinates of the whole element, from the
        outer edge of one terminus to the outer edge of the other.
    strand : str
        Strand as recorded by the extraction step. The extracted sequence is
        always plus-strand regardless.
    length : int
        Element length in base pairs.
    left_uid, right_uid : int
        Uids of the two terminus hits, in genomic order.
    inner_distance : int
        Base pairs between the inner edges of the two termini. Negative when
        the termini overlap.
    """

    pair_id: str
    group_i: int
    element_id: str
    contig_i: int
    start: int
    end: int
    strand: str
    length: int
    left_uid: int
    right_uid: int
    inner_distance: int


@dataclass(frozen=True)
class PairingGroup:
    """
    One pairing procedure, i.e. one row of the pairing map.

    Attributes
    ----------
    group_id : str
        Unique id, ``'<left_model>__<right_model>'``. Does not collapse for
        symmetric pairings, so it stays unique when a run mixes symmetric and
        asymmetric rows sharing a model.
    label : str
        Human-facing label.
    left_model, right_model : str
        Models forming the pairing. Equal for a symmetric pairing.
    orientation : str
        Orientation code the pairing used, e.g. 'F,R'.
    colour : str
        Hex colour used for this group throughout the report on a light
        surface.
    colour_dark : str
        The same hue re-stepped for a dark surface. Selected rather than
        derived: inverting a light-mode colour does not reliably stay
        colour-vision-safe.
    n_pairs, n_elements, n_unpaired : int
        Outcome counts for the group.
    hits_per_model : dict of str to int
        Retained hits per model in this pairing.
    """

    group_id: str
    label: str
    left_model: str
    right_model: str
    orientation: str
    colour: str
    colour_dark: str = ''
    n_pairs: int = 0
    n_elements: int = 0
    n_unpaired: int = 0
    hits_per_model: Dict[str, int] = field(default_factory=dict)


@dataclass(frozen=True)
class ModelInfo:
    """
    Per-model summary.

    Attributes
    ----------
    name : str
        Model name.
    length : int, optional
        Declared model length, when known.
    n_hits, n_contigs, n_paired, n_unpaired : int
        Counts across the whole run.
    median_model_coverage : float, optional
        Median fraction of the model covered, over hits with model
        coordinates.
    frac_full_length : float, optional
        Fraction of hits covering the entire model.
    frac_clipped : float, optional
        Fraction of hits truncated by a contig end rather than by a partial
        model match.
    """

    name: str
    length: Optional[int] = None
    n_hits: int = 0
    n_contigs: int = 0
    n_paired: int = 0
    n_unpaired: int = 0
    median_model_coverage: Optional[float] = None
    frac_full_length: Optional[float] = None
    frac_clipped: Optional[float] = None


@dataclass(frozen=True)
class ContigInfo:
    """
    One contig carrying at least one hit.

    Attributes
    ----------
    name : str
        Sequence name.
    length : int
        Contig length in base pairs.
    length_source : str
        'source' when read from the genome or BLAST database, 'inferred' when
        estimated from the hits because no sequence source was available.
    n_hits, n_pairs, n_elements : int
        Counts on this contig.
    n_rows : int
        Number of stacked annotation rows the track needs.
    hit_lo, hit_hi : int
        Half-open slice of the hit columns belonging to this contig. Hits are
        sorted by contig then start, so the slice is contiguous and the
        browser can binary-search within it.
    element_bp : int
        Base pairs covered by predicted elements.
    """

    name: str
    length: int
    length_source: str
    n_hits: int = 0
    n_pairs: int = 0
    n_elements: int = 0
    n_rows: int = 1
    hit_lo: int = 0
    hit_hi: int = 0
    element_bp: int = 0


@dataclass(frozen=True)
class MsaRow:
    """
    One row of a terminus alignment panel.

    Attributes
    ----------
    uid : int
        Uid of the hit this row shows.
    group_i : int, optional
        Index into the report's `groups` list, or None for an unpaired hit
        outside every group.
    role : str, optional
        'left' or 'right', or None when unpaired.
    seq : str
        Aligned sequence, one character per column of the panel. '-' marks a
        position that does not exist, either an alignment gap or sequence
        beyond the end of the contig.
    pad : list of [int, int, str]
        Runs of padded columns as ``[start_column, length, kind]``. Kind 'm'
        means the model was not covered there although the genome continues,
        and is drawn grey; kind 'g' means the contig ended, and is drawn as
        empty space.
    """

    uid: int
    group_i: Optional[int]
    role: Optional[str]
    seq: str
    pad: List[List[Any]] = field(default_factory=list)


@dataclass(frozen=True)
class MsaPanel:
    """
    Stacked alignment of every hit to one terminus model.

    Attributes
    ----------
    model : str
        Terminus model the panel is for.
    aligner : str
        'mafft' when MAFFT aligned the rows, 'anchor' when they were placed by
        their model coordinates, or 'left' when neither was possible and rows
        are simply left-aligned.
    n_cols : int
        Panel width in columns.
    n_rows_shown, n_rows_total : int
        Rows displayed and rows available, which differ when the row cap
        applied.
    rows : list of MsaRow
        The rows themselves.
    note : str, optional
        Caveat to show in the panel caption.
    """

    model: str
    aligner: str
    n_cols: int
    n_rows_shown: int
    n_rows_total: int
    rows: List[MsaRow] = field(default_factory=list)
    note: Optional[str] = None


@dataclass
class SequenceStore:
    """
    Element sequences embedded for the copy-to-clipboard popup.

    Attributes
    ----------
    embedded : bool
        False when sequence embedding was disabled or unavailable, in which
        case the popup shows coordinates and a FASTA path instead.
    limit_bytes : int
        Byte budget that was applied.
    total_bytes : int
        Bytes actually embedded.
    omitted : int
        Number of elements whose sequence did not fit the budget.
    seq : dict
        Mapping of pair id to plus-strand sequence.
    """

    embedded: bool = False
    limit_bytes: int = 0
    total_bytes: int = 0
    omitted: int = 0
    seq: Dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class FigureSpec:
    """
    A static figure embedded as inline SVG.

    Attributes
    ----------
    id : str
        Slug used as the DOM id and as the prefix that keeps this figure's
        internal SVG ids from colliding with another figure's.
    title : str
        Figure heading.
    caption : str
        Explanatory caption.
    svg : str
        The SVG markup, with its internal ids already namespaced.
    """

    id: str
    title: str
    caption: str
    svg: str = ''


@dataclass
class ReportData:
    """
    Everything one HTML report needs.

    Attributes
    ----------
    kind : str
        'pair' or 'search'.
    schema_version : int
        Version of this structure, checked by the browser.
    tirmite_version : str
        Version of TIRmite that produced the report.
    generated : str
        ISO-8601 UTC timestamp.
    command : str
        The command line that produced the run.
    title : str
        Report heading.
    params : dict
        Run parameters worth surfacing.
    filter_stats : dict
        Hit-filtering statistics.
    models, groups, contigs, elements : list
        The report's entities.
    hits : HitColumns
        Columnar hit data.
    sequences : SequenceStore
        Embedded element sequences.
    msa : list of MsaPanel
        Terminus alignment panels.
    stats : dict
        Pre-computed tables for the statistics section.
    figures : list of FigureSpec
        Static figures.
    warnings : list of str
        Caveats to show prominently, e.g. that a cap truncated the data.
    """

    kind: str = 'pair'
    schema_version: int = SCHEMA_VERSION
    tirmite_version: str = ''
    generated: str = ''
    command: str = ''
    title: str = 'TIRmite report'
    params: Dict[str, Any] = field(default_factory=dict)
    filter_stats: Dict[str, Any] = field(default_factory=dict)
    models: List[ModelInfo] = field(default_factory=list)
    groups: List[PairingGroup] = field(default_factory=list)
    contigs: List[ContigInfo] = field(default_factory=list)
    hits: HitColumns = field(default_factory=HitColumns)
    elements: List[ElementRecord] = field(default_factory=list)
    sequences: SequenceStore = field(default_factory=SequenceStore)
    msa: List[MsaPanel] = field(default_factory=list)
    stats: Dict[str, Any] = field(default_factory=dict)
    figures: List[FigureSpec] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    def as_dict(self) -> Dict[str, Any]:
        """
        Return the report as plain dictionaries and lists.

        Returns
        -------
        dict
            The full structure, with figures' SVG markup removed: the SVG is
            written into the document body, not into the JSON payload, so
            embedding it twice would double the file size.
        """
        data = asdict(self)
        data['figures'] = [
            {'id': f['id'], 'title': f['title'], 'caption': f['caption']}
            for f in data['figures']
        ]
        return data

    def to_json(self) -> str:
        """
        Serialise the report to JSON safe for embedding in a script element.

        Returns
        -------
        str
            Compact JSON in which every character that could terminate the
            enclosing script element, or break a JavaScript string literal, is
            escaped.
        """
        return to_json(self.as_dict())


def to_json(data: Any) -> str:
    """
    Serialise a value to JSON safe for embedding in an HTML script element.

    Parameters
    ----------
    data : any
        JSON-serialisable value.

    Returns
    -------
    str
        Compact JSON.

    Notes
    -----
    Model and contig names come from user-supplied files, so a name containing
    ``</script>`` would otherwise close the element early and inject the rest
    of the payload as markup. Escaping ``<`` prevents that. U+2028 and U+2029
    are valid in JSON strings but terminate a line in older JavaScript parsers,
    so they are escaped too.
    """
    text = json.dumps(data, separators=(',', ':'), allow_nan=False, default=str)
    return (
        text.replace('<', '\\u003c')
        .replace('\u2028', '\\u2028')
        .replace('\u2029', '\\u2029')
    )
