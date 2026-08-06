# HTML Reports

Both `tirmite pair` and `tirmite search` can write a single self-contained HTML
file summarising the run.

The file has no external references — no stylesheets, scripts, fonts or images
are fetched. It works from a local disk, over email, and years later from an
archive. Open it by double-clicking; no server is needed.

| Report | What it answers |
|---|---|
| [`tirmite pair --report`](#tirmite-pair-report) | Which termini paired into elements, how long they are, and how well each hit matches its model |
| [`tirmite search --report`](#tirmite-search-report) | What each query found, where two models claimed the same locus, and what clustering and filtering did |

Both share the same annotation tracks, alignment panels and table behaviour, so
what you learn from one applies to the other.

## `tirmite pair --report`

```bash
tirmite pair \
  --genome $GENOME \
  --nhmmer-file $NHMMERFILE \
  --hmm-file $HMMFILE \
  --pairing-map pairing_map.txt \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --gff \
  --report \
  --outdir MY_TIR_OUTPUT
```

This writes `MY_TIR_OUTPUT/tirmite_pair_report.html` (or
`<prefix>_tirmite_pair_report.html` when `--prefix` is set).

The report is a **single top-level file** covering every pairing group, even
when a `--pairing-map` sends the per-pair FASTA and GFF3 outputs into their own
subdirectories.

!!! tip "Set `--maxdist`"
    Without a distance limit, distant terminus hits can pair into spurious
    multi-hundred-kilobase "elements". The report will show them, and very large
    ones are held out of the embedded sequences with a warning — but the real
    fix is a `--maxdist` appropriate to your element family.

### What the pairing report contains

#### Annotation tracks

One track per sequence carrying at least one hit, drawn as a genome browser.

- **Swipe or drag sideways** to pan, or use the scrollbar under each track.
  **Scroll vertically** to zoom about the cursor. **Double-click** an
  annotation to frame it. With a track focused, the arrow keys pan and
  <kbd>+</kbd>/<kbd>-</kbd> zoom.
- **Hover** any annotation for a summary: model, coordinates, strand, e-value,
  score, model coverage, pairing group and partner.
- **Click** an annotation to open its details popup. For a paired terminus this
  has two tabs — the metadata, and the element sequence with copy, download and
  reverse-complement. The same popup opens from an element name in the tables.

Each annotation encodes several things at once:

| Cue | Meaning |
|-----|---------|
| Arrow direction | Strand of the hit |
| Fill colour | Pairing group (also named in the legend and tables) |
| Faded, dashed outline | Unpaired terminus |
| Thin connecting line | Links the two termini of one element |
| **Torn edge** | The hit did not match the whole model at that end |
| **Grey cap bar** | The hit is short because the contig ends there |

The last two are deliberately different. A partial model match and a hit running
off the end of an assembly are different findings, and the report never shows one
as the other.

Sequences without any terminus pair are hidden by default; the **Show sequences
without pairs** toggle reveals them, and the search box filters by name.

#### Terminus alignments

Every hit for a terminus model, stacked and aligned in model orientation.
Minus-strand hits are reverse-complemented before alignment, so every row reads
in the model's own direction and MAFFT is never asked to align a terminus
against its own reverse complement.

Two kinds of missing sequence are shown differently, for the same reason as
above:

- **Gaps** — model positions the hit did not match. This is the default: the
  panel shows only what the hit actually matched.
- **Grey** — with `--report-pad-model` (or `--padlen`, which implies it in
  `tirmite pair`), those positions instead show the sequence that sits beside
  the hit, greyed to mark it as unclaimed. Blank then means sequence that does
  not exist because the contig ended.

Padding is off by default because the extra bases are not evidence for the
model — they are whatever happens to lie next to the hit, and an alignment
invites reading them as part of the match. Turn it on when you want to see
whether a truncated hit stops at a real boundary or simply falls below the
score threshold part-way through the model: if the model's motif continues in
the grey, the hit was cut short by the search rather than by the sequence.

Hover a row for its hit; click to jump to that hit on its contig track. Each
panel can be downloaded as FASTA, either aligned (gaps kept) or unaligned, with
headers carrying each hit's genomic coordinates and strand.

Alignment uses **MAFFT** when it is on your `PATH`. Without it, hits are placed
by their model coordinates instead, which is exact wherever those coordinates
exist. The panel caption always states which was used.

!!! note "BED input"
    BED input carries no model coordinates, so coverage, torn edges and model
    anchoring are all unavailable. The report says so rather than guessing.

#### Distributions

Predicted element lengths, pairing outcome per group, model coverage, hits per
sequence and hit significance. All are print-ready inline SVG at Nature's
single-column measure.

The model coverage plot marks your `--mincov` setting, so the left edge of the
distribution is legible: everything below the line was removed before pairing.

!!! warning "Two different coverages"
    The plot shows the fraction of the model the *alignment* matched. The
    `--mincov` filter uses a different quantity — the hit's genomic span divided
    by the model length — which can exceed 1 when a hit contains insertions.
    Both are reported, labelled distinctly.

#### Hits and elements

Every terminus hit and every predicted element, filterable and sortable.
Coordinates are links: clicking one reveals that locus on its contig track.
Element names open the same details popup as clicking the annotation. Both
tables download in full as TSV.

Large tables render only their first rows, but that is a display limit alone —
the filter searches every row, the count states the true total, and the
download writes everything.

#### Statistics

Sortable tables covering pairing groups, terminus models, sequences and the
filters applied. Every value shown on a track or in a tooltip also appears
here, so nothing is reachable only by hovering.

### Options

| Option | Description |
|--------|-------------|
| `--report` | Write the HTML report |
| `--report-out` | Path for the report (implies `--report`) |
| `--report-title` | Heading shown at the top of the report |
| `--no-report-sequences` | Do not embed element sequences |
| `--report-max-seq-mb` | Budget for embedded sequences (default: 20 MB) |
| `--report-msa` | `auto`, `mafft`, `anchor` or `off` (default: `auto`) |
| `--report-pad-model` | Fill unmatched model positions with the neighbouring bases, drawn grey (implied by `--padlen`) |
| `--report-msa-max-rows` | Hits per alignment panel (default: 500) |
| `--report-max-hits` | Hits included in the report (default: 200 000) |
| `--report-max-rows` | Stacked annotation rows per sequence (default: 30) |

#### Keeping the report small

Element sequences are what make a report large. They are embedded shortest-first
until the `--report-max-seq-mb` budget is spent, and any single element over
50 kb is skipped regardless. Anything left out is counted in a banner at the top
of the report, and its coordinates and a `samtools faidx` command are still shown
in the element popup.

To make the file much smaller, drop the sequences entirely:

```bash
tirmite pair ... --report --no-report-sequences
```

For a rough sense of scale: a run producing 100 000 terminus hits makes an
18 MB report before any sequences are embedded, and takes about four seconds to
build on top of the pairing run itself.

#### When alignments are slow or unwanted

`--report-msa off` skips the alignment panels. `--report-msa anchor` keeps them
but never calls MAFFT, which is faster and fully deterministic.
`--report-msa mafft` fails the run's argument check if MAFFT is not installed,
rather than silently falling back.

### Interactions with other options

| Option | Effect on the report |
|--------|----------------------|
| `--nopairing` | A hits-only report: tracks, alignments and hit statistics, but no elements or pairs |
| `--no-elements` | Element sequences cannot be extracted, so they are not embedded. TIRmite warns and continues |
| `--prefix` | Applied to the report filename |
| `--keep-temp` | Unrelated; MAFFT's intermediate files live in the run's temporary directory and are cleaned up with it |

## `tirmite search --report`

```bash
tirmite search \
  --hmm MY_TIR.hmm \
  --genome $GENOME \
  --lengths-file model_lengths.txt \
  --cluster-map clusters.tsv \
  --outdir SEARCH_OUTPUT \
  --prefix myrun \
  --report
```

This writes `SEARCH_OUTPUT/myrun_report.html` alongside the usual
`myrun_hits.tab`.

A search report answers a different question from a pairing report — not
"which termini pair up" but "what did each query find, and where did the models
disagree". It carries no elements and no terminus pairs.

### Query termini

One row per query that produced a hit: hits, sequences, median model coverage,
how many matched their model in full, and how many were truncated by a contig
end. When a cluster map was used, "before merging" shows how many hits each
cluster absorbed from its components.

### Cross-matches

Loci where hits from two different groups overlap — two models claiming the
same sequence. Shown twice: a group-by-group matrix of how many loci each pair
contests, and a table of every one with spans and overlap length.

This is worth reading whether or not a later filter resolved the conflict. The
log reports only the first ten such overlaps; the report has them all.

### Shared hit loci

A model-by-model heatmap of how many loci each pair of queries both hit,
counted **before** clustering — the only point at which the model behind each
hit is still known.

Models are ordered by cluster, then by name, and a box is drawn around each
cluster's block. That makes the reading simple: **colour inside a box is
expected** — a cluster asserts that its members describe the same terminus —
while **colour outside a box is two unrelated models claiming the same
sequence**. The diagonal counts two hits of one model at the same locus, which
is redundancy rather than confusion between models.

!!! note "Why the diagonal is usually empty"
    A model's genuine hits land at distinct loci, so most models never overlap
    themselves and the diagonal reads zero. It is deliberately *not* the
    model's total hit count: totals are an order of magnitude larger than any
    overlap, so putting them on the same colour ramp would flatten the
    off-diagonal signal the figure exists to show, and would mix two different
    quantities on one scale. Totals appear in the clustered heatmap's row
    labels instead.

A second heatmap shows the same counts with the axes ordered by
**average-linkage clustering on the overlaps themselves**, with the tree drawn
above. Groups there come from the hits alone, so reading the two together is
the point: a block in the clustered view that crosses a box in the first one is
a cluster map at odds with its own evidence.

Under each leaf of that tree sits a **coloured symbol naming the cluster the
model was assigned to**, keyed to the right of the heatmap. That is what makes
the comparison readable at a glance: a branch whose leaves all carry the same
symbol is a cluster the hits agree with, and a branch carrying mixed symbols is
a grouping the hits support but the cluster map splits — or, where one symbol
turns up on distant branches, the reverse. Models with no cluster get a hollow
grey circle. Shape changes only once the eight palette hues are used up, so no
two clusters ever share both a colour and a shape.

Row labels on the clustered heatmap carry each model's **total hit count** in
brackets, which is what makes an off-diagonal count interpretable — five shared
loci out of six hits is a very different claim from five out of five hundred.

Distance is a Dice dissimilarity — twice the shared loci over the two models'
total hits — so a model that shares most of its few hits ranks as close as one
sharing many of many. Two models that never share a locus sit at distance 1 and
join last.

The exact counts are also in a table below the heatmaps, labelled by whether
each pair is in the same cluster, different clusters, or unclustered.

!!! note "Direct left/right partners are excluded"
    The two termini of one element often share a short stretch of identity, so
    they hit each other's loci. `tirmite search` already resolves that during
    filtering — the locus goes to whichever model scored better — and counting
    those pairs here would light up exactly the relationships the heatmap is
    meant to leave alone, burying the cross-family overlaps that matter. Pairs
    named in `--pairing-map` are therefore not counted.

!!! warning "A model in two clusters"
    If any query belongs to more than one cluster, the report says so in its
    banner. Such a run is ambiguous rather than merely untidy: merging assigns
    that model's hits to whichever cluster claimed them, and the heatmap has to
    place the model in one block or the other — it uses the first cluster
    alphabetically.

### Clustering

What each cluster was assembled from, how many hits merging collapsed, and —
when a pairing map was used — which model won each contested locus, from the
nested-hit and cross-model filter steps.

### Query alignments

Every hit for each query, stacked and aligned — the same panel as the pairing
report, with the same grey-versus-blank distinction between "the model was not
matched here" and "this sequence does not exist". Each panel downloads as FASTA,
aligned or unaligned.

This is often the fastest way to see that a query is matching two different
things: hits that belong together align cleanly, and ones that do not stand out
immediately as a block with a different pattern.

Panels are **off by default** here — they re-read sequence for every hit, which
a search run has no other reason to do. Turn them on with `--report-msa auto`
and supply `--genome` or `--genome-list`.

```bash
tirmite search ... --genome $GENOME --report --report-msa auto
```

!!! tip "Supply a genome even without panels"
    The report indexes it for true sequence lengths on the track axes. Without
    one, every axis is estimated from the hits and the report says so.

With `--genome-list`, panels read from every genome in the list, not just the
first. Sequence names are unique across a run's genomes, so a hit is always
read from the assembly it was found in.

### Filtering

Hits remaining after each pipeline stage, with the change at each step, so a
run that lost everything shows *where*.

### Options

| Option | Description |
|--------|-------------|
| `--report` | Write the HTML report |
| `--report-out` | Path for the report (implies `--report`) |
| `--report-title` | Heading shown at the top |
| `--report-msa` | `auto`, `mafft`, `anchor` or `off` (default: `off`) |
| `--report-pad-model` | Fill unmatched model positions with the neighbouring bases, drawn grey |
| `--report-msa-max-rows` | Hits per alignment panel (default: 500) |
| `--report-max-hits` | Hits included (default: 200 000) |
| `--report-max-rows` | Stacked annotation rows per sequence (default: 30) |

!!! note "Alignment panels are off by default here"
    They re-read sequence for every hit, which a search run has no other reason
    to do, and they need `--genome` or `--genome-list`. Turn them on with
    `--report-msa auto`.

### Model coverage after clustering

Cluster merging renames every hit to its cluster, but a lengths file is keyed
by component. The report resolves a cluster's model length from its components
when they all declare the same one, which is the usual case for a cluster of
variants of one terminus.

When components disagree there is no single denominator — a merged hit's
alignment coordinates come from whichever component scored best, and that
component is not recoverable from the merged row. Those clusters are left
without a length, coverage and jagged edges are unavailable for them, and the
report says so in its banner rather than inventing a number.

## If a report fails

Report generation runs **after** every other output is on disk. If it fails, the
error is logged, the run still exits successfully, and your GFF3, FASTA and
summary files are unaffected. A visualisation is not worth losing a completed
analysis over.

## Accessibility

The report is designed to be readable without relying on colour:

- Every pairing group is named in the legend and in the statistics tables.
- Nucleotide colours in the alignment panels were chosen by testing every
  four-hue subset of the report palette for separation under simulated
  protanopia and deuteranopia. Base letters are drawn whenever the cells are
  large enough.
- Every value shown graphically is also in a table.
- Light and dark modes are both supported and follow your system setting.
- Tracks are keyboard-navigable once focused.
