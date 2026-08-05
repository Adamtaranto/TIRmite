# HTML Reports

`tirmite pair --report` writes a single self-contained HTML file summarising the
whole run: annotation tracks for every sequence carrying a hit, stacked
alignments of each terminus model's hits, distribution plots and statistics
tables.

The file has no external references — no stylesheets, scripts, fonts or images
are fetched. It works from a local disk, over email, and years later from an
archive. Open it by double-clicking; no server is needed.

## Quick start

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

## What the report contains

### Annotation tracks

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

### Terminus alignments

Every hit for a terminus model, stacked and aligned in model orientation
(minus-strand hits are reverse-complemented).

Two kinds of missing sequence are shown differently, for the same reason as
above:

- **Grey** — a model position the hit did not match, although the genome
  continues there. The sequence exists; the alignment did not claim it.
- **Blank** — sequence that does not exist because the contig ended.

Hover a row for its hit; click to jump to that hit on its contig track. Each
panel can be downloaded as FASTA, either aligned (gaps kept) or unaligned, with
headers carrying each hit's genomic coordinates and strand.

Alignment uses **MAFFT** when it is on your `PATH`. Without it, hits are placed
by their model coordinates instead, which is exact wherever those coordinates
exist. The panel caption always states which was used.

!!! note "BED input"
    BED input carries no model coordinates, so coverage, torn edges and model
    anchoring are all unavailable. The report says so rather than guessing.

### Distributions

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

### Hits and elements

Every terminus hit and every predicted element, filterable and sortable.
Coordinates are links: clicking one reveals that locus on its contig track.
Element names open the same details popup as clicking the annotation. Both
tables download in full as TSV.

Large tables render only their first rows, but that is a display limit alone —
the filter searches every row, the count states the true total, and the
download writes everything.

### Statistics

Sortable tables covering pairing groups, terminus models, sequences and the
filters applied. Every value shown on a track or in a tooltip also appears
here, so nothing is reachable only by hovering.

## Options

| Option | Description |
|--------|-------------|
| `--report` | Write the HTML report |
| `--report-out` | Path for the report (implies `--report`) |
| `--report-title` | Heading shown at the top of the report |
| `--no-report-sequences` | Do not embed element sequences |
| `--report-max-seq-mb` | Budget for embedded sequences (default: 20 MB) |
| `--report-msa` | `auto`, `mafft`, `anchor` or `off` (default: `auto`) |
| `--report-msa-max-rows` | Hits per alignment panel (default: 500) |
| `--report-max-hits` | Hits included in the report (default: 200 000) |
| `--report-max-rows` | Stacked annotation rows per sequence (default: 30) |

### Keeping the report small

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

### When alignments are slow or unwanted

`--report-msa off` skips the alignment panels. `--report-msa anchor` keeps them
but never calls MAFFT, which is faster and fully deterministic.
`--report-msa mafft` fails the run's argument check if MAFFT is not installed,
rather than silently falling back.

## Interactions with other options

| Option | Effect on the report |
|--------|----------------------|
| `--nopairing` | A hits-only report: tracks, alignments and hit statistics, but no elements or pairs |
| `--no-elements` | Element sequences cannot be extracted, so they are not embedded. TIRmite warns and continues |
| `--prefix` | Applied to the report filename |
| `--keep-temp` | Unrelated; MAFFT's intermediate files live in the run's temporary directory and are cleaned up with it |

## If the report fails

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
