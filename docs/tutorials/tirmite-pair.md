# Using tirmite pair

`tirmite pair` takes precomputed nhmmer or BLAST search results and applies the TIRmite pairing algorithm to identify valid transposon termini pairs. It outputs paired elements, hit sequences, and GFF3 annotations. When `--flank-len` is set, flanking sequences are extracted. Add `--insertion-site` to enable reconstruction of the original **target site** for each paired element.

## Overview

```mermaid
flowchart TD
    A[nhmmer .tab\nor BLAST .tab] --> B[Load hits]
    B --> C[Filter by e-value\nand model coverage]
    C --> C2{--max-offset\nset?}
    C2 -->|Yes| C3[Anchor filter:\nremove hits not reaching\nouter model edge]
    C2 -->|No| D
    C3 --> D{Pairing map\nprovided?}
    D -->|Yes| E[Split hits by model pair\nRun independent pairing\nfor each combination]
    D -->|No| F[Run single pairing\nprocedure]
    E --> G[Pairing algorithm\nper model pair]
    F --> G
    G --> H[Output paired elements,\nhit sequences, GFF3]
    H --> I{--flank-len\nset?}
    I -->|Yes| J[Extract external flanking\nsequences per terminus]
    I -->|No| K[Done]
    J --> L{--insertion-site\nset?}
    L -->|No| K
    L -->|Yes| M{--tsd-in-model?}
    M -->|No\nTSD is in flank| N[Trim TSD from\nflank boundary]
    M -->|Yes\nTSD is in hit| O[Extract TSD from\nhit inner boundary]
    N --> P[Write target_sites.fasta\nand interleaved_flanks.fasta]
    O --> P
    P --> K
```

## Input Types

`tirmite pair` accepts three types of hit input:

| Input flag | Format | Description |
|-----------|--------|-------------|
| `--nhmmer-file` | nhmmer tabular (`--tblout`) | Results from nhmmer search |
| `--blast-file` | BLAST tabular (`-outfmt 6`) | Results from blastn or megablast |
| Pre-filtered results | Either format | Output from `tirmite search` can be passed directly |

### When is the HMM file required?

Provide `--hmm-file` when using nhmmer input and you want TIRmite to calculate hit **coverage** (the fraction of the model that the hit spans). This requires knowing the model length, which is read directly from the HMM file.

If you don't have the HMM file, you can provide `--lengths-file` with a tab-delimited file mapping model names to their lengths.

### When is the lengths file required?

Provide `--lengths-file` when:

- Your input contains hits from **multiple models** but you don't have HMM files for all of them
- You are using BLAST input and want coverage filtering (BLAST doesn't record query lengths in the hit table)
- You have pre-filtered results from `tirmite search` where the original HMM is no longer accessible

Format of the lengths file (tab-delimited, one model per line):

```
# model_name    model_length
MY_TIR    150
LEFT_TERMINUS    120
RIGHT_TERMINUS    135
```

## The Pairing Algorithm

The pairing algorithm identifies the most parsimonious set of valid terminus pairs. The core logic iterates through unpaired hits and attempts to find a valid partner:

```mermaid
flowchart TD
    A[All filtered hits] --> B[For each hit: find candidate partners]
    B --> C{Candidate on same\nsequence in correct\norientation within\nmaxdist?}
    C -->|No| D[Mark as unpaired]
    C -->|Yes| E[Rank candidates\nby proximity]
    E --> F[Pair reciprocal\ntop candidates]
    F --> G{Unpaired hits\nremaining?}
    G -->|Yes| H[Find nearest unpaired\ncandidate]
    H --> I{Reciprocal?}
    I -->|Yes| J[Pair them]
    I -->|No| K{2nd-order\nreciprocal?}
    K -->|Yes| J
    K -->|No| L[Leave unpaired\nfor this iteration]
    J --> G
    L --> M{Iterations without\nnew pairing\n≥ stable-reps?}
    M -->|No| G
    M -->|Yes| N[End: report paired\nand unpaired hits]
```

### Algorithm steps

1. Use nhmmer (or BLAST) to query genome with termini models/sequences.
2. Import all hits under `--maxeval` threshold.
3. For each significant terminus match, identify candidate partners, where:
   - Hit is on the same sequence
   - Hit is in correct relative orientation
   - Distance is ≤ `--maxdist`
   - Hit length is ≥ (model/query length × `--mincov`)
4. Rank candidate partners by distance downstream of positive-strand hits, and upstream of negative-strand hits.
5. Pair reciprocal top candidate hits.
6. For unpaired hits, find nearest unpaired candidate partner and check for reciprocity.
7. If the first unpaired candidate is non-reciprocal, check for 2nd-order reciprocity (is the outbound top-candidate of the current candidate reciprocal?).
8. Iterate steps 6–7 until all termini hits are paired OR number of iterations without a new pairing exceeds `--stable-reps`.

### Orientation settings

The `--orientation` flag controls which strand combinations are considered valid pairs:

| Setting | Meaning | Use case |
|---------|---------|----------|
| `F,R` | Forward hit paired with Reverse hit | TIR elements (same model) |
| `F,F` | Forward hit paired with Forward hit | LTR retrotransposons |
| `R,F` | Reverse hit paired with Forward hit | Rarely used |
| `R,R` | Reverse hit paired with Reverse hit | Rarely used |

For **asymmetric elements** using a pairing map, the orientation is applied independently for each left/right model pair defined in the pairing map.

```mermaid
flowchart LR
    subgraph "F,R pairing (TIRs)"
        direction LR
        L1[→ LEFT hit\nForward strand] --- E1[Element] --- R1[← RIGHT hit\nReverse strand]
    end
    subgraph "F,F pairing (LTRs)"
        direction LR
        L2[→ LEFT hit\nForward strand] --- E2[Element] --- R2[→ RIGHT hit\nForward strand]
    end
```

## Examples

### Basic usage with nhmmer output

```bash
GENOME="genome.fa"
NHMMERFILE="MY_TIR_nhmmer_hits.tab"
HMMFILE="MY_TIR.hmm"

tirmite pair \
  --genome $GENOME \
  --nhmmer-file $NHMMERFILE \
  --hmm-file $HMMFILE \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --outdir MY_TIR_OUTPUT \
  --gff-out
```

### With BLAST input

```bash
tirmite pair \
  --genome $GENOME \
  --blast-file MY_TIR_blast_hits.tab \
  --query-len 100 \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --outdir MY_TIR_BLAST_OUTPUT \
  --gff-out
```

### With pre-filtered results from `tirmite search`

```bash
tirmite pair \
  --genome $GENOME \
  --blast-file SEARCH_OUTPUT/merged_hits.tab \
  --lengths-file model_lengths.txt \
  --pairing-map pairing_map.txt \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --outdir PAIR_OUTPUT \
  --gff-out
```

### Using a BLAST database for sequence extraction

If your BLAST database was created with `-parse_seqids`, TIRmite can extract sequences directly from the database:

```bash
# Create BLAST database with parsed sequence IDs
makeblastdb -in $GENOME -dbtype nucl -out genome_db -parse_seqids

tirmite pair \
  --blastdb genome_db \
  --blast-file MY_TIR_blast_hits.tab \
  --query-len 100 \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --outdir MY_TIR_OUTPUT
```

### End-to-end: flank extraction and insertion site reconstruction

This example demonstrates a full run that extracts external flanking sequences for confirmed paired termini and reconstructs the original insertion site (target site) for each element:

```bash
GENOME="genome.fa"
NHMMERFILE="MY_TIR_nhmmer_hits.tab"
HMMFILE="MY_TIR.hmm"

tirmite pair \
  --genome $GENOME \
  --nhmmer-file $NHMMERFILE \
  --hmm-file $HMMFILE \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --gff-report all \
  --gff-out \
  --flanks-paired \
  --flank-len 30 \
  --insertion-site \
  --tsd-length 8 \
  --outdir MY_TIR_OUTPUT
```

The `--flanks-paired` flag writes **external** flanks only for confirmed paired termini (upstream of left terminus hit, downstream of right terminus hit). Adding `--insertion-site` then uses those flanks to reconstruct the pre-insertion target site for each element.

Key output files in `MY_TIR_OUTPUT/`:

| File | Description |
|------|-------------|
| `MY_TIR_paired_left_term_hits_N.fasta` | Left terminus hit sequences for all pairs |
| `MY_TIR_paired_right_term_hits_N.fasta` | Right terminus hit sequences for all pairs (reverse complemented) |
| `MY_TIR_elements_N.fasta` | Full element sequences (left terminus to right terminus) |
| `MY_TIR_paired_left_flank_N.fasta` | External flanks upstream of paired left terminus hits |
| `MY_TIR_paired_right_flank_N.fasta` | External flanks downstream of paired right terminus hits |
| `MY_TIR_target_sites.fasta` | Reconstructed pre-insertion target sites with TSD metadata |
| `MY_TIR_interleaved_flanks.fasta` | Left/right flank pairs shown interleaved with TSD highlighted |
| `tirmite_pair_report.gff3` | GFF3 annotation of elements and termini |

To additionally write flanks for **all** hits (including unpaired), add `--flanks`:

```bash
tirmite pair \
  --genome $GENOME \
  --nhmmer-file $NHMMERFILE \
  --hmm-file $HMMFILE \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --gff-report all \
  --gff-out \
  --flanks \
  --flanks-paired \
  --flank-len 30 \
  --insertion-site \
  --tsd-length 8 \
  --outdir MY_TIR_OUTPUT
```

This produces `MY_TIR_left_flank_N.fasta` and `MY_TIR_right_flank_N.fasta` for all hits (paired and unpaired) in addition to the paired-only flank files.

### End-to-end: asymmetric pairing with multiple models

For elements with distinct left and right terminus models, use `--pairing-map` with an asymmetric map and `--tsd-length-map` if TSD lengths differ by model pair:

```bash
# pairing_map.txt
# HEAD_TIR  TAIL_TIR
# (asymmetric: different models for left and right termini)

# tsd_lengths.tsv
# HEAD_TIR  TAIL_TIR  8

tirmite pair \
  --genome $GENOME \
  --nhmmer-file multi_model_hits.tab \
  --lengths-file model_lengths.txt \
  --pairing-map pairing_map.txt \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --gff-report all \
  --gff-out \
  --flanks-paired \
  --flank-len 30 \
  --insertion-site \
  --tsd-length-map tsd_lengths.tsv \
  --outdir MULTI_MODEL_OUTPUT
```

For each model pair defined in the pairing map, TIRmite creates a sub-directory (e.g., `MULTI_MODEL_OUTPUT/HEAD_TIR_TAIL_TIR/`) containing all output files for that pair.

## Extracting Flanking Regions

Use `--flanks` to extract the genomic sequence immediately **outside** each terminus hit – upstream of the left terminus and downstream of the right terminus. These external flanks are useful for:

- Examining Target Site Duplications (TSDs)
- Checking whether element boundaries are correctly identified
- Providing flanking context for downstream analysis

The flank length defaults to 50 bp but can be adjusted with `--flank-len`.

!!! note "Difference between `--padlen` and `--flank-len`"
    `--padlen` pads the hit or element sequence itself (i.e. adds bases to both sides of the hit/element in the output FASTA). `--flank-len` extracts the sequence that lies *outside* the element boundary – the true external flank. Use `--flanks` when you want to study TSDs or the sequence context beyond the termini.

### Flank output files

When `--flanks` is set, TIRmite writes the following FASTA files to the output directory:

| Filename pattern | Contents |
|-----------------|----------|
| `{prefix}{model}_left_flank_{count}.fasta` | External flanks for all left terminus hits |
| `{prefix}{model}_right_flank_{count}.fasta` | External flanks for all right terminus hits |

When `--flanks-paired` is set (with or without `--flanks`), TIRmite also writes:

| Filename pattern | Contents |
|-----------------|----------|
| `{prefix}{model}_paired_left_flank_{count}.fasta` | External flanks for **paired** left terminus hits only (element ID in header) |
| `{prefix}{model}_paired_right_flank_{count}.fasta` | External flanks for **paired** right terminus hits only (element ID in header) |

When `--gff-report all` or `--gff-report paired` is set, TIRmite also writes the terminus hit sequences themselves to two separate files per model:

| Filename pattern | Contents |
|-----------------|----------|
| `{prefix}{model}_paired_left_term_hits_{count}.fasta` | Left terminus hit sequences for all confirmed pairs |
| `{prefix}{model}_paired_right_term_hits_{count}.fasta` | Right terminus hit sequences (reverse complemented) for all confirmed pairs |

The **paired-only** flank files include the element identifier (`Element_N`) in each sequence header, making it easy to link flanks back to specific annotated elements.

Example FASTA header for a paired flank:

```
>Element_1_MY_TIR_1_L  chr1:1000-1019(+)
```

!!! warning "Symmetric same-strand orientations (F,F or R,R)"
    For symmetric same-strand orientations, a single model is present at both ends of an element in the same orientation, so individual unpaired hits cannot be classified as "left" or "right". When `--flanks` is used with these orientations, both left and right flanks are extracted for all hits and a warning is raised advising use of `--flanks-paired` instead.

### Offset correction

The external flank starts at the genomic coordinate corresponding to model position 1 (for a left terminus) or model position `model_len` (for a right terminus). If the hit alignment does not reach that outer model edge, TIRmite corrects the flank coordinate by the gap (offset).

Use `--flank-max-offset` to skip hits where this offset is too large (i.e. the hit does not extend close enough to the terminus edge):

```bash
tirmite pair \
  --genome $GENOME \
  --nhmmer-file $NHMMERFILE \
  --hmm-file $HMMFILE \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --flanks \
  --flank-len 20 \
  --flank-max-offset 5 \
  --outdir MY_TIR_OUTPUT \
  --gff-out
```

### Paired-only flank extraction

Use `--flanks-paired` to write outer flanks only for termini that were assigned to pairs. If `--flanks` is also set, all-hit flanks are written first, then paired flanks are written to separate files (reusing already-extracted sequences to avoid redundant extraction):

```bash
tirmite pair \
  --genome $GENOME \
  --nhmmer-file $NHMMERFILE \
  --hmm-file $HMMFILE \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --flanks \
  --flanks-paired \
  --flank-len 20 \
  --outdir MY_TIR_OUTPUT \
  --gff-out
```

!!! tip
    When both `--flanks` and `--flanks-paired` are used together, the `_left_flank_` and `_right_flank_` files contain all hits (paired + unpaired) while the `_paired_left_flank_` and `_paired_right_flank_` files contain only paired hits with element IDs in the headers.

## Reconstructing Target Sites

When `--flanks` or `--flanks-paired` is set together with `--insertion-site`, TIRmite reconstructs the original **target site** for each paired element. The target site is the genomic sequence that was present at the insertion location before the transposon arrived — a TSD-flanked sequence that appears once in empty sites and duplicated around inserted elements. Use `--tsd-length` (or `--tsd-length-map`) and optionally `--tsd-in-model` to configure how the TSD is handled during reconstruction.

When processing multiple terminus model pairs via `--pairing-map`, one output FASTA is written per pair of models. Target sites are written as single-line non-wrapped FASTA.

### How the reconstruction works

Each transposon insertion results in one copy of the TSD at each end of the element. When the flanks are joined the duplicate TSD must be reduced to a single copy:

```
Inserted element
←[left_flank][TSD]→[LEFT_TERMINUS]…[RIGHT_TERMINUS]←[TSD][right_flank]→
                 ↑ one copy trimmed here
```

#### When `--tsd-in-model` is NOT set (TSD is outside the model)

The TSD sequence appears as the innermost `n` bases of the extracted flanks — immediately adjacent to the terminus hit boundary. TIRmite reads the TSD from the flank boundary, computes the Hamming distance between the two copies, and trims one before joining.

```
Genomic layout:
  ←[upstream]←[TSD_L]←[LEFT_HIT]……[RIGHT_HIT]→[TSD_R]→[downstream]→

  left_flank  = [upstream] + [TSD_L]   (TSD at right/inner end of left flank)
  right_flank = [TSD_R] + [downstream] (TSD at left/inner end of right flank)

Reconstruction: left_flank + right_flank[tsd_length:]
             = [upstream] + [TSD_L] + [downstream]
```

#### When `--tsd-in-model` IS set (TSD is inside the model)

The TSD is encoded within the terminus HMM model itself (at the element-facing end). The extracted flanks are clean sequences with no TSD. TIRmite uses the hit coordinates and HMM alignment positions to locate the TSD within the hit, extracts it directly, and inserts one copy between the clean flanks.

```
Genomic layout:
  ←[upstream]←[LEFT_HIT][TSD_L]……[TSD_R][RIGHT_HIT]→[downstream]→

  left_flank  = [upstream]             (clean, no TSD)
  right_flank = [downstream]           (clean, no TSD)
  TSD_L extracted from LEFT_HIT inner boundary

Reconstruction: left_flank + TSD_L + right_flank
             = [upstream] + [TSD_L] + [downstream]
```

### Flank coordinate and offset corrections

The external flank boundary is aligned to the **outer edge of the terminus model** (model position 1 for left terminus, model position `N` for right terminus). If the hit alignment does not reach the model outer edge, the flank coordinate is shifted inward by the number of uncovered model positions (the offset).

Similarly, for `tsd_in_model=True`, the TSD inner boundary coordinates are computed from the hit start/end and HMM alignment positions (`hmmStart`/`hmmEnd`), accounting for any offset between the alignment end and the model inner edge. This ensures that orientation, strand, and hit offset are all correctly accounted for in the reconstruction.

### TSD length map for multiple model pairs

When processing hits from multiple models via `--pairing-map`, use `--tsd-length-map` instead of a single `--tsd-length`. The map is a three-column TSV file:

```
# tsd_lengths.tsv
# left_model    right_model    tsd_length
LEFT_TIR        RIGHT_TIR      8
ITR_5prime      ITR_3prime     5
```

```bash
tirmite pair \
  --genome $GENOME \
  --nhmmer-file multi_model_hits.tab \
  --lengths-file model_lengths.txt \
  --pairing-map pairing_map.txt \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --flank-len 30 \
  --insertion-site \
  --tsd-length-map tsd_lengths.tsv \
  --outdir OUTPUT \
  --gff-out
```

### Target site output files

| File | Contents |
|------|----------|
| `{prefix}{pair_label}_target_sites.fasta` | One reconstructed target site per paired element |
| `{prefix}{pair_label}_interleaved_flanks.fasta` | Gap-padded left/right flank pairs showing TSD position |

For single-model runs `pair_label` is the model name. For asymmetric pairs it is `{left_model}_{right_model}`.

FASTA headers in `target_sites.fasta` include metadata tags:

| Tag | Description |
|-----|-------------|
| `flank_len` | Extracted flank length (bp) |
| `tsd_len` | User-specified TSD length |
| `tsd_in_model` | Whether TSD is inside the terminus model |
| `left_model` / `right_model` | Names of terminus models for this pair |
| `contig` | Sequence ID in the target genome |
| `left_flank_hit` | Strand and hit coordinates for the left terminus |
| `right_flank_hit` | Strand and hit coordinates for the right terminus |
| `tsd_hamming` | Hamming distance between left and right TSD copies |
| `left_tsd` / `right_tsd` | TSD sequences (when > 0 bp) |

!!! note "Next step: validation"
    Use the `target_sites.fasta` output as input to [`tirmite validate`](tirmite-validate.md) to search a reference genome for matching empty insertion sites and independently confirm the TSD length.

## Hit Anchoring: Filtering by Model Edge Proximity

The `--max-offset` flag filters out hits whose alignment does not reach close enough to the **outer edge** of the query model. This removes spurious interior fragment hits that are unlikely to represent true element termini.

### What is the outer edge?

Each terminus model has an external end that faces away from the element interior:

- **Left terminus**: the outer edge is model position 1 (the very start of the model)
- **Right terminus**: the outer edge is model position `model_len` (the very end of the model)

If a hit's alignment starts at model position 5 (for a left terminus hit on the + strand), the offset is 4. Setting `--max-offset 5` would keep this hit; setting `--max-offset 3` would remove it.

```mermaid
flowchart LR
    subgraph "Left terminus hit on + strand"
        direction LR
        M1[Model pos 1\n← outer edge] -.offset.- A1[alignment\nstart] --> A2[alignment\nend]
    end
    subgraph "Right terminus hit on + strand"
        direction LR
        B1[alignment\nstart] --> B2[alignment\nend] -.offset.- M2[Model pos N\nouter edge →]
    end
```

### Terminus type determination

TIRmite determines whether each hit is a left or right terminus using:

1. **Pairing map** (if provided): model name determines the terminus type
2. **F,R orientation** (strands differ): + strand hits → left terminus; − strand hits → right terminus
3. **F,F or R,R** (same-strand, no pairing map): the hit must be within `--max-offset` of **both** ends of the query model (i.e. the full model must be nearly covered)

### Example

```bash
tirmite pair \
  --genome $GENOME \
  --nhmmer-file $NHMMERFILE \
  --hmm-file $HMMFILE \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --max-offset 5 \
  --outdir MY_TIR_OUTPUT \
  --gff-out
```

!!! note
    `--max-offset` is applied **before** the pairing step, so it reduces the pool of candidate hits available for pairing. It is complementary to `--mincov`, which filters on overall fractional coverage of the model rather than the proximity of the alignment to the outer edge.

## Multiple Models: Using the Pairing Map

When your input file contains hits from **multiple HMM models or BLAST queries**, you must provide a `--pairing-map` file. Without it, TIRmite cannot know which models should be paired together and may produce incorrect results.

### Pairing map format

A tab-delimited file with two columns: `left_feature` and `right_feature`:

```
# pairing_map.txt
# left_feature    right_feature
```

**Symmetric pairing** (same model on both ends):

```
model1    model1
model2    model2
```

**Asymmetric pairing** (different models):

```
LEFT_TIR    RIGHT_TIR
ITR_5prime  ITR_3prime
```

### How pairing map logic works

When a pairing map is provided, TIRmite:

1. Splits hits into groups corresponding to each left/right model pair defined in the map
2. Runs the pairing algorithm **independently** for each pair combination
3. Tracks unpaired hits across all procedures (a hit may be unpaired in one combination but paired in another)

!!! note "Hits in multiple combinations"
    Features can appear in multiple pairing combinations if needed. For example, if `model1` is listed as left in one pair and right in another, TIRmite runs independent pairing procedures for each combination and correctly tracks unpaired hits across all.

### Example with multiple models

```bash
tirmite pair \
  --genome $GENOME \
  --nhmmer-file multi_model_hits.tab \
  --lengths-file model_lengths.txt \
  --pairing-map pairing_map.txt \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --outdir OUTPUT \
  --gff-out
```

## Reporting Options

The `--gff-report` flag controls which hit categories are written to output files:

| Value | Description |
|-------|-------------|
| `all` | Report all hits: paired elements, paired hits, and unpaired hits |
| `paired` | Report only paired elements and their constituent hits |
| `unpaired` | Report only unpaired hits |

### GFF3 output

Enable GFF3 output with `--gff-out`. The GFF3 file includes:

- Predicted element features (paired left+right termini)
- Individual terminus hit features (if `--gff-report all`)
- Attributes including model name, e-value, coverage, and orientation

```bash
tirmite pair \
  --genome $GENOME \
  --nhmmer-file $NHMMERFILE \
  --hmm-file $HMMFILE \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --gff-report all \
  --gff-out \
  --outdir MY_TIR_OUTPUT
```

## Key Options Reference

| Option | Description |
|--------|-------------|
| `--genome` | Path to target genome FASTA |
| `--blastdb` | Path to BLAST database (alternative to `--genome` for sequence extraction) |
| `--nhmmer-file` | nhmmer tabular output file |
| `--blast-file` | BLAST tabular output file (format 6) |
| `--hmm-file` | HMM file (used to determine model lengths for coverage calculation) |
| `--lengths-file` | Tab-delimited file of model lengths (alternative to `--hmm-file`) |
| `--query-len` | Query length for BLAST input (when single-model BLAST results) |
| `--pairing-map` | Tab-delimited file linking left and right model names |
| `--orientation` | Strand orientation for valid pairs (e.g. `F,R`, `F,F`) |
| `--mincov` | Minimum fraction of model covered by hit (0–1) |
| `--maxdist` | Maximum distance between paired terminus hits (bp) |
| `--maxeval` | Maximum e-value for hit filtering |
| `--max-offset` | Maximum unaligned model positions between hit edge and outer terminus edge (anchor filter) |
| `--stable-reps` | Iterations without new pairing before stopping (default: 2) |
| `--padlen` | Flanking bases to pad the extracted hit/element sequences |
| `--flanks` | Enable writing of external flanking sequences for all hits |
| `--flanks-paired` | Write outer flanking sequences for paired termini only |
| `--flank-len` | Length of external flanking region to extract per terminus (default: 50 bp) |
| `--flank-max-offset` | Skip flank extraction for hits where offset from outer model edge exceeds this value |
| `--insertion-site` | Enable insertion site reconstruction and reporting (requires `--flanks` or `--flanks-paired`) |
| `--tsd-length` | TSD/DR length for target site reconstruction (requires `--insertion-site`) |
| `--tsd-length-map` | Per-model-pair TSD length map (TSV: left_model, right_model, tsd_length) |
| `--tsd-in-model` | The TSD is encoded at the inner end of the terminus HMM model |
| `--gff-report` | Reporting mode: `all`, `paired`, or `unpaired` |
| `--gff-out` | Write GFF3 annotation file |
| `--logfile` | Write log to file |
