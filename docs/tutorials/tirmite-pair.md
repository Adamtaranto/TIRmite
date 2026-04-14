# Using tirmite pair

`tirmite pair` takes precomputed nhmmer or BLAST search results and applies the TIRmite pairing algorithm to identify valid transposon termini pairs. It outputs paired elements, hit sequences, and GFF3 annotations.

## Overview

```mermaid
flowchart TD
    A[nhmmer .tab\nor BLAST .tab] --> B[Load hits]
    B --> C[Filter by e-value\nand model coverage]
    C --> D{Pairing map\nprovided?}
    D -->|Yes| E[Split hits by model pair\nRun independent pairing\nfor each combination]
    D -->|No| F[Run single pairing\nprocedure]
    E --> G[Pairing algorithm\nper model pair]
    F --> G
    G --> H[Output paired elements,\nhit sequences, GFF3]
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

## Extracting Flanking Regions

Use `--padlen` to include flanking sequence adjacent to each hit or predicted element in the output FASTA files. This is useful for:

- Checking whether element boundaries are correctly identified
- Examining Target Site Duplications (TSDs)
- Providing context for downstream analysis

```bash
tirmite pair \
  --genome $GENOME \
  --nhmmer-file $NHMMERFILE \
  --hmm-file $HMMFILE \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --padlen 20 \
  --outdir MY_TIR_OUTPUT \
  --gff-out
```

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

The `--report` flag controls which hit categories are written to output files:

| Value | Description |
|-------|-------------|
| `all` | Report all hits: paired elements, paired hits, and unpaired hits |
| `paired` | Report only paired elements and their constituent hits |
| `unpaired` | Report only unpaired hits |

### GFF3 output

Enable GFF3 output with `--gff-out`. The GFF3 file includes:

- Predicted element features (paired left+right termini)
- Individual terminus hit features (if `--report all`)
- Attributes including model name, e-value, coverage, and orientation

```bash
tirmite pair \
  --genome $GENOME \
  --nhmmer-file $NHMMERFILE \
  --hmm-file $HMMFILE \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --report all \
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
| `--stable-reps` | Iterations without new pairing before stopping (default: 2) |
| `--padlen` | Flanking bases to include in extracted sequences |
| `--report` | Reporting mode: `all`, `paired`, or `unpaired` |
| `--gff-out` | Write GFF3 annotation file |
| `--logfile` | Write log to file |
