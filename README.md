[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI version](https://badge.fury.io/py/TIRmite.svg)](https://badge.fury.io/py/TIRmite)
[![codecov](https://codecov.io/gh/Adamtaranto/TIRmite/graph/badge.svg?token=DFEEPKDFZ0)](https://codecov.io/gh/Adamtaranto/TIRmite)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/tirmite/README.html)
[![Downloads](https://pepy.tech/badge/tirmite)](https://pepy.tech/project/tirmite)

<p align="center">
<img src="https://raw.githubusercontent.com/Adamtaranto/TIRmite/main/docs/tirmite_hexlogo.jpg" width="256" height="256" title="tirmite_hex" />
</p>

# TIRmite

Autonomous examples of transposons, belonging to many [distinct super-families, share two common properties](https://doi.org/10.1266/ggs.18-00024 "Yes yes, except when they don't. Don't @ me, nerds."): A gene or genes encoding the mode of transposition; and terminal sequence features that are recognised by these gene products as the element boundaries.

Proper classification of transposons and grouping into families relies on both phylogeny of conserved sequences and conservation of transposition mechanism.

However, not all TE instances are created equal — inhabiting the nulear soup of their host genome, where your brother's transposase is as good as your own, non-autonomous variants (lacking their own functional hardware) proliferate.

MITEs are a classic example of this - derived from autonomous DNA elements with Terminal Inverted Repeats, they are Miniature Inverted-repeat Transposable Elements, sometimes little more than a pair of TIRs.

When non-autonomous structural variants of a TE vastly outnumber their parent element, and include forms that capture novel genes (or other full transposons!), it becomes difficult to correctly cluster related elements based on the limited signal present in terminal sequences (TIRs, LTRs, etc).

**TIRmite** employs profile Hidden Markov Models (HMMs) OR ensembles of sequence queries to model natural variation in transposon termini and recover divergent and degraded hits that are often missed by sequence-based aligners like BLAST with single queries.

An iterative pairing algorithm is then used to annotate cryptic transposon variants with variable internal sequence compositions.

The elements extracted by TIRmite generally represent structural variants derived from an autonomous ancestor and may be further clustered into families.

## Table of contents

* [About TIRmite](#about-tirmite)
* [Options and usage](#options-and-usage)
  * [Installing TIRmite](#installing-tirmite)
  * [Example usage](#example-usage)
* [License](#license)

## About TIRmite

TIRmite can use profile-HMM models of Transposon Terminal Repeats OR BLAST sequence hits for genome-wide annotation of transposon families.

You can search for TE families with symmetrical termini (i.e. TIRs or LTRs) or asymmetrical elements, where there are different conserved features at either end (i.e. Helitrons, Helentrons, and Starship elements).

Three classes of output are produced:

  1. All significant termini hit sequences are written to fasta (per terminus model query).
  2. Candidate elements comprised of paired termini are written to fasta (per terminus model query).
  3. Genomic annotations of candidate elements and, optionally, terminus model hits
  (paired and unpaired) are written as a single GFF3 file.

## Options and usage

### Installing TIRmite

Dependencies:

* Python >= v3.10
* [HMMER3](http://hmmer.org)
* [mafft](https://mafft.cbrc.jp/alignment/software/)
* [BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (Optional)

You can create a Conda environment with these dependencies using the `environment.yml` file in this repo.

```bash
conda env create -f environment.yml

conda activate tirmite
```

Installation options:

1) `pip install` the latest development version directly from this repo.

```bash
pip install git+https://github.com/Adamtaranto/TIRmite.git
```

2) Install latest release (with dependencies) from Bioconda.

```bash
conda install -c bioconda tirmite
```

Test installation.

```bash
# Print version number and exit.
% tirmite --version
tirmite 1.5.0

# Get usage information
% tirmite --help
```

### Example usage

See the [online tutorials](https://adamtaranto.github.io/TIRmite) for detailed walkthroughs of each module:

- [**Building HMMs**](https://adamtaranto.github.io/TIRmite/tutorials/building-hmms/) — build a profile HMM from a seed TIR/LTR sequence using `tirmite seed`.
- [**Using tirmite search**](https://adamtaranto.github.io/TIRmite/tutorials/tirmite-search/) — ensemble BLAST/nhmmer search and hit merging.
- [**Using tirmite pair**](https://adamtaranto.github.io/TIRmite/tutorials/tirmite-pair/) — pairing terminus hits, flank extraction and target site reconstruction.
- [**Reconstructing and validating target sites**](https://adamtaranto.github.io/TIRmite/tutorials/tirmite-validate/) — validate reconstructed target sites with `tirmite validate`.
- [**HTML reports**](https://adamtaranto.github.io/TIRmite/tutorials/html-reports/) — interactive annotation tracks, terminus alignments and summary plots from `tirmite pair --report` and `tirmite search --report`.

#### Quick start

```bash
GENOME="genome.fa"
HMMFILE="MY_TIR.hmm"
NHMMERFILE="MY_TIR_nhmmer_hits.tab"

# 1. Search genome for terminus hits
nhmmer --dna --cpu 8 --tblout $NHMMERFILE $HMMFILE $GENOME

# 2. Pair hits and write elements + GFF3 + an HTML report
tirmite pair \
  --genome $GENOME \
  --nhmmer-file $NHMMERFILE \
  --hmm-file $HMMFILE \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --gff-report all \
  --gff \
  --report \
  --outdir MY_TIR_OUTPUT
```

Note: `--report` writes `MY_TIR_OUTPUT/tirmite_pair_report.html`

### Ensemble search with `tirmite search` (optional pre-processing step)

For complex scenarios with multiple sub-type HMMs or asymmetric element families, use `tirmite search` to merge and filter hits before pairing:

  1. Run BLAST and/or nhmmer with multiple query models simultaneously.
  2. Optionally restrict hits to those anchored near the outer edge of their model (`--max-offset`).
  3. Optionally merge overlapping hits from clustered component features (cluster map). Gap tolerance is set with `--merge-max-gap`.
  4. When a pairing map is provided:
     - **Step 0**: Exclude hits from models not listed in the pairing map.
     - **Step 1**: Remove nested weak hits within each direct left/right pair.
     - **Step 2**: Remove lower-quality cross-model hits at shared genomic loci.
     - Steps 1 and 2 share a decisiveness threshold set with `--min-score-ratio` (default 1.5): the weaker hit is removed only when the better one outscores it by at least this factor.
     - Emit a structured **filter summary report** covering all three steps (per-model exclusion counts, nesting relationships, and per-pair cross-model removal counts).
  5. Output a filtered, merged hit table ready for `tirmite pair`.
  6. Optionally write separate left/right hit files (`--split-paired-output`) for asymmetric elements.

The cluster map names the cluster **first**, followed by its component models:
`cluster_name<TAB>component1<TAB>component2...`. The pairing map is two columns,
`left_feature<TAB>right_feature`; a row naming the same feature twice describes a
symmetric element.

## License

Software provided under GPL-3 license.

## Star History

<a href="https://www.star-history.com/?repos=Adamtaranto%2FTIRmite&type=date&legend=top-left">
 <picture>
   <source media="(prefers-color-scheme: dark)" srcset="https://api.star-history.com/chart?repos=Adamtaranto/TIRmite&type=date&theme=dark&legend=top-left&sealed_token=AgKIg1I-jujfwocJ2QImXxgdy4jbNsaL9HL_xCA_RAW1EUtcTf3lMZ2JkRaB5h_sZqW2zjhQiZwjgR8xJ6tOAx_B0JRGkKpEmu6nSFL9FaFJwYHfsCNPJcVRvmmRXgUaqMhb8Yjyjcohb_17HVreVp7-M2VU82LfaxYwDQSla8ZLN29DvWbA-1TGbXIh" />
   <source media="(prefers-color-scheme: light)" srcset="https://api.star-history.com/chart?repos=Adamtaranto/TIRmite&type=date&legend=top-left&sealed_token=AgKIg1I-jujfwocJ2QImXxgdy4jbNsaL9HL_xCA_RAW1EUtcTf3lMZ2JkRaB5h_sZqW2zjhQiZwjgR8xJ6tOAx_B0JRGkKpEmu6nSFL9FaFJwYHfsCNPJcVRvmmRXgUaqMhb8Yjyjcohb_17HVreVp7-M2VU82LfaxYwDQSla8ZLN29DvWbA-1TGbXIh" />
   <img alt="Star History Chart" src="https://api.star-history.com/chart?repos=Adamtaranto/TIRmite&type=date&legend=top-left&sealed_token=AgKIg1I-jujfwocJ2QImXxgdy4jbNsaL9HL_xCA_RAW1EUtcTf3lMZ2JkRaB5h_sZqW2zjhQiZwjgR8xJ6tOAx_B0JRGkKpEmu6nSFL9FaFJwYHfsCNPJcVRvmmRXgUaqMhb8Yjyjcohb_17HVreVp7-M2VU82LfaxYwDQSla8ZLN29DvWbA-1TGbXIh" />
 </picture>
</a>
