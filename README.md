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

**TIRmite** employs profile Hidden Markov Models (HMMs) to model natural variation in transposon termini and recover divergent and degraded hits that are often missed by sequence-based aligners like BLAST.

An iterative pairing algorithm is then used to annotate cryptic transposon variants with variable internal sequence compositions.

The elements extracted by TIRmite generally represent structuaral variants derived from an autonomous ancestor and may be further clustered into families.

# Table of contents

* [About TIRmite](#about-tirmite)
* [Options and usage](#options-and-usage)
  * [Installing TIRmite](#installing-tirmite)
  * [Example usage](#example-usage)
  * [Standard options](#standard-options)
* [Algorithm overview](#algorithm-overview)
* [Contributing](#contributing)
* [Issues](#issues)
* [License](#license)

## About TIRmite

TIRmite will use profile-HMM models of Transposon Terminal Repeats for genome-wide annotation of transposon families. You can search for TE families with symmetrical termini (i.e. TIRs or LTRs) or asymmetrical elements with different conserved features at either end (i.e. Helitrons, Helentrons, and Starship elements).

Three classes of output are produced:

  1. All significant termini hit sequences are written to fasta (per query HMM).
  2. Candidate elements comprised of paired termini are written to fasta (per query HMM).
  3. Genomic annotations of candidate elements and, optionally, HMM hits
  (paired and unpaired) are written as a single GFF3 file.

## Options and usage

### Installing TIRmite

TIRmite requires Python >= v3.9

Dependencies:

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

2) Install latest release from PyPi.

```bash
pip install tirmite
```

3) Install latest release (with dependencies) from Bioconda.

```bash
conda install -c bioconda tirmite
```

Test installation.

```bash
# Print version number and exit.
% tirmite --version
tirmite 1.3.0

# Get usage information
% tirmite --help
```

### Example usage

See the [online tutorials](https://adamtaranto.github.io/TIRmite) for detailed walkthroughs of each module:

- [**Building HMMs**](https://adamtaranto.github.io/TIRmite/tutorials/building-hmms/) — build a profile HMM from a seed TIR/LTR sequence using `tirmite seed`.
- [**Using tirmite search**](https://adamtaranto.github.io/TIRmite/tutorials/tirmite-search/) — ensemble BLAST/nhmmer search and hit merging.
- [**Using tirmite pair**](https://adamtaranto.github.io/TIRmite/tutorials/tirmite-pair/) — pairing terminus hits, flank extraction and target site reconstruction.
- [**Reconstructing and validating target sites**](https://adamtaranto.github.io/TIRmite/tutorials/tirmite-validate/) — validate reconstructed target sites with `tirmite validate`.

#### Quick start

```bash
GENOME="genome.fa"
HMMFILE="MY_TIR.hmm"
NHMMERFILE="MY_TIR_nhmmer_hits.tab"

# 1. Search genome for terminus hits
nhmmer --dna --cpu 8 --tblout $NHMMERFILE $HMMFILE $GENOME

# 2. Pair hits and write elements + GFF3
tirmite pair \
  --genome $GENOME \
  --nhmmer-file $NHMMERFILE \
  --hmm-file $HMMFILE \
  --orientation F,R \
  --mincov 0.4 \
  --maxdist 20000 \
  --gff-report all \
  --gff \
  --outdir MY_TIR_OUTPUT
```

To also reconstruct target sites (see [tutorial](https://adamtaranto.github.io/TIRmite/tutorials/tirmite-validate/)):

```bash
tirmite pair \
  --genome $GENOME \
  --nhmmer-file $NHMMERFILE \
  --hmm-file $HMMFILE \
  --orientation F,R \
  --mincov 0.6 \
  --maxdist 20000 \
  --flank-len 30 \
  --tsd-length 2 \
  --outdir MY_TIR_OUTPUT \
  --gff
```

### Standard options

Run `tirmite --help` to view available subcommands:

```
tirmite --help
usage: tirmite [-h] [--version] COMMAND ...

TIRmite: Transposon Terminal Repeat detection suite

positional arguments:
  COMMAND     Available subcommands
    legacy    Original TIRmite workflow (HMM search + pairing)
    seed      Build HMM models from seed sequences
    pair      Pair precomputed nhmmer hits
    search    Ensemble search: merge hits from clustered features
    validate  Validate reconstructed target sites

options:
  -h, --help  show this help message and exit
  --version   show program's version number and exit
```

## Algorithm overview

  1. Use nhmmer (or BLAST) to query genome with termini models/sequences.
  2. Import all hits under *--maxeval* threshold.
  3. For each significant terminus match, identify candidate partners, where:
    - Hit is on the same sequence.
    - Hit is in correct relative orientation.
    - Distance is <= *--maxdist*.
    - Hit length is >= (model/query length * *--mincov* prop)
  4. Rank candidate partners by distance downstream of positive-strand hits, and upstream of negative-strand hits.
  5. Pair reciprocal top candidate hits.
  6. For unpaired hits, find nearest unpaired candidate partner and check for reciprocity.
  7. If the first unpaired candidate is non-reciprocal, check for 2nd-order reciprocity (is outbound top-candidate of current candidate reciprocal.)
  8. Iterate steps 6-7 until all termini hits are paired OR number of iterations without new pairing exceeds *--stable-reps*.

### Ensemble search with `tirmite search` (optional pre-processing step)

For complex scenarios with multiple sub-type HMMs or asymmetric element families, use `tirmite search` to merge and filter hits before pairing:

  1. Run BLAST and/or nhmmer with multiple query models simultaneously.
  2. Optionally merge overlapping hits from clustered component features (cluster map).
  3. When a pairing map is provided:
     - **Step 0**: Exclude hits from models not listed in the pairing map.
     - **Step 1**: Remove nested weak hits within each direct left/right pair.
     - **Step 2**: Remove lower-quality cross-model hits at shared genomic loci.
     - Emit a structured **filter summary report** covering all three steps (per-model exclusion counts, nesting relationships, and per-pair cross-model removal counts).
  4. Output a filtered, merged hit table ready for `tirmite pair`.
  5. Optionally write separate left/right hit files (`--split-paired-output`) for asymmetric elements.

## Contributing

If you would like to add a new feature or fix a bug, please see our [contribution guidelines](https://github.com/Adamtaranto/TIRmite?tab=contributing-ov-file#readme).

- Open an issue
- Fork the repo
- Follow the dev env setup instructions below

```bash
# Clone this repo (or your own fork)
git clone https://github.com/Adamtaranto/TIRmite.git && cd TIRmite
# Install custom conda env
conda env create -f environment.yml
# Activate conda env
conda activate tirmite
# Install an editable copy of the package
pip install -e '.[dev]'
# Enable pre-commit checks
pre-commit install
```

## Issues

Submit feedback to the [Issue Tracker](https://github.com/Adamtaranto/TIRmite/issues)

## License

Software provided under GPL-3 license.

## Star History

[![Star History
Chart](https://api.star-history.com/svg?repos=adamtaranto/tirmite&type=Date)](https://star-history.com/#adamtaranto/tirmite&Date)
