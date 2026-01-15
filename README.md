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

TIRmite requires Python >= v3.8

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

First, you will need to build a pHMM of your element's terminal sequence/s.

If you have a draft TE model (i.e. from RepeatModeler or EDTA) and want to identify the TIR's or LTR's to use with TIRmite - I recommend using [*tSplit*](https://github.com/Adamtaranto/TE-splitter/) a tool for extraction of terminal repeats from complete transposons.


1) Extract single TIR from sample element:

```bash
# Uses BLASTn to detect TIRs of min 40% identity and min 10 bp length
tsplit TIR -i TIR_element.fa -d tsplit_results --minid 0.4 --method blastn --minterm 10 --splitmode external
```

2) Build a pHMM from the seed:

```bash
GENOME="genome.fa" # Path to fasta containing one or more genomes to search for matches to seed sequence.

tirmite seed --left-seed tsplit_results/TIR_element_tsplit_output.fasta --model-name MY_TIR --outdir MY_TIR_HMM --genome $GENOME --max-gap 10 --save-blast-hits --threads 8

# Note: Setting `--flank-size 10` will output additional flanking bases outside the TIR, conservation in the flank accross many independent insertions may indicate your seed was truncated. Always check and adjust seed as required.
```

3) Use `nhmmer` to locate hits to the TIR-pHMM in a target genome.

```bash
HMMFILE="MY_TIR_HMM/MY_TIR.hmm"
NHMMERFILE="MY_TIR_nhmmer_hits.tab"
nhmmer --dna --cpu 8 --tblout $NHMMERFILE $HMMFILE $GENOME
```

**Custom DNA Matrices**

Note: nhmmer can be supplied with custom DNA score matrices for assessing hmm match scores.
Standard NCBI-BLAST matrices such as NUC.4.4 are compatible. (See: ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/NUC.4.4)

4) Use `tirmite pair` to identify valid TIR pairs. Outputs hits, elements, and annotations.

```bash
tirmite pair --genome $GENOME  --nhmmerFile $NHMMERFILE --hmmFile $HMMFILE --orientation F,R --mincov 0.4 --report all  --maxdist 20000 --stableReps 2 --outdir MY_TIR_PAIRING_OUTPUT --padlen 20 --maxeval 0.001 --gffOut --logfile
```

#### Legacy mode

TIRmite `legacy` mode will take a TIR-pHMM and target genome fasta as input and run the full standard workflow, reporting all hits, valid pairings, and write GFF3 annotation file.

Note: This usage will be phased out in a later release in favour of custom workflows.

```bash
# Use HMM search to pull more divergent TIR hits from your query genome.
# TIR hits are paired in Fwd/Rev orientation
# Fwd/Rev pairs must be within 20Kbp of each other
# Hits must cover >= 40% of the TIR-pHMM
tirmite legacy --genome $GENOME --hmmFile $HMMFILE--orientation F,R \
--outdir results \
--stableReps 2 \
--report all \
--gffOut --maxdist 20000 --mincov 0.4
```

If you don't have a HMM of your TIR, `tirmite legacy` can create one for you using an aligned sample of your TIR provided with `--alnFile`.

TIRs should always be oriented 5\`- 3\` with the lefthand TIR.

In this example the two TIRs should be oriented to begin with "GA".

5\` **GA\>\>\>\>\>\>\>** ATGC <<<<<<<TC 3\`
3\` CT>>>>>>>>  TACG <<<<<<<AG 5\`

### Standard options

Run `tirmite --help` to view the program's most commonly used options:

```code
tirmite --help
usage: tirmite [-h] [--version] COMMAND ...

TIRmite: Transposon Terminal Repeat detection suite

positional arguments:
  COMMAND     Available subcommands
    legacy    Original TIRmite workflow (HMM search + pairing)
    seed      Build HMM models from seed sequences
    pair      Pair precomputed nhmmer hits

options:
  -h, --help  show this help message and exit
  --version   show program's version number and exit

Available subcommands:
  legacy    Original TIRmite workflow (HMM search + pairing)
  seed      Build HMM models from seed sequences
  pair      Pair precomputed nhmmer hits

Examples:
  tirmite legacy --genome genome.fa --hmmFile model.hmm
  tirmite seed --left-seed left.fa --model-name myTE --genome genome.fa
  tirmite pair --genome genome.fa --nhmmerFile hits.out --hmmFile model.hmm
```

## Algorithm overview

  1. Use nhmmer to query genome with termini HMMs.
  2. Import all hits under *--maxeval* threshold.
  3. For each significant terminus match, identify candidate partners, where:
    - Hit is on the same sequence.
    - Hit is in coreect relative orientation.
    - Distance is <= *--maxdist*.
    - Hit length is >= (model length * *--mincov* prop)
  4. Rank candidate partners by distance downstream of positive-strand hits, and upstream of negative-strand hits.
  5. Pair reciprocal top candidate hits.
  6. For unpaired hits, find nearest unpaired candidate partner and check for reciprocity.
  7. If the first unpaired candidate is non-reciprocal, check for 2nd-order reciprocity (is outbound top-candidate of current candidate reciprocal.)
  8. Iterate steps 6-7 until all termini hits are paired OR number of iterations without new pairing exceeds *--stableReps*.

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
