[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/TIRmite.svg)](https://badge.fury.io/py/TIRmite)
[![codecov](https://codecov.io/gh/Adamtaranto/TIRmite/graph/badge.svg?token=DFEEPKDFZ0)](https://codecov.io/gh/Adamtaranto/TIRmite)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/tirmite/README.html)
[![Downloads](https://pepy.tech/badge/tirmite)](https://pepy.tech/project/tirmite)

<p align="center">
<img src="https://raw.githubusercontent.com/Adamtaranto/TIRmite/main/docs/tirmite_hexlogo.jpg" width="256" height="256" title="tirmite_hex" />
</p>

# TIRmite

Build and map profile Hidden Markov Models (HMMs) of Transposon Terminal Repeat
families to genomic sequences for annotation of cryptic variants with variable internal sequence composition.  

If you have a draft TE model (i.e. from RepeatModeler or EDTA) and want to identify the TIR's or LTR's to use with TIRmite - I recommend using [*tSplit*](https://github.com/Adamtaranto/TE-splitter/) a tool for extraction of terminal repeats from complete transposons.

# Table of contents

* [About TIRmite](#about-tirmite)
* [Algorithm overview](#algorithm-overview)
* [Options and usage](#options-and-usage)
  * [Installing TIRmite](#installing-tirmite)
  * [Example usage](#example-usage)
  * [Standard options](#standard-options)
  * [Custom DNA matrices](#custom-dna-matrices)
* [Issues](#issues)
* [License](#license)

## About TIRmite

TIRmite will use profile-HMM models of Transposon Terminal Repeats for genome-wide annotation of transposon families. You can search for TE families with symmetrical termini (i.e. TIRs or LTRs) or asymmetrical elements with different conserved features at either end (i.e. Starship elements).

Three classes of output are produced:

  1. All significant termini hit sequences are written to fasta (per query HMM).
  2. Candidate elements comprised of paired termini are written to fasta (per query HMM).
  3. Genomic annotations of candidate elements and, optionally, HMM hits
  (paired and unpaired) are written as a single GFF3 file.

## Algorithm overview

  1. Use nhmmer to query genome with termini HMMs.
  2. Import all hits under *--maxeval* threshold.
  3. For each significant terminus match, identify candidate partners, where:  
    - Is on the same sequence.  
    - Hit is in complementary orientation.  
    - Distance is <= *--maxdist*.  
    - Hit length is >= (model length * *--mincov* prop)  
  4. Rank candidate partners by distance downstream of positive-strand hits, and upstream of negative-strand hits.
  5. Pair reciprocal top candidate hits.
  6. For unpaired hits, find first unpaired candidate partner and check for reciprocity.
  7. If the first unpaired candidate is non-reciprocal, check for 2nd-order reciprocity (is outbound top-candidate of current candidate reciprocal.)
  8. Iterate steps 6-7 until all termini hits are paired OR number of iterations without new pairing exceeds *--stableReps*.

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

pip install the latest development version directly from this repo.

```bash
% pip install git+https://github.com/Adamtaranto/TIRmite.git
```

Install latest release from PyPi.

```bash
% pip install tirmite
```

Install from Bioconda.

```bash
% conda install -c bioconda tirmite
```

Clone from this repository and install as a local Python package.

Only do this if you want to edit the code.

```bash
git clone https://github.com/Adamtaranto/TIRmite.git && cd TIRmite && pip install -e '.[dev]'
```

Test installation.

```bash
# Print version number and exit.
% tirmite --version
tirmite 1.2.1

# Get usage information
% tirmite --help
```

### Example usage

Report all hits and valid pairings of TIR_A in target.fasta (interval <= 10000, hits cover > 40% len of hmm model), and write GFF3 annotation file.

```bash
% tirmite legacy --genome target.fasta --hmmFile TIR_A.hmm --gffOut TIR_elements_in_Target.gff3 --maxdist 10000 --mincov 0.4 --orientation F,R
```

If you don't have a HMM of your TIR, TIRmite can create one for you using an aligned sample of your TIR with `--alnFile`.

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

### Custom DNA Matrices

nhmmer can be supplied with custom DNA score matrices for assessing hmm match scores.
Standard NCBI-BLAST matrices such as NUC.4.4 are compatible. (See: ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/NUC.4.4)

## Issues

Submit feedback to the [Issue Tracker](https://github.com/Adamtaranto/TIRmite/issues)

## License

Software provided under MIT license.
