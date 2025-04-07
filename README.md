[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/TIRmite.svg)](https://badge.fury.io/py/TIRmite)
[![codecov](https://codecov.io/gh/Adamtaranto/TIRmite/graph/badge.svg?token=DFEEPKDFZ0)](https://codecov.io/gh/Adamtaranto/TIRmite)

<p align="center">
<img src="https://raw.githubusercontent.com/Adamtaranto/TIRmite/main/docs/tirmite_hexlogo.jpg" width="256" height="256" title="tirmite_hex" />
</p>

# TIRmite

Build and map profile Hidden Markov Models for Terminal Inverted Repeat
families (TIR-pHMMs) to genomic sequences for annotation of MITES and complete
DNA-Transposons with variable internal sequence composition.  

If you have a draft TE model (i.e. from RepeatModeler or EDTA) and want to identify the TIR's to use with TIRmite - we recommend using [*tSplit*](https://github.com/Adamtaranto/TE-splitter/) a tool for extraction of terminal repeats from complete transposons.

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

TIRmite will use profile-HMM models of Terminal Inverted Repeats (TIRs) for 
genome-wide annotation of TIR families. These can be provided by the user or
built from aligned TIRs oriented as 5' outer edge --> 3' inner edge.

Three classes of output are produced:

  1. All significant TIR hit sequences written to fasta (per query HMM).
  2. Candidate elements comprised of paired TIRs are written to fasta (per query HMM).
  3. Genomic annotations of candidate elements and, optionally, TIR hits
  (paired and unpaired) are written as a single GFF3 file.

## Algorithm overview

  1. Use nhmmer genome with TIR-pHMM.
  2. Import all hits below *--maxeval* threshold.
  3. For each significant TIR match identify candidate partners, where:  
    * Is on the same sequence.  
    * Hit is in complementary orientation.  
    * Distance is <= *--maxdist*.  
    * Hit length is >= model length \* *--mincov*.  
  4. Rank candidate partners by distance downstream of positive-strand hits, and upstream of negative-strand hits.
  5. Pair reciprocal top candidate hits.
  6. For unpaired hits, find first unpaired candidate partner and check for reciprocity.
  7. If the first unpaired candidate is non-reciprocal, check for 2nd-order reciprocity (is outbound top-candidate of current candidate reciprocal.)
  8. Iterate steps 6-7 until all TIRs are paired OR number of iterations without new pairing exceeds *--stableReps*.

## Options and usage

### Installing TIRmite

TIRmite requires Python >= v3.8

Dependencies:

* TIR-pHMM build and search
  * [HMMER3](http://hmmer.org)
* Extract terminal repeats from predicted TEs
  * [pymummer](https://github.com/sanger-pathogens/pymummer) version >= 0.10.3 with wrapper for nucmer option *--diagfactor*.
  * [MUMmer](https://github.com/mummer4/mummer)
  * [BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (Optional)

You can create a Conda environment with these dependencies using the YAML files in this repo.

```bash
conda env create -f environment.yml

conda activate tirmite
```

Note: If you are using a Mac with an ARM64 (Apple Silicon) processor, BLAST is not currently available from Bioconda for this architecture. You can instead create a virtual OSX64 env like this:

```bash
conda env create -f env_osx64.yml

conda activate tirmite-osx64
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

Do this if you want to edit the code.

```bash
git clone https://github.com/Adamtaranto/TIRmite.git && cd TIRmite && pip install -e '.[dev]'
```

Test installation.

```bash
# Print version number and exit.
% tirmite --version
tirmite 1.2.0

# Get usage information
% tirmite --help
```

### Example usage

Report all hits and valid pairings of TIR_A in target.fasta (interval <= 10000, hits cover > 40% len of hmm model), and write GFF3 annotation file.

```bash
% tirmite --genome target.fasta --hmmFile TIR_A.hmm --gffOut TIR_elements_in_Target.gff3 --maxdist 10000 --mincov 0.4
```

If you don't have a HMM of your TIR, TIRmite can create one for you using an aligned sample of your TIR with `--alnFile`.

To skip HMM search and run the pairing algorithm on a custom set of TIR hits (i.e. from blastn), you can provide hits in BED format with `--pairbed`.

TIRs should always be oriented 5\`- 3\` with the lefthand TIR.

In this example the two TIRs should be oriented to begin with "GA".

5\` **GA\>\>\>\>\>\>\>** ATGC <<<<<<<TC 3\`  
3\` CT>>>>>>>>  TACG <<<<<<<AG 5\`

### Standard options

Run `tirmite --help` to view the program's most commonly used options:

```code
tirmite [-h] [--version] --genome GENOME [--hmmDir HMMDIR]
               [--hmmFile HMMFILE] [--alnDir ALNDIR] [--alnFile ALNFILE]
               [--alnFormat {clustal,fasta,nexus,phylip,stockholm}]
               [--pairbed PAIRBED] [--stableReps STABLEREPS] [--outdir OUTDIR]
               [--prefix PREFIX] [--nopairing] [--gffOut]
               [--reportTIR {None,all,paired,unpaired}] [--padlen PADLEN]
               [--keeptemp] [-v] [--cores CORES] [--maxeval MAXEVAL]
               [--maxdist MAXDIST] [--nobias] [--matrix MATRIX]
               [--mincov MINCOV] [--hmmpress HMMPRESS] [--nhmmer NHMMER]
               [--hmmbuild HMMBUILD]

Info: 
  -h, --help            Show this help message and exit
  --version             Show program's version number and exit
  
Input options:
  --genome              Path to target genome that will be queried with HMMs.
                          Note: Sequence names must be unique. (required)
  --hmmDir              Directory containing pre-prepared TIR-pHMMs.
  --hmmFile             Path to single TIR-pHMM file. Incompatible with "--hmmDir".
  --alnDir              Path to directory containing only TIR alignments to be
                          converted to HMM.
  --alnFile             Provide a single TIR alignment to be converted to HMM.
                          Incompatible with "--alnDir".
  --alnFormat           Alignments provided with "--alnDir" or "--alnFile" are
                          all in this format.
                          Choices=["clustal","fasta","nexus","phylip", "stockholm"]
  --pairbed             If set TIRmite will preform pairing on TIRs from
                          custom bedfile only.

Pairing heuristics:
  --stableReps          Number of times to iterate pairing procedure when no
                         additional pairs are found AND remaining unpaired hits > 0.
                         (Default = 0)

Output and housekeeping:
  --outdir OUTDIR       All output files will be written to this directory.
  --prefix PREFIX       Add prefix to all TIRs and Paired elements detected in
                          this run. Useful when running same TIR-pHMM against
                          many genomes.
                          (Default = None)
  --nopairing           If set, only report TIR-pHMM hits. Do not attempt
                          pairing.
                          (Default = False)
  --gffOut              If set report features as prefix.gff3. File saved to
                          outdir.
                          (Default = False)
  --reportTIR           Options for reporting TIRs in GFF annotation file.
                          Choices=[None,'all','paired','unpaired']
                          (Default = 'all')
  --padlen              Extract x bases either side of TIR when writing TIRs to fasta.
                          (Default = None)
  --keeptemp            If set do not delete temp file directory.
                          (Default = False)
  -v, --verbose         Set syscall reporting to verbose.
  
HMMER options:
  --cores               Set number of cores available to hmmer software.
  --maxeval             Maximum e-value allowed for valid hit.
                          (Default = 0.001)
  --maxdist             Maximum distance allowed between TIR candidates to
                          consider valid pairing.
                          (Default = None)
  --nobias              Turn OFF bias correction of scores in nhmmer.
                          (Default = False)
  --matrix              Use custom DNA substitution matrix with nhmmer.
  --mincov              Minimum valid hit length as prop of model length.
                          (Default = 0.5)

Non-standard HMMER paths:
  --hmmpress            Set location of hmmpress if not in PATH.
  --nhmmer              Set location of nhmmer if not in PATH.
  --hmmbuild            Set location of hmmbuild if not in PATH.
```

### Custom DNA Matrices

nhmmer can be supplied with custom DNA score matrices for assessing hmm match scores.
Standard NCBI-BLAST matrices such as NUC.4.4 are compatible. (See: ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/NUC.4.4)

## Issues

Submit feedback to the [Issue Tracker](https://github.com/Adamtaranto/TIRmite/issues)

## License

Software provided under MIT license.
