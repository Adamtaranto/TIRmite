# TIRmite

Map TIR-pHMM models to genomic sequences for annotation of MITES and complete DNA-Transposons.  

# Table of contents

* [About TIRmite](#about-tirmite)
* [Algorithm overview](#algorithm-overview)
* [Options and usage](#options-and-usage)
    * [Installing TIRmite](#installing-tirmite)
    * [Example usage](#example-usage)
    * [Standard options](#standard-options)
    * [Custom DNA matrices](#custom-dna-matrices)

* [License](#license)


# About TIRmite

TIRmite will use profile-HMM models of Terminal Inverted Repeats (TIRs) for 
genome-wide annotation of TIR families. These can be provided by the user or
built from aligned TIRs oriented as 5' outer edge --> 3' inner edge.


Three classes of output are produced:
  1. All significant TIR hit sequences written to fasta (per query HMM).
  2. Candidate elements comprised of paired TIRs are written to fasta (per query HMM).
  3. Genomic annotations of candidate elements and, optionally, TIR hits are written as a single GFF3 file.

# Algorithm overview

  1. Use nhmmer genome with TIR-pHMM
  2. Import all hits below *--maxeval* threshold
  3. For each significant TIR match identify candidate partners, where:
    * Is on same sequence
    * Hit is in complementary orientation
    * Distance is <= *--maxdist*
  4. Rank candidate partners by distance downstream of positive-strand hits, and upstream of negative-strand hits.
  5. Pair reciprocal top candidate hits 
  6. For unpaired hits, find first unpaired candidate partner and check for reciprocity.
  7. If first unpaired candidate is non-reciprocal, check for 2nd-order reciprocity (is outbound top-candidate of currect candidate reciprocal.)
  8. Iterate steps 6-7 until all TIRs are paired OR number of iterations without new pairing exceeds *--stableReps*

# Options and usage

## Installing TIRmite

Dependencies:
TIRmite requires [HMMER3](http://hmmer.org) to be installed.

Install from PyPi:
```
pip install tirmite
```

Clone and install from this repository:
```
git clone https://github.com/Adamtaranto/TIRmite.git && cd TIRmite && pip install -e .
```

## Example usage

Example: Report all hits and vaild pairings of TIR_A in target.fasta (interval <= 50000), and write gff annotation file.

```
./tirmite.py --genome target.fasta --hmmFile TIR_A.hmm --gffOut TIR_elements_in_Target.gff3 --maxdist 50000
```

## Standard options

Run `tirmite --help` to view the program's most commonly used options:

```
usage: tirmite [-h] --genome GENOME [--hmmDir HMMDIR] [--hmmFile HMMFILE]
                 [--alnDir ALNDIR] [--alnFile ALNFILE]
                 [--alnFormat {clustal,emboss,fasta,fasta-m10,ig,maf,mauve,nexus,phylip,phylip-sequential,phylip-relaxed,stockholm}]
                 [--stableReps STABLEREPS] [--outdir OUTDIR] [--prefix PREFIX]
                 [--nopairing][--gffOut GFFOUT]
                 [--reportTIR {None,all,paired,unpaired}] [--keeptemp] [-v]
                 [--cores CORES] [--maxeval MAXEVAL] [--maxdist MAXDIST]
                 [--nobias] [--matrix MATRIX] [--hmmpress HMMPRESS]
                 [--nhmmer NHMMER] [--hmmbuild HMMBUILD]

Help:
  -h, --help              Show this help message and exit

Input options:
  --genome                Path to target genome that will be queried with HMMs.
                            Note: Sequence names must be unique.
                            (required)
  --hmmDir                Directory containing pre-prepared TIR-pHMMs.
  --hmmFile               Path to single TIR-pHMM file. 
                            Incompatible with "--hmmDir".
  --alnDir                Path to directory containing only TIR alignments to be converted to HMM.
  --alnFile               Provide a single TIR alignment to be converted to HMM. 
                            Incompatible with "--alnDir".
  --alnFormat             Alignments provided with "--alnDir" or "--alnFile" are all in this format.
                            Choices=["clustal","emboss","fasta","fasta-m10","ig","maf","mauve",
                            "nexus","phylip","phylip-sequential","phylip-relaxed","stockholm"]


Pairing heuristics:
  --stableReps            Number of times to iterate pairing procedure when no additional pairs are 
                            found AND remaining unpaired hits > 0.
                            (Default = 0)


Output and housekeeping:
  --outdir                All output files will be written to this directory.
  --gffOut                GFF3 annotation filename.
  --reportTIR             Options for reporting TIRs in GFF annotation file.
                            Choices=[None,'all','paired','unpaired']
                            (Default = 'all')
  --prefix                Add prefix to all TIRs and Paired elements detected in this run. 
                            Useful when running same TIR-pHMM against many genomes.
                            (Default = None)
  --nopairing             If set, only report TIR-pHMM hits. Do not attempt pairing.
                            (Default = False)
  --keeptemp              If set do not delete temp file directory.
                            (Default = False)
  -v, --verbose           Set syscall reporting to verbose.

HMMER options:
  --cores                 Set number of cores available to hmmer software.
                            (Default = 1)
  --maxeval               Maximum e-value allowed for valid hit.
                            (Default = 0.001)
  --maxdist               Maximum distance allowed between TIR candidates to consider valid pairing.
  --nobias                Turn OFF bias correction of scores in nhmmer.
                            (Default = False)
  --matrix                Use custom DNA substitution matrix with nhmmer.


Non-standard HMMER paths:
  --hmmpress              Set location of hmmpress if not in path.
  --nhmmer                Set location of nhmmer if not in path.
  --hmmbuild              Set location of hmmbuild if not in path.
```

## Custom DNA Matrices

nhmmer can be supplied with custom DNA score matrices for assessing hmm match scores. 
Standard NCBI-BLAST matrices such as NUC.4.4 are compatible. (See: ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/NUC.4.4) 

# License

Software provided under MIT license.
