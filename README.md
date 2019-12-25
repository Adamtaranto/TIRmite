[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<p align="center">
  <img src="docs/tirmite_hexlogo.jpg"  width="350" height="350" title="tirmite_hex">
</p>

# TIRmite

Build and map profile Hidden Markov Models for Terminal Inverted Repeat 
families (TIR-pHMMs) to genomic sequences for annotation of MITES and complete 
DNA-Transposons with variable internal sequence composition.  


TIRmite is packaged with *tSplit* a tool for extraction of terminal repeats 
from complete transposons.

Current version: 1.1.4

# Table of contents

* [About TIRmite](#about-tirmite)
* [Algorithm overview](#algorithm-overview)
* [Options and usage](#options-and-usage)
    * [Installing TIRmite](#installing-tirmite)
    * [Example usage](#example-usage)
    * [Standard options](#standard-options)
    * [Custom DNA matrices](#custom-dna-matrices)
* [Additional tools](additional-tools)
    * [tSplit](tsplit)
    * [tSplit algorithm overview](tsplit-algorithm-overview)
    * [tSplit options and usage](tsplit-options-and-usage)
* [Issues](#issues)
* [License](#license)
* [Logo](#logo)


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

TIRmite requires Python >= v3.6

Dependencies:  
  - TIR-pHMM build and search
    * [HMMER3](http://hmmer.org)
  - Extract terminal repeats from predicted TEs
    * [pymummer](https://pypi.python.org/pypi/pymummer) version >= 0.10.3 with wrapper for nucmer option *--diagfactor*.
    * [MUMmer](http://mummer.sourceforge.net/)
    * [BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (Optional)

Installation options:  

Clone from this repository and install as a local Python package.

```bash
% git clone https://github.com/Adamtaranto/TIRmite.git && cd TIRmite && pip install -e .
```

Install from PyPi.

```bash
% pip install tirmite
```

Install from Bioconda.
```bash
% conda install -c bioconda tirmite
```

Test installation.

```bash
# Print version number and exit.
% tirmite --version
tirmite 1.1.4

# Get usage information
% tirmite --help
```

### Example usage

Report all hits and valid pairings of TIR_A in target.fasta (interval <= 10000, hits cover > 40% len of hmm model), 
and write GFF3 annotation file.

```bash
% tirmite --genome target.fasta --hmmFile TIR_A.hmm --gffOut TIR_elements_in_Target.gff3 --maxdist 10000 --mincov 0.4
```

### Standard options

Run `tirmite --help` to view the program's most commonly used options:

```
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

## Additional tools

### tSplit

Extract Terminal Inverted Repeats (TIRs) DNA transposons.  

### tSplit algorithm overview

tSplit attempts to identify terminal repeats in transposable elements by 
first aligning each element to itself using nucmer, and then applying a set of 
tuneable heuristics to select an alignment pair most likely to represent a TIR.  

  1. Exclude all diagonal/self-matches
  2. If tsplit-TIR: Retain only alignment pairs on opposite strands (inverse repeats)
  3. Retain pairs for which the 5' match begins within x bases of element start
     and whose 3' match ends within x bases of element end
  4. Exclude alignment pairs which overlap (potential SSRs)
  5. If multiple candidates remain select alignment pair with largest internal segment 
  (i.e. closest to element ends)

### tSplit options and usage  

### tSplit example usage  

For each element in *dna-transposons.fasta* split into internal and external (TIR) segments. 
Split segments will be written to *TIR_split_TE-splitter_output.fasta* with suffix "_I" for 
internal or "_TIR" for external segments. TIRs must be at least 10bp in length and share 80% 
identity and occur within 10bp of each end of the input element. Additionally, synthetic 
MITEs will be constructed by concatenation of left and right TIRs, with internal segments 
excised.


```bash
% tsplit-TIR -i dna-transposons.fasta -p TIR_split
```

### tSplit options

Run `tsplit-TIR --help` to view the programs' most commonly used 
options:

```
Usage: tsplit-TIR [-h] -i INFILE [-p PREFIX] [-d OUTDIR]
                        [--splitmode {all,split,internal,external,None}]
                        [--makemites] [--keeptemp] [-v] [-m MAXDIST]
                        [--minid MINID] [--minterm MINTERM] [--minseed MINSEED]
                        [--diagfactor DIAGFACTOR] [--method {blastn,nucmer}]

Help:
  -h, --help         Show this help message and exit.

Input:
  -i, --infile       Multifasta containing complete elements. 
                       (Required)  

Output:
  -p, --prefix       All output files begin with this string.  (Default:[infile basename])  
  -d, --outdir       Write output files to this directory. (Default: cwd)  
  --keeptemp         If set do not remove temp directory on completion.
  -v, --verbose      If set, report progress.

Report settings:
  --splitmode        Options: {all,split,internal,external,None} 
                       all = Report input sequence as well as internal and external segments.  
                       split = Report internal and external segments after splitting.  
                       internal = Report only internal segments.  
                       external = Report only terminal repeat segments.  
                       None = Only report synthetic MITES (when --makemites is also set).  
                       (Default: split)  
  --makemites        Experimental function: Attempt to construct synthetic MITE sequences from TIRs by concatenating 
                       5' and 3' TIRs. Available only in 'tsplit-TIR' mode 

Alignment settings:
  --method          Select alignment tool. Note: blastn may perform better on very short high-identity TRs,
                      while nucmer is more robust to small indels.
                      Options: {blastn,nucmer} 
                      (Default: nucmer)
  --minid           Minimum identity between terminal repeat pairs. As float. 
                      (Default: 80.0)  
  --minterm         Minimum length for a terminal repeat to be considered.  
                      Equivalent to nucmer "--mincluster" 
                      (Default: 10)  
  -m, --maxdist     Terminal repeat candidates must be no more than this many bases from ends of an input element. 
                      Note: Increase this value if you suspect that your element is nested within some flanking sequence. 
                      (Default: 10)
  --minseed         Minimum length of a maximal exact match to be included in final match cluster. 
                      Equivalent to nucmer "--minmatch". 
                      (Default: 5)
  --diagfactor      Maximum diagonal difference factor for clustering of matches within nucmer, 
                      i.e. diagonal difference / match separation 
                      (default 0.20) 
                      Note: Increase value for greater tolerance of indels between terminal repeats.
```

## Issues

Submit feedback to the [Issue Tracker](https://github.com/Adamtaranto/TIRmite/issues)

## License

Software provided under MIT license.

## Logo

Termite hex-sticker was designed by [@Super_Coleider](www.instagram.com/Super_Coleider).