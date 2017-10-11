# TIRmite

Build and map profile Hidden Markov Models for Terminal Inverted repeat 
families (TIR-pHMMs) to genomic sequences for annotation of MITES and complete 
DNA-Transposons with variable internal sequence composition.  


TIRmite is packaged with *tSplit* a tool for extraction of terminal repeats 
from complete transposons.

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
  * [tSplit example usage](tsplit-example-usage)
  * [tsplit-LTR](tsplit-ltr)
  * [tsplit-TIR](tsplit-tir)
  * [tSplit options](tsplit-options)
* [License](#license)

# About TIRmite

TIRmite will use profile-HMM models of Terminal Inverted Repeats (TIRs) for 
genome-wide annotation of TIR families. These can be provided by the user or
built from aligned TIRs oriented as 5' outer edge --> 3' inner edge.


Three classes of output are produced:
  1. All significant TIR hit sequences written to fasta (per query HMM).
  2. Candidate elements comprised of paired TIRs are written to fasta (per query HMM).
  3. Genomic annotations of candidate elements and, optionally, TIR hits 
  (paired and unpaired) are written as a single GFF3 file.

# Algorithm overview

  1. Use nhmmer genome with TIR-pHMM
  2. Import all hits below *--maxeval* threshold
  3. For each significant TIR match identify candidate partners, where:
    * Is on the same sequence
    * Hit is in complementary orientation
    * Distance is <= *--maxdist*
  4. Rank candidate partners by distance downstream of positive-strand hits, and upstream of negative-strand hits.
  5. Pair reciprocal top candidate hits 
  6. For unpaired hits, find first unpaired candidate partner and check for reciprocity.
  7. If the first unpaired candidate is non-reciprocal, check for 2nd-order reciprocity (is outbound top-candidate of current candidate reciprocal.)
  8. Iterate steps 6-7 until all TIRs are paired OR number of iterations without new pairing exceeds *--stableReps*

# Options and usage

## Installing TIRmite

Dependencies:  
  - TIR-pHMM build and search
    * [HMMER3](http://hmmer.org)
  - Optional for experimental Bowtie2 mapping mode:
    - [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
    - [samtools](https://github.com/samtools/samtools)
    - [bedtools](http://bedtools.readthedocs.io/en/latest/)
  - Extract terminal repeats from predicted TEs
    * [pymummer](https://pypi.python.org/pypi/pymummer) version >= 0.10.3 with wrapper for nucmer option *--diagfactor*.
    * [MUMmer](http://mummer.sourceforge.net/)
    * [BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (Optional)

Installation options:

```bash
# Install from PyPi:
pip install tirmite

# Clone and install from this repository:
git clone https://github.com/Adamtaranto/TIRmite.git && cd TIRmite && pip install -e .
```

## Example usage

Report all hits and valid pairings of TIR_A in target.fasta (interval <= 50000), 
and write GFF3 annotation file.

```
tirmite --genome target.fasta --hmmFile TIR_A.hmm --gffOut TIR_elements_in_Target.gff3 --maxdist 50000
```

## Standard options

Run `tirmite --help` to view the program's most commonly used options:

```
Usage: tirmite [-h] --genome GENOME [--hmmDir HMMDIR] [--hmmFile HMMFILE]
               [--alnDir ALNDIR] [--alnFile ALNFILE]
               [--alnFormat {clustal,emboss,fasta,fasta-m10,ig,maf,mauve,nexus,phylip,phylip-sequential,phylip-relaxed,stockholm}]
               [--useBowtie2] [--btTIR BTTIR] [--bowtie2 BOWTIE2]
               [--bt2build BT2BUILD] [--samtools SAMTOOLS]
               [--bedtools BEDTOOLS] [--stableReps STABLEREPS]
               [--outdir OUTDIR] [--prefix PREFIX] [--nopairing]
               [--gffOut GFFOUT] [--reportTIR {None,all,paired,unpaired}]
               [--keeptemp] [-v] [--cores CORES] [--maxeval MAXEVAL]
               [--maxdist MAXDIST] [--nobias] [--matrix MATRIX]
               [--hmmpress HMMPRESS] [--nhmmer NHMMER] [--hmmbuild HMMBUILD]

Help:
  -h, --help              Show this help message and exit.

Input options:
  --genome                Path to target genome that will be queried with HMMs.
                            Note: Sequence names must be unique.
                            (required)
  --hmmDir                Directory containing pre-prepared TIR-pHMMs.
  --hmmFile               Path to single TIR-pHMM file. 
                            Incompatible with "--hmmDir".
  --alnDir                Path to the directory containing only TIR alignments to be converted to HMM.
  --alnFile               Provide a single TIR alignment to be converted to HMM. 
                            Incompatible with "--alnDir".
  --alnFormat             Alignments provided with "--alnDir" or "--alnFile" are all in this format.
                            Choices=["clustal","emboss","fasta","fasta-m10","ig","maf","mauve",
                            "nexus","phylip","phylip-sequential","phylip-relaxed","stockholm"]


Alternative search methods:
  --useBowtie2            If set, map short TIR to genome with bowtie2. 
                            Experimental method, potentially useful for very short though highly conserved TIRs where 
                            TIR-pHMM hits return high e-values.
  --btTIR                 Fasta file containing a single TIR to be mapped with bowtie2.
  --bowtie2               Set location of bowtie2 if not in PATH.
  --bt2build              Set location of bowtie2-build if not in PATH.
  --samtools              Set location of samtools if not in PATH.
  --bedtools              Set location of bedtools if not in PATH.
  

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
  --cores                 Set the number of cores available to hmmer software.
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

# Additional tools

# tSplit

Extract terminal repeats from retrotransposons (LTRs) or DNA transposons (TIRs).  

# tSplit algorithm overview

tSplit attempts to identify terminal repeats in transposable elements by 
first aligning each element to itself using nucmer, and then applying a set of 
tuneable heuristics to select an alignment pair most likely to represent an LTR or TIR.  

  1. Exclude all diagonal/self-matches 
  2. If tsplit-LTR: Retain only alignment pairs on the same strand (tandem repeats)
  3. If tsplit-TIR: Retain only alignment pairs on opposite strands (inverse repeats)
  4. Retain pairs for which the 5' match begins within x bases of element start
     and whose 3' match ends within x bases of element end
  5. Exclude alignment pairs which overlap (potential SSRs)
  6. If multiple candidates remain select alignment pair with largest internal segment 
  (i.e. closest to element ends)

# tSplit options and usage  

## tSplit example usage  

TE-splitter contains two programs: tsplit-LTR and tsplit-TIR, for extracting long terminal 
repeats and terminal inverted repeats, respectively. Options are the same
for each.  

## tsplit-LTR 

For each element in *retroelements.fasta* split into internal and external segments. 
Split segments will be written to *LTR_split_TE-splitter_output.fasta* with suffix "_I" 
for internal or "_LTR" for external segments. LTRs must be at least 10bp in length and 
share 80% identity and occur within 10bp of each end of the input element.

```bash
tsplit-LTR -i retroelements.fasta -p LTR_split
```

## tsplit-TIR

For each element in *dna-transposons.fasta* split into internal and external (TIR) segments. 
Split segments will be written to *TIR_split_TE-splitter_output.fasta* with suffix "_I" for 
internal or "_TIR" for external segments. TIRs must be at least 10bp in length and share 80% 
identity and occur within 10bp of each end of the input element. Additionally, synthetic 
MITEs will be constructed by concatenation of left and right TIRs, with internal segments 
excised.

```bash
tsplit-TIR -i dna-transposons.fasta -p TIR_split --makemites
```

## tSplit options

Run `tsplit-LTR --help` or `tsplit-TIR --help` to view the programs' most commonly used 
options:

```
Usage: tsplit-[LTR or TIR] [-h] -i INFILE [-p PREFIX] [-d OUTDIR]
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

# License

Software provided under MIT license.
