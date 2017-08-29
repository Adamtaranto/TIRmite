# mapmite
Map TIR-pHMM models to genomic sequences for annotation of MITES and complete DNA-Transposons.

## Basic Usage

Example: Report all hits and vaild pairings of TIR_A in target.fasta (interval <= 50000), and write gff annotation file.

```
./mapmite.py --genome target.fasta --hmmFile TIR_A.hmm --gffOut TIR_elements_in_Target.gff3 --maxdist 50000
```

## Options

**Input options**
  - *--genome*
    Path to target genome that will be queried with HMMs.
  - *--hmmDir*
    'Directory containing pre-prepared TIR-pHMMs.
  - *--hmmFile*
    Path to single TIR-pHMM file. Incompatible with "--hmmDir".
  - *--alnDir*
    Path to directory containing only TIR alignments to be converted to HMM.
  - *--alnFile*
    Provide a single TIR alignment to be converted to HMM. Incompatible with "--alnDir".
  - *--alnFormat*
    Alignments provided with "--alnDir" or "--alnFile" are all in this format.
    choices=["clustal","emboss","fasta","fasta-m10","ig","maf","mauve","nexus","phylip","phylip-sequential","phylip-relaxed","stockholm"]


**Pairing heuristics**
  - *--stableReps*
    Number of times to iterate pairing procedure when no additional pairs are found AND remaining unpaired hits > 0.
    default=0


**Output and housekeeping**
  - *--outdir*
    All output files will be written to this directory.
  - *--gffOut*
    GFF3 annotation filename.
  - *--reportTIR*
    Options for reporting TIRs in GFF annotation file.
    choices=[None,'all','paired','unpaired']
    default='all'
  - *--keeptemp*
    If set do not delete temp file directory.
    default=False


**HMMER options**
  - *--cores*
    Set number of cores available to hmmer software.
    default=1
  - *--maxeval*
    Maximum e-value allowed for valid hit. 
    default = 0.001
  - *--maxdist*
    default=None,help='Maximum distance allowed between TIR candidates to consider valid pairing.
  - *--nobias*
    Turn OFF bias correction of scores in nhmmer.
    default=False
  - *--matrix*
    Use custom DNA substitution matrix with nhmmer.


**Non-standard HMMER paths**
  - *--hmmpress*
    Set location of hmmpress if not in path.
  - *--nhmmer*
    Set location of nhmmer if not in path.
  - *--hmmbuild*
    Set location of hmmbuild if not in path.



