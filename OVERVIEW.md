# TIRmite Project Overview

## Purpose

TIRmite is a bioinformatics tool designed for genome-wide annotation of transposable elements (TEs) using profile Hidden Markov Models (HMMs). The tool specializes in identifying and mapping Terminal Inverted Repeats (TIRs) and Long Terminal Repeats (LTRs) of transposon families in genomic sequences. It enables detection of cryptic transposon variants that may have variable internal sequences but conserved terminal repeats.

The key innovation of TIRmite is its ability to work with both symmetrical terminal repeats (TIRs/LTRs) and asymmetrical elements with different conserved features at either end (such as Starship elements).

## Project Structure

```
TIRmite/
├── src/tirmite/          # Main source code
│   ├── cli/              # Command-line interface modules
│   ├── runners/          # External tool wrappers (HMMER, BLAST, MAFFT)
│   ├── utils/            # Utility functions and logging
│   └── tirmitetools.py   # Core algorithm implementation
├── tests/                # Test suite
├── docs/                 # Documentation
├── dist/                 # Distribution packages
├── pyproject.toml        # Project configuration and dependencies
├── environment.yml       # Conda environment specification
└── README.md             # User documentation
```

## Key Modules

### CLI Module (`src/tirmite/cli/`)

The command-line interface provides three main subcommands:

1. **legacy** (`cli.py`, `legacy.py`) - Original TIRmite workflow combining HMM search with pairing
2. **seed** (`hmm_build.py`) - Build HMM models from seed sequences
3. **pair** (`hmm_pair.py`) - Pair precomputed nhmmer hits

### Core Algorithm (`tirmitetools.py`)

The core module implements the pairing algorithm and data processing:

- **Data Import Functions**:
  - `import_nhmmer()` - Parse nhmmer tabular output
  - `import_BED()` - Read BED format hit files
  - `table2dict()` - Convert hit tables to structured dictionaries

- **Filtering Functions**:
  - `filterHitsLen()` - Filter hits by minimum coverage
  - `filterHitsEval()` - Filter hits by E-value threshold

- **Pairing Algorithm**:
  - `parseHits()` - Identify candidate partners for each terminus hit
  - `getPairs()` - Perform reciprocal pairing of termini
  - `iterateGetPairs()` - Iteratively pair unpaired hits
  - `isfirstUnpaired()` - Check reciprocity of candidate pairings

- **Output Generation**:
  - `extractTIRs()` - Extract terminal repeat sequences
  - `fetchElements()` - Retrieve complete element sequences
  - `writeTIRs()` - Write TIR sequences to FASTA
  - `writeElements()` - Write paired elements to FASTA
  - `gffWrite()` - Generate GFF3 format annotations

### External Tool Wrappers (`runners/`)

- **hmmer_wrappers.py** - Interface with HMMER3 suite (nhmmer for searching)
- **wrapping.py** - MAFFT wrapper for multiple sequence alignment
- **runBlastn.py** - BLAST+ integration for sequence comparison

### Utilities (`utils/`)

- **utils.py** - General utility functions (ID cleaning, file handling)
- **logs.py** - Logging configuration and management

## Core Algorithm Overview

TIRmite uses a sophisticated pairing algorithm to identify complete transposon elements:

1. **HMM Search**: Query genome sequences with terminal repeat HMMs using nhmmer
2. **Hit Import**: Load all hits below E-value threshold (--maxeval)
3. **Candidate Identification**: For each terminus hit, find potential partners that are:
   - On the same sequence
   - In complementary orientation (+ paired with -)
   - Within maximum distance (--maxdist)
   - Meeting minimum coverage (--mincov)
4. **Ranking**: Order candidates by proximity (downstream for + strand, upstream for - strand)
5. **Reciprocal Pairing**: Pair hits that are mutual top candidates
6. **Iterative Pairing**: Continue pairing unpaired hits with reciprocity checks until convergence (--stableReps iterations without new pairs)
7. **Output Generation**: Produce FASTA files and GFF3 annotations

## Basic Use Cases

### 1. Standard TIR Annotation

Annotate a genome with an existing TIR HMM model:

```bash
tirmite legacy --genome target.fasta \
               --hmmFile TIR_model.hmm \
               --gffOut annotations.gff3 \
               --maxdist 10000 \
               --mincov 0.4 \
               --orientation F,R
```

### 2. Build HMM from Seed Sequences

Create an HMM model from aligned TIR sequences:

```bash
tirmite seed --left-seed left_TIRs.fa \
             --model-name MyTE \
             --genome genome.fa
```

### 3. Pair Precomputed Hits

Work with existing nhmmer results:

```bash
tirmite pair --genome genome.fa \
             --nhmmerFile precomputed_hits.out \
             --hmmFile model.hmm \
             --maxdist 15000
```

### 4. Asymmetric Elements

Annotate Starship-like elements with different termini:

```bash
tirmite seed --left-seed left_terminus.fa \
             --right-seed right_terminus.fa \
             --model-name Starship1 \
             --genome genome.fa
```

## Dependencies

**Required**:
- Python >= 3.8
- BioPython >= 1.70
- pandas >= 0.23.4
- pyfaidx
- pyhmmer
- rich

**External Tools**:
- HMMER3 (nhmmer)
- MAFFT (for alignment)
- BLAST+ (optional, for additional analysis)

## Output Files

1. **FASTA files**: Terminal repeat sequences and complete elements (per HMM model)
2. **GFF3 file**: Genome annotations of elements and optionally individual HMM hits
3. **Log files**: Detailed processing information and statistics

## Installation Methods

- PyPI: `pip install tirmite`
- Bioconda: `conda install -c bioconda tirmite`
- Development: `pip install git+https://github.com/Adamtaranto/TIRmite.git`

## Key Features

- Profile HMM-based transposon annotation
- Support for symmetric (TIR/LTR) and asymmetric elements
- Flexible pairing algorithm with reciprocity checks
- Multiple output formats (FASTA, GFF3)
- Configurable distance and coverage thresholds
- Integration with standard bioinformatics tools (HMMER, MAFFT, BLAST)
