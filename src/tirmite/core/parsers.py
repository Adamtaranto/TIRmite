"""Readers for the hit formats TIRmite consumes.

Converts nhmmer ``--tblout``, BLAST tabular (outfmt 6) and BED files into one
normalised hit table, and converts alignments between formats.

Notes
-----
Every importer returns a ``pandas`` DataFrame whose values are **strings**,
with the columns ``model, target, hitStart, hitEnd, strand, evalue, score,
bias, hmmStart, hmmEnd``. Reverse-strand coordinates are swapped on import so
that ``hitStart < hitEnd`` always holds, and strand is carried separately.
Downstream consumers rely on both properties.
"""

import glob
import logging
import os
from typing import Optional

from Bio import AlignIO  # type: ignore[import-not-found]
import pandas as pd  # type: ignore[import-untyped]

# Library modules acquire a named logger and attach a NullHandler, so that
# importing TIRmite as a library emits nothing until the host application
# configures logging. Handler setup belongs to the CLI, not here.
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Canonical sort order for a hit table. Coordinates are compared numerically,
# everything else lexicographically.
_HIT_SORT_KEYS = ['model', 'target', 'hitStart', 'hitEnd', 'strand']
_NUMERIC_SORT_KEYS = ('hitStart', 'hitEnd')


def sort_hit_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sort a hit table by model, target, genomic coordinates and strand.

    Parameters
    ----------
    df : pandas.DataFrame
        Hit table with at least the columns ``model``, ``target``,
        ``hitStart``, ``hitEnd`` and ``strand``.

    Returns
    -------
    pandas.DataFrame
        The table sorted in place of the input, with the index reset.

    Notes
    -----
    ``hitStart`` and ``hitEnd`` are stored as **strings**, so sorting on them
    directly compares lexicographically: ``'1000'`` sorts before ``'999'``.
    They are therefore sorted on temporary integer columns, which are dropped
    before returning so that the table's all-string dtype contract is
    preserved -- ``table2dict``, the pairing engine and the GFF writer all
    depend on it.

    Row order is not merely cosmetic. ``filter_best_model_per_locus`` walks the
    group in order and removes greedily, so the ordering decides which of two
    mutually-eliminating hits survives.

    Examples
    --------
    >>> sorted_hits = sort_hit_table(hit_table)
    """
    if df.empty:
        return df.reset_index(drop=True)

    numeric_columns = {}
    for column in _NUMERIC_SORT_KEYS:
        helper = f'__{column}_sortkey'
        # errors='coerce' keeps a malformed coordinate from aborting the sort;
        # such rows sort last rather than raising.
        numeric_columns[helper] = pd.to_numeric(df[column], errors='coerce')

    df = df.assign(**numeric_columns)
    sort_keys = [
        f'__{key}_sortkey' if key in _NUMERIC_SORT_KEYS else key
        for key in _HIT_SORT_KEYS
    ]
    df = df.sort_values(sort_keys, ascending=[True] * len(sort_keys))
    df = df.drop(columns=list(numeric_columns))
    return df.reset_index(drop=True)


def convertAlign(
    alnDir: Optional[str] = None,
    alnFile: Optional[str] = None,
    inFormat: str = 'fasta',
    tempDir: Optional[str] = None,
) -> str:
    """
    Convert input sequence alignments to Stockholm format for HMMER.

    Parameters
    ----------
    alnDir : str, optional
        Glob pattern to match multiple alignment files (e.g., "path/*.fasta").
        Mutually exclusive with alnFile.
    alnFile : str, optional
        Path to a single alignment file to convert.
        Mutually exclusive with alnDir.
    inFormat : str, default 'fasta'
        Format of input alignment files (e.g., 'fasta', 'clustal', 'phylip').
    tempDir : str, optional
        Base directory for creating temporary alignment output directory.
        If None, uses current working directory.

    Returns
    -------
    str
        Path to directory containing converted Stockholm format alignments.

    Notes
    -----
    Creates a 'temp_aln' subdirectory in tempDir (or cwd) to store output files.
    Each alignment file is converted and saved with a .stockholm extension.
    """
    # Construct out model path
    if tempDir:
        alnOutDir = os.path.join(os.path.abspath(tempDir), 'temp_aln')
    else:
        alnOutDir = os.path.join(os.getcwd(), 'temp_aln')
    # Create if does not exist
    if not os.path.isdir(alnOutDir):
        os.makedirs(alnOutDir)
    # Get list of alignment files to process
    if alnFile:
        alignments = [alnFile]
    elif alnDir:
        alignments = glob.glob(alnDir)
    else:
        raise ValueError('Either alnFile or alnDir must be provided')
    # Do conversion on each
    for infile in alignments:
        # Log file being processed
        logger.info(f'Converting alignment file: {infile}')
        # Get basename
        inBase = os.path.splitext(os.path.basename(infile))[0]
        # Make outpath
        outAln = os.path.join(alnOutDir, inBase + '.sto')
        # Open files
        input_handle = open(infile, 'r')
        output_handle = open(outAln, 'w')
        # Read alignment
        alignments = AlignIO.parse(input_handle, inFormat)
        # Log count of sequences in 'alignments' object generated with Align.IO.parse
        logger.info(f'Number of sequences in alignment: {len({alignments})}')
        # Write as stockholm
        logger.info(f'Writing converted alignment to: {outAln}')
        AlignIO.write(alignments, output_handle, 'stockholm')
        # Close handles
        output_handle.close()
        input_handle.close()
    return alnOutDir


def import_nhmmer(
    infile: Optional[str] = None,
    hitTable: Optional[pd.DataFrame] = None,
    prefix: Optional[str] = None,
) -> pd.DataFrame:
    """
    Parse nhmmer tabular output file into a pandas DataFrame.

    Parameters
    ----------
    infile : str
        Path to nhmmer tabular output file containing hit records.
    hitTable : pandas.DataFrame, optional
        Existing DataFrame of hits to concatenate with new hits.
    prefix : str, optional
        Prefix to add to model names (currently unused but reserved for future use).

    Returns
    -------
    pandas.DataFrame
        DataFrame containing parsed hit records with columns:
        model, target, hitStart, hitEnd, strand, evalue, score, bias,
        hmmStart, hmmEnd. Sorted by model, target, hitStart, hitEnd, strand.

    Notes
    -----
    Handles strand orientation: for reverse strand hits (-), hitStart and
    hitEnd are swapped to ensure hitStart < hitEnd in genomic coordinates.
    """
    hitRecords = []
    if not infile:
        raise ValueError('infile parameter is required')
    with open(infile, 'r') as f:
        for line in f.readlines():
            li = line.strip()
            if not li.startswith('#'):
                li_split = li.split()
                if li_split[11] == '+':
                    hitRecords.append(
                        {
                            'target': li_split[0],
                            'model': li_split[2],
                            'hmmStart': li_split[4],
                            'hmmEnd': li_split[5],
                            'hitStart': li_split[6],
                            'hitEnd': li_split[7],
                            'strand': li_split[11],
                            'evalue': li_split[12],
                            'score': li_split[13],
                            'bias': li_split[14],
                        }
                    )
                elif li_split[11] == '-':
                    hitRecords.append(
                        {
                            'target': li_split[0],
                            'model': li_split[2],
                            'hmmStart': li_split[4],
                            'hmmEnd': li_split[5],
                            'hitStart': li_split[7],
                            'hitEnd': li_split[6],
                            'strand': li_split[11],
                            'evalue': li_split[12],
                            'score': li_split[13],
                            'bias': li_split[14],
                        }
                    )
    # Define expected columns
    cols = [
        'model',
        'target',
        'hitStart',
        'hitEnd',
        'strand',
        'evalue',
        'score',
        'bias',
        'hmmStart',
        'hmmEnd',
    ]
    # Convert list of dicts into dataframe
    df = pd.DataFrame(hitRecords)
    # Handle empty results - return DataFrame with correct columns but no rows
    if df.empty:
        df = pd.DataFrame(columns=cols)
    else:
        # Reorder columns
        df = df.loc[:, cols]
    if hitTable is not None:
        # If an existing table was passed, concatenate
        df = pd.concat([df, hitTable], ignore_index=True)
    # Sort hits by HMM, Chromosome, location, and strand (coordinates numerically)
    df = sort_hit_table(df)
    # if prefix:
    #    df['model'] = str(prefix) + '_' + df['model'].astype(str)
    return df


def import_BED(
    infile: Optional[str] = None,
    hitTable: Optional[pd.DataFrame] = None,
    prefix: Optional[str] = None,
) -> pd.DataFrame:
    """
    Parse BED format file of TIR hits into a pandas DataFrame.

    Parameters
    ----------
    infile : str
        Path to BED format file with TIR hit coordinates.
    hitTable : pandas.DataFrame, optional
        Existing DataFrame of hits to concatenate with new hits.
    prefix : str, optional
        Prefix to add to model names (currently unused but reserved for future use).

    Returns
    -------
    pandas.DataFrame
        DataFrame containing parsed hit records with columns:
        model, target, hitStart, hitEnd, strand, evalue, score, bias,
        hmmStart, hmmEnd. Score and bias set to 'NA' for BED format.
        Sorted by model, target, hitStart, hitEnd, strand.

    Notes
    -----
    Expected BED format: chromosome, start, end, name, evalue, strand.
    Since BED format lacks HMM alignment coordinates, hmmStart and hmmEnd
    are set to 'NA'.
    """
    # Format: Chrm, start, end, name, evalue, strand
    hitRecords = []
    if not infile:
        raise ValueError('infile parameter is required')
    with open(infile, 'r') as f:
        for line in f.readlines():
            li = line.strip()
            if not li.startswith('#'):
                li_split = li.split()
                hitRecords.append(
                    {
                        'target': li_split[0],
                        'model': li_split[3],
                        'hmmStart': 'NA',
                        'hmmEnd': 'NA',
                        'hitStart': li_split[1],
                        'hitEnd': li_split[2],
                        'strand': li_split[5],
                        'evalue': li_split[4],
                        'score': 'NA',
                        'bias': 'NA',
                    }
                )
    # Define expected columns
    cols = [
        'model',
        'target',
        'hitStart',
        'hitEnd',
        'strand',
        'evalue',
        'score',
        'bias',
        'hmmStart',
        'hmmEnd',
    ]
    # Convert list of dicts into dataframe
    df = pd.DataFrame(hitRecords)
    # Handle empty results - return DataFrame with correct columns but no rows
    if df.empty:
        df = pd.DataFrame(columns=cols)
    else:
        # Reorder columns
        df = df.loc[:, cols]
    if hitTable is not None:
        # If an existing table was passed, concatenate
        df = pd.concat([df, hitTable], ignore_index=True)
    # Sort hits by HMM, Chromosome, location, and strand (coordinates numerically)
    df = sort_hit_table(df)
    # if prefix:
    #    df['model'] = str(prefix) + '_' + df['model'].astype(str)
    return df


def import_blast(
    infile: Optional[str] = None,
    hitTable: Optional[pd.DataFrame] = None,
    prefix: Optional[str] = None,
    query_name: Optional[str] = None,
) -> pd.DataFrame:
    """
    Parse BLAST tabular output file into a pandas DataFrame.

    Parameters
    ----------
    infile : str
        Path to BLAST tabular output file (outfmt 6 or 7).
    hitTable : pandas.DataFrame, optional
        Existing DataFrame of hits to concatenate with new hits.
    prefix : str, optional
        Prefix to add to model names (currently unused but reserved for future use).
    query_name : str, optional
        Name to use for the query/model. If not provided, uses the query ID from first hit.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing parsed hit records with columns:
        model, target, hitStart, hitEnd, strand, evalue, score, bias,
        hmmStart, hmmEnd. Sorted by model, target, hitStart, hitEnd, strand.

    Notes
    -----
    Expected BLAST tabular format (outfmt 6):
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

    For reverse strand hits (sstart > send), coordinates are swapped and strand set to '-'.
    The query alignment coordinates (qstart, qend) are stored as hmmStart, hmmEnd.
    """
    hitRecords = []
    if not infile:
        raise ValueError('infile parameter is required')

    with open(infile, 'r') as f:
        for line in f.readlines():
            li = line.strip()
            # Skip comments and empty lines
            if not li or li.startswith('#'):
                continue

            li_split = li.split('\t')
            # BLAST tabular format has 12 columns minimum
            if len(li_split) < 12:
                logger.warning(f'Skipping malformed BLAST line: {li}')
                continue

            query_id = li_split[0]
            subject_id = li_split[1]
            qstart = int(li_split[6])
            qend = int(li_split[7])
            sstart = int(li_split[8])
            send = int(li_split[9])
            evalue = float(li_split[10])
            bitscore = float(li_split[11])

            # Determine model name for this hit
            # If query_name parameter was provided, use it for all hits
            # Otherwise, use the query ID from this specific hit
            hit_model_name = query_name if query_name is not None else query_id

            # Determine strand based on subject coordinates
            if sstart <= send:
                # Forward strand
                strand = '+'
                hitStart = sstart
                hitEnd = send
            else:
                # Reverse strand (coordinates are reversed in BLAST output)
                strand = '-'
                hitStart = send
                hitEnd = sstart

            hitRecords.append(
                {
                    'target': subject_id,
                    'model': hit_model_name,
                    'hmmStart': str(qstart),
                    'hmmEnd': str(qend),
                    'hitStart': str(hitStart),
                    'hitEnd': str(hitEnd),
                    'strand': strand,
                    'evalue': str(evalue),
                    'score': str(bitscore),
                    'bias': 'NA',
                }
            )

    # Define expected columns
    cols = [
        'model',
        'target',
        'hitStart',
        'hitEnd',
        'strand',
        'evalue',
        'score',
        'bias',
        'hmmStart',
        'hmmEnd',
    ]

    # Convert list of dicts into dataframe
    df = pd.DataFrame(hitRecords)

    # Handle empty results - return DataFrame with correct columns but no rows
    if df.empty:
        df = pd.DataFrame(columns=cols)
    else:
        # Reorder columns
        df = df.loc[:, cols]

    if hitTable is not None:
        # If an existing table was passed, concatenate
        # Order: new data (df) first, existing (hitTable) second (matches import_nhmmer)
        df = pd.concat([df, hitTable], ignore_index=True)

    # Sort hits by model, target, location, and strand (coordinates numerically)
    df = sort_hit_table(df)

    return df


def detect_input_format(infile: str) -> str:
    """
    Detect whether input file is BLAST or nhmmer format.

    Parameters
    ----------
    infile : str
        Path to input file to analyze.

    Returns
    -------
    str
        Format type: 'blast', 'nhmmer', or 'unknown'.

    Notes
    -----
    Detection heuristics, applied to the first non-comment line:

    - BLAST (outfmt 6): tab-delimited, 12+ columns, with numeric subject
      coordinates in columns 8 and 9 and a numeric e-value in column 10.
    - nhmmer (``--tblout``): whitespace-delimited, 15+ columns, with the
      strand in column 11 and a numeric e-value in column 12.

    The nhmmer column indices deliberately match those used by
    :func:`import_nhmmer`. They previously did not: this function looked for
    the strand in column 9 and the e-value in column 15, which are "env to"
    and the description field respectively, so a genuine nhmmer file was
    reported as ``'unknown'``.
    """
    try:
        with open(infile, 'r') as f:
            for line in f:
                li = line.strip()
                # Skip comments and empty lines
                if not li or li.startswith('#'):
                    continue

                # Found first data line - analyze it
                # Try tab-delimited first (BLAST)
                tab_split = li.split('\t')
                if len(tab_split) >= 12:
                    # Likely BLAST format - check if columns 8-9 are numeric (subject coords)
                    try:
                        int(tab_split[8])
                        int(tab_split[9])
                        float(tab_split[10])  # evalue
                        return 'blast'
                    except (ValueError, IndexError):
                        pass

                # Try whitespace-delimited (nhmmer --tblout). Column order is
                # target, accession, query, accession, hmmfrom, hmm to,
                # alifrom, ali to, envfrom, env to, sq len, strand, E-value,
                # score, bias, description. The description is optional and
                # may itself contain spaces, so require only 15 fields.
                space_split = li.split()
                if len(space_split) >= 15:
                    try:
                        if space_split[11] in ('+', '-'):
                            float(space_split[12])  # E-value
                            return 'nhmmer'
                    except (ValueError, IndexError):
                        pass

                # If we got here, couldn't determine format from first line
                logger.warning('Could not determine format from first data line')
                return 'unknown'

    except Exception as e:
        logger.error(f'Error reading file {infile}: {e}')
        return 'unknown'

    # Empty file or no data lines
    return 'unknown'
