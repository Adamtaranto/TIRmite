"""
TIRmite core tools for transposon terminal repeat detection and pairing.

This module implements the core algorithms for:
- Parsing nhmmer search results
- Filtering hits by coverage and e-value
- Pairing terminal inverted repeats
- Extracting and writing transposon sequences
- Supporting both symmetric and asymmetric pairing modes

The pairing algorithms use reciprocal best-match approaches with
configurable strand orientations and distance constraints.
"""

from collections import namedtuple
import glob
import logging
from operator import attrgetter
import os
from typing import Any, Dict, List, NamedTuple, Optional, Set, Tuple, Union

from Bio import AlignIO, Seq, SeqIO  # type: ignore[import-not-found]
from Bio.SeqRecord import SeqRecord  # type: ignore[import-not-found]
import pandas as pd  # type: ignore[import-untyped]

from tirmite.utils.extract import (
    SequenceSource,
    annotate,
    clamp_region,
    fetch_region_padded,
    fetch_sequence,
    make_source,
)
from tirmite.utils.utils import cleanID


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
        logging.info(f'Converting alignment file: {infile}')
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
        logging.info(f'Number of sequences in alignment: {len({alignments})}')
        # Write as stockholm
        logging.info(f'Writing converted alignment to: {outAln}')
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
    # Sort hits by HMM, Chromosome, location, and strand
    df = df.sort_values(
        ['model', 'target', 'hitStart', 'hitEnd', 'strand'],
        ascending=[True, True, True, True, True],
    )
    # Reindex
    df = df.reset_index(drop=True)
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
    # Sort hits by HMM, Chromosome, location, and strand
    df = df.sort_values(
        ['model', 'target', 'hitStart', 'hitEnd', 'strand'],
        ascending=[True, True, True, True, True],
    )
    # Reindex
    df = df.reset_index(drop=True)
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
                logging.warning(f'Skipping malformed BLAST line: {li}')
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

    # Sort hits by model, target, location, and strand
    df = df.sort_values(
        ['model', 'target', 'hitStart', 'hitEnd', 'strand'],
        ascending=[True, True, True, True, True],
    )

    # Reindex
    df = df.reset_index(drop=True)

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
    Detection heuristics:
    - BLAST: Tab-delimited, typically 12+ columns, numeric subject coordinates
    - nhmmer: Space-delimited, has specific column patterns with model names
    - Checks first non-comment line for format signature
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

                # Try space-delimited (nhmmer)
                space_split = li.split()
                if len(space_split) >= 16:
                    # nhmmer has specific pattern: column 9 is strand (+/-)
                    # and column 15 is evalue (scientific notation)
                    try:
                        if space_split[9] in ['+', '-']:
                            float(space_split[15])  # evalue
                            return 'nhmmer'
                    except (ValueError, IndexError):
                        pass

                # If we got here, couldn't determine format from first line
                logging.warning('Could not determine format from first data line')
                return 'unknown'

    except Exception as e:
        logging.error(f'Error reading file {infile}: {e}')
        return 'unknown'

    # Empty file or no data lines
    return 'unknown'


def filterHitsLen(
    hmmDB: Optional[str] = None,
    mincov: Optional[float] = None,
    hitTable: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """
    Filter hit table to remove hits with insufficient model coverage.

    Parameters
    ----------
    hmmDB : str
        Path to directory containing HMM files (.hmm extension).
    mincov : float
        Minimum coverage threshold as a fraction of model length (0.0 to 1.0).
        Hits shorter than (model_length * mincov) will be removed.
    hitTable : pandas.DataFrame
        DataFrame of hits to filter, must contain 'model', 'hitStart', 'hitEnd' columns.

    Returns
    -------
    pandas.DataFrame
        Filtered DataFrame with short hits removed.

    Notes
    -----
    Extracts model lengths from HMM files by parsing LENG and NAME fields.
    For each model, calculates minimum acceptable hit length as model_length * mincov.
    """
    if not hmmDB or not mincov or hitTable is None:
        raise ValueError('hmmDB, mincov, and hitTable are required parameters')

    modelLens = {}
    for hmm in glob.glob(os.path.join(hmmDB, '*.hmm')):
        hmmLen = None
        hmmName = None
        with open(hmm, 'r') as f:
            for line in f.readlines():
                li = line.strip()
                if li.startswith('LENG'):
                    hmmLen = int(li.split()[1])
                if li.startswith('NAME'):
                    hmmName = str(li.split()[1])
            if hmmLen and hmmName:
                modelLens[hmmName] = hmmLen
    for model in modelLens.keys():
        minlen = modelLens[model] * mincov
        hitTable = hitTable.loc[  # type: ignore[assignment]
            ~(
                (hitTable['model'] == model)
                & (
                    (hitTable['hitEnd'].astype(int) - hitTable['hitStart'].astype(int))
                    + 1
                    < minlen  # type: ignore[operator]
                )
            )
        ]
    return hitTable


def filterHitsEval(
    maxeval: Optional[float] = None, hitTable: Optional[pd.DataFrame] = None
) -> pd.DataFrame:
    """
    Filter hit table to remove hits with e-values exceeding threshold.

    Parameters
    ----------
    maxeval : float
        Maximum acceptable e-value threshold. Hits with e-values greater than
        this value will be removed.
    hitTable : pandas.DataFrame
        DataFrame of hits to filter, must contain 'evalue' column.

    Returns
    -------
    pandas.DataFrame
        Filtered DataFrame containing only hits with e-value < maxeval.
    """
    if maxeval is None or hitTable is None:
        raise ValueError('maxeval and hitTable are required parameters')
    hitTable = hitTable.loc[((hitTable['evalue'].astype(float)) < float(maxeval))]
    return hitTable


def table2dict(
    hitTable: pd.DataFrame,
) -> Tuple[Dict[str, Dict[str, List[Any]]], Dict[str, Dict[int, Dict[str, Any]]]]:
    """
    Convert pandas DataFrame of hits into nested dictionaries for pairing analysis.

    Parameters
    ----------
    hitTable : pandas.DataFrame
        DataFrame containing hit records with columns: model, target, hitStart,
        hitEnd, strand, evalue, hmmStart, hmmEnd.

    Returns
    -------
    hitsDict : dict
        Nested dictionary structure: hitsDict[model][chromosome] = [list of hit records]
        where hit records are namedtuples containing hit information.
    hitIndex : dict
        Nested dictionary structure: hitIndex[model][row_index] = {rec, partner, candidates}
        where 'rec' is the hit record namedtuple, 'partner' tracks pairing status,
        and 'candidates' is list of potential pairing partners.

    Notes
    -----
    Creates namedtuple 'Elem' with fields: model, target, hitStart, hitEnd, strand,
    idx (DataFrame row index), evalue. The idx field links back to the original DataFrame.
    """
    # Set up empty dict
    hitsDict: Dict[str, Dict[str, List[Any]]] = {}
    hitIndex: Dict[str, Dict[int, Dict[str, Any]]] = {}
    # Populate keys from dataframe
    for hmm in hitTable.model.unique():
        hitsDict[hmm] = {}
        hitIndex[hmm] = {}
        for chr in hitTable[hitTable['model'] == hmm].target.unique():
            hitsDict[hmm][chr] = []
    # Set up named tuple
    hitTup = namedtuple(
        'hitTup',
        ['model', 'target', 'hitStart', 'hitEnd', 'strand', 'idx', 'evalue'],  # type: ignore[name-match]
    )
    # Add each record to dicts
    for row in hitTable.iterrows():
        record = hitTup(
            row[1].model,
            row[1].target,
            int(row[1].hitStart),
            int(row[1].hitEnd),
            row[1].strand,
            row[0],
            row[1].evalue,
        )
        # Log hit for model on chromosome
        hitsDict[row[1].model][row[1].target].append(record)
        # Populate tracker - FIX: use row[1].model not hmm
        hitIndex[row[1].model][row[0]] = {  # type: ignore[index]
            'rec': record,
            'partner': None,
            'candidates': [],
        }
    # Return master rec object and pairing tracker
    return hitsDict, hitIndex


def parseHits(
    hitsDict: Optional[Dict[str, Dict[str, List[Any]]]] = None,
    hitIndex: Optional[Dict[str, Dict[int, Dict[str, Any]]]] = None,
    maxDist: Optional[int] = None,
) -> Dict[str, Dict[int, Dict[str, Any]]]:
    """
    Identify potential pairing partners for each hit based on strand and distance.

    Parameters
    ----------
    hitsDict : dict
        Nested dictionary of hits: hitsDict[model][chromosome] = [hit_records].
    hitIndex : dict
        Nested dictionary tracking pairing: hitIndex[model][idx] = {rec, partner, candidates}.
    maxDist : int, optional
        Maximum distance in base pairs between paired TIR elements.
        If None, uses infinite distance (no distance constraint).

    Returns
    -------
    dict
        Updated hitIndex with populated 'candidates' lists for each hit.

    Notes
    -----
    For forward strand (+) hits, searches for downstream reverse strand (-) partners.
    For reverse strand (-) hits, searches for upstream forward strand (+) partners.
    Candidates are sorted by proximity to the reference hit.
    """
    assert hitsDict is not None, 'hitsDict cannot be None'
    assert hitIndex is not None, 'hitIndex cannot be None'
    if not maxDist:
        maxDist = float('inf')  # type: ignore[assignment]
    for hmm in hitIndex.keys():
        for UID in hitIndex[hmm].keys():
            ref = hitIndex[hmm][UID]['rec']
            if ref.strand == '+':
                for localhit in hitsDict[ref.model][ref.target]:
                    if (
                        localhit.strand == '-'
                        and localhit.hitStart >= ref.hitEnd
                        and localhit.hitStart - ref.hitEnd <= maxDist
                    ):
                        hitIndex[hmm][UID]['candidates'].append(localhit)
                # Sort candidate hit records from low to high on hitStart vals
                hitIndex[hmm][UID]['candidates'] = sorted(
                    hitIndex[hmm][UID]['candidates'],
                    key=attrgetter('hitStart', 'hitEnd'),
                )
            if ref.strand == '-':
                for localhit in hitsDict[ref.model][ref.target]:
                    if (
                        localhit.strand == '+'
                        and localhit.hitEnd <= ref.hitStart
                        and ref.hitStart - localhit.hitEnd <= maxDist
                    ):
                        hitIndex[hmm][UID]['candidates'].append(localhit)
                # Sort candidate hit records from high to low on hitEnd values
                hitIndex[hmm][UID]['candidates'] = sorted(
                    hitIndex[hmm][UID]['candidates'],
                    key=attrgetter('hitEnd', 'hitStart'),
                    reverse=True,
                )
    # hitIndex[model][idx].keys() == [rec,candidates,partner]
    return hitIndex


def isfirstUnpaired(
    ref: Optional[int] = None,
    mate: Optional[int] = None,
    model: Optional[str] = None,
    index: Optional[Dict[str, Dict[int, Dict[str, Any]]]] = None,
) -> Tuple[Optional[Set[int]], Dict[str, Dict[int, Dict[str, Any]]], Optional[int]]:
    """
    Check for reciprocal best unpaired partner relationship between two hits.

    Parameters
    ----------
    ref : int
        Index (DataFrame row number) of reference hit to check.
    mate : int
        Index of potential partner hit.
    model : str
        Name of HMM model for these hits.
    index : dict
        Hit index dictionary: index[model][idx] = {rec, partner, candidates}.

    Returns
    -------
    found : set or None
        If reciprocal match found, returns set {ref, mate}. Otherwise None.
    index : dict
        Updated index with partner assignments if match found.
    mateFUP : int or None
        If no reciprocal match, returns the index of mate's first unpaired
        candidate (for second-degree reciprocity checking).

    Notes
    -----
    Searches mate's candidate list for ref. If ref is mate's first unpaired
    candidate, they form a reciprocal pair and are marked as partners.
    """
    assert ref is not None, 'ref cannot be None'
    assert mate is not None, 'mate cannot be None'
    assert model is not None, 'model cannot be None'
    assert index is not None, 'index cannot be None'
    # Init result trackers
    found = None
    mateFUP = None

    # Scan candidate partners of 'mate' looking for ref
    for matePartner in index[model][mate]['candidates']:
        # Get the model that this candidate belongs to
        candidate_model = matePartner.model
        candidate_idx = matePartner.idx

        # Check if this candidate is unpaired and matches our ref
        if (
            candidate_model in index
            and candidate_idx in index[candidate_model]
            and index[candidate_model][candidate_idx]['partner'] is None
            and candidate_idx == ref
        ):
            found = {candidate_idx, mate}
            index[model][ref]['partner'] = mate
            index[candidate_model][mate]['partner'] = ref
            return found, index, mateFUP

        # If first unpaired candidate partner for mate is not ref
        elif (
            candidate_model in index
            and candidate_idx in index[candidate_model]
            and index[candidate_model][candidate_idx]['partner'] is None
        ):
            mateFUP = candidate_idx
            return found, index, mateFUP
        else:
            continue

    # If mate candidates include no unpaired reps, return unchanged index
    return found, index, mateFUP


def getPairs(
    hitIndex: Optional[Dict[str, Dict[int, Dict[str, Any]]]] = None,
    paired: Optional[Dict[str, List[Set[int]]]] = None,
) -> Tuple[Dict[str, Dict[int, Dict[str, Any]]], Dict[str, List[Set[int]]]]:
    """
    Identify reciprocal pairs using two-degree candidate matching.

    Parameters
    ----------
    hitIndex : dict
        Hit index dictionary: hitIndex[model][idx] = {rec, partner, candidates}.
    paired : dict, optional
        Existing dictionary of paired hits: paired[model] = [list of pair sets].
        If None, creates new empty dictionary.

    Returns
    -------
    hitIndex : dict
        Updated index with new partner assignments.
    paired : dict
        Updated dictionary with newly identified pairs added to each model's list.

    Notes
    -----
    Implements multi-degree reciprocity checking:
    1. First-degree: Check if hit A and candidate B are each other's best match.
    2. Second-degree: If not, check if B's best match C is reciprocal with B.
    This allows pairing in cases where direct reciprocity is blocked by
    competing candidates.
    """
    assert hitIndex is not None, 'hitIndex cannot be None'
    # If pair tracker not given
    if not paired:
        # Create dict of empty lists, keyed by model name
        paired_dict: Dict[str, List[Set[int]]] = {}
        for model in hitIndex.keys():
            paired_dict[model] = []
    else:
        paired_dict = paired
    # For each HMM model
    for model in hitIndex.keys():
        # Ask each hit in genome
        for refID in hitIndex[model].keys():
            # If it has been asigned a partner
            if hitIndex[model][refID]['partner'] is None:
                # If not partnered, start checking candidate partners
                for Can1 in hitIndex[model][refID]['candidates']:
                    # For a candidate that is also unpartnered
                    if hitIndex[model][Can1.idx]['partner'] is None:
                        # Check if unpartnered candidate is a reciprocal
                        # match for our hit
                        found, hitIndex, mateFUP = isfirstUnpaired(
                            ref=refID, mate=Can1.idx, model=model, index=hitIndex
                        )
                        if found:
                            # If current hit is also the best return match of
                            # our candidate, store as pair.
                            paired_dict[model].append(found)
                        elif mateFUP:
                            # Else if not a return match, check candidate's
                            # first outbound match for reciprocity.
                            found, hitIndex, mateFUP = isfirstUnpaired(
                                ref=Can1.idx, mate=mateFUP, model=model, index=hitIndex
                            )
                            if found:
                                # Store if found.
                                paired_dict[model].append(found)
    return hitIndex, paired_dict


def countUnpaired(hitIndex: Dict[str, Dict[int, Dict[str, Any]]]) -> int:
    """
    Count the total number of unpaired hits across all models.

    Parameters
    ----------
    hitIndex : dict
        Hit index dictionary: hitIndex[model][idx] = {rec, partner, candidates}.

    Returns
    -------
    int
        Total count of hits without assigned partners.
    """
    count = 0
    for model in hitIndex.keys():
        for hitID in hitIndex[model].keys():
            if hitIndex[model][hitID]['partner'] is None:
                count += 1
    return count


def listunpaired(hitIndex: Dict[str, Dict[int, Dict[str, Any]]]) -> List[int]:
    """
    Collect indices of all unpaired hits across all models.

    Parameters
    ----------
    hitIndex : dict
        Hit index dictionary: hitIndex[model][idx] = {rec, partner, candidates}.

    Returns
    -------
    list of int
        List of DataFrame indices for all hits without assigned partners.
    """
    unpaired = []
    for model in hitIndex.keys():
        for hitID in hitIndex[model].keys():
            if hitIndex[model][hitID]['partner'] is None:
                unpaired.append(hitID)
    return unpaired


def iterateGetPairs(
    hitIndex: Dict[str, Dict[int, Dict[str, Any]]], stableReps: int = 0
) -> Tuple[Dict[str, Dict[int, Dict[str, Any]]], Dict[str, List[Set[int]]], List[int]]:
    """
    Repeatedly apply pairing algorithm until convergence or iteration limit.

    Parameters
    ----------
    hitIndex : dict
        Hit index dictionary: hitIndex[model][idx] = {rec, partner, candidates}.
    stableReps : int, default 0
        Maximum number of iterations to continue after no new pairs are found.
        If 0, stops immediately when no new pairs are found.

    Returns
    -------
    hitIndex : dict
        Updated index with final partner assignments.
    paired : dict
        Dictionary of all identified pairs: paired[model] = [list of pair sets].
    unpaired : list of int
        List of DataFrame indices for hits that remain unpaired.

    Notes
    -----
    Iterates pairing until either all hits are paired or the unpaired count
    remains stable for 'stableReps' consecutive iterations. This allows
    pairing to progress through complex candidate competition scenarios.
    """
    # Init stable repeat counter
    reps = 0
    # Run initial pairing
    hitIndex, paired = getPairs(hitIndex=hitIndex)
    # Count remaining unpaired hits
    countUP = countUnpaired(hitIndex)
    # Iterate pairing procedure until either no unpaired remain
    # OR max number of interations without new pairing is reached
    while countUP > 0 and reps < stableReps:
        # Re-run pairing procedure
        hitIndex, paired = getPairs(hitIndex=hitIndex, paired=paired)
        # Store previous unpaired hit count
        lastCountUP = countUP
        # Update unpaired hit count
        countUP = countUnpaired(hitIndex)
        # If no change in upaired hit count, iterate stable rep counter
        if lastCountUP == countUP:
            reps += 1
    # Get IDs of remaining unpaired hits
    unpaired = listunpaired(hitIndex)
    # Return results
    return hitIndex, paired, unpaired


def fetch_padded_hit(
    source: SequenceSource,
    seqid: str,
    start: int,
    end: int,
    strand: str = '+',
    padlen: Optional[int] = None,
    pad: bool = True,
) -> Optional[str]:
    """
    Fetch a hit region, optionally with flanking sequence marked in lowercase.

    Parameters
    ----------
    source : FastaSource or BlastDBSource
        Sequence source to read from.
    seqid : str
        Sequence identifier.
    start, end : int
        1-based inclusive hit coordinates on the plus strand.
    strand : str, default '+'
        Strand of the hit. The result is reverse-complemented when '-'.
    padlen : int, optional
        Number of flanking bases to include either side of the hit. Flanking
        bases are lowercased; the hit itself stays uppercase.
    pad : bool, default True
        When the padded window runs past a contig boundary, pad it with N so
        the record is always ``(end - start + 1) + 2 * padlen`` bases. Padding
        is lowercased along with the rest of the flank, keeping the "not part
        of the hit" convention intact.

    Returns
    -------
    str or None
        Extracted sequence, or None if extraction failed.

    Notes
    -----
    Case marking is applied to the plus-strand sequence and the reverse
    complement is taken afterwards, so the lowercase markers stay attached to
    the flanks they came from. Marking after reverse-complementing would swap
    the 5' and 3' markers.

    Only the pad region can extend past a contig; the hit itself came from the
    contig and is always in bounds.
    """
    if not padlen:
        seq = fetch_sequence(source, seqid, start, end)
        if seq is None:
            return None
        return seq if strand != '-' else str(Seq.Seq(seq).reverse_complement())

    # An inverted window has no hit region to mark. Padding would still produce
    # a valid-looking span, so reject it here rather than emit a nonsense
    # arrangement of lowercase flanks around nothing.
    if start > end:
        logging.debug(f'Inverted hit window for {seqid}: {start}-{end}, skipping')
        return None

    pad_start = start - padlen
    pad_end = end + padlen

    if pad:
        region = fetch_region_padded(source, seqid, pad_start, pad_end)
        if region is None:
            return None
        full_seq = region.seq
    else:
        # Clamp the padded window so the case-marking offsets below are
        # computed against the region actually returned.
        clamped = clamp_region(source, seqid, pad_start, pad_end)
        if clamped is None:
            return None
        pad_start, pad_end = clamped
        maybe_seq = fetch_sequence(source, seqid, pad_start, pad_end)
        if maybe_seq is None:
            return None
        full_seq = maybe_seq

    # Offsets of the hit within the padded window, as 0-based half-open slice
    # bounds. Both coordinate systems here are 1-based inclusive.
    lead = max(0, start - pad_start)
    tail = lead + (min(end, pad_end) - max(start, pad_start) + 1)

    marked = full_seq[:lead].lower() + full_seq[lead:tail] + full_seq[tail:].lower()

    if strand == '-':
        return str(Seq.Seq(marked).reverse_complement())
    return marked


def _model_extension(
    row: Any,
    model: str,
    model_lengths: Optional[Dict[str, int]],
) -> Optional[Tuple[int, int]]:
    """
    Compute bases to add either side of a hit to span the full model.

    Parameters
    ----------
    row : pandas.Series
        Hit row carrying hmmStart, hmmEnd and strand.
    model : str
        Model name, used to look up the length and for the warning.
    model_lengths : dict or None
        Mapping of model name to model length.

    Returns
    -------
    tuple of (int, int) or None
        ``(lower, upper)`` bases to extend at the lower and higher genomic
        coordinate ends, or None if the model length is unavailable.

    Notes
    -----
    The deficits are the same quantities used for flank projection: how many
    model positions the alignment failed to cover at each end. On the plus
    strand the ``hmmStart`` deficit sits at the lower-coordinate end; on the
    minus strand the mapping is reversed.
    """
    model_len = model_lengths.get(model) if model_lengths else None
    if model_len is None:
        logging.warning(
            f'Model length not found for {model}; cannot extend hits to model '
            'length, extracting the aligned region only'
        )
        return None

    try:
        hmm_start = int(row['hmmStart'])
        hmm_end = int(row['hmmEnd'])
    except (ValueError, TypeError, KeyError):
        logging.debug(f'HMM coordinates unavailable for a {model} hit, not extending')
        return None

    start_deficit = _model_deficit(
        hmm_start - 1, 'hmmStart', hmm_start, hmm_end, model_len
    )
    end_deficit = _model_deficit(
        model_len - hmm_end, 'hmmEnd', hmm_start, hmm_end, model_len
    )

    if row['strand'] == '-':
        # hmmStart aligns to the higher genomic coordinate on the minus strand.
        return end_deficit, start_deficit
    return start_deficit, end_deficit


def extractTIRs(
    model: Optional[str] = None,
    hitTable: Optional[pd.DataFrame] = None,
    maxeval: float = 0.001,
    source: Optional[SequenceSource] = None,
    padlen: Optional[int] = None,
    genome_descriptions: Optional[Dict[str, str]] = None,
    model_lengths: Optional[Dict[str, int]] = None,
    extend_to_model: bool = False,
    pad: bool = True,
) -> Tuple[List[SeqRecord], int]:
    """
    Extract TIR sequences for hits of a specific model.

    Parameters
    ----------
    model : str
        Name of HMM model to extract hits for.
    hitTable : pandas.DataFrame
        DataFrame containing all hits with columns: model, target, hitStart, hitEnd,
        strand, evalue, hmmStart, hmmEnd.
    maxeval : float, default 0.001
        Maximum e-value threshold for extracting hits.
    source : FastaSource or BlastDBSource
        Sequence source, from :func:`tirmite.utils.extract.make_source`.
    padlen : int, optional
        Number of flanking bases to extract on each side of hit (shown in lowercase).
    genome_descriptions : dict, optional
        Dictionary mapping sequence IDs to their descriptions.
    model_lengths : dict, optional
        Mapping of model name to model length. Required when
        ``extend_to_model`` is set.
    extend_to_model : bool, default False
        Extend each hit outward by the number of model positions it failed to
        cover, so partial hits are emitted at full model length.
    pad : bool, default True
        Pad with N where the extracted window runs past a contig boundary.

    Returns
    -------
    seqList : list of Bio.SeqRecord.SeqRecord
        List of SeqRecord objects containing extracted TIR sequences.
    hitcount : int
        Number of hits extracted that passed e-value threshold.

    Notes
    -----
    Reverse complement is applied for hits on the negative strand, and such
    records get an '_rc' ID suffix. Padded flanking sequence is shown in
    lowercase letters. Coordinates are 1-based inclusive throughout; see
    :mod:`tirmite.utils.extract` for the extraction contract.

    With ``extend_to_model``, the extension is added on both sides using the
    same model-coverage deficits that :func:`compute_flank_coordinates` uses for
    the external edge and :func:`compute_inner_tsd_coordinates` for the inner
    edge. The extended bases were not matched by the model, so they are a
    projection of where the terminus should end, not evidence that it does.
    """
    assert model is not None, 'model cannot be None'
    assert hitTable is not None, 'hitTable cannot be None'
    assert source is not None, 'source cannot be None'
    hitcount = 0
    seqList = []

    if genome_descriptions is None:
        genome_descriptions = getattr(source, 'descriptions', None)

    eligible_hits = hitTable[
        (hitTable['model'] == model) & (hitTable['evalue'].astype(float) <= maxeval)
    ]
    total_eligible = len(eligible_hits)
    logging.info(
        f'Extracting {total_eligible} hits for model "{model}" '
        f'from {source.describe()}...'
    )
    _log_step = max(1, min(100, total_eligible // 10)) if total_eligible > 0 else 1

    for index, row in eligible_hits.iterrows():
        hitcount += 1
        if hitcount % _log_step == 0 or hitcount == total_eligible:
            logging.info(
                f'  Extracted {hitcount}/{total_eligible} hits for model "{model}"'
            )

        hit_start = int(row['hitStart'])
        hit_end = int(row['hitEnd'])
        extension = ''

        if extend_to_model:
            span = _model_extension(row, model, model_lengths)
            if span is not None:
                lower_ext, upper_ext = span
                hit_start -= lower_ext
                hit_end += upper_ext
                if lower_ext or upper_ext:
                    extension = f' extended:{lower_ext},{upper_ext}'

        hit_seq_str = fetch_padded_hit(
            source,
            row['target'],
            hit_start,
            hit_end,
            row['strand'],
            padlen,
            pad=pad,
        )

        if hit_seq_str is None:
            logging.warning(f'Failed to extract sequence for {model}_{index}, skipping')
            continue

        # Create SeqRecord
        hitrecord = SeqRecord(Seq.Seq(hit_seq_str))
        hitrecord.id = model + '_' + str(index)

        if row['strand'] == '-':
            # Sequence is already reverse-complemented by fetch_padded_hit;
            # only the ID suffix is applied here.
            hitrecord.id = hitrecord.id + '_rc'

        hitrecord.name = hitrecord.id

        # Build description with genome description. Report the extracted span
        # so the header always describes the sequence actually emitted.
        coord_info = '_'.join(
            [
                '[' + str(row['target']) + ':' + str(row['strand']),
                str(hit_start),
                str(hit_end) + ' modelAlignment:' + str(row['hmmStart']),
                str(row['hmmEnd']) + ' E-value:' + str(row['evalue']) + extension + ']',
            ]
        )

        # Add genome description if available
        hitrecord.description = annotate(
            source, row['target'], coord_info, genome_descriptions
        )

        seqList.append(hitrecord)

    return seqList, hitcount


# Fix: Do not load fasta into genome!
def writeTIRs(
    outDir: Optional[str] = None,
    hitTable: Optional[pd.DataFrame] = None,
    maxeval: float = 0.001,
    genome: Any = None,
    prefix: Optional[str] = None,
    padlen: Optional[int] = None,
    genome_descriptions: Optional[Dict[str, str]] = None,
    blastdb: Optional[str] = None,
    model_lengths: Optional[Dict[str, int]] = None,
    extend_hits_to_model: bool = False,
    pad: bool = True,
) -> None:
    """
    Write extracted TIR sequences to FASTA files organized by model.

    Parameters
    ----------
    outDir : str, optional
        Output directory for FASTA files. If None, uses current directory.
    hitTable : pandas.DataFrame
        DataFrame containing all hits with columns: model, target, hitStart, hitEnd,
        strand, evalue, hmmStart, hmmEnd.
    maxeval : float, default 0.001
        Maximum e-value threshold for extracting hits.
    genome : pyfaidx.Fasta, optional
        Indexed genome object for sequence extraction. Required if blastdb is None.
    prefix : str, optional
        Prefix to add to output filenames and sequence IDs.
    padlen : int, optional
        Number of flanking bases to extract on each side of hit (shown in lowercase).
    genome_descriptions : dict, optional
        Dictionary mapping sequence IDs to their descriptions.
    blastdb : str, optional
        Path to BLAST database for sequence extraction. Alternative to genome.
    model_lengths : dict, optional
        Mapping of model name to model length; required for
        ``extend_hits_to_model``.
    extend_hits_to_model : bool, default False
        Emit partial hits at full model length by extending them outward by the
        uncovered model positions.
    pad : bool, default True
        Pad with N where an extracted window runs past a contig boundary.

    Returns
    -------
    None
        Writes FASTA files to disk but returns nothing.

    Notes
    -----
    Creates one FASTA file per model with filename format:
    {prefix}{model}_hits_{count}.fasta

    Either genome or blastdb must be provided for sequence extraction.

    Set ``extend_hits_to_model`` to emit partial hits at full model length;
    ``model_lengths`` is required for that and the option is ignored with a
    warning without it.
    """
    assert hitTable is not None, 'hitTable cannot be None'

    if prefix:
        prefix = cleanID(prefix) + '_'
    else:
        prefix = ''
    if outDir:
        outDir = os.path.abspath(outDir)
        if not os.path.isdir(outDir):
            os.makedirs(outDir)
    else:
        outDir = os.getcwd()

    source = make_source(genome=genome, blastdb=blastdb)

    for model in hitTable['model'].unique():
        # List of TIR seqrecords, and count of hits
        seqList, hitcount = extractTIRs(
            model=model,
            hitTable=hitTable,
            maxeval=maxeval,
            source=source,
            padlen=padlen,
            genome_descriptions=genome_descriptions,
            model_lengths=model_lengths,
            extend_to_model=extend_hits_to_model,
            pad=pad,
        )
        outfile = os.path.join(
            outDir, prefix + model + '_hits_' + str(hitcount) + '.fasta'
        )
        # Write extracted hits to model outfile
        with open(outfile, 'w') as handle:
            for seq in seqList:
                seq.id = prefix + str(seq.id)
                SeqIO.write(seq, handle, 'fasta')


# CS10_Chromosome_02_+_88294_88353_modelAlignment:1_60


def flipTIRs(x: Any, y: Any) -> Tuple[Any, Any]:
    """
    Order two hits by genomic position to determine left and right TIRs.

    Parameters
    ----------
    x : namedtuple
        First hit record with hitStart and hitEnd attributes.
    y : namedtuple
        Second hit record with hitStart and hitEnd attributes.

    Returns
    -------
    tuple
        (left_hit, right_hit) ordered by genomic coordinates (hitStart, then hitEnd).
    """
    left2right = sorted([x, y], key=attrgetter('hitStart', 'hitEnd'))
    return (left2right[0], left2right[1])


def fetchElements(
    paired: Optional[Dict[str, List[Set[int]]]] = None,
    hitIndex: Optional[Dict[str, Dict[int, Dict[str, Any]]]] = None,
    genome: Any = None,
    genome_descriptions: Optional[Dict[str, str]] = None,
    blastdb: Optional[str] = None,
) -> Dict[str, List[Any]]:
    """
    Extract full-length transposon element sequences from paired TIR hits.

    Parameters
    ----------
    paired : dict
        Dictionary of paired hits: paired[model] = [list of pair sets {id1, id2}].
    hitIndex : dict
        Hit index dictionary: hitIndex[model][idx] = {rec, partner, candidates}.
        Supports both nested (model-keyed) and flat structures.
    genome : pyfaidx.Fasta, optional
        Indexed genome object for sequence extraction. Required if blastdb is None.
    genome_descriptions : dict, optional
        Dictionary mapping sequence IDs to their descriptions.
    blastdb : str, optional
        Path to BLAST database for sequence extraction. Alternative to genome.

    Returns
    -------
    dict
        Dictionary of element records keyed by model: TIRelements[model] = [element_records].
        Each element record is a namedtuple with fields: model, chromosome, start, end,
        strand, type, id, leftHit, rightHit, seq, evalue.

    Notes
    -----
    Extracts sequence from leftHit.hitStart to rightHit.hitEnd.
    Handles both symmetric (same model) and asymmetric (different models) pairing.
    Element IDs have format: Element_{counter}.

    Either genome or blastdb must be provided for sequence extraction.
    """
    assert paired is not None, 'paired cannot be None'
    assert hitIndex is not None, 'hitIndex cannot be None'
    source = make_source(
        genome=genome, blastdb=blastdb, descriptions=genome_descriptions
    )
    # Check if hitIndex is nested or flat
    is_nested = isinstance(next(iter(hitIndex.values())), dict)

    def get_hit_record(hit_id: int) -> Any:
        """
        Retrieve hit record from either nested or flat hitIndex structure.

        Parameters
        ----------
        hit_id : int
            Index of the hit record to retrieve.

        Returns
        -------
        namedtuple
            Hit record namedtuple with fields: model, target, hitStart, hitEnd,
            strand, idx, evalue.

        Raises
        ------
        KeyError
            If hit_id is not found in any model within the hitIndex.

        Notes
        -----
        Handles both nested (model-keyed) and flat hitIndex structures for
        backward compatibility.
        """
        if is_nested:
            # Search through all models to find the hit
            for _model_name, model_hits in hitIndex.items():
                if hit_id in model_hits:
                    return model_hits[hit_id]['rec']
            raise KeyError(f'Hit ID {hit_id} not found in any model')
        else:
            # Flat structure - direct access
            return hitIndex[hit_id]['rec']  # type: ignore[index]

    TIRelements: Dict[str, List[Any]] = {}
    gffTup = namedtuple(
        'gffTup',  # type: ignore[name-match]
        [
            'model',
            'chromosome',
            'start',
            'end',
            'strand',
            'type',
            'id',
            'leftHit',
            'rightHit',
            'seq',
            'evalue',
        ],
    )

    # Only process models that actually have pairs
    for model in paired.keys():
        if len(paired[model]) > 0:
            TIRelements[model] = []
            model_counter = 0
            total_pairs = len(paired[model])
            logging.info(
                f'Extracting sequences for {total_pairs} paired elements (model "{model}")...'
            )
            _log_step = max(1, min(100, total_pairs // 10)) if total_pairs > 0 else 1

            for pair in paired[model]:
                model_counter += 1
                if model_counter % _log_step == 0 or model_counter == total_pairs:
                    logging.info(
                        f'  Extracted {model_counter}/{total_pairs} elements for model "{model}"'
                    )
                # Convert set to list for indexing
                hit_ids = list(pair)
                x_id, y_id = hit_ids[0], hit_ids[1]

                # Get hit records using helper function
                x = get_hit_record(x_id)
                y = get_hit_record(y_id)

                leftHit, rightHit = flipTIRs(x, y)

                # Create element ID with counter only (avoids model name duplication)
                eleID = f'Element_{model_counter}'

                # Extract element sequence. Always on the + strand so that the
                # sequence matches the genomic coordinates reported in the GFF,
                # consistent with the flank and target-site extractors.
                ele_seq_str = fetch_sequence(
                    source,
                    leftHit.target,
                    int(leftHit.hitStart),
                    int(rightHit.hitEnd),
                )
                if ele_seq_str is None:
                    logging.warning(f'Failed to extract element {eleID}, skipping')
                    continue

                eleSeq = SeqRecord(Seq.Seq(ele_seq_str))
                eleSeq.id = eleID
                eleSeq.name = eleID

                # Build description with genome description
                coord_info = (
                    '_'.join(
                        [
                            '[' + leftHit.target + ':' + str(leftHit.hitStart),
                            str(rightHit.hitEnd),
                        ]
                    )
                    + ' len='
                    + str(rightHit.hitEnd - leftHit.hitStart)
                    + ']'
                )

                # Add genome description if available
                eleSeq.description = annotate(
                    source, leftHit.target, coord_info, genome_descriptions
                )

                TIRelement = gffTup(
                    model,  # This is the pairing model (left model for asymmetric)
                    leftHit.target,
                    leftHit.hitStart,
                    rightHit.hitEnd,
                    leftHit.strand,
                    'Element',
                    eleID,
                    leftHit,
                    rightHit,
                    eleSeq,
                    'NA',
                )
                TIRelements[model].append(TIRelement)

    return TIRelements


def writeElements(
    outDir: str,
    eleDict: Optional[Dict[str, List[Any]]] = None,
    prefix: Optional[str] = None,
) -> None:
    """
    Write extracted element sequences to FASTA files organized by model.

    Parameters
    ----------
    outDir : str
        Output directory for element FASTA files.
    eleDict : dict, optional
        Dictionary of element records keyed by model: eleDict[model] = [element_records].
    prefix : str, optional
        Prefix to add to output filenames and sequence IDs.

    Returns
    -------
    None
        Writes FASTA files to disk but returns nothing.

    Notes
    -----
    Only creates files for models that have at least one element.
    Output filename format: {prefix}{model}_elements_{count}.fasta
    """
    assert eleDict is not None, 'eleDict cannot be None'
    if prefix:
        prefix = cleanID(prefix) + '_'
    else:
        prefix = ''

    for model in eleDict.keys():
        if len(eleDict[model]) > 0:  # Only write files for models with actual elements
            count = len(eleDict[model])
            outfile = os.path.join(
                outDir,
                prefix + model + '_elements_' + str(count) + '.fasta',
            )
            with open(outfile, 'w') as handle:
                for element in eleDict[model]:
                    element.seq.id = prefix + str(element.seq.id)
                    SeqIO.write(element.seq, handle, 'fasta')


def writePairedTIRs(
    outDir: Optional[str] = None,
    paired: Optional[Dict[str, List[Set[int]]]] = None,
    hitIndex: Optional[Dict[str, Dict[int, Dict[str, Any]]]] = None,
    genome: Any = None,
    prefix: Optional[str] = None,
    padlen: Optional[int] = None,
    genome_descriptions: Optional[Dict[str, str]] = None,
    blastdb: Optional[str] = None,
    config: Any = None,
) -> None:
    """
    Extract and write left and right TIR sequences from paired hits to FASTA.

    Parameters
    ----------
    outDir : str, optional
        Output directory for TIR FASTA files.
    paired : dict
        Dictionary of paired hits: paired[model] = [list of pair sets {id1, id2}].
    hitIndex : dict
        Hit index dictionary: hitIndex[model][idx] = {rec, partner, candidates}.
        Supports both nested (model-keyed) and flat structures.
    genome : pyfaidx.Fasta, optional
        Indexed genome object for sequence extraction. Required if blastdb is None.
    prefix : str, optional
        Prefix to add to output filenames and sequence IDs.
    padlen : int, optional
        Number of flanking bases to extract on each side (shown in lowercase).
    genome_descriptions : dict, optional
        Dictionary mapping sequence IDs to their descriptions.
    blastdb : str, optional
        Path to BLAST database for sequence extraction. Alternative to genome.

    config : PairingConfig, optional
        Pairing configuration. When it describes an asymmetric element, the
        _L/_R labels and output files follow the terminus role rather than
        genomic order, so each model's termini stay together.

    Returns
    -------
    None
        Writes FASTA files to disk but returns nothing.

    Notes
    -----
    Right TIRs are reverse complemented. Outputs are split into two files per model:
    {model}_{counter}_L (left TIR) → {prefix}{model}_paired_left_term_hits_{count}.fasta
    {model}_{counter}_R (right TIR) → {prefix}{model}_paired_right_term_hits_{count}.fasta

    Either genome or blastdb must be provided for sequence extraction.
    """
    assert outDir is not None, 'outDir cannot be None'
    assert paired is not None, 'paired cannot be None'
    assert hitIndex is not None, 'hitIndex cannot be None'
    source = make_source(
        genome=genome, blastdb=blastdb, descriptions=genome_descriptions
    )
    if prefix:
        prefix = cleanID(prefix) + '_'
    else:
        prefix = ''

    # Check if hitIndex is nested (new format) or flat (old format)
    is_nested = isinstance(next(iter(hitIndex.values())), dict)

    def get_hit_record(hit_id: int) -> Any:
        """
        Retrieve hit record from either nested or flat hitIndex structure.

        Parameters
        ----------
        hit_id : int
            Index of the hit record to retrieve.

        Returns
        -------
        namedtuple
            Hit record namedtuple with fields: model, target, hitStart, hitEnd,
            strand, idx, evalue.

        Raises
        ------
        KeyError
            If hit_id is not found in any model within the hitIndex.

        Notes
        -----
        Handles both nested (model-keyed) and flat hitIndex structures for
        backward compatibility.
        """
        if is_nested:
            # Search through all models to find the hit
            for _model_name, model_hits in hitIndex.items():
                if hit_id in model_hits:
                    return model_hits[hit_id]['rec']
            raise KeyError(f'Hit ID {hit_id} not found in any model')
        else:
            # Flat structure - direct access
            return hitIndex[hit_id]['rec']  # type: ignore[index]

    # Only process models that actually have pairs
    for model in paired.keys():
        if len(paired[model]) > 0:  # Only write files for models with actual pairs
            model_counter = 0
            left_seqList: List[Any] = []  # Left terminus hit sequences
            right_seqList: List[Any] = []  # Right terminus hit sequences
            total_pairs = len(paired[model])
            logging.info(
                f'Extracting TIR sequences for {total_pairs} pairs (model "{model}")...'
            )
            _log_step = max(1, min(100, total_pairs // 10)) if total_pairs > 0 else 1

            for pair in paired[model]:
                model_counter += 1
                if model_counter % _log_step == 0 or model_counter == total_pairs:
                    logging.info(
                        f'  Extracted {model_counter}/{total_pairs} TIR pairs for model "{model}"'
                    )
                # Convert set to list for indexing
                hit_ids = list(pair)
                x_id, y_id = hit_ids[0], hit_ids[1]

                # Get hit records using helper function
                x = get_hit_record(x_id)
                y = get_hit_record(y_id)

                leftHit, rightHit = flipTIRs(x, y)
                eleID = f'{model}_{model_counter}'

                # Extract TIR sequences. Both are taken on the + strand; the
                # right TIR is reverse complemented below so that the pair reads
                # 5'->3' inward.
                left_seq_str = fetch_padded_hit(
                    source,
                    leftHit.target,
                    int(leftHit.hitStart),
                    int(leftHit.hitEnd),
                    '+',
                    padlen,
                )
                right_seq_str = fetch_padded_hit(
                    source,
                    rightHit.target,
                    int(rightHit.hitStart),
                    int(rightHit.hitEnd),
                    '+',
                    padlen,
                )

                if left_seq_str is None or right_seq_str is None:
                    logging.warning(f'Failed to extract TIRs for {eleID}, skipping')
                    continue

                # Create SeqRecords for FASTA output only
                eleSeqLeft = SeqRecord(Seq.Seq(left_seq_str))
                eleSeqRight = SeqRecord(Seq.Seq(right_seq_str))
                eleSeqRight = eleSeqRight.reverse_complement()

                # Sequence orientation is positional (the higher-coordinate
                # terminus is reverse complemented so both read inward), but
                # the _L/_R label and the output file follow terminus role, so
                # each model's termini stay together for downstream alignment.
                left_role, right_role = _pair_roles(leftHit, rightHit, config)
                left_label = 'L' if left_role == 'left' else 'R'
                right_label = 'L' if right_role == 'left' else 'R'

                eleSeqLeft.id = eleID + '_' + left_label
                eleSeqLeft.name = eleSeqLeft.id

                # Build left description with genome description
                left_coord = (
                    '_'.join(
                        [
                            '[' + leftHit.target + ':' + str(leftHit.hitStart),
                            str(leftHit.hitEnd),
                        ]
                    )
                    + ']'
                )

                eleSeqLeft.description = annotate(
                    source, leftHit.target, left_coord, genome_descriptions
                )

                eleSeqRight.id = eleID + '_' + right_label
                eleSeqRight.name = eleSeqRight.id

                # Build right description with genome description
                right_coord = (
                    '_'.join(
                        [
                            '[' + leftHit.target + ':' + str(rightHit.hitEnd),
                            str(rightHit.hitStart),
                        ]
                    )
                    + ']'
                )

                eleSeqRight.description = annotate(
                    source, leftHit.target, right_coord, genome_descriptions
                )

                # File routing follows the same role labels.
                for rec, label in (
                    (eleSeqLeft, left_label),
                    (eleSeqRight, right_label),
                ):
                    if label == 'L':
                        left_seqList.append(rec)
                    else:
                        right_seqList.append(rec)

            # Write separate FASTA files for left and right terminus hits
            if left_seqList:
                left_outfile = os.path.join(
                    outDir,
                    prefix
                    + model
                    + '_paired_left_term_hits_'
                    + str(len(left_seqList))
                    + '.fasta',
                )
                with open(left_outfile, 'w') as handle:
                    for seq in left_seqList:
                        seq.id = prefix + str(seq.id)
                        SeqIO.write(seq, handle, 'fasta')

            if right_seqList:
                right_outfile = os.path.join(
                    outDir,
                    prefix
                    + model
                    + '_paired_right_term_hits_'
                    + str(len(right_seqList))
                    + '.fasta',
                )
                with open(right_outfile, 'w') as handle:
                    for seq in right_seqList:
                        seq.id = prefix + str(seq.id)
                        SeqIO.write(seq, handle, 'fasta')


def _model_deficit(
    raw_offset: int,
    which: str,
    hmm_start: int,
    hmm_end: int,
    model_len: int,
) -> int:
    """
    Clamp a model-coverage deficit to a non-negative value.

    Parameters
    ----------
    raw_offset : int
        Computed number of uncovered model positions, which may be negative if
        the model length is wrong.
    which : str
        Name of the coordinate the deficit was derived from, for the warning.
    hmm_start, hmm_end : int
        Alignment coordinates on the model, for the warning.
    model_len : int
        Declared model length.

    Returns
    -------
    int
        ``raw_offset`` if it is non-negative, otherwise 0.

    Notes
    -----
    A negative deficit means the alignment ran past the declared end of the
    model, which is impossible for a correct model length. It indicates a
    mismatched ``--lengths-file``, an HMM that does not correspond to the hit
    table, or ``--query-len`` applied to a table containing several queries.
    Left unclamped it shifts the flank window *into* the hit, silently
    extracting element sequence and labelling it as flanking. The
    ``flank_max_offset`` guard cannot catch it either, since ``offset > max``
    is trivially false for a negative value.
    """
    if raw_offset >= 0:
        return raw_offset

    logging.warning(
        f'Model length {model_len} is inconsistent with alignment coordinates '
        f'{hmm_start}-{hmm_end} ({which} implies a deficit of {raw_offset}). '
        'Check that the HMM or --lengths-file matches this hit table. '
        'Treating the offset as 0.'
    )
    return 0


def compute_flank_coordinates(
    hit_start: int,
    hit_end: int,
    strand: str,
    is_left_terminus: bool,
    hmm_start: int,
    hmm_end: int,
    model_len: int,
    flank_len: int,
) -> Tuple[int, int, int]:
    """
    Compute genomic coordinates for the external flanking region of a terminus hit.

    The "external" end of a terminus hit is the side that faces away from the TE
    body. For the left terminus this is the side at lower genomic coordinates; for
    the right terminus it is the side at higher genomic coordinates.

    When a hit does not cover position 1 of the model (hmmStart > 1) the external
    boundary must be shifted by the number of uncovered model positions so that the
    reported flank begins at the correct genomic position.

    Parameters
    ----------
    hit_start : int
        1-based start coordinate of the hit in genomic coordinates (always < hit_end).
    hit_end : int
        1-based end coordinate of the hit in genomic coordinates.
    strand : str
        Strand of the hit: '+' or '-'.
    is_left_terminus : bool
        True if the hit represents the left (5') terminus of the element.
    hmm_start : int
        Alignment start position on the HMM model (1-based).
        For + strand hits this aligns to hit_start; for - strand hits it aligns to hit_end.
    hmm_end : int
        Alignment end position on the HMM model (1-based).
        For + strand hits this aligns to hit_end; for - strand hits it aligns to hit_start.
    model_len : int
        Total length of the HMM model in positions.
    flank_len : int
        Number of bases to extract in the flanking region.

    Returns
    -------
    flank_start : int
        1-based start coordinate of the flank region.
    flank_end : int
        1-based end coordinate of the flank region (inclusive).
    offset : int
        Number of model positions between the hit alignment and the external
        end of the model (0 means the alignment reaches the model end).

    Notes
    -----
    Coordinate system:
      - For + strand: hmmStart aligns to hit_start, hmmEnd aligns to hit_end.
      - For - strand: hmmStart aligns to hit_end (higher coord), hmmEnd aligns to hit_start.

    Left terminus external boundary:
      - + strand: external_pos = hit_start - (hmm_start - 1)
      - - strand: external_pos = hit_start - (model_len - hmm_end)

    Right terminus external boundary:
      - + strand: external_pos = hit_end + (model_len - hmm_end)
      - - strand: external_pos = hit_end + (hmm_start - 1)
    """
    if is_left_terminus:
        # External end faces LEFT (lower genomic coordinates)
        if strand == '+':
            offset = _model_deficit(
                hmm_start - 1, 'hmmStart', hmm_start, hmm_end, model_len
            )
            external_pos = hit_start - offset
        else:  # '-'
            # hmmStart aligns to hit_end (higher coord); hmmEnd aligns to hit_start
            offset = _model_deficit(
                model_len - hmm_end, 'hmmEnd', hmm_start, hmm_end, model_len
            )
            external_pos = hit_start - offset
        flank_start = external_pos - flank_len
        flank_end = external_pos - 1
    else:
        # External end faces RIGHT (higher genomic coordinates)
        if strand == '+':
            offset = _model_deficit(
                model_len - hmm_end, 'hmmEnd', hmm_start, hmm_end, model_len
            )
            external_pos = hit_end + offset
        else:  # '-'
            # hmmStart aligns to hit_end; external end is at hit_end side
            offset = _model_deficit(
                hmm_start - 1, 'hmmStart', hmm_start, hmm_end, model_len
            )
            external_pos = hit_end + offset
        flank_start = external_pos + 1
        flank_end = external_pos + flank_len

    return flank_start, flank_end, offset


def compute_inner_tsd_coordinates(
    hit_start: int,
    hit_end: int,
    strand: str,
    is_left_terminus: bool,
    hmm_start: int,
    hmm_end: int,
    model_len: int,
    tsd_length: int,
) -> Tuple[int, int]:
    """
    Compute genomic coordinates of the TSD at the inner boundary of a terminus hit.

    The "inner boundary" of a terminus hit is the side that faces the element body.
    For a left terminus this is the RIGHT (higher coordinate) end; for a right
    terminus it is the LEFT (lower coordinate) end.  When the TSD is modelled as
    part of the terminus HMM it occupies the innermost ``tsd_length`` model
    positions of the terminus model.

    Parameters
    ----------
    hit_start : int
        1-based start coordinate of the hit (always ≤ hit_end).
    hit_end : int
        1-based end coordinate of the hit (always ≥ hit_start).
    strand : str
        Strand of the hit: '+' or '-'.
    is_left_terminus : bool
        True if the hit is the left (5') terminus of the element.
    hmm_start : int
        1-based alignment start on the HMM model.
        For + strand hits this aligns with ``hit_start``; for - strand
        hits it aligns with ``hit_end``.
    hmm_end : int
        1-based alignment end on the HMM model.
        For + strand hits this aligns with ``hit_end``; for - strand
        hits it aligns with ``hit_start``.
    model_len : int
        Total number of positions in the HMM model.
    tsd_length : int
        Number of bases in the TSD feature.

    Returns
    -------
    tsd_start : int
        1-based genomic start coordinate of the TSD (always ≤ tsd_end).
    tsd_end : int
        1-based genomic end coordinate of the TSD.

    Notes
    -----
    Coordinate conventions (same as :func:`compute_flank_coordinates`):

    For + strand hits, model positions increase in the same direction as
    genomic coordinates (hmmStart aligns to hit_start).

    For - strand hits, model positions increase in the *opposite* direction
    (hmmStart aligns to hit_end; hmmEnd aligns to hit_start).

    Left terminus inner boundary
        - ``+`` strand: model position ``model_len`` is at the RIGHT (higher)
          end of the hit.  ``inner_pos = hit_end + (model_len - hmm_end)``.
          TSD occupies ``[inner_pos - tsd_length + 1, inner_pos]``.
        - ``-`` strand: model position 1 is at the RIGHT (higher) end of the
          hit.  ``inner_pos = hit_end + (hmm_start - 1)``.
          TSD occupies ``[inner_pos - tsd_length + 1, inner_pos]``.

    Right terminus inner boundary
        - ``+`` strand: model position 1 is at the LEFT (lower) end of the
          hit.  ``inner_pos = hit_start - (hmm_start - 1)``.
          TSD occupies ``[inner_pos, inner_pos + tsd_length - 1]``.
        - ``-`` strand: model position ``model_len`` is at the LEFT (lower)
          end of the hit.  ``inner_pos = hit_start - (model_len - hmm_end)``.
          TSD occupies ``[inner_pos, inner_pos + tsd_length - 1]``.
    """
    # Deficits are clamped to >= 0; see _model_deficit for why a negative value
    # is possible and what it means.
    end_deficit = _model_deficit(
        model_len - hmm_end, 'hmmEnd', hmm_start, hmm_end, model_len
    )
    start_deficit = _model_deficit(
        hmm_start - 1, 'hmmStart', hmm_start, hmm_end, model_len
    )

    if is_left_terminus:
        # Inner end faces RIGHT (higher genomic coords)
        if strand == '+':
            # Model pos model_len aligns to right of hit
            inner_pos = hit_end + end_deficit
        else:  # '-'
            # Model pos 1 (hmmStart) aligns to hit_end (rightmost genomic coord)
            inner_pos = hit_end + start_deficit
        tsd_start = inner_pos - tsd_length + 1
        tsd_end = inner_pos
    else:
        # Right terminus: inner end faces LEFT (lower genomic coords)
        if strand == '+':
            # Model pos 1 (hmmStart) aligns to hit_start (leftmost genomic coord)
            inner_pos = hit_start - start_deficit
        else:  # '-'
            # Model pos model_len aligns to hit_start (leftmost genomic coord)
            inner_pos = hit_start - end_deficit
        tsd_start = inner_pos
        tsd_end = inner_pos + tsd_length - 1

    return tsd_start, tsd_end


class FlankResult(NamedTuple):
    """
    External flanking sequence for one terminus hit.

    Attributes
    ----------
    seq : str
        Flank sequence on the plus strand, always ``flank_len`` bases when
        ``pad`` was requested.
    start, end : int
        1-based inclusive coordinates of the real (non-padded) sequence.
    offset : int
        Uncovered model positions between the hit alignment and the projected
        external edge of the element.
    left_pad, right_pad : int
        Bases synthesised because the flank ran past a contig boundary.
    """

    seq: str
    start: int
    end: int
    offset: int
    left_pad: int
    right_pad: int

    @property
    def is_padded(self) -> bool:
        """
        Report whether any part of this flank was synthesised.

        Returns
        -------
        bool
            True if either pad count is non-zero.
        """
        return bool(self.left_pad or self.right_pad)


def extract_terminus_flank(
    source: SequenceSource,
    hit: Any,
    is_left: bool,
    model_len: Optional[int],
    hmm_start: Optional[int],
    hmm_end: Optional[int],
    flank_len: int,
    flank_max_offset: Optional[int] = None,
    pad: bool = True,
) -> Optional[FlankResult]:
    """
    Extract the external flanking region of a single terminus hit.

    Parameters
    ----------
    source : FastaSource or BlastDBSource
        Sequence source to read from.
    hit : namedtuple
        Hit record with target, hitStart, hitEnd, strand, idx.
    is_left : bool
        True if the hit is a left terminus; False for a right terminus.
    model_len : int or None
        Length of the terminus model. None skips the hit with a warning.
    hmm_start, hmm_end : int or None
        Alignment coordinates on the model. None skips the hit.
    flank_len : int
        Number of flanking bases to extract.
    flank_max_offset : int, optional
        Reject hits whose model-coverage deficit exceeds this.
    pad : bool, default True
        Pad with N when the flank runs past a contig boundary, so every flank
        is ``flank_len`` bases. When False the flank is truncated instead.

    Returns
    -------
    FlankResult or None
        None when the flank cannot be extracted: missing model length or HMM
        coordinates, offset above ``flank_max_offset``, or a region that does
        not overlap the contig at all.

    Notes
    -----
    Flanks are always taken on the plus strand for genomic orientation; the
    hit's own strand only affects which side of the hit the flank is on, via
    :func:`compute_flank_coordinates`.

    This is the single implementation used by both ``writeFlanks`` and
    ``writeTargetSites``. They previously carried near-identical private copies
    that had already drifted in their logging.
    """
    if model_len is None:
        logging.warning(
            f'Model length not found for {hit.model}, skipping flank for hit {hit.idx}'
        )
        return None

    if hmm_start is None or hmm_end is None:
        logging.warning(
            f'HMM coordinates unavailable for hit {hit.idx}, skipping flank'
        )
        return None

    flank_start, flank_end, offset = compute_flank_coordinates(
        hit_start=int(hit.hitStart),
        hit_end=int(hit.hitEnd),
        strand=hit.strand,
        is_left_terminus=is_left,
        hmm_start=hmm_start,
        hmm_end=hmm_end,
        model_len=model_len,
        flank_len=flank_len,
    )

    if flank_max_offset is not None and offset > flank_max_offset:
        logging.debug(
            f'Skipping flank for hit {hit.idx}: offset {offset} > max {flank_max_offset}'
        )
        return None

    # A flank entirely before the contig start is reported distinctly from one
    # that merely overhangs, since it means the element sits at the very edge of
    # the assembly and no flanking context exists at all.
    if flank_end < 1:
        logging.warning(
            f'Flank for hit {hit.idx} on {hit.target} falls entirely before '
            f'contig start (computed coords {flank_start}–{flank_end}), skipping'
        )
        return None

    if pad:
        region = fetch_region_padded(source, hit.target, flank_start, flank_end)
        if region is None:
            logging.debug(f'Empty flank region for hit {hit.idx}, skipping')
            return None
        return FlankResult(
            seq=region.seq,
            start=region.start,
            end=region.end,
            offset=offset,
            left_pad=region.left_pad,
            right_pad=region.right_pad,
        )

    clamped = clamp_region(source, hit.target, flank_start, flank_end)
    if clamped is None:
        logging.debug(f'Empty flank region for hit {hit.idx}, skipping')
        return None
    flank_start, flank_end = clamped

    seq_str = fetch_sequence(source, hit.target, flank_start, flank_end)
    if seq_str is None:
        logging.warning(f'Failed to extract flank sequence for hit {hit.idx}, skipping')
        return None

    if len(seq_str) < flank_len:
        logging.warning(
            f'Flank for hit {hit.idx} on {hit.target} is truncated at '
            f'contig boundary: expected {flank_len}bp, '
            f'extracted {len(seq_str)}bp (coords {flank_start}–{flank_end})'
        )

    return FlankResult(
        seq=seq_str,
        start=flank_start,
        end=flank_end,
        offset=offset,
        left_pad=0,
        right_pad=0,
    )


def _determine_terminus_type(hit: Any, config: Any) -> Optional[str]:
    """
    Determine whether a hit is a left or right terminus based on pairing config.

    Parameters
    ----------
    hit : namedtuple
        Hit record with model and strand attributes.
    config : PairingConfig
        Configuration specifying orientation and model assignments.

    Returns
    -------
    str or None
        'left' if the hit is a left terminus, 'right' if right terminus,
        or None if the terminus type cannot be determined (e.g. same-strand
        symmetric pairing without paired context).
    """
    if config.is_asymmetric:
        if hit.model == config.left_model:
            return 'left'
        elif hit.model == config.right_model:
            return 'right'
        return None
    else:
        # Symmetric model – distinguish by strand when strands differ
        if config.left_strand != config.right_strand:
            if hit.strand == config.left_strand:
                return 'left'
            elif hit.strand == config.right_strand:
                return 'right'
        # Same-strand symmetric pairing (F,F or R,R) – can't determine without pair
        return None


class TerminusAssignment(NamedTuple):
    """
    Which terminus of an element a hit represents, and which way it faces.

    Attributes
    ----------
    role : str
        'left' or 'right' - which end of the *element* this hit is, in the
        element's own 5'->3' frame. For asymmetric elements this follows the
        model that produced the hit; for symmetric elements it follows strand.
    is_lower : bool
        True when this terminus' external edge faces the lower genomic
        coordinate. This is what :func:`compute_flank_coordinates` means by
        ``is_left_terminus``.

    Notes
    -----
    ``role`` and ``is_lower`` coincide only for forward-inserted elements. When
    an element is inserted in reverse, the model that defines its left terminus
    sits at the higher genomic coordinate, so ``role='left'`` pairs with
    ``is_lower=False``. Conflating the two is what caused flanks to be taken
    from inside reverse-inserted elements.
    """

    role: str
    is_lower: bool


def _pair_roles(left_hit: Any, right_hit: Any, config: Any) -> Tuple[str, str]:
    """
    Assign terminus roles to the two hits of a pair.

    Parameters
    ----------
    left_hit, right_hit : namedtuple
        The pair's hits in genomic order, as returned by :func:`flipTIRs`.
    config : PairingConfig or None
        Pairing configuration. Roles follow model identity when it describes an
        asymmetric element.

    Returns
    -------
    tuple of (str, str)
        ``(role_of_left_hit, role_of_right_hit)``, each 'left' or 'right'.

    Notes
    -----
    For a symmetric element, or when no asymmetric config is available, role
    follows genomic order - there is nothing else to distinguish the ends.

    For an asymmetric element the role follows the model that produced each
    hit, so a reverse-inserted element correctly reports its left-terminus
    model as the left terminus even though that hit lies at the higher
    coordinate. Only labelling and output routing use this; the extracted
    sequence is governed by genomic position.
    """
    if config is None or not config.is_asymmetric:
        return 'left', 'right'

    if left_hit.model == config.left_model and right_hit.model == config.right_model:
        return 'left', 'right'
    if left_hit.model == config.right_model and right_hit.model == config.left_model:
        # Reverse-inserted element: genomic order is the mirror of model order.
        return 'right', 'left'

    # Models do not match the config (e.g. a pair carried over from another
    # pairing-map row). Fall back to genomic order rather than guessing.
    return 'left', 'right'


def resolve_terminus(hit: Any, config: Any) -> Optional[TerminusAssignment]:
    """
    Determine a hit's terminus role and which genomic side its outer edge faces.

    Parameters
    ----------
    hit : namedtuple
        Hit record with at least ``model`` and ``strand``.
    config : PairingConfig
        Pairing configuration supplying the models and expected strands.

    Returns
    -------
    TerminusAssignment or None
        None when the terminus cannot be determined from a single hit: an
        unrecognised model, or a symmetric same-strand orientation (F,F / R,R)
        where one hit alone carries no information about which end it is.

    Notes
    -----
    For asymmetric configurations the role comes from the model name, and the
    insertion direction is inferred by comparing the hit's strand with the
    strand the configuration expects for that role. A match means a forward
    insertion, where role and genomic side agree; a mismatch means the element
    is inserted in reverse and the sides swap.

    For symmetric configurations with differing strands (F,R / R,F) the strand
    alone gives both role and side, and is invariant to insertion direction.
    """
    if config is None:
        return None

    if config.is_asymmetric:
        if hit.model == config.left_model:
            role = 'left'
            expected_strand = config.left_strand
        elif hit.model == config.right_model:
            role = 'right'
            expected_strand = config.right_strand
        else:
            return None

        forward_insertion = hit.strand == expected_strand
        if forward_insertion:
            is_lower = role == 'left'
        else:
            # Element inserted in reverse: the left terminus is now at the
            # higher coordinate and vice versa.
            is_lower = role == 'right'
        return TerminusAssignment(role=role, is_lower=is_lower)

    # Symmetric model - distinguish by strand when the strands differ
    if config.left_strand != config.right_strand:
        if hit.strand == config.left_strand:
            return TerminusAssignment(role='left', is_lower=True)
        elif hit.strand == config.right_strand:
            return TerminusAssignment(role='right', is_lower=False)
        return None

    # Same-strand symmetric pairing (F,F or R,R) - can't determine without pair
    return None


def writeFlanks(
    outDir: Optional[str] = None,
    hitTable: Optional[pd.DataFrame] = None,
    model_lengths: Optional[Dict[str, int]] = None,
    paired: Optional[Dict[str, List[Set[int]]]] = None,
    hitIndex: Optional[Dict[str, Dict[int, Dict[str, Any]]]] = None,
    config: Any = None,
    genome: Any = None,
    prefix: Optional[str] = None,
    flank_len: int = 50,
    flank_max_offset: Optional[int] = None,
    write_all: bool = False,
    write_paired: bool = False,
    genome_descriptions: Optional[Dict[str, str]] = None,
    blastdb: Optional[str] = None,
    pad_flanks: bool = True,
) -> None:
    """
    Extract and write external flanking sequences for terminus hits to FASTA files.

    The "external flank" is the genomic sequence immediately outside each terminus
    hit, i.e. upstream of the left terminus and downstream of the right terminus.
    Flank coordinates are corrected for any gap between the hit alignment and the
    external end of the model (offset correction).

    Parameters
    ----------
    outDir : str, optional
        Output directory for flank FASTA files.
    hitTable : pandas.DataFrame
        DataFrame with columns model, target, hitStart, hitEnd, strand, evalue,
        hmmStart, hmmEnd. Used to look up HMM alignment coordinates by hit index.
    model_lengths : dict
        Dictionary mapping model name to model length in positions.
    paired : dict
        Dictionary of paired hits: paired[model] = [list of pair sets {id1, id2}].
    hitIndex : dict
        Hit index: hitIndex[model][idx] = {rec, partner, candidates}.
    config : PairingConfig
        Configuration with orientation and model assignments.
    genome : pyfaidx.Fasta, optional
        Indexed genome for sequence extraction. Required if blastdb is None.
    prefix : str, optional
        Prefix for output filenames and sequence IDs.
    flank_len : int, default 50
        Number of bases to extract in each flanking region.
    flank_max_offset : int, optional
        Maximum allowed offset (uncovered model positions) between hit alignment
        and model end. Hits with offset > this value are skipped.
    write_all : bool, default False
        If True, write flanks for all hits (paired and unpaired) to output files.
        For symmetric same-strand orientations (F,F or R,R), both left and right
        flanks are written to separate files with a warning.
    write_paired : bool, default False
        If True, write outer flanks for termini assigned to pairs to separate files.
    genome_descriptions : dict, optional
        Dictionary mapping sequence IDs to their descriptions.
    blastdb : str, optional
        Path to BLAST database. Alternative to genome.
    pad_flanks : bool, default True
        Pad flanks with N where they run past a contig boundary, so every
        record is ``flank_len`` bases. Padded records carry a ``padded:L,R``
        field in their description. When False, such flanks are truncated.

    Returns
    -------
    None
        Writes FASTA files to disk.

    Notes
    -----
    Output files are named:
      {prefix}{model}_left_flank_{count}.fasta  – flanks for left terminus hits
      {prefix}{model}_right_flank_{count}.fasta – flanks for right terminus hits
      {prefix}{model}_paired_left_flank_{count}.fasta – paired left flanks
      {prefix}{model}_paired_right_flank_{count}.fasta – paired right flanks

    For asymmetric pairings the left and right model names may differ, so each
    model gets its own pair of files. For symmetric pairings the same model name
    is used for both files (distinguished by the _left_/_right_ suffix).

    For symmetric same-strand orientations (F,F or R,R) when write_all is True,
    both left and right flanks are written for all hits to separate files and a
    warning is raised advising the user to use --flanks-paired instead.

    Flanking sequences are always reported in the forward (+) genomic strand
    orientation. Coordinates are 1-based.
    """
    assert outDir is not None, 'outDir cannot be None'
    assert hitTable is not None, 'hitTable cannot be None'
    assert model_lengths is not None, 'model_lengths cannot be None'
    assert paired is not None, 'paired cannot be None'
    assert hitIndex is not None, 'hitIndex cannot be None'
    source = make_source(
        genome=genome, blastdb=blastdb, descriptions=genome_descriptions
    )

    # When config is None (e.g. pairing-map mode with multiple configs) unpaired
    # hits cannot be attributed to a terminus; fall back to paired-only processing.
    if config is None and write_all:
        logging.debug(
            'config is None: cannot determine terminus type for unpaired hits; '
            'processing paired hits only.'
        )
        write_all = False
        write_paired = True

    if prefix:
        prefix = cleanID(prefix) + '_'
    else:
        prefix = ''

    # left_flanks[model] and right_flanks[model] accumulate SeqRecords
    left_flanks: Dict[str, List[Any]] = {}
    right_flanks: Dict[str, List[Any]] = {}

    # ------------------------------------------------------------------
    # Helper: look up hmmStart / hmmEnd for a hit by its DataFrame index
    # ------------------------------------------------------------------
    def get_hmm_coords(hit_idx: int) -> Tuple[Optional[int], Optional[int]]:
        """
        Retrieve HMM alignment coordinates for a hit.

        Parameters
        ----------
        hit_idx : int
            DataFrame row index for the hit.

        Returns
        -------
        tuple of (int or None, int or None)
            (hmmStart, hmmEnd) as integers, or (None, None) if unavailable.
        """
        if hit_idx not in hitTable.index:  # type: ignore[union-attr]
            return None, None
        row = hitTable.loc[hit_idx]  # type: ignore[union-attr]
        try:
            h_start = int(row['hmmStart'])
            h_end = int(row['hmmEnd'])
        except (ValueError, TypeError):
            return None, None
        return h_start, h_end

    # ------------------------------------------------------------------
    # Helper: retrieve hit record from nested or flat hitIndex
    # ------------------------------------------------------------------
    is_nested = bool(hitIndex) and isinstance(next(iter(hitIndex.values())), dict)

    def get_hit_record(hit_id: int) -> Any:
        """
        Retrieve hit record from hitIndex (handles nested or flat structure).

        Parameters
        ----------
        hit_id : int
            Index of the hit record.

        Returns
        -------
        namedtuple
            Hit record with model, target, hitStart, hitEnd, strand, idx, evalue.
        """
        if is_nested:
            for _m, model_hits in hitIndex.items():
                if hit_id in model_hits:
                    return model_hits[hit_id]['rec']
            raise KeyError(f'Hit ID {hit_id} not found in hitIndex')
        return hitIndex[hit_id]['rec']  # type: ignore[index]

    # ------------------------------------------------------------------
    # Helper: extract a flank SeqRecord for one hit
    # ------------------------------------------------------------------
    def build_flank_record(
        hit: Any,
        is_left: bool,
        record_id: str,
        role: Optional[str] = None,
    ) -> Optional[SeqRecord]:
        """
        Build a SeqRecord for the external flank of a single terminus hit.

        Parameters
        ----------
        hit : namedtuple
            Hit record with model, target, hitStart, hitEnd, strand, idx.
        is_left : bool
            True if the hit is a left terminus; False for right terminus.
        record_id : str
            Identifier to assign to the resulting SeqRecord.

        role : str, optional
            Terminus role, 'left' or 'right', used for the _L/_R suffix.
            Defaults to the genomic side given by ``is_left``, which is the
            same thing except for reverse-inserted asymmetric elements.

        Returns
        -------
        SeqRecord or None
            The flank SeqRecord, or None if the flank cannot be extracted
            (missing model length, missing HMM coords, offset exceeds max,
            or a region that does not overlap the contig).
        """
        hmm_start, hmm_end = get_hmm_coords(hit.idx)
        flank = extract_terminus_flank(
            source=source,
            hit=hit,
            is_left=is_left,
            model_len=model_lengths.get(hit.model) if model_lengths else None,  # type: ignore[union-attr]
            hmm_start=hmm_start,
            hmm_end=hmm_end,
            flank_len=flank_len,
            flank_max_offset=flank_max_offset,
            pad=pad_flanks,
        )
        if flank is None:
            return None

        record = SeqRecord(Seq.Seq(flank.seq))
        # The _L/_R suffix names the terminus ROLE. It defaults to the genomic
        # side, which is the same thing for symmetric elements and for forward
        # insertions; callers pass an explicit role where the two can differ.
        effective_role = role if role is not None else ('left' if is_left else 'right')
        side = 'L' if effective_role == 'left' else 'R'
        record.id = f'{record_id}_{side}'
        record.name = record.id

        coord_info = (
            f'[{hit.target}:+ {flank.start}_{flank.end}'
            f' hit:{hit.strand}:{hit.hitStart}_{hit.hitEnd}'
            f' offset:{flank.offset}'
        )
        if flank.is_padded:
            coord_info += f' padded:{flank.left_pad},{flank.right_pad}'
        coord_info += ']'

        record.description = annotate(
            source, hit.target, coord_info, genome_descriptions
        )

        return record

    # ------------------------------------------------------------------
    # Process paired hits
    # ------------------------------------------------------------------
    paired_hit_ids: Set[int] = set()

    # Separate accumulators for paired-only flanks (with element IDs)
    paired_left_flanks: Dict[str, List[Any]] = {}
    paired_right_flanks: Dict[str, List[Any]] = {}

    def _make_paired_flank_record(
        source_rec: SeqRecord, element_id: str, pair_id: str, suffix: str
    ) -> SeqRecord:
        """
        Create a paired-only flank record with element ID in the header.

        Parameters
        ----------
        source_rec : SeqRecord
            Record to copy the sequence and description from.
        element_id : str
            Element identifier for the new record ID.
        pair_id : str
            Pair identifier for the new record ID.
        suffix : str
            Trailing component of the new record ID, e.g. 'L' or 'R'.

        Returns
        -------
        SeqRecord
            A new record carrying the element ID in its header.
        """
        rec = SeqRecord(source_rec.seq)
        rec.id = f'{element_id}_{pair_id}_{suffix}'
        rec.name = rec.id
        rec.description = source_rec.description
        return rec

    for model in paired.keys():
        model_counter = 0
        for pair in paired[model]:
            model_counter += 1
            hit_ids = list(pair)
            x_id, y_id = hit_ids[0], hit_ids[1]
            x = get_hit_record(x_id)
            y = get_hit_record(y_id)
            leftHit, rightHit = flipTIRs(x, y)
            pair_id = f'{model}_{model_counter}'
            element_id = f'Element_{model_counter}'

            # Geometry is positional: the lower-coordinate hit's outer edge
            # faces lower coordinates. That stays true whichever way the
            # element is inserted, so the extracted sequence is unaffected by
            # the routing decision below.
            left_role, right_role = _pair_roles(leftHit, rightHit, config)
            left_rec = build_flank_record(
                leftHit, is_left=True, record_id=pair_id, role=left_role
            )
            right_rec = build_flank_record(
                rightHit, is_left=False, record_id=pair_id, role=right_role
            )

            # Routing is by terminus ROLE, so each model's termini always land
            # in that model's file. For a reverse-inserted asymmetric element
            # the roles are swapped relative to genomic order; grouping by
            # position would put the right model's terminus in the left
            # model's file, mixing the two models' sequences in one output.
            for rec, hit, role in (
                (left_rec, leftHit, left_role),
                (right_rec, rightHit, right_role),
            ):
                if not rec:
                    continue
                model_key = (
                    (config.left_model if role == 'left' else config.right_model)
                    if config is not None and config.is_asymmetric
                    else hit.model
                )
                suffix = 'L' if role == 'left' else 'R'
                if role == 'left':
                    left_flanks.setdefault(model_key, []).append(rec)
                    paired_left_flanks.setdefault(model_key, []).append(
                        _make_paired_flank_record(rec, element_id, pair_id, suffix)
                    )
                else:
                    right_flanks.setdefault(model_key, []).append(rec)
                    paired_right_flanks.setdefault(model_key, []).append(
                        _make_paired_flank_record(rec, element_id, pair_id, suffix)
                    )

            paired_hit_ids.add(leftHit.idx)
            paired_hit_ids.add(rightHit.idx)

    # ------------------------------------------------------------------
    # Process unpaired hits (only when write_all=True)
    # ------------------------------------------------------------------
    # Determine if this is a symmetric same-strand pairing (F,F or R,R)
    is_symmetric_same_strand = (
        config is not None
        and not config.is_asymmetric
        and config.left_strand == config.right_strand
    )

    if write_all and is_symmetric_same_strand:
        logging.warning(
            'Symmetric same-strand orientation detected (%s). '
            'Cannot determine the outer edge for unpaired hits from a single model. '
            'Unpaired hits will be skipped in --flanks output. '
            'Use --flanks-paired for external flanks of confirmed paired termini.',
            ','.join(config.orientation),
        )

    # For asymmetric pairings the hitIndex may contain hits for models that
    # belong to *other* pairs in a multi-pair pairing-map run.  These foreign
    # models should be silently skipped here because their terminus type is
    # determined correctly in their own pair's writeFlanks call.
    config_models: Optional[Set[str]] = None
    if config is not None and config.is_asymmetric:
        # is_asymmetric guarantees both left_model and right_model are non-None,
        # but we filter defensively in case of any future config construction path.
        config_models = {
            m for m in (config.left_model, config.right_model) if m is not None
        }

    if write_all:
        for model in hitIndex.keys():
            # Skip models that do not belong to the current asymmetric pair.
            if config_models is not None and model not in config_models:
                continue

            for hit_id, hit_data in hitIndex[model].items():
                if hit_data['partner'] is not None:
                    continue  # already handled above
                hit = hit_data['rec']

                if is_symmetric_same_strand:
                    # For F,F or R,R symmetric pairings, we cannot determine which
                    # side of an unpaired hit is the external (outer) flank without
                    # knowing whether it is the left or right terminus.  Writing both
                    # flanks would include an internal flank, so we skip unpaired hits
                    # entirely and advise the user to use --flanks-paired instead.
                    logging.debug(
                        f'Skipping unpaired hit {hit_id} (model={hit.model}): '
                        'cannot determine external flank in symmetric same-strand mode'
                    )
                    continue
                else:
                    terminus = resolve_terminus(hit, config)
                    if terminus is None:
                        logging.debug(
                            f'Cannot determine terminus type for unpaired hit {hit.idx} '
                            f'(model={hit.model}, strand={hit.strand}), skipping'
                        )
                        continue

                    # Geometry follows the genomic side the outer edge faces;
                    # grouping follows the terminus role. These differ for a
                    # reverse-inserted asymmetric element, where the left
                    # model's hit lies at the higher coordinate.
                    record_id = f'{model}_{hit_id}_unpaired'
                    rec = build_flank_record(
                        hit, is_left=terminus.is_lower, record_id=record_id
                    )
                    if rec:
                        if terminus.role == 'left':
                            left_flanks.setdefault(hit.model, []).append(rec)
                        else:
                            right_flanks.setdefault(hit.model, []).append(rec)

    # ------------------------------------------------------------------
    # Write output files (all flanks) - only when write_all=True
    # ------------------------------------------------------------------
    if write_all:
        for model, flanks in left_flanks.items():
            if flanks:
                outfile = os.path.join(
                    outDir,
                    prefix + model + '_left_flank_' + str(len(flanks)) + '.fasta',
                )
                with open(outfile, 'w') as handle:
                    for seq in flanks:
                        seq.id = prefix + str(seq.id)
                        SeqIO.write(seq, handle, 'fasta')

        for model, flanks in right_flanks.items():
            if flanks:
                outfile = os.path.join(
                    outDir,
                    prefix + model + '_right_flank_' + str(len(flanks)) + '.fasta',
                )
                with open(outfile, 'w') as handle:
                    for seq in flanks:
                        seq.id = prefix + str(seq.id)
                        SeqIO.write(seq, handle, 'fasta')

    # ------------------------------------------------------------------
    # Write paired-only flank files (with element IDs in headers)
    # Only when write_paired=True
    # ------------------------------------------------------------------
    if write_paired:
        for model, flanks in paired_left_flanks.items():
            if flanks:
                outfile = os.path.join(
                    outDir,
                    prefix
                    + model
                    + '_paired_left_flank_'
                    + str(len(flanks))
                    + '.fasta',
                )
                with open(outfile, 'w') as handle:
                    for seq in flanks:
                        seq.id = prefix + str(seq.id)
                        SeqIO.write(seq, handle, 'fasta')

        for model, flanks in paired_right_flanks.items():
            if flanks:
                outfile = os.path.join(
                    outDir,
                    prefix
                    + model
                    + '_paired_right_flank_'
                    + str(len(flanks))
                    + '.fasta',
                )
                with open(outfile, 'w') as handle:
                    for seq in flanks:
                        seq.id = prefix + str(seq.id)
                        SeqIO.write(seq, handle, 'fasta')


def hamming_distance(seq1: str, seq2: str) -> int:
    """
    Compute the Hamming distance between two equal-length strings.

    Parameters
    ----------
    seq1 : str
        First sequence.
    seq2 : str
        Second sequence.

    Returns
    -------
    int
        Number of positions at which the corresponding characters differ.

    Raises
    ------
    ValueError
        If the sequences are not of equal length.
    """
    if len(seq1) != len(seq2):
        raise ValueError(
            f'Sequences must be equal length, got {len(seq1)} and {len(seq2)}'
        )
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def load_tsd_length_map(tsd_map_file: str) -> Dict[str, int]:
    """
    Load TSD (Target Site Duplication) lengths from a tab-delimited file.

    Parameters
    ----------
    tsd_map_file : str
        Path to tab-delimited file mapping model pair keys to TSD lengths.
        Format: left_model<TAB>right_model<TAB>tsd_length.

    Returns
    -------
    dict
        Dictionary mapping 'left_model<TAB>right_model' keys to integer TSD lengths.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file format is invalid.
    """
    tsd_lengths: Dict[str, int] = {}

    try:
        with open(tsd_map_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) != 3:
                    raise ValueError(
                        f'Invalid format on line {line_num}: expected 3 tab-delimited '
                        f'columns (left_model, right_model, tsd_length), got {len(parts)}'
                    )

                left_model, right_model, tsd_len_str = (
                    parts[0].strip(),
                    parts[1].strip(),
                    parts[2].strip(),
                )
                try:
                    tsd_len = int(tsd_len_str)
                except ValueError:
                    raise ValueError(
                        f'Invalid TSD length on line {line_num}: {tsd_len_str}'
                    ) from None

                if tsd_len < 0:
                    raise ValueError(
                        f'TSD length must be non-negative on line {line_num}: {tsd_len}'
                    )

                key = f'{left_model}\t{right_model}'
                tsd_lengths[key] = tsd_len

    except FileNotFoundError:
        raise FileNotFoundError(
            f'TSD length map file not found: {tsd_map_file}'
        ) from None

    if not tsd_lengths:
        raise ValueError(f'No valid TSD lengths found in {tsd_map_file}')

    logging.info(
        f'Loaded TSD lengths for {len(tsd_lengths)} model pair(s) from {tsd_map_file}'
    )
    return tsd_lengths


def reconstruct_target_site(
    left_flank_seq: str,
    right_flank_seq: str,
    tsd_length: int = 0,
    tsd_in_model: bool = False,
) -> Tuple[str, str, str, Optional[int]]:
    """
    Reconstruct a target site by joining left and right flanking sequences.

    When a TSD (Target Site Duplication) is present, it is de-duplicated
    so only one copy appears in the reconstructed target site.

    This function handles the ``tsd_in_model=False`` case, where the TSD
    is encoded as the innermost ``tsd_length`` bases of the flanking
    sequences (last ``tsd_length`` bases of the left flank and first
    ``tsd_length`` bases of the right flank).  One copy is trimmed before
    joining.

    .. note::

       When ``tsd_in_model=True`` the TSD lies *inside* the terminus hit
       (part of the HMM model), not in the external flanking sequence.
       :func:`writeTargetSites` handles this case by calling
       :func:`compute_inner_tsd_coordinates` to extract the TSD from
       genomic coordinates within the hit before joining the clean flanks.
       Passing ``tsd_in_model=True`` to this function is therefore a
       no-op — it is treated identically to ``tsd_in_model=False``, both
       operating on the boundary of the supplied flank sequences.

    Parameters
    ----------
    left_flank_seq : str
        External flanking sequence upstream of the left terminus.
    right_flank_seq : str
        External flanking sequence downstream of the right terminus.
    tsd_length : int, default 0
        Length of the terminal duplication feature.
    tsd_in_model : bool, default False
        Retained for API compatibility.  Has no effect on the output
        because this function always treats the innermost bases of the
        flanks as the TSD.  When the TSD is inside the model, the correct
        approach is to use :func:`writeTargetSites` which extracts the TSD
        from the hit genomic coordinates via
        :func:`compute_inner_tsd_coordinates`.

    Returns
    -------
    target_site : str
        Reconstructed target site sequence.
    left_tsd : str
        TSD sequence extracted from the left side (empty string if
        ``tsd_length`` is 0).
    right_tsd : str
        TSD sequence extracted from the right side (empty string if
        ``tsd_length`` is 0).
    tsd_hamming : int or None
        Hamming distance between the left and right TSD sequences over
        informative (non-padded) positions. ``0`` when ``tsd_length`` is 0.
        ``None`` when the TSDs cannot be compared — unequal lengths, or every
        position padded at a contig boundary — which callers must render as
        "unverified" rather than as a distance of 0.

    Notes
    -----
    TSD is at the inner boundary of each flank:
        - ``left_tsd``  = last ``tsd_length`` bases of ``left_flank_seq``
        - ``right_tsd`` = first ``tsd_length`` bases of ``right_flank_seq``

    One copy (the right TSD) is trimmed before joining, yielding:
        ``target_site = left_flank_seq + right_flank_seq[tsd_length:]``
    """
    left_tsd = ''
    right_tsd = ''
    tsd_hamming = 0

    if tsd_length > 0:
        # For both tsd_in_model modes, the TSD appears at the inner boundary
        # of the flanks: the tail of the left flank and the head of the right
        # flank. The distinction (tsd_in_model vs not) affects how the user
        # interprets the duplication relative to the termini model, but the
        # trimming logic is the same: remove one copy from the right flank.
        left_tsd = (
            left_flank_seq[-tsd_length:]
            if len(left_flank_seq) >= tsd_length
            else left_flank_seq
        )
        right_tsd = (
            right_flank_seq[:tsd_length]
            if len(right_flank_seq) >= tsd_length
            else right_flank_seq
        )
        # Trim one TSD copy from the right flank to de-duplicate. The min() is
        # for clarity only - slicing past the end of a shorter flank already
        # yields '' - but it makes the intent explicit: never trim more than is
        # there.
        trimmed_right = right_flank_seq[min(tsd_length, len(right_flank_seq)) :]
        target_site = left_flank_seq + trimmed_right

        # None means the TSDs could not be compared (unequal length, or every
        # position padded). It must stay None rather than collapse to 0, which
        # would report an unverifiable TSD as a perfect duplication.
        tsd_hamming, _compared = compare_tsds(left_tsd, right_tsd)
        if tsd_hamming is None:
            logging.warning(
                f'TSD could not be verified: left={left_tsd or "-"} '
                f'right={right_tsd or "-"}'
            )
        elif tsd_hamming > 0:
            logging.warning(
                f'TSD mismatch (hamming={tsd_hamming}): '
                f'left={left_tsd} right={right_tsd}'
            )
    else:
        target_site = left_flank_seq + right_flank_seq

    return target_site, left_tsd, right_tsd, tsd_hamming


def compare_tsds(
    left_tsd: str,
    right_tsd: str,
    pad_char: str = 'N',
) -> Tuple[Optional[int], int]:
    """
    Compare two TSD sequences, ignoring positions that were padded.

    Parameters
    ----------
    left_tsd, right_tsd : str
        TSD sequences of equal length, possibly containing pad characters where
        the region ran past a contig boundary.
    pad_char : str, default 'N'
        Character marking synthesised positions.

    Returns
    -------
    hamming : int or None
        Hamming distance over informative positions, or None if the TSDs cannot
        be compared at all (different lengths, or no informative position).
    compared : int
        Number of positions actually compared.

    Notes
    -----
    A padded position is not a mismatch and is not a match — no base was
    observed there. Counting it either way misreports the evidence: scoring it
    as a match inflates confidence in a duplication that was never seen, and
    scoring it as a mismatch invents a discrepancy. Both are excluded, and a
    ``None`` result means the comparison could not be made and must not be
    rendered as a distance.
    """
    if not left_tsd or not right_tsd or len(left_tsd) != len(right_tsd):
        return None, 0

    pairs = [
        (a, b)
        for a, b in zip(left_tsd.upper(), right_tsd.upper())
        if a != pad_char and b != pad_char
    ]
    if not pairs:
        return None, 0

    return sum(1 for a, b in pairs if a != b), len(pairs)


def format_interleaved_flanks(
    left_flank_seq: str,
    right_flank_seq: str,
    tsd_length: int = 0,
) -> Tuple[str, str]:
    """
    Format left and right flanks as interleaved gap-padded strings.

    Creates two rows where TSD regions overlap visually:
    - Left flank is right-padded with gaps by the length of the right
      flank minus the TSD overlap.
    - Right flank is left-padded by the length of the left flank minus
      the TSD overlap.

    Parameters
    ----------
    left_flank_seq : str
        Left flanking sequence.
    right_flank_seq : str
        Right flanking sequence.
    tsd_length : int, default 0
        Length of the TSD overlap region.

    Returns
    -------
    left_row : str
        Left flank with gap padding on the right.
    right_row : str
        Right flank with gap padding on the left.

    Examples
    --------
    >>> format_interleaved_flanks('AAAAATSD', 'TSDCCCCC', 3)
    ('AAAAATSD-----', '-----TSDCCCCC')
    """
    overlap = min(tsd_length, len(left_flank_seq), len(right_flank_seq))
    right_pad = len(right_flank_seq) - overlap
    left_pad = len(left_flank_seq) - overlap

    left_row = left_flank_seq + '-' * right_pad
    right_row = '-' * left_pad + right_flank_seq

    return left_row, right_row


def writeTargetSites(
    outDir: Optional[str] = None,
    hitTable: Optional[pd.DataFrame] = None,
    model_lengths: Optional[Dict[str, int]] = None,
    paired: Optional[Dict[str, List[Set[int]]]] = None,
    hitIndex: Optional[Dict[str, Dict[int, Dict[str, Any]]]] = None,
    config: Any = None,
    genome: Any = None,
    prefix: Optional[str] = None,
    flank_len: int = 0,
    flank_max_offset: Optional[int] = None,
    tsd_length: int = 0,
    tsd_in_model: bool = False,
    tsd_length_map: Optional[Dict[str, int]] = None,
    genome_descriptions: Optional[Dict[str, str]] = None,
    blastdb: Optional[str] = None,
    pad_flanks: bool = True,
) -> None:
    """
    Reconstruct and write target sites for paired terminus hits.

    Extracts external flanking sequences for each pair, de-duplicates the
    TSD feature, and writes reconstructed target sites, interleaved flanks,
    and a summary report.

    Parameters
    ----------
    outDir : str, optional
        Output directory for target site FASTA files.
    hitTable : pandas.DataFrame
        DataFrame with columns model, target, hitStart, hitEnd, strand, evalue,
        hmmStart, hmmEnd.
    model_lengths : dict
        Dictionary mapping model name to model length.
    paired : dict
        Paired hits: paired[model] = [list of pair sets {id1, id2}].
    hitIndex : dict
        Hit index: hitIndex[model][idx] = {rec, partner, candidates}.
    config : PairingConfig, optional
        Configuration with orientation and model assignments. May be None
        when using pairing-map mode.
    genome : pyfaidx.Fasta, optional
        Indexed genome for sequence extraction.
    prefix : str, optional
        Prefix for output filenames.
    flank_len : int, default 0
        Number of bases to extract in each flanking region.
    flank_max_offset : int, optional
        Maximum allowed offset from model end.
    tsd_length : int, default 0
        Default TSD length (used when tsd_length_map is not provided).
    tsd_in_model : bool, default False
        Whether the TSD is inside the termini model.
    tsd_length_map : dict, optional
        Map of 'left_model\\tright_model' to TSD length.
    genome_descriptions : dict, optional
        Dictionary mapping sequence IDs to descriptions.
    blastdb : str, optional
        Path to BLAST database.
    pad_flanks : bool, default True
        Pad flanks and in-model TSDs with N where they run past a contig
        boundary. Reconstructed sites built from padded sequence are marked
        with a ``padded=`` field, and a TSD that cannot be verified reports
        ``tsd_hamming=NA`` rather than 0.

    Returns
    -------
    None
        Writes FASTA files to disk.

    Notes
    -----
    Output files, written per model pair:
      {prefix}{pair_label}_target_sites_{N}.fasta – reconstructed target sites
      {prefix}{pair_label}_interleaved_flanks_{N}.fasta – interleaved flanks
    """
    assert outDir is not None, 'outDir cannot be None'
    assert hitTable is not None, 'hitTable cannot be None'
    assert model_lengths is not None, 'model_lengths cannot be None'
    assert paired is not None, 'paired cannot be None'
    assert hitIndex is not None, 'hitIndex cannot be None'
    source = make_source(
        genome=genome, blastdb=blastdb, descriptions=genome_descriptions
    )

    if prefix:
        prefix_str = cleanID(prefix) + '_'
    else:
        prefix_str = ''

    target_site_records: List[Any] = []
    interleaved_records: List[Any] = []

    # ------------------------------------------------------------------
    # Helper: look up hmmStart / hmmEnd for a hit by its DataFrame index
    # ------------------------------------------------------------------
    def get_hmm_coords(hit_idx: int) -> Tuple[Optional[int], Optional[int]]:
        """
        Retrieve HMM alignment coordinates for a hit.

        Parameters
        ----------
        hit_idx : int
            DataFrame row index for the hit.

        Returns
        -------
        tuple of (int or None, int or None)
            (hmmStart, hmmEnd), or (None, None) if unavailable.
        """
        if hit_idx not in hitTable.index:  # type: ignore[union-attr]
            return None, None
        row = hitTable.loc[hit_idx]  # type: ignore[union-attr]
        try:
            h_start = int(row['hmmStart'])
            h_end = int(row['hmmEnd'])
        except (ValueError, TypeError):
            return None, None
        return h_start, h_end

    # ------------------------------------------------------------------
    # Helper: retrieve hit record from nested or flat hitIndex
    # ------------------------------------------------------------------
    is_nested = bool(hitIndex) and isinstance(next(iter(hitIndex.values())), dict)

    def get_hit_record(hit_id: int) -> Any:
        """
        Retrieve a hit record from a nested or flat hitIndex.

        Parameters
        ----------
        hit_id : int
            Index of the hit record.

        Returns
        -------
        namedtuple
            Hit record with model, target, hitStart, hitEnd, strand, idx, evalue.

        Raises
        ------
        KeyError
            If hit_id is not present in a nested hitIndex.
        """
        if is_nested:
            for _m, model_hits in hitIndex.items():
                if hit_id in model_hits:
                    return model_hits[hit_id]['rec']
            raise KeyError(f'Hit ID {hit_id} not found in hitIndex')
        return hitIndex[hit_id]['rec']  # type: ignore[index]

    # ------------------------------------------------------------------
    # Helper: extract flank sequence for one hit
    # ------------------------------------------------------------------
    def extract_flank(hit: Any, is_left: bool) -> Optional[FlankResult]:
        """
        Extract the external flank for one terminus hit of a pair.

        Parameters
        ----------
        hit : namedtuple
            Hit record with target, hitStart, hitEnd, strand, idx.
        is_left : bool
            True if the hit is the left terminus of the pair.

        Returns
        -------
        FlankResult or None
            The flank, or None if it could not be extracted.
        """
        hmm_start, hmm_end = get_hmm_coords(hit.idx)
        return extract_terminus_flank(
            source=source,
            hit=hit,
            is_left=is_left,
            model_len=model_lengths.get(hit.model) if model_lengths else None,  # type: ignore[union-attr]
            hmm_start=hmm_start,
            hmm_end=hmm_end,
            flank_len=flank_len,
            flank_max_offset=flank_max_offset,
            pad=pad_flanks,
        )

    # ------------------------------------------------------------------
    # Helper: resolve TSD length for a pair of models
    # ------------------------------------------------------------------
    def get_tsd_length_for_pair(left_model: str, right_model: str) -> int:
        """
        Resolve the TSD length to use for a pair of models.

        Parameters
        ----------
        left_model, right_model : str
            Model names of the paired termini.

        Returns
        -------
        int
            TSD length from the map, falling back to the default tsd_length.
        """
        if tsd_length_map:
            key = f'{left_model}\t{right_model}'
            if key in tsd_length_map:
                return tsd_length_map[key]
            # Try symmetric key
            key_sym = f'{right_model}\t{left_model}'
            if key_sym in tsd_length_map:
                return tsd_length_map[key_sym]
            logging.warning(
                f'No TSD length found for model pair ({left_model}, {right_model}) '
                f'in TSD length map, using default tsd_length={tsd_length}'
            )
        return tsd_length

    # ------------------------------------------------------------------
    # Helper: extract TSD sequence from the inner boundary of a terminus hit
    # Used when tsd_in_model=True (TSD is part of the termini model, not in flank)
    # ------------------------------------------------------------------
    def extract_inner_tsd(hit: Any, is_left: bool, tsd_len: int) -> Optional[str]:
        """
        Extract the TSD sequence from the inner (element-facing) boundary of a hit.

        When ``tsd_in_model=True`` the TSD is encoded at the element-facing end
        of the terminus HMM and therefore lies *within* the hit coordinates
        rather than in the external flanking sequence.  This helper computes the
        genomic coordinates of the TSD and extracts the bases from the genome or
        BLAST database.

        Parameters
        ----------
        hit : namedtuple
            Hit record (model, target, hitStart, hitEnd, strand, idx).
        is_left : bool
            True if the hit is the left (upstream) terminus of the element.
        tsd_len : int
            Length of the TSD to extract.

        Returns
        -------
        PaddedRegion or None
            The TSD region (always on the forward/+ strand), or None if the
            coordinates cannot be determined or the region does not overlap the
            contig.

        Notes
        -----
        The region is padded to ``tsd_len`` when it runs past a contig
        boundary. A padded TSD is not evidence of a duplication and callers
        must exclude the padded positions from any comparison — see
        ``_compare_tsds``.
        """
        model = hit.model
        model_len = model_lengths.get(model) if model_lengths else None  # type: ignore[union-attr]
        if model_len is None:
            logging.warning(
                f'Model length not found for {model}, skipping TSD extraction'
            )
            return None

        hmm_start, hmm_end = get_hmm_coords(hit.idx)
        if hmm_start is None or hmm_end is None:
            logging.debug(
                f'HMM coordinates unavailable for hit {hit.idx}, skipping TSD'
            )
            return None

        tsd_start, tsd_end = compute_inner_tsd_coordinates(
            hit_start=int(hit.hitStart),
            hit_end=int(hit.hitEnd),
            strand=hit.strand,
            is_left_terminus=is_left,
            hmm_start=hmm_start,
            hmm_end=hmm_end,
            model_len=model_len,
            tsd_length=tsd_len,
        )

        if tsd_end < tsd_start:
            logging.debug(
                f'Invalid TSD coordinates for hit {hit.idx}: {tsd_start}-{tsd_end}'
            )
            return None

        # Pad rather than clamp. A TSD silently returned short at a contig
        # boundary used to be compared against a full-length partner, and the
        # length mismatch made the comparison be skipped entirely and reported
        # as hamming 0 - a truncated, unverifiable TSD was indistinguishable
        # from a perfect match.
        return fetch_region_padded(source, hit.target, tsd_start, tsd_end)

    # ------------------------------------------------------------------
    # Determine canonical pair key for file naming
    # ------------------------------------------------------------------
    # Use config model assignments when available so that all pairs for the
    # same model combination are grouped into one output file regardless of
    # which model happens to be at lower genomic coordinates.
    if config is not None and config.left_model is not None:
        _canonical_pair_key = f'{config.left_model}\t{config.right_model}'
    else:
        _canonical_pair_key = None  # will be derived per-pair below

    # ------------------------------------------------------------------
    # Process paired hits – group by model pair
    # ------------------------------------------------------------------
    # Records grouped by model-pair key for per-pair output files
    pair_key_records: Dict[str, List[Any]] = {}
    pair_key_interleaved: Dict[str, List[Any]] = {}

    for model in paired.keys():
        model_counter = 0
        for pair in paired[model]:
            model_counter += 1
            hit_ids = list(pair)
            x_id, y_id = hit_ids[0], hit_ids[1]
            x = get_hit_record(x_id)
            y = get_hit_record(y_id)
            leftHit, rightHit = flipTIRs(x, y)
            pair_id = f'{prefix_str}{model}_{model_counter}'

            left_flank = extract_flank(leftHit, is_left=True)
            right_flank = extract_flank(rightHit, is_left=False)

            if left_flank is None or right_flank is None:
                logging.debug(
                    f'Could not extract flanks for pair {pair_id}, skipping target site'
                )
                continue

            left_seq = left_flank.seq
            right_seq = right_flank.seq
            # Any padded base in the flanks, or in an in-model TSD below, makes
            # the reconstructed target site partly synthetic.
            padded = left_flank.is_padded or right_flank.is_padded

            pair_tsd_len = get_tsd_length_for_pair(leftHit.model, rightHit.model)

            if tsd_in_model and pair_tsd_len > 0:
                # TSD is inside the terminus model – extract from the inner boundary
                # of each hit, not from the external flank.
                # Reconstruction: left_flank + TSD + right_flank
                left_region = extract_inner_tsd(
                    leftHit, is_left=True, tsd_len=pair_tsd_len
                )
                right_region = extract_inner_tsd(
                    rightHit, is_left=False, tsd_len=pair_tsd_len
                )
                left_tsd = left_region.seq if left_region else ''
                right_tsd = right_region.seq if right_region else ''
                tsd_padded = bool(
                    (left_region and left_region.is_padded)
                    or (right_region and right_region.is_padded)
                )
                padded = padded or tsd_padded

                # Use left TSD (arbitrarily) to fill the gap; warn if mismatched
                tsd_seq = left_tsd if left_tsd else right_tsd
                target_site = left_seq + tsd_seq + right_seq

                tsd_hamming, tsd_compared = compare_tsds(left_tsd, right_tsd)
                if tsd_hamming is None:
                    logging.warning(
                        f'TSD for pair {pair_id} could not be verified '
                        f'(left={left_tsd or "-"} right={right_tsd or "-"}); '
                        'reporting as unverified rather than as a match'
                    )
                elif tsd_hamming > 0:
                    logging.warning(
                        f'TSD mismatch for pair {pair_id} '
                        f'(hamming={tsd_hamming} over {tsd_compared}bp): '
                        f'left={left_tsd} right={right_tsd}'
                    )
            else:
                # TSD is outside the model – it is the innermost n bases of each flank.
                # reconstruct_target_site() trims one copy and joins.
                target_site, left_tsd, right_tsd, tsd_hamming = reconstruct_target_site(
                    left_flank_seq=left_seq,
                    right_flank_seq=right_seq,
                    tsd_length=pair_tsd_len,
                    tsd_in_model=False,
                )

            # Report the header in terms of terminus ROLE, so left_model and
            # left_flank_hit always describe the same terminus. For a
            # reverse-inserted asymmetric element that is the hit at the higher
            # coordinate, which element_orientation records explicitly.
            left_role, _right_role = _pair_roles(leftHit, rightHit, config)
            reversed_insertion = left_role == 'right'
            role_left_hit = rightHit if reversed_insertion else leftHit
            role_right_hit = leftHit if reversed_insertion else rightHit

            # Build metadata for FASTA header. A TSD that could not be compared
            # is reported as 'NA', never as 0 - a padded or absent TSD must not
            # read as a confirmed perfect duplication.
            meta_parts = [
                f'flank_len={flank_len}',
                f'tsd_len={pair_tsd_len}',
                f'tsd_in_model={tsd_in_model}',
                f'left_model={role_left_hit.model}',
                f'right_model={role_right_hit.model}',
                f'contig={leftHit.target}',
                f'element_orientation={"reverse" if reversed_insertion else "forward"}',
                f'left_flank_hit={role_left_hit.strand}:{role_left_hit.hitStart}_{role_left_hit.hitEnd}',
                f'right_flank_hit={role_right_hit.strand}:{role_right_hit.hitStart}_{role_right_hit.hitEnd}',
                f'tsd_hamming={"NA" if tsd_hamming is None else tsd_hamming}',
            ]
            if padded:
                meta_parts.append(
                    f'padded=left:{left_flank.left_pad},{left_flank.right_pad};'
                    f'right:{right_flank.left_pad},{right_flank.right_pad}'
                )
            if left_tsd:
                meta_parts.append(f'left_tsd={left_tsd}')
            if right_tsd:
                meta_parts.append(f'right_tsd={right_tsd}')

            description = ' '.join(meta_parts)

            ts_record = SeqRecord(
                Seq.Seq(target_site),
                id=pair_id,
                name=pair_id,
                description=description,
            )
            target_site_records.append(ts_record)

            # Group by model pair for per-pair output — use canonical key
            pair_key = (
                _canonical_pair_key
                if _canonical_pair_key is not None
                else f'{leftHit.model}\t{rightHit.model}'
            )
            pair_key_records.setdefault(pair_key, []).append(ts_record)

            # Build interleaved flanks
            left_row, right_row = format_interleaved_flanks(
                left_flank_seq=left_seq,
                right_flank_seq=right_seq,
                tsd_length=pair_tsd_len,
            )

            il_left = SeqRecord(
                Seq.Seq(left_row),
                id=f'{pair_id}_left',
                name=f'{pair_id}_left',
                description=description,
            )
            il_right = SeqRecord(
                Seq.Seq(right_row),
                id=f'{pair_id}_right',
                name=f'{pair_id}_right',
                description=description,
            )
            interleaved_records.append(il_left)
            interleaved_records.append(il_right)

            # Group by model pair for per-pair output — use same canonical key
            pair_key_interleaved.setdefault(pair_key, []).extend([il_left, il_right])

    # ------------------------------------------------------------------
    # Helper: write single-line non-wrapped FASTA
    # ------------------------------------------------------------------
    def _write_single_line_fasta(records: List[Any], filepath: str) -> None:
        """
        Write SeqRecords as single-line non-wrapped FASTA.

        Parameters
        ----------
        records : list
            SeqRecords to write.
        filepath : str
            Destination path.

        Returns
        -------
        None
            Writes to disk.
        """
        with open(filepath, 'w') as handle:
            for rec in records:
                handle.write(f'>{rec.id} {rec.description}\n')
                handle.write(f'{str(rec.seq)}\n')

    # ------------------------------------------------------------------
    # Write output files – one per model pair
    # ------------------------------------------------------------------
    if target_site_records:
        # Write per-pair target site files
        for pair_key, records in pair_key_records.items():
            left_model, right_model = pair_key.split('\t')
            pair_label = (
                f'{left_model}_{right_model}'
                if left_model != right_model
                else left_model
            )
            ts_outfile = os.path.join(
                outDir, f'{prefix_str}{pair_label}_target_sites_{len(records)}.fasta'
            )
            _write_single_line_fasta(records, ts_outfile)
            logging.info(
                f'Wrote {len(records)} reconstructed target sites to {ts_outfile}'
            )
    else:
        logging.warning('No target sites could be reconstructed')

    if interleaved_records:
        # Write per-pair interleaved flank files
        for pair_key, records in pair_key_interleaved.items():
            left_model, right_model = pair_key.split('\t')
            pair_label = (
                f'{left_model}_{right_model}'
                if left_model != right_model
                else left_model
            )
            il_outfile = os.path.join(
                outDir,
                f'{prefix_str}{pair_label}_interleaved_flanks_{len(records)}.fasta',
            )
            _write_single_line_fasta(records, il_outfile)
            logging.info(f'Wrote interleaved flanks to {il_outfile}')


def fetchUnpaired(
    hitIndex: Optional[Dict[str, Dict[int, Dict[str, Any]]]] = None,
) -> List[Any]:
    """
    Create GFF3-formatted records for all unpaired (orphan) TIR hits.

    Parameters
    ----------
    hitIndex : dict
        Hit index dictionary: hitIndex[model][idx] = {rec, partner, candidates}.

    Returns
    -------
    list
        List of namedtuples representing orphan TIR features. Each tuple has fields:
        model, chromosome, start, end, strand, type ('orphan_term'), id, leftHit,
        rightHit, seq, evalue.

    Notes
    -----
    Unpaired hits are marked with type 'orphan_term' to distinguish them from
    paired terminal repeats in GFF3 output.
    """
    assert hitIndex is not None, 'hitIndex cannot be None'
    orphans = []
    gffTup = namedtuple(
        'gffTup',  # type: ignore[name-match]
        [
            'model',
            'chromosome',
            'start',
            'end',
            'strand',
            'type',
            'id',
            'leftHit',
            'rightHit',
            'seq',
            'evalue',
        ],
    )
    for model in hitIndex.keys():
        for recID in hitIndex[model].keys():
            if hitIndex[model][recID]['partner'] is None:
                x = hitIndex[model][recID]['rec']
                orphan = gffTup(
                    x.model,
                    x.target,
                    x.hitStart,
                    x.hitEnd,
                    x.strand,
                    'orphan_term',
                    x.idx,
                    None,
                    None,
                    None,
                    x.evalue,
                )
                orphans.append(orphan)
    return orphans


def gffWrite(
    outpath: Optional[str] = None,
    featureList: Optional[Dict[str, List[Any]]] = None,
    writeTIRs: Union[bool, str] = True,
    unpaired: Optional[List[Any]] = None,
    suppressMeta: bool = False,
    prefix: Optional[str] = None,
) -> None:
    """
    Write transposon elements and terminal repeats to GFF3 format file.

    Parameters
    ----------
    outpath : str, optional
        Path to output GFF3 file. If None, uses 'tirmite_features.gff3' in cwd.
    featureList : dict, optional
        Dictionary of element records keyed by model.
    writeTIRs : bool or str, default True
        Controls TIR output: True/'all' (all TIRs), 'paired' (only from elements),
        'unpaired' (only orphans), False (no TIRs).
    unpaired : list, optional
        List of orphan TIR records from fetchUnpaired().
    suppressMeta : bool, default False
        If True, suppresses metadata headers (currently unused).
    prefix : str, optional
        Prefix to add to feature IDs.

    Returns
    -------
    None
        Writes GFF3 file to disk but returns nothing.

    Notes
    -----
    Output includes element features with optional child TIR features.
    Elements have type 'Element', paired TIRs 'paired_term', orphans 'orphan_term'.
    Features are sorted by model, chromosome, start, end.
    """
    if featureList is None:
        featureList = {}  # type: ignore[assignment]
    if prefix:
        prefix = cleanID(prefix) + '_'
    else:
        prefix = ''
    # If path to output gff3 file not provided, set default location.
    if not outpath:
        outpath = os.path.join(os.getcwd(), 'tirmite_features.gff3')
    # Unpack element dict to list
    all_features = []
    for model in featureList.keys():
        for record in featureList[model]:
            all_features.append(record)
    # Add list of unpaired TIRs to main featureList if provided.
    if unpaired:
        all_features = all_features + unpaired
    # Sort features
    sorted_features = sorted(
        all_features, key=attrgetter('model', 'chromosome', 'start', 'end')
    )
    # Open GFF handle
    with open(outpath, 'w') as file:
        # Write headers
        file.write('##gff-version 3' + '\n')
        file.write(
            '\t'.join(
                [
                    '#seqid',
                    'source',
                    'type',
                    'start',
                    'end',
                    'score',
                    'strand',
                    'phase',
                    'attributes',
                ]
            )
            + '\n'
        )
        # Format features for GFF3
        for feature in sorted_features:
            if feature.type == 'orphan_term' and writeTIRs in ['all', 'unpaired']:
                file.write(
                    f'{feature.chromosome}\t'
                    f'tirmite\t'
                    f'{feature.type}\t'
                    f'{feature.start}\t'
                    f'{feature.end}\t'
                    f'.\t'
                    f'{feature.strand}\t'
                    f'.\t'
                    f'ID={prefix}{feature.model}_{feature.id};'
                    f'model={feature.model};'
                    f'evalue={feature.evalue};\n'
                )
            if feature.type == 'Element':
                # Write Element line
                # Fix: Create proper element ID with prefix
                element_id = f'{prefix}{feature.model}_{feature.id}'
                file.write(
                    f'{feature.chromosome}\t'
                    f'tirmite\t'
                    f'{feature.type}\t'
                    f'{feature.start}\t'
                    f'{feature.end}\t'
                    f'.\t'
                    f'{feature.strand}\t'
                    f'.\t'
                    f'ID={prefix}{element_id};model={feature.model};\n'
                )
                if writeTIRs in ['all', 'paired']:
                    # Write left TIR line as child
                    left_hit = feature.leftHit
                    # Fix: Use the actual hit's model name, not the element's model
                    left_model = left_hit.model
                    file.write(
                        (
                            f'{left_hit.target}\t'
                            f'tirmite\t'
                            f'paired_term\t'
                            f'{left_hit.hitStart}\t'
                            f'{left_hit.hitEnd}\t'
                            f'.\t'
                            f'{left_hit.strand}\t'
                            f'.\t'
                            f'ID={prefix}{left_model}_{left_hit.idx};'
                            f'model={left_model};'
                            f'Parent={element_id};'
                            f'evalue={left_hit.evalue};\n'
                        )
                    )
                    # Write right TIR line as child
                    right_hit = feature.rightHit
                    # Fix: Use the actual hit's model name, not the element's model
                    right_model = right_hit.model
                    file.write(
                        (
                            f'{right_hit.target}\t'
                            f'tirmite\t'
                            f'paired_term\t'
                            f'{right_hit.hitStart}\t'
                            f'{right_hit.hitEnd}\t'
                            f'.\t'
                            f'{right_hit.strand}\t'
                            f'.\t'
                            f'ID={prefix}{right_model}_{right_hit.idx};'
                            f'model={right_model};'
                            f'Parent={element_id};'
                            f'evalue={right_hit.evalue};\n'
                        )
                    )


# gffTup fields: 'model', 'chromosome', 'start', 'end', 'strand', 'type', 'id', 'score','bias', 'evalue', 'leftHit' , 'rightHit', 'eleSeq'
# Types: "Element", "orphan_term"


"""
# Useful attributes of pymummer objects:
[x.ref_start for x in alignments]
[x.ref_end for x in alignments]
[x.qry_start for x in alignments]
[x.qry_end for x in alignments]
[x.hit_length_ref for x in alignments]
[x.hit_length_qry for x in alignments]
[x.percent_identity for x in alignments]
[x.ref_length for x in alignments]
[x.qry_length for x in alignments]
[x.frame for x in alignments]
[x.ref_name for x in alignments]
[x.qry_name for x in alignments]

## Don't use these, bizzaro format. Not indexed to 0. Cannot sort as ints.
#coord.reverse_query()
#coord.reverse_reference()
#coord.qry_coords()
#coord.ref_coords()
"""


VALID_ORIENTATION_CODES = ('F', 'R')


def parse_orientation(orientation: str) -> List[str]:
    """
    Parse and validate an orientation string such as 'F,R'.

    Parameters
    ----------
    orientation : str
        Two comma-separated codes, each 'F' (forward, + strand) or 'R'
        (reverse, - strand). Case-insensitive and whitespace-tolerant.

    Returns
    -------
    list of str
        Upper-cased ``[left_code, right_code]``.

    Raises
    ------
    ValueError
        If the value is not exactly two valid codes.

    Notes
    -----
    Validating here means a malformed orientation fails loudly at
    configuration time. Previously the string was split without upper-casing,
    so 'f,r' silently produced two '-' strands, and any unrecognised value left
    the strand-combination table empty in
    :func:`parseHitsGeneral`, yielding zero pairs with no explanation.
    """
    if not isinstance(orientation, str):
        raise ValueError(
            f'Orientation must be a string like "F,R", got {type(orientation).__name__}'
        )

    codes = [part.strip().upper() for part in orientation.split(',')]

    if len(codes) != 2 or any(c not in VALID_ORIENTATION_CODES for c in codes):
        raise ValueError(
            f'Invalid orientation {orientation!r}. Expected two comma-separated '
            "codes, each 'F' (forward/+) or 'R' (reverse/-), "
            "e.g. 'F,R' for TIRs or 'F,F' for LTRs."
        )

    return codes


# New configuration class to manage pairing rules
class PairingConfig:
    """
    Configuration for terminal repeat element pairing rules.

    Manages orientation and model selection for symmetric or asymmetric
    transposon terminal repeat pairing.

    Parameters
    ----------
    orientation : str, default 'F,R'
        Comma-separated pair of orientation codes: F=Forward(+), R=Reverse(-).
        Examples: 'F,R' (forward-reverse), 'F,F' (both forward), 'R,R' (both reverse).
    left_model : str, optional
        HMM model name for left terminus. Used for asymmetric pairing.
    right_model : str, optional
        HMM model name for right terminus. Used for asymmetric pairing.
    single_model : str, optional
        HMM model name when using same model for both termini (symmetric pairing).

    Attributes
    ----------
    orientation : list of str
        Parsed orientation codes as list [left_orient, right_orient].
    left_strand : str
        Strand symbol for left terminus: '+' or '-'.
    right_strand : str
        Strand symbol for right terminus: '+' or '-'.
    is_asymmetric : bool
        True if using different models for left and right termini.
    left_model : str
        Model name for left terminus.
    right_model : str
        Model name for right terminus.

    Methods
    -------
    get_model_pairs()
        Returns list of (left_model, right_model) tuples for pairing analysis.
    """

    def __init__(
        self,
        orientation: str = 'F,R',
        left_model: Optional[str] = None,
        right_model: Optional[str] = None,
        single_model: Optional[str] = None,
    ) -> None:
        """
        Configure pairing rules for terminal repeat elements.

        Args:
            orientation: String like 'F,R', 'F,F', 'R,R', 'R,F'
                        F=Forward(+), R=Reverse(-)
            left_model: Model name for left terminus (None for symmetric)
            right_model: Model name for right terminus (None for symmetric)
            single_model: Model name when using same model for both ends
        """
        self.orientation = parse_orientation(orientation)
        self.left_strand = '+' if self.orientation[0] == 'F' else '-'
        self.right_strand = '+' if self.orientation[1] == 'F' else '-'

        # Determine if using symmetric (same model) or asymmetric (different models)
        self.is_asymmetric = left_model is not None and right_model is not None
        self.left_model = left_model if self.is_asymmetric else single_model
        self.right_model = right_model if self.is_asymmetric else single_model

    def get_model_pairs(self) -> List[Tuple[Optional[str], Optional[str]]]:
        """
        Get model pairs for pairing analysis.

        Returns
        -------
        list of tuple
            List containing (left_model, right_model) tuples.
            For symmetric pairing, returns [(model, model)].
            For asymmetric pairing, returns [(left_model, right_model)].
        """
        if self.is_asymmetric:
            return [(self.left_model, self.right_model)]
        else:
            # For symmetric models, pair with themselves
            return [(self.left_model, self.left_model)]


def parseHitsGeneral(
    hitsDict: Optional[Dict[str, Dict[str, List[Any]]]] = None,
    hitIndex: Optional[Dict[str, Dict[int, Dict[str, Any]]]] = None,
    maxDist: Optional[int] = None,
    config: Any = None,
) -> Dict[str, Dict[int, Dict[str, Any]]]:
    """
    Populate candidate partners using configurable strand orientations.

    Parameters
    ----------
    hitsDict : dict
        Nested dictionary of hits: hitsDict[model][chromosome] = [hit_records].
    hitIndex : dict
        Nested dictionary tracking pairing: hitIndex[model][idx] = {rec, partner, candidates}.
    maxDist : int, optional
        Maximum distance in base pairs between paired elements. If None, uses infinity.
    config : PairingConfig
        Configuration object specifying orientation and model pairing rules.

    Returns
    -------
    dict
        Updated hitIndex with populated candidate lists respecting config orientation.

    Notes
    -----
    Handles both symmetric (same model) and asymmetric (different models) pairing.
    For each orientation configuration, searches for valid partners on correct strands
    within the specified distance constraint. Supports all strand combinations:
    F,R (canonical), R,F, F,F, and R,R orientations.
    """
    assert hitsDict is not None, 'hitsDict cannot be None'
    assert hitIndex is not None, 'hitIndex cannot be None'
    assert config is not None, 'config cannot be None'
    logging.debug('=== ENTERING parseHitsGeneral ===')
    logging.debug(
        f'Config: orientation={config.orientation}, left_strand={config.left_strand}, right_strand={config.right_strand}'
    )
    logging.debug(f'Is asymmetric: {config.is_asymmetric}')

    if not maxDist:
        maxDist_value: float = float('inf')
        logging.debug('Using infinite maxDist')
    else:
        maxDist_value = float(maxDist)
        logging.debug(f'Using maxDist: {maxDist_value}')

    model_pairs = config.get_model_pairs()
    logging.debug(f'Model pairs to process: {model_pairs}')

    for left_model, right_model in model_pairs:
        logging.debug(f'Processing pair: {left_model} -> {right_model}')

        if left_model == right_model:
            # Symmetric pairing - enhanced logic for custom orientations
            logging.debug(f'=== SYMMETRIC PAIRING for {left_model} ===')

            if left_model in hitIndex:
                logging.debug(
                    f'Found {len(hitIndex[left_model])} hits for model {left_model}'
                )

                for UID in hitIndex[left_model].keys():
                    ref = hitIndex[left_model][UID]['rec']
                    logging.debug(
                        f'Processing hit {UID}: {ref.strand}:{ref.hitStart}-{ref.hitEnd}'
                    )

                    # For symmetric pairing, a hit can act as either left or right terminus
                    # depending on its strand and the config orientation

                    # A hit may satisfy the left role, the right role, or - when
                    # the orientation is same-strand (F,F / R,R) - both. These
                    # must be independent tests: an if/elif would let only the
                    # left search run for F,F and R,R, so no hit would ever
                    # collect an upstream candidate and reciprocity could never
                    # be established.
                    can_be_left = ref.strand == config.left_strand
                    can_be_right = ref.strand == config.right_strand

                    # The left terminus searches away from its own 5' end: on
                    # '+' that is downstream, on '-' upstream. The right
                    # terminus searches the opposite way.
                    left_direction = (
                        'left_to_right'
                        if config.left_strand == '+'
                        else 'right_to_left'
                    )
                    right_direction = (
                        'right_to_left'
                        if left_direction == 'left_to_right'
                        else 'left_to_right'
                    )

                    if can_be_left:
                        logging.debug(
                            f'Hit {UID} acting as LEFT terminus, searching {left_direction} for RIGHT partners on {config.right_strand} strand'
                        )
                        _find_candidates(
                            ref,
                            right_model,
                            config.right_strand,
                            hitsDict,
                            hitIndex,
                            maxDist_value,
                            left_direction,
                        )

                    if can_be_right:
                        logging.debug(
                            f'Hit {UID} acting as RIGHT terminus, searching {right_direction} for LEFT partners on {config.left_strand} strand'
                        )
                        _find_candidates(
                            ref,
                            left_model,
                            config.left_strand,
                            hitsDict,
                            hitIndex,
                            maxDist_value,
                            right_direction,
                        )

                    if not (can_be_left or can_be_right):
                        logging.debug(
                            f'Hit {UID} on strand {ref.strand} does not match required orientations ({config.left_strand}, {config.right_strand})'
                        )
        else:
            # FIXED: Asymmetric pairing with strand combination handling
            logging.debug(f'=== ASYMMETRIC PAIRING: {left_model} + {right_model} ===')

            # Get all valid strand combinations for this orientation
            strand_combinations = []

            if config.orientation == ['F', 'R']:
                # F,R can appear as: (+,-) on pos strand OR (-,+) on neg strand
                strand_combinations = [('+', '-'), ('-', '+')]
            elif config.orientation == ['R', 'F']:
                # R,F can appear as: (-,+) on pos strand OR (+,-) on neg strand
                strand_combinations = [('-', '+'), ('+', '-')]
            elif config.orientation == ['F', 'F']:
                # F,F can appear as: (+,+) on pos strand OR (-,-) on neg strand
                strand_combinations = [('+', '+'), ('-', '-')]
            elif config.orientation == ['R', 'R']:
                # R,R can appear as: (-,-) on pos strand OR (+,+) on neg strand
                strand_combinations = [('-', '-'), ('+', '+')]

            logging.debug(f'Processing strand combinations: {strand_combinations}')

            for left_strand, right_strand in strand_combinations:
                logging.debug(
                    f'Processing strand combination: left={left_strand}, right={right_strand}'
                )

                # Process hits for the left model with this strand combination
                if left_model in hitIndex:
                    for UID in hitIndex[left_model].keys():
                        ref = hitIndex[left_model][UID]['rec']

                        if ref.strand == left_strand:
                            # Determine search direction based on strand
                            if left_strand == '+':
                                search_direction = 'left_to_right'  # Search downstream
                            else:
                                search_direction = (
                                    'right_to_left'  # Search upstream (neg strand)
                                )

                            logging.debug(
                                f'Left model hit {UID} ({left_strand}) searching {search_direction} for right model on {right_strand}'
                            )

                            _find_candidates(
                                ref,
                                right_model,
                                right_strand,
                                hitsDict,
                                hitIndex,
                                maxDist_value,
                                search_direction,
                            )

                # Process hits for the right model with this strand combination
                if right_model in hitIndex:
                    for UID in hitIndex[right_model].keys():
                        ref = hitIndex[right_model][UID]['rec']

                        if ref.strand == right_strand:
                            # The right model looks for the left model in the opposite
                            # genomic direction from where the left model looks for it.
                            # The direction is determined by the left_strand (where the
                            # left terminus sits relative to the right):
                            #   left(+) is at lower coords  → right looks right_to_left
                            #   left(-) is at higher coords → right looks left_to_right
                            if left_strand == '+':
                                search_direction = 'right_to_left'
                            else:
                                search_direction = 'left_to_right'

                            logging.debug(
                                f'Right model hit {UID} ({right_strand}) searching {search_direction} for left model on {left_strand}'
                            )

                            _find_candidates(
                                ref,
                                left_model,
                                left_strand,
                                hitsDict,
                                hitIndex,
                                maxDist_value,
                                search_direction,
                            )

    logging.debug('=== COMPLETED parseHitsGeneral ===')
    return hitIndex


def inter_hit_distance(ref_hit: Any, candidate: Any, direction: str) -> int:
    """
    Genomic distance between the facing inner edges of two terminus hits.

    Parameters
    ----------
    ref_hit : namedtuple
        Reference hit with hitStart and hitEnd attributes.
    candidate : namedtuple
        Candidate partner hit with hitStart and hitEnd attributes.
    direction : str
        'left_to_right' when the candidate lies downstream of the reference,
        'right_to_left' when it lies upstream.

    Returns
    -------
    int
        Separation between the two hits, i.e. the span of the element interior
        between them. Negative when the hits overlap or are in the wrong order
        for ``direction``.

    Notes
    -----
    This is the quantity ``--maxdist`` limits: the gap between the termini, not
    including the termini themselves. It is measured between the inner edge of
    the upstream hit (its ``hitEnd``) and the inner edge of the downstream hit
    (its ``hitStart``), which makes it independent of strand - the importers
    already normalise every hit so that ``hitStart < hitEnd``.

    The previous implementation measured from the upstream hit's inner edge to
    the *downstream* hit's strand-relative 5' end, which for a minus-strand
    partner is its far edge. That added the whole length of the partner
    terminus to the measured distance, so the same element needed a larger
    ``--maxdist`` purely because its terminus model was longer.
    """
    if direction == 'left_to_right':
        upstream, downstream = ref_hit, candidate
    else:  # 'right_to_left'
        upstream, downstream = candidate, ref_hit

    return int(downstream.hitStart) - int(upstream.hitEnd)


def _check_distance(
    ref_hit: Any, candidate: Any, direction: str, maxDist: float
) -> bool:
    """
    Validate that candidate hit is within distance threshold of reference hit.

    Parameters
    ----------
    ref_hit : namedtuple
        Reference hit with model, strand, hitStart, hitEnd attributes.
    candidate : namedtuple
        Candidate partner hit with model, strand, hitStart, hitEnd attributes.
    direction : str
        Search direction: 'left_to_right' (left looking for right terminus downstream)
        or 'right_to_left' (right looking for left terminus upstream).
    maxDist : float
        Maximum allowed distance in base pairs. Can be infinity for no constraint.

    Returns
    -------
    bool
        True if candidate is valid (positive distance within maxDist), False otherwise.

    Notes
    -----
    Distance is the gap between the two hits, as computed by
    :func:`inter_hit_distance`. A negative value means the candidate is on the
    wrong side of the reference for the search direction, or the two hits
    overlap, and the candidate is rejected.
    """
    distance = inter_hit_distance(ref_hit, candidate, direction)

    logging.debug('=== DISTANCE CHECK DEBUG ===')
    logging.debug(f'Direction: {direction}')
    logging.debug(
        f'Ref hit: {ref_hit.model} strand={ref_hit.strand} coords=({ref_hit.hitStart}, {ref_hit.hitEnd})'
    )
    logging.debug(
        f'Candidate: {candidate.model} strand={candidate.strand} coords=({candidate.hitStart}, {candidate.hitEnd})'
    )
    logging.debug(f'Calculated distance: {distance}')

    # Check for negative distances (invalid pairing)
    if distance < 0:
        logging.debug(
            f'Negative distance ({distance}) between {ref_hit.model} and {candidate.model} '
            f'on {ref_hit.target}. Ref: {ref_hit.strand}:{ref_hit.hitStart}-{ref_hit.hitEnd}, '
            f'Candidate: {candidate.strand}:{candidate.hitStart}-{candidate.hitEnd}'
        )
        return False

    # Check if within max distance
    valid = distance >= 0 and distance <= maxDist
    logging.debug(
        f'Valid distance check: {valid} (distance: {distance}, maxDist: {maxDist})'
    )

    return valid


def _find_candidates(
    ref_hit: Any,
    target_model: str,
    target_strand: str,
    hitsDict: Dict[str, Dict[str, List[Any]]],
    hitIndex: Dict[str, Dict[int, Dict[str, Any]]],
    maxDist: float,
    direction: str,
) -> None:
    """
    Find and store valid candidate partners for a reference hit.

    Parameters
    ----------
    ref_hit : namedtuple
        Reference hit record with model, target, strand, idx, hitStart, hitEnd.
    target_model : str
        Model name to search for candidate partners.
    target_strand : str
        Required strand orientation ('+' or '-') for candidates.
    hitsDict : dict
        Dictionary of all hits: hitsDict[model][chromosome] = [hit_records].
    hitIndex : dict
        Index for storing candidates: hitIndex[model][idx]['candidates'].
    maxDist : float
        Maximum distance constraint for valid partners.
    direction : str
        Search direction: 'left_to_right' or 'right_to_left'.

    Returns
    -------
    None
        Modifies hitIndex in place by appending valid candidates and sorting by distance.

    Notes
    -----
    Candidates are sorted by calculated biological distance with closest partners first.
    Only hits on target chromosome matching target_strand and within maxDist are added.

    The reference hit is never its own candidate. For symmetric same-strand
    orientations (F,F and R,R) the reference and the candidates come from the
    same model on the same strand, so without this the hit would be offered as
    its own partner and could be "paired" with itself.
    """
    import logging

    logging.debug('=== _find_candidates DEBUG ===')
    logging.debug(
        f'Ref hit: {ref_hit.model} {ref_hit.strand}:{ref_hit.hitStart}-{ref_hit.hitEnd}'
    )
    logging.debug(
        f'Looking for target_model: {target_model}, target_strand: {target_strand}'
    )
    logging.debug(f'Direction: {direction}, maxDist: {maxDist}')

    if target_model not in hitsDict or ref_hit.target not in hitsDict[target_model]:
        logging.debug(
            f'No hits found for target_model {target_model} on {ref_hit.target}'
        )
        return

    # Store candidates under the reference hit's model and UID
    model_key = ref_hit.model
    uid_key = ref_hit.idx

    candidates_found = 0

    for candidate in hitsDict[target_model][ref_hit.target]:
        # A hit can never be its own partner. This only bites for symmetric
        # same-strand orientations, where ref and candidate share a model and a
        # strand; on '-' the 5'/3' swap makes the self-distance positive, so the
        # hit would otherwise pass the distance test against itself.
        if candidate.model == ref_hit.model and candidate.idx == ref_hit.idx:
            continue

        if candidate.strand == target_strand:
            logging.debug(
                f'Checking candidate: {candidate.model} {candidate.strand}:{candidate.hitStart}-{candidate.hitEnd}'
            )

            # Calculate distance based on direction and orientation
            valid_distance = _check_distance(ref_hit, candidate, direction, maxDist)

            if valid_distance:
                hitIndex[model_key][uid_key]['candidates'].append(candidate)
                candidates_found += 1
                logging.debug(
                    f'Added valid candidate: {candidate.model}_{candidate.idx}'
                )

    logging.debug(
        f'Found {candidates_found} valid candidates for {ref_hit.model}_{ref_hit.idx}'
    )

    # Sort candidates by distance using the same logic as _check_distance
    # This ensures closest valid partners are prioritized
    if hitIndex[model_key][uid_key]['candidates']:
        # Sort by the same measure the distance filter uses, so "closest" means
        # the same thing when ranking candidates as when accepting them.
        hitIndex[model_key][uid_key]['candidates'] = sorted(
            hitIndex[model_key][uid_key]['candidates'],
            key=lambda x: inter_hit_distance(ref_hit, x, direction),
        )

        logging.debug(
            f'Sorted {len(hitIndex[model_key][uid_key]["candidates"])} candidates by distance'
        )


def iterateGetPairsAsymmetric(
    hitIndex: Dict[str, Dict[int, Dict[str, Any]]], config: Any, stableReps: int = 0
) -> Tuple[Dict[str, Dict[int, Dict[str, Any]]], Dict[str, List[Set[int]]], List[int]]:
    """
    Iterate asymmetric pairing with different left and right HMM models.

    Parameters
    ----------
    hitIndex : dict
        Nested hit index: hitIndex[model][idx] = {rec, partner, candidates}.
    config : PairingConfig
        Configuration with left_model and right_model specified.
    stableReps : int, default 0
        Maximum iterations to continue after no new pairs found.

    Returns
    -------
    hitIndex : dict
        Updated index with partner assignments.
    paired : dict
        Dictionary of pairs: paired[left_model] = [list of pair sets].
    unpaired : list
        List of hit indices that remain unpaired.

    Notes
    -----
    Pairs hits from different HMM models representing left and right termini.
    Handles multiple strand combinations for each orientation configuration.
    Iterates until convergence or stable iteration limit reached.
    """
    import logging

    logging.debug('=== ENTERING iterateGetPairsAsymmetric ===')
    logging.debug(
        f'Config: {config.left_model} ({config.left_strand}) + {config.right_model} ({config.right_strand})'
    )

    # Init stable repeat counter
    reps = 0

    # Initialize paired dict with left model name (convention for asymmetric)
    paired: Dict[str, List[Set[int]]] = {config.left_model: []}

    # Run initial pairing
    hitIndex, paired = getPairsAsymmetric(
        hitIndex=hitIndex, config=config, paired=paired
    )

    # Count remaining unpaired hits
    countUP = countUnpairedAsymmetric(hitIndex, config)

    logging.debug(f'Initial unpaired count: {countUP}')

    # Iterate pairing procedure until either no unpaired remain
    # OR max number of iterations without new pairing is reached
    while countUP > 0 and reps < stableReps:
        # Re-run pairing procedure
        hitIndex, paired = getPairsAsymmetric(
            hitIndex=hitIndex, config=config, paired=paired
        )

        # Store previous unpaired hit count
        lastCountUP = countUP
        # Update unpaired hit count
        countUP = countUnpairedAsymmetric(hitIndex, config)

        logging.debug(
            f'Iteration {reps + 1}: unpaired count {lastCountUP} -> {countUP}'
        )

        # If no change in unpaired hit count, iterate stable rep counter
        if lastCountUP == countUP:
            reps += 1

    # Get IDs of remaining unpaired hits
    unpaired = listunpairedAsymmetric(hitIndex, config)

    total_pairs = sum(len(pairs) for pairs in paired.values())
    logging.debug(
        f'Asymmetric pairing completed: {total_pairs} pairs, {len(unpaired)} unpaired'
    )

    return hitIndex, paired, unpaired


def getPairsAsymmetric(
    hitIndex: Optional[Dict[str, Dict[int, Dict[str, Any]]]] = None,
    config: Any = None,
    paired: Optional[Dict[str, List[Set[int]]]] = None,
) -> Tuple[Dict[str, Dict[int, Dict[str, Any]]], Dict[str, List[Set[int]]]]:
    """
    Perform one round of asymmetric pairing between different models.

    Parameters
    ----------
    hitIndex : dict
        Nested hit index: hitIndex[model][idx] = {rec, partner, candidates}.
    config : PairingConfig
        Configuration specifying left_model, right_model, and strand orientations.
    paired : dict, optional
        Existing pairs dictionary. If None, creates new dictionary.

    Returns
    -------
    hitIndex : dict
        Updated index with new partner assignments.
    paired : dict
        Updated pairs dictionary: paired[left_model] = [list of pair sets].

    Notes
    -----
    Checks reciprocal best-match relationship between left and right model hits.
    Only pairs hits that are each other's closest valid unpaired partners.
    """
    assert hitIndex is not None, 'hitIndex cannot be None'
    assert config is not None, 'config cannot be None'
    import logging

    if not paired:
        paired_dict: Dict[str, List[Set[int]]] = {config.left_model: []}
    else:
        paired_dict = paired

    pairs_found = 0

    # Get hits from left model looking for right model partners
    if config.left_model in hitIndex:
        for leftID in hitIndex[config.left_model].keys():
            if hitIndex[config.left_model][leftID]['partner'] is None:
                left_hit = hitIndex[config.left_model][leftID]['rec']

                # REMOVED: Strand restriction check
                # The parseHitsGeneral already populated candidates with valid combinations
                logging.debug(
                    f'Processing left hit {leftID}: {left_hit.strand}:{left_hit.hitStart}-{left_hit.hitEnd}'
                )

                # Look through candidates (which should be from right model)
                for candidate in hitIndex[config.left_model][leftID]['candidates']:
                    logging.debug(
                        f'Checking candidate: {candidate.model}_{candidate.idx} {candidate.strand}:{candidate.hitStart}-{candidate.hitEnd}'
                    )

                    # Check if candidate is from right model and unpaired
                    if (
                        candidate.model == config.right_model
                        and candidate.idx in hitIndex[config.right_model]
                        and hitIndex[config.right_model][candidate.idx]['partner']
                        is None
                    ):
                        # Check if this left hit is also the best candidate for the right hit
                        found = checkAsymmetricReciprocity(
                            left_model=config.left_model,
                            left_id=leftID,
                            right_model=config.right_model,
                            right_id=candidate.idx,
                            hitIndex=hitIndex,
                            config=config,
                        )

                        if found:
                            # Mark as paired
                            hitIndex[config.left_model][leftID]['partner'] = (
                                candidate.idx
                            )
                            hitIndex[config.right_model][candidate.idx]['partner'] = (
                                leftID
                            )

                            # Add to paired list (store under left model)
                            paired_dict[config.left_model].append(
                                {leftID, candidate.idx}
                            )
                            pairs_found += 1

                            logging.debug(
                                f'Paired: {config.left_model}_{leftID} + {config.right_model}_{candidate.idx}'
                            )
                            break

    logging.debug(f'Found {pairs_found} new asymmetric pairs')
    return hitIndex, paired_dict


def checkAsymmetricReciprocity(
    left_model: str,
    left_id: int,
    right_model: str,
    right_id: int,
    hitIndex: Dict[str, Dict[int, Dict[str, Any]]],
    config: Any,
) -> bool:
    """
    Check if asymmetric pair has reciprocal best-match relationship.

    Parameters
    ----------
    left_model : str
        Model name for left terminus hit.
    left_id : int
        Index of left hit in hitIndex.
    right_model : str
        Model name for right terminus hit.
    right_id : int
        Index of right hit in hitIndex.
    hitIndex : dict
        Hit index dictionary with candidate lists.
    config : PairingConfig
        Configuration with strand requirements.

    Returns
    -------
    bool
        True if left hit is right hit's best unpaired candidate, False otherwise.

    Notes
    -----
    Asymmetric reciprocity requires left hit to be the first valid unpaired
    candidate in right hit's candidate list, accounting for strand compatibility.
    """
    import logging

    # Check if left hit is the best candidate for the right hit
    right_candidates = hitIndex[right_model][right_id]['candidates']

    for candidate in right_candidates:
        if (
            candidate.model == left_model
            # REMOVED: strand compatibility check - already filtered
            and candidate.idx in hitIndex[left_model]
            and hitIndex[left_model][candidate.idx]['partner'] is None
        ):
            if candidate.idx == left_id:
                logging.debug(
                    f'Reciprocal match: {left_model}_{left_id} <-> {right_model}_{right_id}'
                )
                return True  # Reciprocal match found
            else:
                logging.debug(
                    f'Better candidate exists for {right_model}_{right_id}: {candidate.idx}'
                )
                return False  # A better candidate exists

    logging.debug(
        f'No reciprocal match for {left_model}_{left_id} -> {right_model}_{right_id}'
    )
    return False  # No valid candidates


def iterateGetPairsCustom(
    hitIndex: Dict[str, Dict[int, Dict[str, Any]]], config: Any, stableReps: int = 0
) -> Tuple[Dict[str, Dict[int, Dict[str, Any]]], Dict[str, List[Set[int]]], List[int]]:
    """
    Iterate symmetric pairing with custom strand orientations.

    Parameters
    ----------
    hitIndex : dict
        Nested hit index: hitIndex[model][idx] = {rec, partner, candidates}.
    config : PairingConfig
        Configuration with single_model and custom orientation specified.
    stableReps : int, default 0
        Maximum iterations to continue after no new pairs found.

    Returns
    -------
    hitIndex : dict
        Updated index with partner assignments.
    paired : dict
        Dictionary of pairs: paired[model] = [list of pair sets].
    unpaired : list
        List of hit indices that remain unpaired.

    Notes
    -----
    Handles non-standard orientations (F,F or R,R) for symmetric model pairing.
    Useful for elements with non-canonical TIR orientations or inverted structures.
    """
    import logging

    logging.debug('=== ENTERING iterateGetPairsCustom ===')
    logging.debug(
        f'Config: {config.left_model} orientation {config.left_strand},{config.right_strand}'
    )

    model_name = config.left_model

    if model_name not in hitIndex:
        logging.error(f'Model {model_name} not found in hitIndex')
        return hitIndex, {model_name: []}, []

    # Initialize pairing structures
    paired: Dict[str, List[Set[int]]] = {model_name: []}
    reps = 0

    # Run initial pairing
    hitIndex, paired = getPairsSymmetric(
        hitIndex=hitIndex, model_name=model_name, config=config, paired=paired
    )

    # Count remaining unpaired hits
    countUP = countUnpairedSymmetric(hitIndex, model_name, config)

    logging.debug(f'Initial unpaired count: {countUP}')

    # Iterate pairing procedure
    while countUP > 0 and reps < stableReps:
        hitIndex, paired = getPairsSymmetric(
            hitIndex=hitIndex, model_name=model_name, config=config, paired=paired
        )

        lastCountUP = countUP
        countUP = countUnpairedSymmetric(hitIndex, model_name, config)

        logging.debug(
            f'Iteration {reps + 1}: unpaired count {lastCountUP} -> {countUP}'
        )

        if lastCountUP == countUP:
            reps += 1

    # Get unpaired list
    unpaired = listunpairedSymmetric(hitIndex, model_name, config)

    total_pairs = len(paired[model_name])
    logging.debug(
        f'Symmetric pairing completed: {total_pairs} pairs, {len(unpaired)} unpaired'
    )

    return hitIndex, paired, unpaired


def getPairsSymmetric(
    hitIndex: Optional[Dict[str, Dict[int, Dict[str, Any]]]] = None,
    model_name: Optional[str] = None,
    config: Any = None,
    paired: Optional[Dict[str, List[Set[int]]]] = None,
) -> Tuple[Dict[str, Dict[int, Dict[str, Any]]], Dict[str, List[Set[int]]]]:
    """
    Perform one round of symmetric pairing within a single model.

    Parameters
    ----------
    hitIndex : dict
        Nested hit index: hitIndex[model][idx] = {rec, partner, candidates}.
    model_name : str
        Name of HMM model for symmetric pairing.
    config : PairingConfig
        Configuration specifying orientation constraints.
    paired : dict, optional
        Existing pairs dictionary. If None, creates new dictionary.

    Returns
    -------
    hitIndex : dict
        Updated index with new partner assignments.
    paired : dict
        Updated pairs dictionary: paired[model] = [list of pair sets].

    Notes
    -----
    Pairs hits from the same HMM model that meet orientation requirements.
    Each hit must have complementary role (left or right) based on strand
    to form a valid symmetric pair.
    """
    assert hitIndex is not None, 'hitIndex cannot be None'
    assert model_name is not None, 'model_name cannot be None'
    assert config is not None, 'config cannot be None'
    import logging

    if model_name not in hitIndex:
        if paired is None:
            paired = {}
        return hitIndex, paired

    if not paired:
        paired_dict: Dict[str, List[Set[int]]] = {model_name: []}
    else:
        paired_dict = paired

    pairs_found = 0

    for refID in hitIndex[model_name].keys():
        if hitIndex[model_name][refID]['partner'] is None:
            ref_hit = hitIndex[model_name][refID]['rec']

            # Check if this hit can act as a left or right terminus based on strand
            can_be_left = ref_hit.strand == config.left_strand
            can_be_right = ref_hit.strand == config.right_strand

            if not (can_be_left or can_be_right):
                logging.debug(
                    f"Hit {refID} on strand {ref_hit.strand} doesn't match orientation {config.left_strand},{config.right_strand}"
                )
                continue

            logging.debug(
                f'Processing hit {refID}: {ref_hit.strand}:{ref_hit.hitStart}-{ref_hit.hitEnd} (can_be_left: {can_be_left}, can_be_right: {can_be_right})'
            )

            # Check candidates for this hit
            for candidate in hitIndex[model_name][refID]['candidates']:
                logging.debug(
                    f'Checking candidate: {candidate.model}_{candidate.idx} {candidate.strand}:{candidate.hitStart}-{candidate.hitEnd}'
                )

                # Candidate should be from the same model for symmetric pairing
                if (
                    candidate.model == model_name
                    and candidate.idx in hitIndex[model_name]
                    and hitIndex[model_name][candidate.idx]['partner'] is None
                ):
                    # Check strand compatibility for symmetric pairing
                    candidate_can_be_left = candidate.strand == config.left_strand
                    candidate_can_be_right = candidate.strand == config.right_strand

                    # For symmetric pairing, we need complementary roles
                    compatible = False
                    if can_be_left and candidate_can_be_right:
                        compatible = True
                        logging.debug(
                            f'Compatible: {refID} (left) + {candidate.idx} (right)'
                        )
                    elif can_be_right and candidate_can_be_left:
                        compatible = True
                        logging.debug(
                            f'Compatible: {refID} (right) + {candidate.idx} (left)'
                        )

                    if compatible:
                        # Check reciprocity
                        if checkSymmetricReciprocity(
                            model_name, refID, candidate.idx, hitIndex, config
                        ):
                            # Mark as paired
                            hitIndex[model_name][refID]['partner'] = candidate.idx
                            hitIndex[model_name][candidate.idx]['partner'] = refID

                            # Add to paired list
                            paired_dict[model_name].append({refID, candidate.idx})
                            pairs_found += 1

                            logging.debug(
                                f'Paired: {model_name}_{refID} + {model_name}_{candidate.idx}'
                            )
                            break

    logging.debug(f'Found {pairs_found} new symmetric pairs')
    return hitIndex, paired_dict


def checkSymmetricReciprocity(
    model_name: str,
    ref_id: int,
    candidate_id: int,
    hitIndex: Dict[str, Dict[int, Dict[str, Any]]],
    config: Any,
) -> bool:
    """
    Check reciprocal best-match for symmetric pairing with orientation constraints.

    Parameters
    ----------
    model_name : str
        Name of HMM model.
    ref_id : int
        Index of reference hit.
    candidate_id : int
        Index of candidate partner hit.
    hitIndex : dict
        Hit index dictionary with candidate lists.
    config : PairingConfig
        Configuration with strand orientation requirements.

    Returns
    -------
    bool
        True if ref and candidate are reciprocal best unpaired matches, False otherwise.

    Notes
    -----
    Verifies ref_id appears as first valid unpaired candidate in candidate_id's
    candidate list, with both hits having complementary strand roles.
    """
    import logging

    # Check if ref_id is the best unpaired candidate for candidate_id
    for mate_candidate in hitIndex[model_name][candidate_id]['candidates']:
        if (
            mate_candidate.idx in hitIndex[model_name]
            and hitIndex[model_name][mate_candidate.idx]['partner'] is None
        ):
            # Check strand compatibility
            mate_hit = hitIndex[model_name][mate_candidate.idx]['rec']
            candidate_hit = hitIndex[model_name][candidate_id]['rec']

            # Determine if they can form a valid pair
            mate_can_be_left = mate_hit.strand == config.left_strand
            mate_can_be_right = mate_hit.strand == config.right_strand
            candidate_can_be_left = candidate_hit.strand == config.left_strand
            candidate_can_be_right = candidate_hit.strand == config.right_strand

            # Check if this candidate pair is strand-compatible
            strand_compatible = (mate_can_be_left and candidate_can_be_right) or (
                mate_can_be_right and candidate_can_be_left
            )

            if strand_compatible:
                reciprocal = bool(mate_candidate.idx == ref_id)
                logging.debug(
                    f'Reciprocal check: {candidate_id} -> {mate_candidate.idx} == {ref_id}? {reciprocal}'
                )
                return reciprocal

    logging.debug(f'No reciprocal match found for {ref_id} -> {candidate_id}')
    return False


# Update helper functions to include config parameter
def countUnpairedAsymmetric(
    hitIndex: Dict[str, Dict[int, Dict[str, Any]]], config: Any
) -> int:
    """
    Count unpaired hits across both left and right asymmetric models.

    Parameters
    ----------
    hitIndex : dict
        Nested hit index: hitIndex[model][idx] = {rec, partner, candidates}.
    config : PairingConfig
        Configuration specifying left_model and right_model.

    Returns
    -------
    int
        Total number of unpaired hits across both models.
    """
    count = 0
    for model in [config.left_model, config.right_model]:
        if model in hitIndex:
            for hitID in hitIndex[model].keys():
                if hitIndex[model][hitID]['partner'] is None:
                    count += 1  # Count all unpaired hits regardless of strand
    return count


def listunpairedAsymmetric(
    hitIndex: Dict[str, Dict[int, Dict[str, Any]]], config: Any
) -> List[int]:
    """
    List all unpaired hit indices for asymmetric models.

    Parameters
    ----------
    hitIndex : dict
        Nested hit index: hitIndex[model][idx] = {rec, partner, candidates}.
    config : PairingConfig
        Configuration specifying left_model and right_model.

    Returns
    -------
    list of int
        List of hit indices without assigned partners from both models.
    """
    unpaired = []
    for model in [config.left_model, config.right_model]:
        if model in hitIndex:
            for hitID in hitIndex[model].keys():
                if hitIndex[model][hitID]['partner'] is None:
                    unpaired.append(
                        hitID
                    )  # Include all unpaired hits regardless of strand
    return unpaired


def countUnpairedSymmetric(
    hitIndex: Dict[str, Dict[int, Dict[str, Any]]], model_name: str, config: Any
) -> int:
    """
    Count unpaired hits for symmetric model considering orientation constraints.

    Parameters
    ----------
    hitIndex : dict
        Nested hit index: hitIndex[model][idx] = {rec, partner, candidates}.
    model_name : str
        Name of HMM model.
    config : PairingConfig
        Configuration with strand requirements.

    Returns
    -------
    int
        Number of unpaired hits on valid strands for this model.

    Notes
    -----
    Only counts hits whose strand matches either left_strand or right_strand
    in the configuration, as other hits cannot participate in pairing.
    """
    if model_name not in hitIndex:
        return 0

    count = 0
    for hitID in hitIndex[model_name].keys():
        if hitIndex[model_name][hitID]['partner'] is None:
            hit = hitIndex[model_name][hitID]['rec']
            # Only count hits that can participate in pairing
            if hit.strand in [config.left_strand, config.right_strand]:
                count += 1
    return count


def listunpairedSymmetric(
    hitIndex: Dict[str, Dict[int, Dict[str, Any]]], model_name: str, config: Any
) -> List[int]:
    """
    List unpaired hit indices for symmetric model with orientation constraints.

    Parameters
    ----------
    hitIndex : dict
        Nested hit index: hitIndex[model][idx] = {rec, partner, candidates}.
    model_name : str
        Name of HMM model.
    config : PairingConfig
        Configuration with strand requirements.

    Returns
    -------
    list of int
        List of unpaired hit indices on valid strands.

    Notes
    -----
    Only includes hits whose strand can participate in pairing based on
    left_strand and right_strand in configuration.
    """
    if model_name not in hitIndex:
        return []

    unpaired = []
    for hitID in hitIndex[model_name].keys():
        if hitIndex[model_name][hitID]['partner'] is None:
            hit = hitIndex[model_name][hitID]['rec']
            # Only include hits that can participate in pairing
            if hit.strand in [config.left_strand, config.right_strand]:
                unpaired.append(hitID)
    return unpaired
