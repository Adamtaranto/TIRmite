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
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from Bio import AlignIO, Seq, SeqIO  # type: ignore[import-not-found]
from Bio.SeqRecord import SeqRecord  # type: ignore[import-not-found]
import pandas as pd  # type: ignore[import-untyped]

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


def extract_from_blastdb(
    blastdb: str, seqid: str, start: int, end: int, strand: str = '+'
) -> Optional[str]:
    """
    Extract sequence from BLAST database using blastdbcmd.

    Parameters
    ----------
    blastdb : str
        Path to BLAST database (without file extension).
    seqid : str
        Sequence identifier to extract from.
    start : int
        Start position (1-based, inclusive).
    end : int
        End position (1-based, inclusive).
    strand : str, default '+'
        Strand to extract: '+' for forward, '-' for reverse complement.

    Returns
    -------
    str or None
        Extracted sequence string, or None if extraction failed.

    Notes
    -----
    Uses blastdbcmd with -range parameter for sequence extraction.
    Coordinates are 1-based as expected by blastdbcmd.
    """
    import subprocess

    try:
        cmd = [
            'blastdbcmd',
            '-db',
            blastdb,
            '-entry',
            seqid,
            '-range',
            f'{start}-{end}',
        ]

        if strand == '-':
            cmd.append('-strand')
            cmd.append('minus')

        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, timeout=30
        )

        # Parse FASTA output - validate format and skip header line
        lines = result.stdout.strip().split('\n')
        if len(lines) < 2 or not lines[0].startswith('>'):
            logging.error(
                f'Invalid FASTA output from blastdbcmd for {seqid}:{start}-{end}'
            )
            return None

        # Concatenate sequence lines (skip header at index 0)
        seq = ''.join(lines[1:])
        return seq

    except subprocess.CalledProcessError as e:
        logging.error(f'blastdbcmd failed: {e.stderr}')
        return None
    except subprocess.TimeoutExpired:
        logging.error(f'blastdbcmd timed out for {seqid}:{start}-{end}')
        return None
    except FileNotFoundError:
        logging.error('blastdbcmd command not found. Is BLAST+ installed?')
        return None


def extractTIRs_blastdb(
    model: Optional[str] = None,
    hitTable: Optional[pd.DataFrame] = None,
    maxeval: float = 0.001,
    blastdb: Optional[str] = None,
    padlen: Optional[int] = None,
) -> Tuple[List[SeqRecord], int]:
    """
    Extract TIR sequences from BLAST database for hits of a specific model.

    Parameters
    ----------
    model : str
        Name of model to extract hits for.
    hitTable : pandas.DataFrame
        DataFrame containing all hits with columns: model, target, hitStart, hitEnd,
        strand, evalue, hmmStart, hmmEnd.
    maxeval : float, default 0.001
        Maximum e-value threshold for extracting hits.
    blastdb : str
        Path to BLAST database (without extension).
    padlen : int, optional
        Number of flanking bases to extract on each side of hit (shown in lowercase).

    Returns
    -------
    seqList : list of Bio.SeqRecord.SeqRecord
        List of SeqRecord objects containing extracted TIR sequences.
    hitcount : int
        Number of hits extracted that passed e-value threshold.

    Notes
    -----
    Uses blastdbcmd for sequence extraction from BLAST database.
    Coordinates are adjusted for blastdbcmd (1-based, inclusive).
    Padded flanking sequence is shown in lowercase letters.
    """
    assert model is not None, 'model cannot be None'
    assert hitTable is not None, 'hitTable cannot be None'
    assert blastdb is not None, 'blastdb cannot be None'

    hitcount = 0
    seqList = []

    for index, row in hitTable[hitTable['model'] == model].iterrows():
        if float(row['evalue']) <= maxeval:
            hitcount += 1

            # blastdbcmd uses 1-based coordinates
            start = int(row['hitStart'])
            end = int(row['hitEnd'])

            # Extract sequence (possibly with padding)
            if padlen:
                # Extract with padding
                pad_start = max(1, start - padlen)
                pad_end = end + padlen

                # Extract the full padded region
                full_seq = extract_from_blastdb(
                    blastdb, row['target'], pad_start, pad_end, row['strand']
                )

                if full_seq is None:
                    logging.warning(
                        f'Failed to extract sequence for {model}_{index}, skipping'
                    )
                    continue

                # Calculate positions within extracted sequence to apply case formatting
                # blastdbcmd uses 1-based inclusive coordinates
                # Convert to 0-based for Python string slicing
                hit_start_in_seq = start - pad_start
                hit_end_in_seq = (
                    end - pad_start + 1
                )  # +1 for inclusive end in blastdbcmd

                # Convert to lowercase/uppercase
                hit_seq_str = (
                    full_seq[:hit_start_in_seq].lower()
                    + full_seq[hit_start_in_seq:hit_end_in_seq]
                    + full_seq[hit_end_in_seq:].lower()
                )
            else:
                # Extract without padding
                hit_seq_str = extract_from_blastdb(
                    blastdb, row['target'], start, end, row['strand']
                )

                if hit_seq_str is None:
                    logging.warning(
                        f'Failed to extract sequence for {model}_{index}, skipping'
                    )
                    continue

            # Create SeqRecord
            hitrecord = SeqRecord(Seq.Seq(hit_seq_str))
            hitrecord.id = model + '_' + str(index)
            hitrecord.name = hitrecord.id

            # Build description
            hitrecord.description = '_'.join(
                [
                    '[' + str(row['target']) + ':' + str(row['strand']),
                    str(row['hitStart']),
                    str(row['hitEnd']) + ' modelAlignment:' + str(row['hmmStart']),
                    str(row['hmmEnd']) + ' E-value:' + str(row['evalue']) + ']',
                ]
            )

            seqList.append(hitrecord)

    return seqList, hitcount


# Fix: Do not load fasta into genome!
def extractTIRs(
    model: Optional[str] = None,
    hitTable: Optional[pd.DataFrame] = None,
    maxeval: float = 0.001,
    genome: Any = None,
    padlen: Optional[int] = None,
    genome_descriptions: Optional[Dict[str, str]] = None,
) -> Tuple[List[SeqRecord], int]:
    """
    Extract TIR sequences from genome for hits of a specific model.

    Parameters
    ----------
    model : str
        Name of HMM model to extract hits for.
    hitTable : pandas.DataFrame
        DataFrame containing all hits with columns: model, target, hitStart, hitEnd,
        strand, evalue, hmmStart, hmmEnd.
    maxeval : float, default 0.001
        Maximum e-value threshold for extracting hits.
    genome : pyfaidx.Fasta
        Indexed genome object for sequence extraction.
    padlen : int, optional
        Number of flanking bases to extract on each side of hit (shown in lowercase).
    genome_descriptions : dict, optional
        Dictionary mapping sequence IDs to their descriptions.

    Returns
    -------
    seqList : list of Bio.SeqRecord.SeqRecord
        List of SeqRecord objects containing extracted TIR sequences.
    hitcount : int
        Number of hits extracted that passed e-value threshold.

    Notes
    -----
    Reverse complement is applied for hits on the negative strand.
    Padded flanking sequence is shown in lowercase letters.
    Uses 0-based indexing for pyfaidx extraction.
    """
    assert model is not None, 'model cannot be None'
    assert hitTable is not None, 'hitTable cannot be None'
    assert genome is not None, 'genome cannot be None'
    hitcount = 0
    seqList = []

    for index, row in hitTable[hitTable['model'] == model].iterrows():
        if float(row['evalue']) <= maxeval:
            hitcount += 1

            # Extract sequence using pyfaidx (0-based indexing)
            chrom = genome[row['target']]
            start = int(row['hitStart']) - 1  # Convert to 0-based
            end = int(row['hitEnd'])  # End is exclusive in slicing

            if padlen:
                # Extract with padding
                pad_start = max(0, start - padlen)
                pad_end = min(len(chrom), end + padlen)

                # Build sequence with padding in lowercase
                seq_parts = []
                if start > pad_start:
                    seq_parts.append(str(chrom[pad_start:start]).lower())
                seq_parts.append(str(chrom[start:end]))
                if end < pad_end:
                    seq_parts.append(str(chrom[end:pad_end]).lower())

                hit_seq_str = ''.join(seq_parts)
            else:
                hit_seq_str = str(chrom[start:end])

            # Create SeqRecord
            hitrecord = SeqRecord(Seq.Seq(hit_seq_str))
            hitrecord.id = model + '_' + str(index)

            if row['strand'] == '-':
                hitrecord = hitrecord.reverse_complement(id=hitrecord.id + '_rc')

            hitrecord.name = hitrecord.id

            # Build description with genome description
            coord_info = '_'.join(
                [
                    '[' + str(row['target']) + ':' + str(row['strand']),
                    str(row['hitStart']),
                    str(row['hitEnd']) + ' modelAlignment:' + str(row['hmmStart']),
                    str(row['hmmEnd']) + ' E-value:' + str(row['evalue']) + ']',
                ]
            )

            # Add genome description if available
            if genome_descriptions and row['target'] in genome_descriptions:
                genome_desc = genome_descriptions[row['target']]
                hitrecord.description = f'{coord_info} {genome_desc}'
            else:
                hitrecord.description = coord_info

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

    Returns
    -------
    None
        Writes FASTA files to disk but returns nothing.

    Notes
    -----
    Creates one FASTA file per model with filename format:
    {prefix}{model}_hits_{count}.fasta

    Either genome or blastdb must be provided for sequence extraction.
    """
    assert hitTable is not None, 'hitTable cannot be None'
    assert genome is not None or blastdb is not None, (
        'Either genome or blastdb must be provided'
    )

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

    for model in hitTable['model'].unique():
        # List of TIR seqrecords, and count of hits
        if blastdb:
            seqList, hitcount = extractTIRs_blastdb(
                model=model,
                hitTable=hitTable,
                maxeval=maxeval,
                blastdb=blastdb,
                padlen=padlen,
            )
        else:
            seqList, hitcount = extractTIRs(
                model=model,
                hitTable=hitTable,
                maxeval=maxeval,
                genome=genome,
                padlen=padlen,
                genome_descriptions=genome_descriptions,
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
    assert genome is not None or blastdb is not None, (
        'Either genome or blastdb must be provided'
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

            for pair in paired[model]:
                model_counter += 1
                # Convert set to list for indexing
                hit_ids = list(pair)
                x_id, y_id = hit_ids[0], hit_ids[1]

                # Get hit records using helper function
                x = get_hit_record(x_id)
                y = get_hit_record(y_id)

                leftHit, rightHit = flipTIRs(x, y)

                # Create element ID with counter only (avoids model name duplication)
                eleID = f'Element_{model_counter}'

                # Extract element sequence
                if blastdb:
                    # Extract from BLAST database
                    ele_seq_str = extract_from_blastdb(
                        blastdb,
                        leftHit.target,
                        int(leftHit.hitStart),
                        int(rightHit.hitEnd),
                        leftHit.strand,
                    )
                    if ele_seq_str is None:
                        logging.warning(f'Failed to extract element {eleID}, skipping')
                        continue
                else:
                    # Extract using pyfaidx (0-based indexing)
                    chrom = genome[leftHit.target]
                    start = int(leftHit.hitStart) - 1  # Convert to 0-based
                    end = int(rightHit.hitEnd)  # End is exclusive in slicing
                    ele_seq_str = str(chrom[start:end])

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
                if genome_descriptions and leftHit.target in genome_descriptions:
                    genome_desc = genome_descriptions[leftHit.target]
                    eleSeq.description = f'{coord_info} {genome_desc}'
                else:
                    eleSeq.description = coord_info

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
    Output filename format: {prefix}{model}_elements.fasta
    """
    assert eleDict is not None, 'eleDict cannot be None'
    if prefix:
        prefix = cleanID(prefix) + '_'
    else:
        prefix = ''

    for model in eleDict.keys():
        if len(eleDict[model]) > 0:  # Only write files for models with actual elements
            outfile = os.path.join(outDir, prefix + model + '_elements.fasta')
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

    Returns
    -------
    None
        Writes FASTA files to disk but returns nothing.

    Notes
    -----
    Right TIRs are reverse complemented. Outputs two sequences per pair:
    {model}_{counter}_L (left TIR) and {model}_{counter}_R (right TIR).
    Filename format: {prefix}{model}_paired_term_hits_{count}.fasta

    Either genome or blastdb must be provided for sequence extraction.
    """
    assert outDir is not None, 'outDir cannot be None'
    assert paired is not None, 'paired cannot be None'
    assert hitIndex is not None, 'hitIndex cannot be None'
    assert genome is not None or blastdb is not None, (
        'Either genome or blastdb must be provided'
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
            seqList = []  # Just collect sequences for FASTA output

            for pair in paired[model]:
                model_counter += 1
                # Convert set to list for indexing
                hit_ids = list(pair)
                x_id, y_id = hit_ids[0], hit_ids[1]

                # Get hit records using helper function
                x = get_hit_record(x_id)
                y = get_hit_record(y_id)

                leftHit, rightHit = flipTIRs(x, y)
                eleID = f'{model}_{model_counter}'

                # Extract TIR sequences
                if blastdb:
                    # Extract from BLAST database
                    if padlen:
                        # Extract left TIR with padding
                        left_start = max(1, int(leftHit.hitStart) - padlen)
                        left_end = int(leftHit.hitEnd) + padlen
                        left_full = extract_from_blastdb(
                            blastdb,
                            leftHit.target,
                            left_start,
                            left_end,
                            leftHit.strand,
                        )
                        if left_full is None:
                            logging.warning(f'Failed to extract left TIR for {eleID}')
                            continue

                        # Apply case formatting
                        # blastdbcmd uses 1-based inclusive coordinates
                        hit_start_in_seq = int(leftHit.hitStart) - left_start
                        hit_end_in_seq = (
                            int(leftHit.hitEnd) - left_start + 1
                        )  # +1 for inclusive end
                        left_seq_str = (
                            left_full[:hit_start_in_seq].lower()
                            + left_full[hit_start_in_seq:hit_end_in_seq]
                            + left_full[hit_end_in_seq:].lower()
                        )

                        # Extract right TIR with padding
                        right_start = max(1, int(rightHit.hitStart) - padlen)
                        right_end = int(rightHit.hitEnd) + padlen
                        right_full = extract_from_blastdb(
                            blastdb, rightHit.target, right_start, right_end, '+'
                        )
                        if right_full is None:
                            logging.warning(f'Failed to extract right TIR for {eleID}')
                            continue

                        # Apply case formatting (before reverse complement)
                        # blastdbcmd uses 1-based inclusive coordinates
                        hit_start_in_seq = int(rightHit.hitStart) - right_start
                        hit_end_in_seq = (
                            int(rightHit.hitEnd) - right_start + 1
                        )  # +1 for inclusive end
                        right_seq_str = (
                            right_full[:hit_start_in_seq].lower()
                            + right_full[hit_start_in_seq:hit_end_in_seq]
                            + right_full[hit_end_in_seq:].lower()
                        )
                    else:
                        # Extract without padding
                        left_seq_str = extract_from_blastdb(
                            blastdb,
                            leftHit.target,
                            int(leftHit.hitStart),
                            int(leftHit.hitEnd),
                            leftHit.strand,
                        )
                        right_seq_str = extract_from_blastdb(
                            blastdb,
                            rightHit.target,
                            int(rightHit.hitStart),
                            int(rightHit.hitEnd),
                            '+',  # Don't reverse complement yet
                        )

                        if left_seq_str is None or right_seq_str is None:
                            logging.warning(f'Failed to extract TIRs for {eleID}')
                            continue
                else:
                    # Extract using pyfaidx
                    chrom = genome[leftHit.target]
                    left_start = int(leftHit.hitStart) - 1  # Convert to 0-based
                    left_end = int(leftHit.hitEnd)

                    # Extract right TIR sequence using pyfaidx
                    right_start = int(rightHit.hitStart) - 1  # Convert to 0-based
                    right_end = int(rightHit.hitEnd)

                    if padlen:
                        # Extract with padding for left TIR
                        left_pad_start = max(0, left_start - padlen)
                        left_pad_end = min(len(chrom), left_end + padlen)

                        left_seq_parts = []
                        if left_start > left_pad_start:
                            left_seq_parts.append(
                                str(chrom[left_pad_start:left_start]).lower()
                            )
                        left_seq_parts.append(str(chrom[left_start:left_end]))
                        if left_end < left_pad_end:
                            left_seq_parts.append(
                                str(chrom[left_end:left_pad_end]).lower()
                            )

                        left_seq_str = ''.join(left_seq_parts)

                        # Extract with padding for right TIR
                        right_pad_start = max(0, right_start - padlen)
                        right_pad_end = min(len(chrom), right_end + padlen)

                        right_seq_parts = []
                        if right_start > right_pad_start:
                            right_seq_parts.append(
                                str(chrom[right_pad_start:right_start]).lower()
                            )
                        right_seq_parts.append(str(chrom[right_start:right_end]))
                        if right_end < right_pad_end:
                            right_seq_parts.append(
                                str(chrom[right_end:right_pad_end]).lower()
                            )

                        right_seq_str = ''.join(right_seq_parts)
                    else:
                        left_seq_str = str(chrom[left_start:left_end])
                        right_seq_str = str(chrom[right_start:right_end])

                # Create SeqRecords for FASTA output only
                eleSeqLeft = SeqRecord(Seq.Seq(left_seq_str))
                eleSeqRight = SeqRecord(Seq.Seq(right_seq_str))
                eleSeqRight = eleSeqRight.reverse_complement()

                eleSeqLeft.id = eleID + '_L'
                eleSeqLeft.name = eleID + '_L'

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

                if genome_descriptions and leftHit.target in genome_descriptions:
                    genome_desc = genome_descriptions[leftHit.target]
                    eleSeqLeft.description = f'{left_coord} {genome_desc}'
                else:
                    eleSeqLeft.description = left_coord

                eleSeqRight.id = eleID + '_R'
                eleSeqRight.name = eleID + '_R'

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

                if genome_descriptions and leftHit.target in genome_descriptions:
                    genome_desc = genome_descriptions[leftHit.target]
                    eleSeqRight.description = f'{right_coord} {genome_desc}'
                else:
                    eleSeqRight.description = right_coord

                # Add to sequence list for FASTA output
                seqList.append(eleSeqLeft)
                seqList.append(eleSeqRight)

            # Write FASTA file for this model only if we have sequences
            if seqList:
                outfile = os.path.join(
                    outDir,
                    prefix
                    + model
                    + '_paired_term_hits_'
                    + str(len(seqList))
                    + '.fasta',
                )
                with open(outfile, 'w') as handle:
                    for seq in seqList:
                        seq.id = prefix + str(seq.id)
                        SeqIO.write(seq, handle, 'fasta')


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
            offset = hmm_start - 1
            external_pos = hit_start - offset
        else:  # '-'
            # hmmStart aligns to hit_end (higher coord); hmmEnd aligns to hit_start
            offset = model_len - hmm_end
            external_pos = hit_start - offset
        flank_start = external_pos - flank_len
        flank_end = external_pos - 1
    else:
        # External end faces RIGHT (higher genomic coordinates)
        if strand == '+':
            offset = model_len - hmm_end
            external_pos = hit_end + offset
        else:  # '-'
            # hmmStart aligns to hit_end; external end is at hit_end side
            offset = hmm_start - 1
            external_pos = hit_end + offset
        flank_start = external_pos + 1
        flank_end = external_pos + flank_len

    return flank_start, flank_end, offset


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


def writeFlanks(
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
    paired_only: bool = False,
    genome_descriptions: Optional[Dict[str, str]] = None,
    blastdb: Optional[str] = None,
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
    flank_len : int, default 0
        Number of bases to extract in each flanking region.
    flank_max_offset : int, optional
        Maximum allowed offset (uncovered model positions) between hit alignment
        and model end. Hits with offset > this value are skipped.
    paired_only : bool, default False
        If True, only extract flanks for hits that are part of a pair.
        If False, also extract flanks for unpaired hits.
    genome_descriptions : dict, optional
        Dictionary mapping sequence IDs to their descriptions.
    blastdb : str, optional
        Path to BLAST database. Alternative to genome.

    Returns
    -------
    None
        Writes FASTA files to disk.

    Notes
    -----
    Output files are named:
      {prefix}{model}_left_flank_{count}.fasta  – flanks for left terminus hits
      {prefix}{model}_right_flank_{count}.fasta – flanks for right terminus hits

    For asymmetric pairings the left and right model names may differ, so each
    model gets its own pair of files. For symmetric pairings the same model name
    is used for both files (distinguished by the _left_/_right_ suffix).

    Flanking sequences are always reported in the forward (+) genomic strand
    orientation. Coordinates are 1-based.
    """
    assert outDir is not None, 'outDir cannot be None'
    assert hitTable is not None, 'hitTable cannot be None'
    assert model_lengths is not None, 'model_lengths cannot be None'
    assert paired is not None, 'paired cannot be None'
    assert hitIndex is not None, 'hitIndex cannot be None'
    assert genome is not None or blastdb is not None, (
        'Either genome or blastdb must be provided'
    )

    # When config is None (e.g. pairing-map mode with multiple configs) unpaired
    # hits cannot be attributed to a terminus; fall back to paired-only processing.
    if config is None and not paired_only:
        logging.debug(
            'config is None: cannot determine terminus type for unpaired hits; '
            'processing paired hits only.'
        )
        paired_only = True

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

        Returns
        -------
        SeqRecord or None
            The flank SeqRecord, or None if the flank cannot be extracted
            (missing model length, missing HMM coords, offset exceeds max,
            or sequence extraction failure).
        """
        model = hit.model
        model_len = model_lengths.get(model) if model_lengths else None  # type: ignore[union-attr]
        if model_len is None:
            logging.warning(f'Model length not found for {model}, skipping flank')
            return None

        hmm_start, hmm_end = get_hmm_coords(hit.idx)
        if hmm_start is None or hmm_end is None:
            logging.debug(f'HMM coordinates unavailable for hit {hit.idx}, skipping')
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

        # Clamp to valid 1-based range
        flank_start = max(1, flank_start)
        flank_end = max(1, flank_end)

        if flank_start > flank_end:
            logging.debug(f'Empty flank region for hit {hit.idx}, skipping')
            return None

        # Extract sequence (always on + strand for genomic orientation)
        if blastdb:
            seq_str = extract_from_blastdb(
                blastdb, hit.target, flank_start, flank_end, '+'
            )
        else:
            chrom = genome[hit.target]
            # Clamp end to chromosome length
            chrom_len = len(chrom)
            flank_end = min(flank_end, chrom_len)
            if flank_start > flank_end:
                return None
            seq_str = str(chrom[flank_start - 1 : flank_end])

        if seq_str is None:
            logging.warning(
                f'Failed to extract flank sequence for hit {hit.idx}, skipping'
            )
            return None

        record = SeqRecord(Seq.Seq(seq_str))
        side = 'L' if is_left else 'R'
        record.id = f'{record_id}_{side}'
        record.name = record.id

        coord_info = (
            f'[{hit.target}:+ {flank_start}_{flank_end}'
            f' hit:{hit.strand}:{hit.hitStart}_{hit.hitEnd}'
            f' offset:{offset}]'
        )
        if genome_descriptions and hit.target in genome_descriptions:
            record.description = f'{coord_info} {genome_descriptions[hit.target]}'
        else:
            record.description = coord_info

        return record

    # ------------------------------------------------------------------
    # Process paired hits
    # ------------------------------------------------------------------
    paired_hit_ids: Set[int] = set()

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

            left_rec = build_flank_record(leftHit, is_left=True, record_id=pair_id)
            right_rec = build_flank_record(rightHit, is_left=False, record_id=pair_id)

            if left_rec:
                left_flanks.setdefault(leftHit.model, []).append(left_rec)
            if right_rec:
                right_flanks.setdefault(rightHit.model, []).append(right_rec)

            paired_hit_ids.add(leftHit.idx)
            paired_hit_ids.add(rightHit.idx)

    # ------------------------------------------------------------------
    # Process unpaired hits (only when paired_only=False)
    # ------------------------------------------------------------------
    if not paired_only:
        for model in hitIndex.keys():
            for hit_id, hit_data in hitIndex[model].items():
                if hit_data['partner'] is not None:
                    continue  # already handled above
                hit = hit_data['rec']
                terminus_type = _determine_terminus_type(hit, config)
                if terminus_type is None:
                    logging.debug(
                        f'Cannot determine terminus type for unpaired hit {hit.idx} '
                        f'(model={hit.model}, strand={hit.strand}), skipping'
                    )
                    continue

                is_left = terminus_type == 'left'
                record_id = f'{model}_{hit_id}_unpaired'
                rec = build_flank_record(hit, is_left=is_left, record_id=record_id)
                if rec:
                    if is_left:
                        left_flanks.setdefault(hit.model, []).append(rec)
                    else:
                        right_flanks.setdefault(hit.model, []).append(rec)

    # ------------------------------------------------------------------
    # Write output files
    # ------------------------------------------------------------------
    for model, flanks in left_flanks.items():
        if flanks:
            outfile = os.path.join(
                outDir, prefix + model + '_left_flank_' + str(len(flanks)) + '.fasta'
            )
            with open(outfile, 'w') as handle:
                for seq in flanks:
                    seq.id = prefix + str(seq.id)
                    SeqIO.write(seq, handle, 'fasta')

    for model, flanks in right_flanks.items():
        if flanks:
            outfile = os.path.join(
                outDir, prefix + model + '_right_flank_' + str(len(flanks)) + '.fasta'
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
        Format: left_model<TAB>right_model<TAB>tsd_length

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
) -> Tuple[str, str, str, int]:
    """
    Reconstruct a target site by joining left and right flanking sequences.

    When a TSD (Target Site Duplication) is present, it is de-duplicated
    so only one copy appears in the reconstructed target site.

    Parameters
    ----------
    left_flank_seq : str
        External flanking sequence upstream of the left terminus.
    right_flank_seq : str
        External flanking sequence downstream of the right terminus.
    tsd_length : int, default 0
        Length of the terminal duplication feature.
    tsd_in_model : bool, default False
        If True, the TSD is part of the termini model and appears at the
        edges of the hit (last ``tsd_length`` bases of the left flank's
        terminus end and the first ``tsd_length`` bases of the right
        flank's terminus end). If False, the TSD occurs as the first
        ``tsd_length`` bases of the flanking region immediately adjacent
        to the terminus hits.

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
    tsd_hamming : int
        Hamming distance between the left and right TSD sequences
        (0 if ``tsd_length`` is 0).

    Notes
    -----
    When ``tsd_in_model`` is True:
        The TSD is located at the inner boundary of each flank (tail of
        left flank, head of right flank). These are trimmed from one side
        before joining.

    When ``tsd_in_model`` is False:
        The TSD is the first ``tsd_length`` bases of the flanking region
        next to the terminus. For the left flank, it's at the right end;
        for the right flank, it's at the left end. One copy is trimmed.
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
        # Trim TSD from right flank to de-duplicate
        trimmed_right = right_flank_seq[tsd_length:]
        target_site = left_flank_seq + trimmed_right

        if len(left_tsd) == len(right_tsd) and len(left_tsd) > 0:
            tsd_hamming = hamming_distance(left_tsd.upper(), right_tsd.upper())
        if tsd_hamming > 0:
            logging.warning(
                f'TSD mismatch (hamming={tsd_hamming}): '
                f'left={left_tsd} right={right_tsd}'
            )
    else:
        target_site = left_flank_seq + right_flank_seq

    return target_site, left_tsd, right_tsd, tsd_hamming


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

    Returns
    -------
    None
        Writes FASTA files to disk.

    Notes
    -----
    Output files:
      {prefix}target_sites.fasta – reconstructed target sites
      {prefix}interleaved_flanks.fasta – interleaved left/right flanks
    """
    assert outDir is not None, 'outDir cannot be None'
    assert hitTable is not None, 'hitTable cannot be None'
    assert model_lengths is not None, 'model_lengths cannot be None'
    assert paired is not None, 'paired cannot be None'
    assert hitIndex is not None, 'hitIndex cannot be None'
    assert genome is not None or blastdb is not None, (
        'Either genome or blastdb must be provided'
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
        if is_nested:
            for _m, model_hits in hitIndex.items():
                if hit_id in model_hits:
                    return model_hits[hit_id]['rec']
            raise KeyError(f'Hit ID {hit_id} not found in hitIndex')
        return hitIndex[hit_id]['rec']  # type: ignore[index]

    # ------------------------------------------------------------------
    # Helper: extract flank sequence for one hit
    # ------------------------------------------------------------------
    def extract_flank(hit: Any, is_left: bool) -> Optional[str]:
        model = hit.model
        model_len = model_lengths.get(model) if model_lengths else None  # type: ignore[union-attr]
        if model_len is None:
            return None

        hmm_start, hmm_end = get_hmm_coords(hit.idx)
        if hmm_start is None or hmm_end is None:
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
            return None

        flank_start = max(1, flank_start)
        flank_end = max(1, flank_end)

        if flank_start > flank_end:
            return None

        if blastdb:
            seq_str = extract_from_blastdb(
                blastdb, hit.target, flank_start, flank_end, '+'
            )
        else:
            chrom = genome[hit.target]
            chrom_len = len(chrom)
            flank_end = min(flank_end, chrom_len)
            if flank_start > flank_end:
                return None
            seq_str = str(chrom[flank_start - 1 : flank_end])

        return seq_str

    # ------------------------------------------------------------------
    # Helper: resolve TSD length for a pair of models
    # ------------------------------------------------------------------
    def get_tsd_length_for_pair(left_model: str, right_model: str) -> int:
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
    # Process paired hits
    # ------------------------------------------------------------------
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

            left_seq = extract_flank(leftHit, is_left=True)
            right_seq = extract_flank(rightHit, is_left=False)

            if left_seq is None or right_seq is None:
                logging.debug(
                    f'Could not extract flanks for pair {pair_id}, skipping target site'
                )
                continue

            pair_tsd_len = get_tsd_length_for_pair(leftHit.model, rightHit.model)

            target_site, left_tsd, right_tsd, tsd_hamming = reconstruct_target_site(
                left_flank_seq=left_seq,
                right_flank_seq=right_seq,
                tsd_length=pair_tsd_len,
                tsd_in_model=tsd_in_model,
            )

            # Build metadata for FASTA header
            meta_parts = [
                f'flank_len={flank_len}',
                f'tsd_len={pair_tsd_len}',
                f'tsd_in_model={tsd_in_model}',
                f'left_model={leftHit.model}',
                f'right_model={rightHit.model}',
                f'contig={leftHit.target}',
                f'left_flank_hit={leftHit.strand}:{leftHit.hitStart}_{leftHit.hitEnd}',
                f'right_flank_hit={rightHit.strand}:{rightHit.hitStart}_{rightHit.hitEnd}',
                f'tsd_hamming={tsd_hamming}',
            ]
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

    # ------------------------------------------------------------------
    # Write output files
    # ------------------------------------------------------------------
    if target_site_records:
        ts_outfile = os.path.join(outDir, f'{prefix_str}target_sites.fasta')
        with open(ts_outfile, 'w') as handle:
            for rec in target_site_records:
                SeqIO.write(rec, handle, 'fasta')
        logging.info(
            f'Wrote {len(target_site_records)} reconstructed target sites to {ts_outfile}'
        )
    else:
        logging.warning('No target sites could be reconstructed')

    if interleaved_records:
        il_outfile = os.path.join(outDir, f'{prefix_str}interleaved_flanks.fasta')
        with open(il_outfile, 'w') as handle:
            for rec in interleaved_records:
                SeqIO.write(rec, handle, 'fasta')
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
        self.orientation = orientation.split(',')
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

                    if ref.strand == config.left_strand:
                        # This hit matches the left terminus strand requirement
                        # It should look for right terminus partners

                        # Determine search direction based on orientation and strand
                        if config.left_strand == '+' and config.right_strand == '+':
                            # F,F on positive strand: left looks downstream (higher coords)
                            search_direction = 'left_to_right'
                        elif config.left_strand == '+' and config.right_strand == '-':
                            # F,R: left(+) looks for right(-), still downstream
                            search_direction = 'left_to_right'
                        elif config.left_strand == '-' and config.right_strand == '+':
                            # R,F: left(-) looks for right(+), upstream in genomic coords
                            search_direction = 'right_to_left'
                        elif config.left_strand == '-' and config.right_strand == '-':
                            # R,R on negative strand: left(-) looks upstream (lower coords)
                            search_direction = 'right_to_left'

                        logging.debug(
                            f'Hit {UID} acting as LEFT terminus, searching {search_direction} for RIGHT partners on {config.right_strand} strand'
                        )

                        _find_candidates(
                            ref,
                            right_model,
                            config.right_strand,
                            hitsDict,
                            hitIndex,
                            maxDist_value,
                            search_direction,
                        )

                    elif ref.strand == config.right_strand:
                        # This hit matches the right terminus strand requirement
                        # It should look for left terminus partners

                        # Determine search direction (opposite of left terminus search)
                        if config.left_strand == '+' and config.right_strand == '+':
                            # F,F: right(+) looks upstream for left(+)
                            search_direction = 'right_to_left'
                        elif config.left_strand == '+' and config.right_strand == '-':
                            # F,R: right(-) looks upstream for left(+)
                            search_direction = 'right_to_left'
                        elif config.left_strand == '-' and config.right_strand == '+':
                            # R,F: right(+) looks downstream for left(-)
                            search_direction = 'left_to_right'
                        elif config.left_strand == '-' and config.right_strand == '-':
                            # R,R: right(-) looks downstream for left(-)
                            search_direction = 'left_to_right'

                        logging.debug(
                            f'Hit {UID} acting as RIGHT terminus, searching {search_direction} for LEFT partners on {config.left_strand} strand'
                        )

                        _find_candidates(
                            ref,
                            left_model,
                            config.left_strand,
                            hitsDict,
                            hitIndex,
                            maxDist_value,
                            search_direction,
                        )
                    else:
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
    Handles strand-aware distance calculation. For negative strand hits, coordinates
    are inverted in nhmmer output (hitStart > hitEnd). Calculates biological distance
    from terminus 3' end to partner 5' end based on strand orientations.
    """

    def get_terminus_position(hit: Any) -> Dict[str, int]:
        """
        Extract terminus positions based on strand orientation and biological direction.

        Parameters
        ----------
        hit : namedtuple
            Hit record with strand, hitStart, and hitEnd attributes.

        Returns
        -------
        dict
            Dictionary with keys: 'start', 'end', '5_prime', '3_prime' containing
            genomic coordinates adjusted for strand orientation.

        Notes
        -----
        For positive strand (+): hitStart < hitEnd, 5' = hitStart, 3' = hitEnd.
        For negative strand (-): hitStart > hitEnd (flipped), 5' = hitEnd, 3' = hitStart.
        The 5' and 3' positions represent biological directionality.
        """
        if hit.strand == '+':
            # Positive strand: hitStart < hitEnd
            # For left terminus: use 3' end (hitEnd)
            # For right terminus: use 5' end (hitStart)
            return {
                'start': hit.hitStart,
                'end': hit.hitEnd,
                '5_prime': hit.hitStart,
                '3_prime': hit.hitEnd,
            }
        else:
            # Negative strand: hitStart > hitEnd (coordinates flipped)
            # For left terminus: use 3' end (hitStart)
            # For right terminus: use 5' end (hitEnd)
            return {
                'start': hit.hitEnd,
                'end': hit.hitStart,
                '5_prime': hit.hitEnd,
                '3_prime': hit.hitStart,
            }

    # Get terminus positions for both hits
    ref_pos = get_terminus_position(ref_hit)
    cand_pos = get_terminus_position(candidate)

    # Debug logging
    logging.debug('=== DISTANCE CHECK DEBUG ===')
    logging.debug(f'Direction: {direction}')
    logging.debug(
        f'Ref hit: {ref_hit.model} strand={ref_hit.strand} coords=({ref_hit.hitStart}, {ref_hit.hitEnd})'
    )
    logging.debug(
        f'Candidate: {candidate.model} strand={candidate.strand} coords=({candidate.hitStart}, {candidate.hitEnd})'
    )
    logging.debug(f'Ref positions: {ref_pos}')
    logging.debug(f'Candidate positions: {cand_pos}')

    # Calculate distance based on biological relationship
    if direction == 'left_to_right':
        # Left terminus looking for right terminus downstream
        # Distance from left terminus 3' end to right terminus 5' end
        if ref_hit.strand == '+' and candidate.strand == '+':
            # F,F orientation on + strand: left_end -> right_start
            distance = cand_pos['5_prime'] - ref_pos['3_prime']
        elif ref_hit.strand == '+' and candidate.strand == '-':
            # F,R orientation: left_end -> right_start (right is on - strand)
            distance = cand_pos['5_prime'] - ref_pos['3_prime']
        elif ref_hit.strand == '-' and candidate.strand == '+':
            # R,F orientation: left_end -> right_start (left is on - strand)
            distance = cand_pos['5_prime'] - ref_pos['3_prime']
        elif ref_hit.strand == '-' and candidate.strand == '-':
            # R,R orientation on - strand: left_end -> right_start
            distance = cand_pos['5_prime'] - ref_pos['3_prime']
        else:
            logging.error(
                f'Unexpected strand combination: {ref_hit.strand}, {candidate.strand}'
            )
            return False

    else:  # 'right_to_left'
        # Right terminus looking for left terminus upstream
        # Distance from right terminus 5' end to left terminus 3' end (should be negative, so we flip)
        if ref_hit.strand == '+' and candidate.strand == '+':
            # F,F orientation: right_start -> left_end (upstream)
            distance = ref_pos['5_prime'] - cand_pos['3_prime']
        elif ref_hit.strand == '-' and candidate.strand == '+':
            # R,F orientation: right_start -> left_end
            distance = ref_pos['5_prime'] - cand_pos['3_prime']
        elif ref_hit.strand == '+' and candidate.strand == '-':
            # F,R orientation: right_start -> left_end
            distance = ref_pos['5_prime'] - cand_pos['3_prime']
        elif ref_hit.strand == '-' and candidate.strand == '-':
            # R,R orientation: right_start -> left_end
            distance = ref_pos['5_prime'] - cand_pos['3_prime']
        else:
            logging.error(
                f'Unexpected strand combination: {ref_hit.strand}, {candidate.strand}'
            )
            return False

    logging.debug(f'Calculated distance: {distance}')

    # Check for negative distances (invalid pairing)
    if distance < 0:
        logging.warning(
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
    def get_distance_for_sorting(ref: Any, cand: Any, direction: str) -> int:
        """
        Calculate biological distance between hits for sorting candidates.

        Parameters
        ----------
        ref : namedtuple
            Reference hit with strand, hitStart, hitEnd attributes.
        cand : namedtuple
            Candidate partner hit with strand, hitStart, hitEnd attributes.
        direction : str
            Search direction: 'left_to_right' or 'right_to_left'.

        Returns
        -------
        int
            Calculated distance in base pairs from reference to candidate based
            on terminus positions and biological orientation.

        Notes
        -----
        Matches the distance calculation logic in _check_distance function.
        For left_to_right: distance from left terminus 3' end to right terminus 5' end.
        For right_to_left: distance from right terminus 5' end to left terminus 3' end.
        """

        def get_terminus_position(hit: Any) -> Dict[str, int]:
            """
            Extract 5' and 3' terminus positions accounting for strand orientation.

            Parameters
            ----------
            hit : namedtuple
                Hit record with strand, hitStart, and hitEnd attributes.

            Returns
            -------
            dict
                Dictionary with '5_prime' and '3_prime' keys containing genomic coordinates
                adjusted for strand directionality.
            """
            if hit.strand == '+':
                return {'5_prime': hit.hitStart, '3_prime': hit.hitEnd}
            else:
                return {'5_prime': hit.hitEnd, '3_prime': hit.hitStart}

        ref_pos = get_terminus_position(ref)
        cand_pos = get_terminus_position(cand)

        if direction == 'left_to_right':
            # Distance from left terminus 3' end to right terminus 5' end
            return int(cand_pos['5_prime'] - ref_pos['3_prime'])
        else:
            # Distance from right terminus 5' end to left terminus 3' end
            return int(ref_pos['5_prime'] - cand_pos['3_prime'])

    if hitIndex[model_key][uid_key]['candidates']:
        # Sort by calculated distance (closest first)
        hitIndex[model_key][uid_key]['candidates'] = sorted(
            hitIndex[model_key][uid_key]['candidates'],
            key=lambda x: get_distance_for_sorting(ref_hit, x, direction),
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
