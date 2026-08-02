"""The terminus pairing engine.

Finds pairs of terminus hits that plausibly bound a single element, using
reciprocal best-match search under a configurable strand orientation and
distance constraint.

Notes
-----
Two generations of the algorithm live here. ``getPairs`` and its helpers are
the original symmetric-only path, kept because ``tirmite legacy`` still uses
it. ``PairingConfig`` and the ``*Asymmetric`` / ``*Symmetric`` functions are
the general engine used by ``tirmite pair``, which handles distinct left and
right models and all four strand orientations.
"""

from collections import namedtuple
import logging
from operator import attrgetter
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd  # type: ignore[import-untyped]

# Library modules acquire a named logger and attach a NullHandler, so that
# importing TIRmite as a library emits nothing until the host application
# configures logging. Handler setup belongs to the CLI, not here.
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


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
    logger.debug('=== ENTERING parseHitsGeneral ===')
    logger.debug(
        f'Config: orientation={config.orientation}, left_strand={config.left_strand}, right_strand={config.right_strand}'
    )
    logger.debug(f'Is asymmetric: {config.is_asymmetric}')

    if not maxDist:
        maxDist_value: float = float('inf')
        logger.debug('Using infinite maxDist')
    else:
        maxDist_value = float(maxDist)
        logger.debug(f'Using maxDist: {maxDist_value}')

    model_pairs = config.get_model_pairs()
    logger.debug(f'Model pairs to process: {model_pairs}')

    for left_model, right_model in model_pairs:
        logger.debug(f'Processing pair: {left_model} -> {right_model}')

        if left_model == right_model:
            # Symmetric pairing - enhanced logic for custom orientations
            logger.debug(f'=== SYMMETRIC PAIRING for {left_model} ===')

            if left_model in hitIndex:
                logger.debug(
                    f'Found {len(hitIndex[left_model])} hits for model {left_model}'
                )

                for UID in hitIndex[left_model].keys():
                    ref = hitIndex[left_model][UID]['rec']
                    logger.debug(
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
                        logger.debug(
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
                        logger.debug(
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
                        logger.debug(
                            f'Hit {UID} on strand {ref.strand} does not match required orientations ({config.left_strand}, {config.right_strand})'
                        )
        else:
            # FIXED: Asymmetric pairing with strand combination handling
            logger.debug(f'=== ASYMMETRIC PAIRING: {left_model} + {right_model} ===')

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

            logger.debug(f'Processing strand combinations: {strand_combinations}')

            for left_strand, right_strand in strand_combinations:
                logger.debug(
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

                            logger.debug(
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

                            logger.debug(
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

    logger.debug('=== COMPLETED parseHitsGeneral ===')
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

    logger.debug('=== DISTANCE CHECK DEBUG ===')
    logger.debug(f'Direction: {direction}')
    logger.debug(
        f'Ref hit: {ref_hit.model} strand={ref_hit.strand} coords=({ref_hit.hitStart}, {ref_hit.hitEnd})'
    )
    logger.debug(
        f'Candidate: {candidate.model} strand={candidate.strand} coords=({candidate.hitStart}, {candidate.hitEnd})'
    )
    logger.debug(f'Calculated distance: {distance}')

    # Check for negative distances (invalid pairing)
    if distance < 0:
        logger.debug(
            f'Negative distance ({distance}) between {ref_hit.model} and {candidate.model} '
            f'on {ref_hit.target}. Ref: {ref_hit.strand}:{ref_hit.hitStart}-{ref_hit.hitEnd}, '
            f'Candidate: {candidate.strand}:{candidate.hitStart}-{candidate.hitEnd}'
        )
        return False

    # Check if within max distance
    valid = distance >= 0 and distance <= maxDist
    logger.debug(
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

    logger.debug('=== _find_candidates DEBUG ===')
    logger.debug(
        f'Ref hit: {ref_hit.model} {ref_hit.strand}:{ref_hit.hitStart}-{ref_hit.hitEnd}'
    )
    logger.debug(
        f'Looking for target_model: {target_model}, target_strand: {target_strand}'
    )
    logger.debug(f'Direction: {direction}, maxDist: {maxDist}')

    if target_model not in hitsDict or ref_hit.target not in hitsDict[target_model]:
        logger.debug(
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
            logger.debug(
                f'Checking candidate: {candidate.model} {candidate.strand}:{candidate.hitStart}-{candidate.hitEnd}'
            )

            # Calculate distance based on direction and orientation
            valid_distance = _check_distance(ref_hit, candidate, direction, maxDist)

            if valid_distance:
                hitIndex[model_key][uid_key]['candidates'].append(candidate)
                candidates_found += 1
                logger.debug(
                    f'Added valid candidate: {candidate.model}_{candidate.idx}'
                )

    logger.debug(
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

        logger.debug(
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

    logger.debug('=== ENTERING iterateGetPairsAsymmetric ===')
    logger.debug(
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

    logger.debug(f'Initial unpaired count: {countUP}')

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

        logger.debug(f'Iteration {reps + 1}: unpaired count {lastCountUP} -> {countUP}')

        # If no change in unpaired hit count, iterate stable rep counter
        if lastCountUP == countUP:
            reps += 1

    # Get IDs of remaining unpaired hits
    unpaired = listunpairedAsymmetric(hitIndex, config)

    total_pairs = sum(len(pairs) for pairs in paired.values())
    logger.debug(
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
                logger.debug(
                    f'Processing left hit {leftID}: {left_hit.strand}:{left_hit.hitStart}-{left_hit.hitEnd}'
                )

                # Look through candidates (which should be from right model)
                for candidate in hitIndex[config.left_model][leftID]['candidates']:
                    logger.debug(
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

                            logger.debug(
                                f'Paired: {config.left_model}_{leftID} + {config.right_model}_{candidate.idx}'
                            )
                            break

    logger.debug(f'Found {pairs_found} new asymmetric pairs')
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
                logger.debug(
                    f'Reciprocal match: {left_model}_{left_id} <-> {right_model}_{right_id}'
                )
                return True  # Reciprocal match found
            else:
                logger.debug(
                    f'Better candidate exists for {right_model}_{right_id}: {candidate.idx}'
                )
                return False  # A better candidate exists

    logger.debug(
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

    logger.debug('=== ENTERING iterateGetPairsCustom ===')
    logger.debug(
        f'Config: {config.left_model} orientation {config.left_strand},{config.right_strand}'
    )

    model_name = config.left_model

    if model_name not in hitIndex:
        logger.error(f'Model {model_name} not found in hitIndex')
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

    logger.debug(f'Initial unpaired count: {countUP}')

    # Iterate pairing procedure
    while countUP > 0 and reps < stableReps:
        hitIndex, paired = getPairsSymmetric(
            hitIndex=hitIndex, model_name=model_name, config=config, paired=paired
        )

        lastCountUP = countUP
        countUP = countUnpairedSymmetric(hitIndex, model_name, config)

        logger.debug(f'Iteration {reps + 1}: unpaired count {lastCountUP} -> {countUP}')

        if lastCountUP == countUP:
            reps += 1

    # Get unpaired list
    unpaired = listunpairedSymmetric(hitIndex, model_name, config)

    total_pairs = len(paired[model_name])
    logger.debug(
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
                logger.debug(
                    f"Hit {refID} on strand {ref_hit.strand} doesn't match orientation {config.left_strand},{config.right_strand}"
                )
                continue

            logger.debug(
                f'Processing hit {refID}: {ref_hit.strand}:{ref_hit.hitStart}-{ref_hit.hitEnd} (can_be_left: {can_be_left}, can_be_right: {can_be_right})'
            )

            # Check candidates for this hit
            for candidate in hitIndex[model_name][refID]['candidates']:
                logger.debug(
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
                        logger.debug(
                            f'Compatible: {refID} (left) + {candidate.idx} (right)'
                        )
                    elif can_be_right and candidate_can_be_left:
                        compatible = True
                        logger.debug(
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

                            logger.debug(
                                f'Paired: {model_name}_{refID} + {model_name}_{candidate.idx}'
                            )
                            break

    logger.debug(f'Found {pairs_found} new symmetric pairs')
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
                logger.debug(
                    f'Reciprocal check: {candidate_id} -> {mate_candidate.idx} == {ref_id}? {reciprocal}'
                )
                return reciprocal

    logger.debug(f'No reciprocal match found for {ref_id} -> {candidate_id}')
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
