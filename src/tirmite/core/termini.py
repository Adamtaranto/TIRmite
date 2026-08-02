"""Shared terminus geometry: which end of an element a hit represents.

Several parts of TIRmite need to answer the same question -- given a hit and a
pairing configuration, is this the left or the right terminus of the element,
and which model edge faces outward? The helpers here answer it once.

Notes
-----
This module exists to keep the dependency graph acyclic. Sequence extraction,
flank extraction and target-site reconstruction all need these helpers, and
they in turn need nothing from those modules, so ``termini`` sits below all
three.
"""

import logging
from operator import attrgetter
from typing import Any, NamedTuple, Optional, Tuple

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
