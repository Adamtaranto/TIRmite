"""Backwards-compatible façade over :mod:`tirmite.core`.

This module used to hold every core algorithm in a single 5,200-line file. Its
contents now live in focused modules under :mod:`tirmite.core`:

===========================  ================================================
Module                       Contents
===========================  ================================================
:mod:`tirmite.core.parsers`      nhmmer / BLAST / BED readers, format detection
:mod:`tirmite.core.filters`      length and e-value thresholds
:mod:`tirmite.core.termini`      shared left/right terminus geometry
:mod:`tirmite.core.extraction`   TIR and element sequence extraction
:mod:`tirmite.core.flanks`       external flank extraction
:mod:`tirmite.core.tsd`          target site reconstruction and comparison
:mod:`tirmite.core.pairing`      the terminus pairing engine
:mod:`tirmite.core.output`       GFF3 and unpaired-hit output
===========================  ================================================

Every public name is re-exported here so that both established import styles
keep working::

    from tirmite.tirmitetools import extractTIRs
    import tirmite.tirmitetools as tirmite; tirmite.getPairs(...)

New code should import from the specific :mod:`tirmite.core` module instead.
Names are listed explicitly rather than star-imported, so that ``__all__`` and
the type checker both see the real surface.
"""

from tirmite.core.extraction import (  # noqa: F401
    _model_extension,
    extractTIRs,
    fetch_padded_hit,
    fetchElements,
    writeElements,
    writePairedTIRs,
    writeTIRs,
)
from tirmite.core.filters import (  # noqa: F401
    filterHitsEval,
    filterHitsLen,
)
from tirmite.core.flanks import (  # noqa: F401
    FlankResult,
    compute_flank_coordinates,
    compute_inner_tsd_coordinates,
    extract_terminus_flank,
    writeFlanks,
)
from tirmite.core.output import (  # noqa: F401
    fetchUnpaired,
    gffWrite,
)
from tirmite.core.pairing import (  # noqa: F401
    VALID_ORIENTATION_CODES,
    PairingConfig,
    _check_distance,
    _find_candidates,
    candidate_separation,
    checkAsymmetricReciprocity,
    checkSymmetricReciprocity,
    countUnpaired,
    countUnpairedAsymmetric,
    countUnpairedSymmetric,
    getPairs,
    getPairsAsymmetric,
    getPairsSymmetric,
    inter_hit_distance,
    isfirstUnpaired,
    iterateGetPairs,
    iterateGetPairsAsymmetric,
    iterateGetPairsCustom,
    listunpaired,
    listunpairedAsymmetric,
    listunpairedSymmetric,
    parse_orientation,
    parseHits,
    parseHitsGeneral,
    table2dict,
)
from tirmite.core.parsers import (  # noqa: F401
    convertAlign,
    detect_input_format,
    import_BED,
    import_blast,
    import_nhmmer,
)
from tirmite.core.termini import (  # noqa: F401
    TerminusAssignment,
    _determine_terminus_type,
    _model_deficit,
    _pair_roles,
    flipTIRs,
    resolve_terminus,
)
from tirmite.core.tsd import (  # noqa: F401
    compare_tsds,
    format_interleaved_flanks,
    hamming_distance,
    load_tsd_length_map,
    reconstruct_target_site,
    writeTargetSites,
)

__all__ = [
    'candidate_separation',
    'FlankResult',
    'PairingConfig',
    'TerminusAssignment',
    'VALID_ORIENTATION_CODES',
    '_check_distance',
    '_determine_terminus_type',
    '_find_candidates',
    '_model_deficit',
    '_model_extension',
    '_pair_roles',
    'checkAsymmetricReciprocity',
    'checkSymmetricReciprocity',
    'compare_tsds',
    'compute_flank_coordinates',
    'compute_inner_tsd_coordinates',
    'convertAlign',
    'countUnpaired',
    'countUnpairedAsymmetric',
    'countUnpairedSymmetric',
    'detect_input_format',
    'extractTIRs',
    'extract_terminus_flank',
    'fetchElements',
    'fetchUnpaired',
    'fetch_padded_hit',
    'filterHitsEval',
    'filterHitsLen',
    'flipTIRs',
    'format_interleaved_flanks',
    'getPairs',
    'getPairsAsymmetric',
    'getPairsSymmetric',
    'gffWrite',
    'hamming_distance',
    'import_BED',
    'import_blast',
    'import_nhmmer',
    'inter_hit_distance',
    'isfirstUnpaired',
    'iterateGetPairs',
    'iterateGetPairsAsymmetric',
    'iterateGetPairsCustom',
    'listunpaired',
    'listunpairedAsymmetric',
    'listunpairedSymmetric',
    'load_tsd_length_map',
    'parseHits',
    'parseHitsGeneral',
    'parse_orientation',
    'reconstruct_target_site',
    'resolve_terminus',
    'table2dict',
    'writeElements',
    'writeFlanks',
    'writePairedTIRs',
    'writeTIRs',
    'writeTargetSites',
]
