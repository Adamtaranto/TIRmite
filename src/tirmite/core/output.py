"""Annotation output: GFF3 and unpaired-hit reporting."""

from collections import namedtuple
from operator import attrgetter
import os
from typing import Any, Dict, List, Optional, Union

from tirmite.utils.utils import cleanID


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
        'gffTup',
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
        featureList = {}
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
