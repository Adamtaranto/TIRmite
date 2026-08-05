"""Threshold filters applied to a hit table.

Coverage and e-value filtering, applied before pairing so that the pairing
engine only ever sees hits worth pairing.
"""

import glob
import os
from typing import Optional

import pandas as pd


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
        hitTable = hitTable.loc[
            ~(
                (hitTable['model'] == model)
                & (
                    (hitTable['hitEnd'].astype(int) - hitTable['hitStart'].astype(int))
                    + 1
                    < minlen
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
