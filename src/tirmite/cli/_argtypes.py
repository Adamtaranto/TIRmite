"""Shared ``argparse`` type validators for the TIRmite subcommands.

These are passed as ``type=`` callables to :meth:`argparse.ArgumentParser.
add_argument`.  argparse calls them with the raw string from the command line
and expects either the converted value or an
:class:`argparse.ArgumentTypeError`, which it renders as a usage error rather
than a traceback.

Every validator here previously existed as two independent copies, one in
:mod:`tirmite.cli.hmm_build` and one in :mod:`tirmite.cli.ensemble_search`.
The copies applied identical bounds but worded their errors differently, so
the same bad input produced different messages depending on the subcommand.
Messages are now uniform: ``<option> must be <constraint>, got <value>``.
"""

import argparse


def validate_evalue(value: str) -> float:
    """
    Validate an e-value threshold.

    Parameters
    ----------
    value : str
        Raw argument string from argparse.

    Returns
    -------
    float
        The parsed e-value.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is not a strictly positive number.

    Notes
    -----
    Zero is rejected as well as negatives: an e-value of 0 would reject every
    hit, which is never what a user means and is a common typo for a very
    small threshold such as ``1e-10``.
    """
    try:
        fvalue = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"e-value must be a number, got '{value}'"
        ) from None

    if fvalue <= 0:
        raise argparse.ArgumentTypeError(f'e-value must be positive, got {fvalue}')
    return fvalue


def validate_identity(value: str) -> float:
    """
    Validate a percent-identity threshold.

    Parameters
    ----------
    value : str
        Raw argument string from argparse.

    Returns
    -------
    float
        The parsed identity, as a percentage in the range 0-100.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is not a number between 0 and 100 inclusive.

    See Also
    --------
    validate_coverage : Takes a 0-1 *fraction*, not a percentage.
    """
    try:
        fvalue = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"identity must be a number, got '{value}'"
        ) from None

    if not 0.0 <= fvalue <= 100.0:
        raise argparse.ArgumentTypeError(
            f'identity must be between 0 and 100, got {fvalue}'
        )
    return fvalue


def validate_coverage(value: str) -> float:
    """
    Validate a coverage threshold expressed as a fraction.

    Parameters
    ----------
    value : str
        Raw argument string from argparse.

    Returns
    -------
    float
        The parsed coverage, as a fraction in the range 0-1.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is not a number between 0.0 and 1.0 inclusive.

    Notes
    -----
    Coverage is a fraction here while identity is a percentage. That asymmetry
    is inherited from the original CLI and is preserved for compatibility, so
    ``--min-coverage 0.8`` and ``--min-identity 80`` express the same
    stringency.
    """
    try:
        fvalue = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"coverage must be a number, got '{value}'"
        ) from None

    if not 0.0 <= fvalue <= 1.0:
        raise argparse.ArgumentTypeError(
            f'coverage must be between 0.0 and 1.0, got {fvalue}'
        )
    return fvalue


def validate_threads(value: str) -> int:
    """
    Validate a thread count.

    Parameters
    ----------
    value : str
        Raw argument string from argparse.

    Returns
    -------
    int
        The parsed thread count.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is not an integer of at least 1.
    """
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"threads must be an integer, got '{value}'"
        ) from None

    if ivalue < 1:
        raise argparse.ArgumentTypeError(f'threads must be at least 1, got {ivalue}')
    return ivalue


def validate_word_size(value: str) -> int:
    """
    Validate a BLAST word size.

    Parameters
    ----------
    value : str
        Raw argument string from argparse.

    Returns
    -------
    int
        The parsed word size.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is not an integer of at least 4.

    Notes
    -----
    The lower bound of 4 is imposed by ``blastn`` itself, which refuses any
    smaller word size for a nucleotide search.
    """
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"word size must be an integer, got '{value}'"
        ) from None

    if ivalue < 4:
        raise argparse.ArgumentTypeError(f'word size must be at least 4, got {ivalue}')
    return ivalue
