#!/usr/bin/env python3
"""
TIRmite command-line interface with subcommands.
"""

import argparse
import logging
import sys
from typing import Dict

from tirmite._version import __version__

logger = logging.getLogger(__name__)


def create_parser() -> argparse.ArgumentParser:
    """
    Create the main argument parser with subcommand structure.

    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser with legacy, seed, pair, and search subcommands.

    Notes
    -----
    Subcommands available:
    - legacy: Original TIRmite workflow (HMM search + pairing)
    - seed: Build HMM models from seed sequences
    - pair: Pair precomputed nhmmer hits
    - search: Ensemble search with hit merging from clustered features
    """
    parser = argparse.ArgumentParser(
        prog='tirmite',
        description='TIRmite: Transposon Terminal Repeat detection suite',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Available subcommands:
  legacy    Original TIRmite workflow (HMM search + pairing)
  seed      Build HMM models from seed sequences
  pair      Pair precomputed nhmmer hits
  search    Ensemble search: merge hits from clustered features
  validate  Validate reconstructed target sites

Examples:
  tirmite legacy --genome genome.fa --hmm-file model.hmm
  tirmite seed --left-seed left.fa --model-name myTE --genome genome.fa
  tirmite pair --genome genome.fa --nhmmer-file hits.out --hmm-file model.hmm
  tirmite search --blast-results hits.tab --cluster-map clusters.tsv
  tirmite validate --target-sites targets.fa --blastdb validation_db
        """,
    )

    parser.add_argument('--version', action='version', version=f'tirmite {__version__}')

    # Create subparsers
    subparsers = parser.add_subparsers(
        dest='command', help='Available subcommands', metavar='COMMAND'
    )

    # Import and add subcommands
    from tirmite.cli.ensemble_search import add_search_parser
    from tirmite.cli.hmm_build import add_seed_parser
    from tirmite.cli.hmm_pair import add_pair_parser
    from tirmite.cli.legacy import add_legacy_parser
    from tirmite.cli.validate import add_validate_parser

    add_legacy_parser(subparsers)
    add_seed_parser(subparsers)
    add_pair_parser(subparsers)
    add_search_parser(subparsers)
    add_validate_parser(subparsers)

    return parser


def _subcommand_parsers(
    parser: argparse.ArgumentParser,
) -> Dict[str, argparse.ArgumentParser]:
    """
    Return the subcommand parsers registered on ``parser``.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The top-level parser built by :func:`create_parser`.

    Returns
    -------
    dict of str to argparse.ArgumentParser
        Mapping of subcommand name to its parser. Empty if the parser has no
        subcommands.

    Notes
    -----
    argparse does not expose a public accessor for this, so the subparsers
    action is located by type and its ``choices`` mapping read. That mapping
    has been stable across every Python version TIRmite supports.
    """
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return dict(action.choices)
    return {}


def main() -> int:
    """
    Main CLI entry point for TIRmite.

    Returns
    -------
    int
        Exit code from subcommand execution.

    Notes
    -----
    Parses command-line arguments and dispatches to the appropriate subcommand
    handler. Running ``tirmite`` with no arguments prints the top-level help;
    running a subcommand with no arguments prints *that subcommand's* help.
    Both exit with status 2, argparse's conventional code for a usage error.
    """
    parser = create_parser()

    # A bare `tirmite` gets the top-level help.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(2)

    # A bare `tirmite <subcommand>` gets that subcommand's help rather than a
    # late failure from inside its main(). Previously only the top-level case
    # was handled, so e.g. `tirmite search` parsed successfully, created its
    # output directory and initialised logging before reporting that it had no
    # inputs.
    if len(sys.argv) == 2:
        subcommands = _subcommand_parsers(parser)
        subcommand = subcommands.get(sys.argv[1])
        if subcommand is not None:
            subcommand.print_help()
            sys.exit(2)

    args = parser.parse_args()

    # Dispatch to appropriate subcommand
    if args.command == 'legacy':
        from tirmite.cli.legacy import main as legacy_main

        result = legacy_main(args)
        return int(result) if result is not None else 0
    elif args.command == 'seed':
        from tirmite.cli.hmm_build import main as seed_main

        result = seed_main(args)
        return int(result) if result is not None else 0
    elif args.command == 'pair':
        from tirmite.cli.hmm_pair import main as pair_main

        result = pair_main(args)
        return int(result) if result is not None else 0
    elif args.command == 'search':
        from tirmite.cli.ensemble_search import main as search_main

        result = search_main(args)
        return int(result) if result is not None else 0
    elif args.command == 'validate':
        from tirmite.cli.validate import main as validate_main

        result = validate_main(args)
        return int(result) if result is not None else 0
    else:
        # Reached when argparse accepted the arguments but dispatch found no
        # matching handler, i.e. a subcommand was registered on the parser but
        # never wired up here. Logged so the mismatch is visible rather than
        # looking like an ordinary no-arguments invocation.
        logger.debug(f'No dispatch handler for command: {args.command!r}')
        parser.print_help()
        sys.exit(2)


if __name__ == '__main__':
    sys.exit(main())
