#!/usr/bin/env python3
"""
TIRmite command-line interface with subcommands.
"""

import argparse
import sys

from tirmite._version import __version__  # type: ignore[import-not-found]


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

Examples:
  tirmite legacy --genome genome.fa --hmm-file model.hmm
  tirmite seed --left-seed left.fa --model-name myTE --genome genome.fa
  tirmite pair --genome genome.fa --nhmmer-file hits.out --hmm-file model.hmm
  tirmite search --blast-results hits.tab --cluster-map clusters.tsv
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

    add_legacy_parser(subparsers)
    add_seed_parser(subparsers)
    add_pair_parser(subparsers)
    add_search_parser(subparsers)

    return parser


def main() -> int:
    """
    Main CLI entry point for TIRmite.

    Returns
    -------
    int
        Exit code from subcommand execution or 1 if no subcommand specified.

    Notes
    -----
    Parses command-line arguments and dispatches to appropriate subcommand handler.
    Prints help message if run without arguments.
    """
    parser = create_parser()

    # Parse arguments
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

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
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    sys.exit(main())
