#!/usr/bin/env python3
"""TIRmite-pair: Pair precomputed nhmmer hits for transposon terminal repeat detection.

This module takes precomputed nhmmer results and runs the pairing workflow:
1. Import nhmmer hits from tabular files
2. Filter hits by model coverage and e-value
3. Run pairing algorithms (symmetric or asymmetric)
4. Output paired elements, individual hits, and GFF annotations
"""

import logging
import os
from pathlib import Path
import sys

from tirmite._version import __version__
import tirmite.tirmitetools as tirmite
from tirmite.utils.logs import init_logging
from tirmite.utils.utils import (
    cleanup_temp_directory,
    indexGenome,
    setup_directories,
)


def get_hmm_model_length(hmm_file_path):
    """Extract model length from HMM file.

    Args:
        hmm_file_path: Path to HMM file

    Returns:
        dict: {model_name: model_length}

    """
    model_lengths = {}

    try:
        with open(hmm_file_path, 'r') as f:
            current_model = None
            for line in f:
                line = line.strip()
                if line.startswith('NAME  '):
                    current_model = line.split()[1]
                elif line.startswith('LENG  '):
                    if current_model:
                        length = int(line.split()[1])
                        model_lengths[current_model] = length
                elif line == '//':
                    current_model = None

    except Exception as e:
        logging.error(f'Error reading HMM file {hmm_file_path}: {e}')

    return model_lengths


def load_model_lengths_file(lengths_file):
    """Load model lengths from tab-delimited file.

    Format: model_name<TAB>model_length

    Args:
        lengths_file: Path to tab-delimited file

    Returns:
        dict: {model_name: model_length}

    """
    model_lengths = {}

    try:
        with open(lengths_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) != 2:
                    logging.warning(
                        f'Skipping malformed line {line_num} in {lengths_file}: {line}'
                    )
                    continue

                model_name, model_length_str = parts
                try:
                    model_length = int(model_length_str)
                    model_lengths[model_name] = model_length
                except ValueError:
                    logging.warning(
                        f'Invalid model length on line {line_num}: {model_length_str}'
                    )

    except Exception as e:
        logging.error(f'Error reading model lengths file {lengths_file}: {e}')

    return model_lengths


def calculate_hit_coverage(hitTable, model_lengths):
    """Calculate coverage for hits based on model lengths.

    Args:
        hitTable: DataFrame with hits
        model_lengths: Dict of {model_name: model_length}

    Returns:
        DataFrame with coverage column added

    """
    hitTable = hitTable.copy()
    coverage_values = []

    for _idx, row in hitTable.iterrows():
        model = row['model']
        hit_length = abs(int(row['hitEnd']) - int(row['hitStart'])) + 1

        if model in model_lengths:
            model_length = model_lengths[model]
            coverage = hit_length / model_length
        else:
            logging.warning(f'Model length not found for {model}, using coverage = 0')
            coverage = 0.0

        coverage_values.append(coverage)

    hitTable['coverage'] = coverage_values
    return hitTable


def filter_hits_coverage(hitTable, mincov):
    """Filter hits by coverage threshold.

    Args:
        hitTable: DataFrame with coverage column
        mincov: Minimum coverage threshold

    Returns:
        Filtered DataFrame

    """
    return hitTable[hitTable['coverage'] >= mincov]


def extract_model_name_from_path(model_path):
    """Extract model name from HMM file path by reading the HMM file."""
    if not model_path:
        return None

    try:
        with open(model_path, 'r') as f:
            for line in f:
                if line.startswith('NAME  '):
                    return line.split()[1].strip()
    except (FileNotFoundError, IOError):
        # Fallback to basename without extension
        return Path(model_path).stem

    return Path(model_path).stem


def add_pair_parser(subparsers):
    """Add pair subcommand parser."""
    parser = subparsers.add_parser(
        'pair',
        help='Pair precomputed nhmmer hits',
        description='Pair precomputed nhmmer hits for transposon detection',
    )

    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__),
    )

    parser.add_argument(
        '--loglevel',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        help="Set logging level. Default: 'INFO'",
    )

    parser.add_argument(
        '--logfile',
        action='store_true',
        default=False,
        help='Write log messages to file in output directory.',
    )

    # Required inputs
    parser.add_argument(
        '--genome',
        type=str,
        required=True,
        help='Path to target genome FASTA file.',
    )

    # nhmmer input files (mutually exclusive groups)
    nhmmer_group = parser.add_mutually_exclusive_group(required=True)
    nhmmer_group.add_argument(
        '--nhmmerFile',
        type=str,
        help='Path to single nhmmer output file (requires --hmmFile or --lengthsFile).',
    )
    nhmmer_group.add_argument(
        '--leftNhmmer',
        type=str,
        help='Path to nhmmer output for left model (use with --rightNhmmer).',
    )

    parser.add_argument(
        '--rightNhmmer',
        type=str,
        help='Path to nhmmer output for right model (use with --leftNhmmer).',
    )

    # Model length sources (mutually exclusive)
    length_group = parser.add_mutually_exclusive_group()
    length_group.add_argument(
        '--hmmFile',
        type=str,
        help='Path to HMM file for extracting model lengths (for single model pairing).',
    )
    length_group.add_argument(
        '--leftModel',
        type=str,
        help='Path to left HMM model file (for asymmetric pairing).',
    )
    length_group.add_argument(
        '--lengthsFile',
        type=str,
        help='Path to tab-delimited file with model_name and model_length columns.',
    )

    parser.add_argument(
        '--rightModel',
        type=str,
        help='Path to right HMM model file (for asymmetric pairing).',
    )

    # Filtering parameters
    parser.add_argument(
        '--maxeval',
        type=float,
        default=0.001,
        help='Maximum e-value allowed for valid hit. Default: 0.001',
    )

    parser.add_argument(
        '--mincov',
        type=float,
        default=0.5,
        help='Minimum hit coverage as proportion of model length. Default: 0.5',
    )

    parser.add_argument(
        '--maxdist',
        type=int,
        default=None,
        help='Maximum distance allowed between termini for pairing.',
    )

    # Pairing configuration
    parser.add_argument(
        '--orientation',
        type=str,
        default='F,R',
        help='Orientation pattern for pairing. F=Forward(+), R=Reverse(-). Options: F,R (TIR), F,F (LTR), R,R, R,F. Default: F,R',
    )

    parser.add_argument(
        '--stableReps',
        type=int,
        default=0,
        help='Number of iterations when no new pairs found. Default: 0',
    )

    # Output options
    parser.add_argument(
        '--outdir',
        type=str,
        default=None,
        help='Output directory. Default: current directory',
    )

    parser.add_argument(
        '--prefix',
        type=str,
        default=None,
        help='Prefix for output files.',
    )

    parser.add_argument(
        '--nopairing',
        action='store_true',
        default=False,
        help='Only report individual hits, skip pairing.',
    )

    parser.add_argument(
        '--gffOut',
        action='store_true',
        default=False,
        help='Generate GFF3 output file.',
    )

    parser.add_argument(
        '--report',
        default='all',
        choices=['all', 'paired', 'unpaired'],
        help='Types of hits to include in GFF output.',
    )

    parser.add_argument(
        '--padlen',
        type=int,
        default=None,
        help='Extract N bases flanking each hit in FASTA output.',
    )

    # Utility options
    parser.add_argument(
        '--tempdir',
        type=str,
        default=None,
        help='Base directory for temporary files.',
    )

    parser.add_argument(
        '--keep-temp',
        action='store_true',
        default=False,
        help='Preserve temporary directory.',
    )

    return parser


def validate_arguments(args):
    """Validate argument combinations and file existence."""
    # Check asymmetric pairing requirements
    if args.leftNhmmer and not args.rightNhmmer:
        raise ValueError('--leftNhmmer requires --rightNhmmer')
    if args.rightNhmmer and not args.leftNhmmer:
        raise ValueError('--rightNhmmer requires --leftNhmmer')
    if args.leftModel and not args.rightModel:
        raise ValueError('--leftModel requires --rightModel')
    if args.rightModel and not args.leftModel:
        raise ValueError('--rightModel requires --leftModel')

    # Check model length source requirements
    if args.nhmmerFile:
        if not (args.hmmFile or args.lengthsFile):
            raise ValueError('--nhmmerFile requires either --hmmFile or --lengthsFile')

    if args.leftNhmmer and args.rightNhmmer:
        if not (args.leftModel and args.rightModel) and not args.lengthsFile:
            raise ValueError(
                'Asymmetric pairing requires --leftModel/--rightModel or --lengthsFile'
            )

    # Check file existence
    required_files = [args.genome]

    if args.nhmmerFile:
        required_files.append(args.nhmmerFile)
    if args.leftNhmmer:
        required_files.extend([args.leftNhmmer, args.rightNhmmer])
    if args.hmmFile:
        required_files.append(args.hmmFile)
    if args.leftModel:
        required_files.extend([args.leftModel, args.rightModel])
    if args.lengthsFile:
        required_files.append(args.lengthsFile)

    for file_path in required_files:
        if not Path(file_path).exists():
            raise FileNotFoundError(f'Required file not found: {file_path}')


def main(args=None):
    """Main entry point for tirmite-pair."""
    try:
        # Validate arguments
        try:
            validate_arguments(args)
        except (ValueError, FileNotFoundError) as e:
            logging.error(f'Argument validation failed: {e}')
            sys.exit(1)

        # Set up directories first (we need output dir for default logfile location)
        try:
            outDir, tempDir = setup_directories(args)
        except (OSError, FileNotFoundError, PermissionError) as e:
            # If directory setup fails, we can't even set up logging properly
            print(f'Directory setup failed: {e}', file=sys.stderr)
            sys.exit(1)

        # Determine logfile path
        logfile_path = None
        if args.logfile:
            # Create logfile in output directory
            if args.prefix:
                logfile_name = f'{args.prefix}_tirmite_pair.log'
            else:
                logfile_name = 'tirmite_pair.log'
            logfile_path = outDir / logfile_name

        # Set up logging with or without file output
        init_logging(loglevel=args.loglevel, logfile=logfile_path)

        logging.info(f'TIRmite-pair version {__version__}')
        logging.info(f'Output directory: {outDir}')
        logging.debug(f'Temporary directory: {tempDir}')

        # Index genome
        logging.info(f'Indexing genome: {args.genome}')
        genome, genome_descriptions = indexGenome(args.genome)

        # Load model lengths
        logging.info('Loading model lengths...')
        model_lengths = {}

        if args.lengthsFile:
            model_lengths = load_model_lengths_file(args.lengthsFile)
        elif args.hmmFile:
            model_lengths = get_hmm_model_length(args.hmmFile)
        elif args.leftModel and args.rightModel:
            left_lengths = get_hmm_model_length(args.leftModel)
            right_lengths = get_hmm_model_length(args.rightModel)
            model_lengths = {**left_lengths, **right_lengths}

        if not model_lengths:
            logging.error('No model lengths could be loaded')
            cleanup_temp_directory(tempDir, args.keep_temp)
            sys.exit(1)

        logging.info(f'Loaded lengths for models: {", ".join(model_lengths.keys())}')

        # Import nhmmer hits
        logging.info('Importing nhmmer hits...')
        hitTable = None

        if args.nhmmerFile:
            # Single file mode
            hitTable = tirmite.import_nhmmer(
                infile=args.nhmmerFile, hitTable=None, prefix=args.prefix
            )
        else:
            # Asymmetric mode - import from both files
            left_model_name = extract_model_name_from_path(args.leftModel)
            right_model_name = extract_model_name_from_path(args.rightModel)

            hitTable = tirmite.import_nhmmer(
                infile=args.leftNhmmer,
                hitTable=None,
                prefix=args.prefix,
            )
            hitTable = tirmite.import_nhmmer(
                infile=args.rightNhmmer,
                hitTable=hitTable,
                prefix=args.prefix,
            )

        logging.info(f'Imported {len(hitTable)} total hits')

        # Calculate hit coverage
        logging.info('Calculating hit coverage...')
        hitTable = calculate_hit_coverage(hitTable, model_lengths)

        # Apply coverage filtering
        logging.info(f'Filtering hits with coverage < {args.mincov}')
        hitCount = len(hitTable)

        # Log pre-filtering counts
        model_counts_before = hitTable['model'].value_counts().to_dict()
        for model, count in model_counts_before.items():
            logging.info(f'Model {model}: {count} hits before coverage filtering')

        hitTable = filter_hits_coverage(hitTable, args.mincov)

        # Log post-filtering counts
        model_counts_after = hitTable['model'].value_counts().to_dict()
        total_excluded = hitCount - len(hitTable)
        logging.info(f'Excluded {total_excluded} hits on coverage criteria')

        for model in model_counts_before:
            before = model_counts_before[model]
            after = model_counts_after.get(model, 0)
            excluded = before - after
            logging.info(f'Model {model}: {excluded} excluded, {after} remaining')

        # Apply e-value filtering
        logging.info(f'Filtering hits with e-value > {args.maxeval}')
        hitCount = len(hitTable)

        model_counts_before = hitTable['model'].value_counts().to_dict()
        hitTable = tirmite.filterHitsEval(maxeval=args.maxeval, hitTable=hitTable)

        model_counts_after = hitTable['model'].value_counts().to_dict()
        total_excluded = hitCount - len(hitTable)
        logging.info(f'Excluded {total_excluded} hits on e-value criteria')

        for model in model_counts_before:
            before = model_counts_before[model]
            after = model_counts_after.get(model, 0)
            excluded = before - after
            logging.info(f'Model {model}: {excluded} excluded, {after} remaining')

        logging.info(f'Total remaining hits: {len(hitTable)}')

        # Convert to dict structure
        hitsDict, hitIndex = tirmite.table2dict(hitTable)

        # Create pairing configuration
        if args.leftNhmmer and args.rightNhmmer:
            # Asymmetric pairing
            left_model_name = extract_model_name_from_path(args.leftModel)
            right_model_name = extract_model_name_from_path(args.rightModel)

            config = tirmite.PairingConfig(
                orientation=args.orientation,
                left_model=left_model_name,
                right_model=right_model_name,
            )
        else:
            # Symmetric pairing
            available_models = list(hitsDict.keys())
            if not available_models:
                logging.error('No models found in hits')
                cleanup_temp_directory(tempDir, args.keep_temp)
                sys.exit(1)

            single_model = available_models[0]
            config = tirmite.PairingConfig(
                orientation=args.orientation, single_model=single_model
            )

        # Write individual hits
        logging.info('Writing individual hits to FASTA...')
        tirmite.writeTIRs(
            outDir=outDir,
            hitTable=hitTable,
            maxeval=args.maxeval,
            genome=genome,
            prefix=args.prefix,
            padlen=args.padlen,
            genome_descriptions=genome_descriptions,
        )

        # Skip pairing if requested
        if args.nopairing:
            logging.info('Pairing disabled. Analysis complete.')
            cleanup_temp_directory(tempDir, args.keep_temp)
            return

        # Run pairing
        logging.info('Searching for candidate pairings...')
        hitIndex = tirmite.parseHitsGeneral(
            hitsDict=hitsDict, hitIndex=hitIndex, maxDist=args.maxdist, config=config
        )

        logging.info('Performing iterative pairing...')
        if config.is_asymmetric:
            logging.info(
                f'Using asymmetric pairing: {config.left_model} + {config.right_model}'
            )
            hitIndex, paired, unpaired = tirmite.iterateGetPairsAsymmetric(
                hitIndex, config, stableReps=args.stableReps
            )
        else:
            logging.info(
                f'Using symmetric pairing with orientation {config.orientation}'
            )
            hitIndex, paired, unpaired = tirmite.iterateGetPairsCustom(
                hitIndex, config, stableReps=args.stableReps
            )

        # Log pairing results
        total_pairs = sum(len(pairs) for pairs in paired.values())
        total_unpaired = len(unpaired)
        logging.info(
            f'Pairing completed: {total_pairs} pairs, {total_unpaired} unpaired'
        )

        # Write paired TIRs
        if args.report in ['all', 'paired']:
            logging.info('Writing paired termini to FASTA...')
            tirmite.writePairedTIRs(
                outDir=outDir,
                paired=paired,
                hitIndex=hitIndex,
                genome=genome,
                prefix=args.prefix,
                padlen=args.padlen,
                genome_descriptions=genome_descriptions,
            )

        # Extract and write elements
        logging.info('Extracting paired elements...')
        pairedEles = tirmite.fetchElements(
            paired=paired,
            hitIndex=hitIndex,
            genome=genome,
            genome_descriptions=genome_descriptions,
        )

        logging.info('Writing paired elements to FASTA...')
        tirmite.writeElements(outDir, eleDict=pairedEles, prefix=args.prefix)

        # Write GFF if requested
        if args.gffOut:
            # Get unpaired TIRs if needed
            unpairedTIRs = None
            if args.report in ['all', 'unpaired']:
                unpairedTIRs = tirmite.fetchUnpaired(hitIndex=hitIndex)
                logging.info(f'Found {len(unpairedTIRs)} unpaired termini')

            # Set GFF output path
            if args.prefix:
                gffOutPath = os.path.join(outDir, f'{args.prefix}.gff3')
            else:
                gffOutPath = os.path.join(outDir, 'tirmite_pair_report.gff3')

            logging.info(f'Writing GFF3 output: {gffOutPath}')
            tirmite.gffWrite(
                outpath=gffOutPath,
                featureList=pairedEles,
                writeTIRs=args.report,
                unpaired=unpairedTIRs,
                prefix=args.prefix,
            )

        logging.info('TIRmite-pair analysis completed successfully')

    except KeyboardInterrupt:
        logging.info('Analysis interrupted by user')
        sys.exit(130)
    except Exception as e:
        logging.error(f'Unexpected error: {e}')
        logging.exception('Full traceback:')  # This will log the full stack trace
        sys.exit(1)
    finally:
        if 'tempDir' in locals():
            cleanup_temp_directory(tempDir, args.keep_temp)

        # Log completion message with logfile location if enabled
        if 'logfile_path' in locals() and logfile_path and args.logfile:
            logging.info(f'Log file saved to: {logfile_path}')


if __name__ == '__main__':
    main()
