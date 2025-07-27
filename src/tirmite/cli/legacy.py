import argparse
import glob
import logging
import os
from pathlib import Path
import shutil
import sys

from tirmite._version import __version__
from tirmite.runners.hmmer_wrappers import process_hmmer_workflow
import tirmite.tirmitetools as tirmite
from tirmite.utils.logs import init_logging
from tirmite.utils.utils import (
    cleanup_temp_directory,
    indexGenome,
    setup_directories,
)


def mainArgs():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Map HMM models of transposon termini to genomic sequences for annotation \
        of variable non-autonomous and complete transposons.',
        prog='tirmite',
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
        help="Set logging level. Default: 'DEBUG'",
    )
    # Input
    parser.add_argument(
        '--genome',
        type=str,
        required=True,
        help='Path to target genome that will be queried with HMMs.',
    )
    parser.add_argument(
        '--hmmDir',
        type=str,
        default=None,
        help='Directory containing pre-prepared TIR-pHMMs.',
    )
    parser.add_argument(
        '--hmmFile',
        type=str,
        default=None,
        help='Path to single HMM file. Incompatible with "--hmmDir".',
    )
    parser.add_argument(
        '--alnDir',
        type=str,
        default=None,
        help='Path to directory containing only TIR alignments in FASTA format to be converted to HMM.',
    )
    parser.add_argument(
        '--alnFile',
        type=str,
        default=None,
        help='Provide a single TIR alignment in FASTA format to be converted to HMM. Incompatible with "--alnDir".',
    )
    parser.add_argument(
        '--pairbed',
        type=str,
        default=None,
        help='If set TIRmite will preform pairing on TIRs from custom bedfile only.',
    )

    parser.add_argument(
        '--stableReps',
        type=int,
        default=0,
        help='Number of times to iterate pairing procedure when no additional pairs are found AND remaining unpaired hits > 0.',
    )
    # Output and housekeeping
    parser.add_argument(
        '--outdir',
        type=str,
        default=None,
        help='All output files will be written to this directory.',
    )
    parser.add_argument(
        '--prefix',
        type=str,
        default=None,
        help='Add prefix to all hits and paired elements detected in this run. (Default = None)',
    )
    parser.add_argument(
        '--nopairing',
        action='store_true',
        default=False,
        help='If set, only report HMM hits. Do not attempt pairing.',
    )
    parser.add_argument(
        '--gffOut',
        action='store_true',
        default=False,
        help='If set report features as prefix.gff3. File saved to outdir. Default: False',
    )
    parser.add_argument(
        '--report',
        default='all',
        choices=[None, 'all', 'paired', 'unpaired'],
        help='Options for reporting model hits in GFF annotation file.',
    )
    parser.add_argument(
        '--padlen',
        type=int,
        default=None,
        help='Extract x bases either side of model hit when writing hits to fasta.',
    )
    parser.add_argument(
        '--keep-temp',
        action='store_true',
        default=False,
        help='If set do not delete temp file directory.',
    )
    # HMMER options
    parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help='Set number of threads available to hmmer software.',
    )
    parser.add_argument(
        '--maxeval',
        type=float,
        default=0.001,
        help='Maximum e-value allowed for valid hit. Default = 0.001',
    )
    parser.add_argument(
        '--maxdist',
        type=int,
        default=None,
        help='Maximum distance allowed between termini candidates to consider valid pairing.',
    )
    parser.add_argument(
        '--nobias',
        action='store_true',
        default=False,
        help='Turn OFF bias correction of scores in nhmmer.',
    )
    parser.add_argument(
        '--matrix',
        type=str,
        default=None,
        help='Use custom DNA substitution matrix with nhmmer.',
    )
    parser.add_argument(
        '--mincov',
        type=float,
        default=0.5,
        help='Minimum valid hit length as prop of model length. Defaults to 0.5',
    )
    # Non-standard HMMER paths
    parser.add_argument(
        '--hmmpress',
        type=str,
        default='hmmpress',
        help='Set location of hmmpress if not in PATH.',
    )
    parser.add_argument(
        '--nhmmer',
        type=str,
        default='nhmmer',
        help='Set location of nhmmer if not in PATH.',
    )
    parser.add_argument(
        '--hmmbuild',
        type=str,
        default='hmmbuild',
        help='Set location of hmmbuild if not in PATH.',
    )
    parser.add_argument(
        '--orientation',
        type=str,
        default='F,R',
        help='Orientation pattern for pairing hits. F=Forward, R=Reverse. Options: F,R (TIR), F,F (LTR), R,R, R,F',
    )

    parser.add_argument(
        '--leftModel',
        type=str,
        default=None,
        help='HMM model for left terminus. Use with --rightModel for asymmetric elements.',
    )

    parser.add_argument(
        '--rightModel',
        type=str,
        default=None,
        help='HMM model for right terminus. Use with --leftModel for asymmetric elements.',
    )
    # Add new temp directory option
    parser.add_argument(
        '--tempdir',
        type=str,
        default=None,
        help='Base directory for temporary files. Uses system temp if not specified.',
    )
    parser.add_argument(
        '--nhmmerFile',
        type=str,
        default=None,
        help='Path to precomputed nhmmer output file. Requires --hmmFile for model length calculation.',
    )
    parser.add_argument(
        '--leftNhmmer',
        type=str,
        default=None,
        help='Path to precomputed nhmmer output for left model. Use with --rightNhmmer and --leftModel/--rightModel for asymmetric elements.',
    )
    parser.add_argument(
        '--rightNhmmer',
        type=str,
        default=None,
        help='Path to precomputed nhmmer output for right model. Use with --leftNhmmer and --leftModel/--rightModel for asymmetric elements.',
    )
    args = parser.parse_args()
    return args


def missing_tool(tool_name):
    path = shutil.which(tool_name)
    if path is None:
        return [tool_name]
    else:
        return []


def validate_pairbed_compatibility(hitTable, config, args):
    """
    Validate that BED file contents are compatible with pairing configuration.

    Args:
        hitTable: DataFrame of hits from BED file
        config: PairingConfig object with orientation and model settings
        args: Command line arguments

    Returns:
        bool: True if compatible, False otherwise
    """
    logging.info('Validating BED file compatibility with pairing configuration...')

    # Check if any hits were loaded
    if len(hitTable) == 0:
        logging.error('No valid hits found in BED file')
        return False

    # Get available models from BED file
    available_models = set(hitTable['model'].unique())
    logging.info(f'Models found in BED file: {", ".join(sorted(available_models))}')

    # Check required models are present
    if config.is_asymmetric:
        required_models = {config.left_model, config.right_model}
        missing_models = required_models - available_models
        if missing_models:
            logging.error(
                f'Required models missing from BED file: {", ".join(missing_models)}'
            )
            logging.error(f'Available models: {", ".join(sorted(available_models))}')
            return False
        logging.info(
            f'Asymmetric pairing: Left={config.left_model}, Right={config.right_model}'
        )
    else:
        # For symmetric pairing, need at least one model
        if config.left_model not in available_models:
            logging.error(f'Required model {config.left_model} not found in BED file')
            logging.error(f'Available models: {", ".join(sorted(available_models))}')
            return False
        logging.info(f'Symmetric pairing using model: {config.left_model}')

    # Check strand orientations are compatible
    available_strands = set(hitTable['strand'].unique())
    required_strands = {config.left_strand, config.right_strand}

    if config.left_strand == config.right_strand:
        # Same orientation pairing (e.g., F,F for LTR)
        if config.left_strand not in available_strands:
            logging.error(
                f'Required strand orientation {config.left_strand} not found in BED file'
            )
            logging.error(f'Available strands: {", ".join(sorted(available_strands))}')
            return False
        logging.info(
            f'Same orientation pairing: {config.orientation[0]},{config.orientation[1]}'
        )
    else:
        # Different orientation pairing (e.g., F,R for TIR)
        missing_strands = required_strands - available_strands
        if missing_strands:
            logging.error(
                f'Required strand orientations missing from BED file: {", ".join(missing_strands)}'
            )
            logging.error(f'Available strands: {", ".join(sorted(available_strands))}')
            return False
        logging.info(
            f'Inverted orientation pairing: {config.orientation[0]},{config.orientation[1]}'
        )

    # Count hits per model and strand combination
    hit_counts = {}
    for model in available_models:
        model_hits = hitTable[hitTable['model'] == model]
        hit_counts[model] = {}
        for strand in ['+', '-']:
            count = len(model_hits[model_hits['strand'] == strand])
            hit_counts[model][strand] = count
            if count > 0:
                logging.info(f'Model {model} strand {strand}: {count} hits')

    # Warn about potential pairing issues
    if config.is_asymmetric:
        left_hits = hit_counts.get(config.left_model, {}).get(config.left_strand, 0)
        right_hits = hit_counts.get(config.right_model, {}).get(config.right_strand, 0)
        if left_hits == 0 or right_hits == 0:
            logging.warning(
                f'Potential pairing issue: Left model has {left_hits} hits, Right model has {right_hits} hits'
            )
    else:
        left_hits = hit_counts.get(config.left_model, {}).get(config.left_strand, 0)
        right_hits = hit_counts.get(config.left_model, {}).get(config.right_strand, 0)
        if left_hits == 0 or right_hits == 0:
            logging.warning(
                f'Potential pairing issue: {left_hits} left-strand hits, {right_hits} right-strand hits for model {config.left_model}'
            )

    return True


# Extract model names from file paths for validation
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


def add_legacy_parser(subparsers):
    """Add legacy subcommand parser."""
    parser = subparsers.add_parser(
        'legacy',
        help='Original TIRmite workflow (HMM search + pairing)',
        description='Map HMM models of transposon termini to genomic sequences',
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
        help="Set logging level. Default: 'DEBUG'",
    )
    # Input
    parser.add_argument(
        '--genome',
        type=str,
        required=True,
        help='Path to target genome that will be queried with HMMs.',
    )
    parser.add_argument(
        '--hmmDir',
        type=str,
        default=None,
        help='Directory containing pre-prepared TIR-pHMMs.',
    )
    parser.add_argument(
        '--hmmFile',
        type=str,
        default=None,
        help='Path to single HMM file. Incompatible with "--hmmDir".',
    )
    parser.add_argument(
        '--alnDir',
        type=str,
        default=None,
        help='Path to directory containing only TIR alignments in FASTA format to be converted to HMM.',
    )
    parser.add_argument(
        '--alnFile',
        type=str,
        default=None,
        help='Provide a single TIR alignment in FASTA format to be converted to HMM. Incompatible with "--alnDir".',
    )
    parser.add_argument(
        '--pairbed',
        type=str,
        default=None,
        help='If set TIRmite will preform pairing on TIRs from custom bedfile only.',
    )

    parser.add_argument(
        '--stableReps',
        type=int,
        default=0,
        help='Number of times to iterate pairing procedure when no additional pairs are found AND remaining unpaired hits > 0.',
    )
    # Output and housekeeping
    parser.add_argument(
        '--outdir',
        type=str,
        default=None,
        help='All output files will be written to this directory.',
    )
    parser.add_argument(
        '--prefix',
        type=str,
        default=None,
        help='Add prefix to all hits and paired elements detected in this run. (Default = None)',
    )
    parser.add_argument(
        '--nopairing',
        action='store_true',
        default=False,
        help='If set, only report HMM hits. Do not attempt pairing.',
    )
    parser.add_argument(
        '--gffOut',
        action='store_true',
        default=False,
        help='If set report features as prefix.gff3. File saved to outdir. Default: False',
    )
    parser.add_argument(
        '--report',
        default='all',
        choices=[None, 'all', 'paired', 'unpaired'],
        help='Options for reporting model hits in GFF annotation file.',
    )
    parser.add_argument(
        '--padlen',
        type=int,
        default=None,
        help='Extract x bases either side of model hit when writing hits to fasta.',
    )
    parser.add_argument(
        '--keep-temp',
        action='store_true',
        default=False,
        help='If set do not delete temp file directory.',
    )
    # HMMER options
    parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help='Set number of threads available to hmmer software.',
    )
    parser.add_argument(
        '--maxeval',
        type=float,
        default=0.001,
        help='Maximum e-value allowed for valid hit. Default = 0.001',
    )
    parser.add_argument(
        '--maxdist',
        type=int,
        default=None,
        help='Maximum distance allowed between termini candidates to consider valid pairing.',
    )
    parser.add_argument(
        '--nobias',
        action='store_true',
        default=False,
        help='Turn OFF bias correction of scores in nhmmer.',
    )
    parser.add_argument(
        '--matrix',
        type=str,
        default=None,
        help='Use custom DNA substitution matrix with nhmmer.',
    )
    parser.add_argument(
        '--mincov',
        type=float,
        default=0.5,
        help='Minimum valid hit length as prop of model length. Defaults to 0.5',
    )
    # Non-standard HMMER paths
    parser.add_argument(
        '--hmmpress',
        type=str,
        default='hmmpress',
        help='Set location of hmmpress if not in PATH.',
    )
    parser.add_argument(
        '--nhmmer',
        type=str,
        default='nhmmer',
        help='Set location of nhmmer if not in PATH.',
    )
    parser.add_argument(
        '--hmmbuild',
        type=str,
        default='hmmbuild',
        help='Set location of hmmbuild if not in PATH.',
    )
    parser.add_argument(
        '--orientation',
        type=str,
        default='F,R',
        help='Orientation pattern for pairing hits. F=Forward, R=Reverse. Options: F,R (TIR), F,F (LTR), R,R, R,F',
    )

    parser.add_argument(
        '--leftModel',
        type=str,
        default=None,
        help='HMM model for left terminus. Use with --rightModel for asymmetric elements.',
    )

    parser.add_argument(
        '--rightModel',
        type=str,
        default=None,
        help='HMM model for right terminus. Use with --leftModel for asymmetric elements.',
    )
    # Add new temp directory option
    parser.add_argument(
        '--tempdir',
        type=str,
        default=None,
        help='Base directory for temporary files. Uses system temp if not specified.',
    )
    parser.add_argument(
        '--nhmmerFile',
        type=str,
        default=None,
        help='Path to precomputed nhmmer output file. Requires --hmmFile for model length calculation.',
    )
    parser.add_argument(
        '--leftNhmmer',
        type=str,
        default=None,
        help='Path to precomputed nhmmer output for left model. Use with --rightNhmmer and --leftModel/--rightModel for asymmetric elements.',
    )
    parser.add_argument(
        '--rightNhmmer',
        type=str,
        default=None,
        help='Path to precomputed nhmmer output for right model. Use with --leftNhmmer and --leftModel/--rightModel for asymmetric elements.',
    )

    return parser


def main(args=None):
    """Do the work."""
    # Get cmd line args
    if args is None:
        args = mainArgs()

    # TODO: Remove use of verbose option
    # Manually add args.verbose and set to True
    args.verbose = True

    # Set up logging
    init_logging(loglevel=args.loglevel)

    # Check for required programs.
    tools = [args.hmmpress, args.nhmmer, args.hmmbuild]

    missing_tools = []
    for tool in tools:
        missing_tools += missing_tool(tool)
    if missing_tools:
        logging.warning(
            'Some tools required by tirmite could not be found: '
            + ', '.join(missing_tools)
        )
        logging.warning('You may need to install them to use all features.')

    # Set up directories with proper validation and error handling
    try:
        outDir, tempDir = setup_directories(args)
        logging.info(f'Output directory: {outDir}')
        logging.debug(f'Temporary directory: {tempDir}')
    except (OSError, FileNotFoundError, PermissionError) as e:
        logging.error(f'Directory setup failed: {e}')
        sys.exit(1)

    # Use context manager for better temp directory management
    try:
        # Index reference genome using pyfaidx for efficient access
        logging.info('Indexing genome from: %s ' % args.genome)
        genome, genome_descriptions = indexGenome(args.genome)

        # Import custom TIR hits from BEDfile.
        if args.pairbed:
            # Die if no input file
            if not glob.glob(os.path.abspath(args.pairbed)):
                logging.error('BED file %s not found. Quitting.' % args.pairbed)
                cleanup_temp_directory(tempDir, args.keep_temp)
                sys.exit(1)

            logging.info(
                'Skipping HMM search. Using custom TIRs from file: %s' % args.pairbed
            )

            # Import hits from BED file
            hitTable = None
            logging.info('Loading custom termini hits from: %s' % str(args.pairbed))
            hitTable = tirmite.import_BED(
                infile=args.pairbed, hitTable=hitTable, prefix=args.prefix
            )

            logging.info(f'Loaded {len(hitTable)} hits from BED file')

            # Create pairing configuration BEFORE validation
            if args.leftModel and args.rightModel:
                # Asymmetric pairing with different models
                config = tirmite.PairingConfig(
                    orientation=args.orientation,
                    left_model=args.leftModel,
                    right_model=args.rightModel,
                )
            else:
                # For symmetric pairing with BED file, need to determine single model
                available_models = (
                    list(hitTable['model'].unique()) if len(hitTable) > 0 else []
                )

                if not available_models:
                    logging.error('No models found in BED file')
                    cleanup_temp_directory(tempDir, args.keep_temp)
                    sys.exit(1)

                # Use first available model or user-specified model
                if args.leftModel and args.leftModel in available_models:
                    single_model = args.leftModel
                elif args.rightModel and args.rightModel in available_models:
                    single_model = args.rightModel
                else:
                    single_model = available_models[0]
                    logging.info(f'No specific model provided, using: {single_model}')

                config = tirmite.PairingConfig(
                    orientation=args.orientation, single_model=single_model
                )

            # Validate BED file compatibility with pairing configuration
            if not validate_pairbed_compatibility(hitTable, config, args):
                logging.error(
                    'BED file incompatible with pairing configuration. Quitting.'
                )
                cleanup_temp_directory(tempDir, args.keep_temp)
                sys.exit(1)

            # Apply hit e-value filters
            logging.info('Filtering hits with e-value > %s' % str(args.maxeval))
            hitCount = len(hitTable.index)
            hitTable = tirmite.filterHitsEval(maxeval=args.maxeval, hitTable=hitTable)
            logging.info(
                'Excluded %s hits on e-value criteria.'
                % str(hitCount - len(hitTable.index))
            )
            logging.info('Remaining hits: %s ' % str(len(hitTable.index)))

            # Group hits by model and chromosome (hitsDict), and initiate hit tracker hitIndex to manage pair-searching
            hitsDict, hitIndex = tirmite.table2dict(hitTable)

            # If pairing is off, just report the hits.
            if args.nopairing:
                tirmite.writeTIRs(
                    outDir=outDir,
                    hitTable=hitTable,
                    maxeval=args.maxeval,
                    genome=genome,
                    padlen=args.padlen,
                )
                logging.info('Pairing is off. Reporting hits only.')
                cleanup_temp_directory(tempDir, args.keep_temp)
                sys.exit(0)

        # Else run nhmmer and load TIR hits.
        else:
            # If raw alignments provided, convert to stockholm format.
            if args.alnDir or args.alnFile:
                stockholmDir = tirmite.convertAlign(
                    alnDir=args.alnDir,
                    alnFile=args.alnFile,
                    inFormat='fasta',
                    tempDir=tempDir,
                )
            else:
                stockholmDir = None

            # If pre-built HMM provided, check correct format.
            if args.hmmFile:
                if (
                    os.path.splitext(os.path.basename(args.hmmFile))[1].lstrip('.')
                    != 'hmm'
                ):
                    logging.error('--hmmFile has non-hmm extension. Exiting.')
                    cleanup_temp_directory(tempDir, args.keep_temp)
                    sys.exit(1)

            # Use new HMMER workflow
            try:
                executable_paths = {
                    'hmmbuild': args.hmmbuild,
                    'hmmpress': args.hmmpress,
                    'nhmmer': args.nhmmer,
                }

                search_params = {
                    'evalue': args.maxeval,
                    'cores': args.threads,
                    'nobias': args.nobias,
                    'matrix': args.matrix,
                }

                resultDir, hmmDB = process_hmmer_workflow(
                    hmm_dir=args.hmmDir,
                    hmm_file=args.hmmFile,
                    alignment_dir=stockholmDir,
                    left_model=getattr(args, 'leftModel', None),  # Add this line
                    right_model=getattr(args, 'rightModel', None),  # Add this line
                    temp_dir=tempDir,
                    genome_path=args.genome,
                    executable_paths=executable_paths,
                    search_params=search_params,
                    verbose=args.verbose,
                )

            except Exception as e:
                logging.error(f'HMMER workflow failed: {e}')
                cleanup_temp_directory(tempDir, args.keep_temp)
                sys.exit(1)

            # Die if no hits found
            if not glob.glob(os.path.join(os.path.abspath(resultDir), '*.out')):
                logging.error('No hits found in %s . Quitting.' % resultDir)
                cleanup_temp_directory(tempDir, args.keep_temp)
                sys.exit(1)

            # TODO: Modify this code to allow import of hits from external nhmmer results.
            # Import hits from nhmmer result files
            hitTable = None
            modelCount = 0
            for resultfile in glob.glob(
                os.path.join(os.path.abspath(resultDir), '*.out')
            ):
                logging.info('Loading nhmmer hits from: %s ' % resultfile)
                hitTable = tirmite.import_nhmmer(
                    infile=resultfile, hitTable=hitTable, prefix=args.prefix
                )
                modelCount += 1

            logging.info(
                'Imported %s hits from %s models. '
                % (str(len(hitTable.index)), str(modelCount))
            )

            # Apply hit length filters with detailed logging per model
            logging.info('Filtering hits with < %s model coverage. ' % str(args.mincov))
            hitCount = len(hitTable.index)

            # Log pre-filtering counts per model
            model_counts_before = hitTable['model'].value_counts().to_dict()
            for model, count in model_counts_before.items():
                logging.info(f'Model {model}: {count} hits before coverage filtering')

            hitTable = tirmite.filterHitsLen(
                hmmDB=hmmDB, mincov=args.mincov, hitTable=hitTable
            )

            # Log post-filtering counts per model
            model_counts_after = hitTable['model'].value_counts().to_dict()
            total_excluded = hitCount - len(hitTable.index)
            logging.info(f'Excluded {total_excluded} hits on coverage criteria.')

            for model in model_counts_before:
                before = model_counts_before[model]
                after = model_counts_after.get(model, 0)
                excluded = before - after
                logging.info(
                    f'Model {model}: {excluded} hits excluded, {after} hits remaining'
                )

            logging.info('Total remaining hits: %s ' % str(len(hitTable.index)))

            # Apply hit e-value filters with detailed logging per model
            logging.info('Filtering hits with e-value > %s' % str(args.maxeval))
            hitCount = len(hitTable.index)

            # Log pre-filtering counts per model
            model_counts_before = hitTable['model'].value_counts().to_dict()
            for model, count in model_counts_before.items():
                logging.info(f'Model {model}: {count} hits before e-value filtering')

            hitTable = tirmite.filterHitsEval(maxeval=args.maxeval, hitTable=hitTable)

            # Log post-filtering counts per model
            model_counts_after = hitTable['model'].value_counts().to_dict()
            total_excluded = hitCount - len(hitTable.index)
            logging.info(f'Excluded {total_excluded} hits on e-value criteria.')

            for model in model_counts_before:
                before = model_counts_before[model]
                after = model_counts_after.get(model, 0)
                excluded = before - after
                logging.info(
                    f'Model {model}: {excluded} hits excluded, {after} hits remaining'
                )

            logging.info('Total remaining hits: %s ' % str(len(hitTable.index)))

            # Group hits by model and chromosome (hitsDict), and initiate hit tracker hitIndex to manage pair-searching
            hitsDict, hitIndex = tirmite.table2dict(hitTable)

            # Create pairing configuration for HMM workflow
            if args.leftModel and args.rightModel:
                # Extract actual model names from HMM files
                left_model_name = extract_model_name_from_path(args.leftModel)
                right_model_name = extract_model_name_from_path(args.rightModel)

                # Asymmetric pairing with different models
                config = tirmite.PairingConfig(
                    orientation=args.orientation,
                    left_model=left_model_name,
                    right_model=right_model_name,
                )
                required_models = [left_model_name, right_model_name]
            else:
                # Symmetric pairing with single model (backward compatible)
                single_model = list(hitsDict.keys())[0] if hitsDict else None
                config = tirmite.PairingConfig(
                    orientation=args.orientation, single_model=single_model
                )
                required_models = [single_model]

            # Validate that required models are available
            for model in required_models:
                if model not in hitsDict:
                    available_models = list(hitsDict.keys())
                    logging.error(f'Required model {model} not found in search results')
                    logging.error(f'Available models: {", ".join(available_models)}')
                    cleanup_temp_directory(tempDir, args.keep_temp)
                    sys.exit(1)

            # If pairing is off, just report the hits
            if args.nopairing:
                tirmite.writeTIRs(
                    outDir=outDir,
                    hitTable=hitTable,
                    maxeval=args.maxeval,
                    genome=genome,
                    padlen=args.padlen,
                )
                logging.info('Pairing is off. Reporting hits only.')
                cleanup_temp_directory(tempDir, args.keep_temp)
                sys.exit(0)

        logging.debug('=== HIT TABLE DEBUG ===')
        logging.debug(f'Hit table columns: {list(hitTable.columns)}')
        logging.debug(f'Hit table shape: {hitTable.shape}')
        logging.debug('Hit table contents:')
        for idx, row in hitTable.iterrows():
            logging.debug(f'  Row {idx}: {dict(row)}')

        logging.debug('=== HITS DICT DEBUG ===')
        logging.debug(f'hitsDict keys: {list(hitsDict.keys())}')
        for model, model_hits in hitsDict.items():
            logging.debug(f'  Model {model}:')
            for chrom, chrom_hits in model_hits.items():
                logging.debug(f'    Chrom {chrom}: {len(chrom_hits)} hits')
                for hit in chrom_hits:
                    logging.debug(f'      Hit: {hit}')

        logging.debug('=== HIT INDEX DEBUG ===')
        logging.debug(f'hitIndex type: {type(hitIndex)}')
        logging.debug(
            f'hitIndex keys: {list(hitIndex.keys()) if hasattr(hitIndex, "keys") else "Not a dict"}'
        )
        if hasattr(hitIndex, 'keys'):
            for key, value in hitIndex.items():
                logging.debug(f'  hitIndex[{key}]: {value}')

        # Use generalized parsing instead of original parseHits
        logging.info('Searching for candidate pairings...')
        hitIndex = tirmite.parseHitsGeneral(
            hitsDict=hitsDict, hitIndex=hitIndex, maxDist=args.maxdist, config=config
        )

        # Choose appropriate pairing algorithm based on configuration
        logging.info('Performing iterative pairing...')
        try:
            logging.debug(f'Before pairing - hitIndex type: {type(hitIndex)}')
            logging.debug(
                f'Before pairing - hitIndex keys: {list(hitIndex.keys()) if hasattr(hitIndex, "keys") else "Not a dict"}'
            )

            if config.is_asymmetric:
                # Use asymmetric pairing for different left/right models
                logging.info(
                    f'Using asymmetric pairing: {config.left_model} + {config.right_model}'
                )
                hitIndex, paired, unpaired = tirmite.iterateGetPairsAsymmetric(
                    hitIndex, config, stableReps=args.stableReps
                )
            else:
                # Use enhanced symmetric pairing that supports custom orientations
                logging.info(
                    f'Using symmetric pairing with orientation {config.orientation}'
                )
                hitIndex, paired, unpaired = tirmite.iterateGetPairsCustom(
                    hitIndex, config, stableReps=args.stableReps
                )

            logging.debug(f'After pairing - paired type: {type(paired)}')
            logging.debug(f'After pairing - unpaired type: {type(unpaired)}')

        except Exception as e:
            logging.error(f'Error in pairing: {e} (type: {type(e)})')
            raise

        # Log pairing results
        try:
            total_pairs = sum(len(pairs) for pairs in paired.values())
            total_unpaired = len(unpaired)
            logging.info(
                f'Pairing completed: {total_pairs} pairs formed, {total_unpaired} hits remain unpaired'
            )
        except Exception as e:
            logging.error(f'Error calculating pairing results: {e} (type: {type(e)})')
            logging.debug(f'paired: {paired}')
            logging.debug(f'unpaired: {unpaired}')
            raise

        # Write TIR hits to fasta for each pHMM
        logging.info('Writing all valid terminus hits to fasta.')
        try:
            tirmite.writeTIRs(
                outDir=outDir,
                hitTable=hitTable,
                maxeval=args.maxeval,
                genome=genome,
                prefix=args.prefix,
                padlen=args.padlen,
                genome_descriptions=genome_descriptions,  # Add this
            )
        except Exception as e:
            logging.error(f'Error in writeTIRs: {e} (type: {type(e)})')
            raise

        # Write paired TIR hits to fasta
        if args.report in ['all', 'paired']:
            logging.info('Writing successfully paired termini to fasta.')
            try:
                tirmite.writePairedTIRs(
                    outDir=outDir,
                    paired=paired,
                    hitIndex=hitIndex,
                    genome=genome,
                    prefix=args.prefix,
                    padlen=args.padlen,
                    genome_descriptions=genome_descriptions,  # Add this
                )
            except Exception as e:
                logging.error(f'Error in writePairedTIRs: {e} (type: {type(e)})')
                raise

        # Extract paired hit regions
        try:
            pairedEles = tirmite.fetchElements(
                paired=paired,
                hitIndex=hitIndex,
                genome=genome,
                genome_descriptions=genome_descriptions,  # Add this
            )
        except Exception as e:
            logging.error(f'Error in fetchElements: {e} (type: {type(e)})')
            raise

        # Debug logging
        total_elements = sum(len(elements) for elements in pairedEles.values())
        logging.info(
            f'Created {total_elements} paired elements from {sum(len(pairs) for pairs in paired.values())} pairs'
        )
        for model, elements in pairedEles.items():
            logging.info(f'Model {model}: {len(elements)} elements')

        # Write paired-TIR features to fasta
        logging.info('Writing paired-termini elements to fasta.')
        tirmite.writeElements(outDir, eleDict=pairedEles, prefix=args.prefix)

        # Write paired features to gff3, optionally also report paired/unpaired TIRs
        if args.gffOut:
            # Get unpaired TIRs
            if args.report in ['all', 'unpaired']:
                unpairedTIRs = tirmite.fetchUnpaired(hitIndex=hitIndex)
                logging.info(f'Found {len(unpairedTIRs)} unpaired termini')
            else:
                unpairedTIRs = None

            # Set gff file path
            if args.prefix:
                gffOutPath = os.path.join(outDir, args.prefix + '.gff3')
            else:
                gffOutPath = os.path.join(outDir, 'tirmite_report.gff3')

            # Write gff3
            logging.info('Writing features to gff: %s ' % gffOutPath)
            tirmite.gffWrite(
                outpath=gffOutPath,
                featureList=pairedEles,
                writeTIRs=args.report,
                unpaired=unpairedTIRs,
                prefix=args.prefix,
            )

        logging.info('TIRmite analysis completed successfully')

    except KeyboardInterrupt:
        logging.info('Analysis interrupted by user')
        sys.exit(130)  # Standard exit code for SIGINT

    except Exception as e:
        logging.error(f'Unexpected error during analysis: {e}')
        sys.exit(1)

    finally:
        # Always clean up temporary directory
        cleanup_temp_directory(tempDir, args.keep_temp)


if __name__ == '__main__':
    main()
