import argparse
import glob
import logging
import os
import shutil
import sys

from tirmite._version import __version__
from tirmite.hmmer_wrappers import cmdScriptHMMER
from tirmite.logs import init_logging
import tirmite.tirmitetools as tirmite
from tirmite.utils import dochecks, importFasta
from tirmite.wrapping import run_cmd


def mainArgs():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Map TIR-pHMM models to genomic sequences for annotation \
        of MITES and complete DNA-Transposons.",
        prog="tirmite",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )
    parser.add_argument(
        "--loglevel",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging level. Default: 'DEBUG'",
    )
    # Input
    parser.add_argument(
        "--genome",
        type=str,
        required=True,
        help="Path to target genome that will be queried with HMMs.",
    )
    parser.add_argument(
        "--hmmDir",
        type=str,
        default=None,
        help="Directory containing pre-prepared TIR-pHMMs.",
    )
    parser.add_argument(
        "--hmmFile",
        type=str,
        default=None,
        help='Path to single TIR-pHMM file. Incompatible with "--hmmDir".',
    )
    parser.add_argument(
        "--alnDir",
        type=str,
        default=None,
        help="Path to directory containing only TIR alignments to be converted to HMM.",
    )
    parser.add_argument(
        "--alnFile",
        type=str,
        default=None,
        help='Provide a single TIR alignment to be converted to HMM. Incompatible with "--alnDir".',
    )
    parser.add_argument(
        "--alnFormat",
        default="fasta",
        choices=["clustal", "fasta", "nexus", "phylip", "stockholm"],
        help='Alignments provided with "--alnDir" or "--alnFile" are all in this format.',
    )
    parser.add_argument(
        "--pairbed",
        type=str,
        default=None,
        help="If set TIRmite will preform pairing on TIRs from custom bedfile only.",
    )

    parser.add_argument(
        "--stableReps",
        type=int,
        default=0,
        help="Number of times to iterate pairing procedure when no additional pairs are found AND remaining unpaired hits > 0.",
    )
    # Output and housekeeping
    parser.add_argument(
        "--outdir",
        type=str,
        default=None,
        help="All output files will be written to this directory.",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default=None,
        help="Add prefix to all TIRs and Paired elements detected in this run. Useful when running same TIR-pHMM against many genomes.(Default = None)",
    )
    parser.add_argument(
        "--nopairing",
        action="store_true",
        default=False,
        help="If set, only report TIR-pHMM hits. Do not attempt pairing.",
    )
    parser.add_argument(
        "--gffOut",
        action="store_true",
        default=False,
        help="If set report features as prefix.gff3. File saved to outdir. Default: False",
    )
    parser.add_argument(
        "--reportTIR",
        default="all",
        choices=[None, "all", "paired", "unpaired"],
        help="Options for reporting TIRs in GFF annotation file.",
    )
    parser.add_argument(
        "--padlen",
        type=int,
        default=None,
        help="Extract x bases either side of TIR when writing TIRs to fasta.",
    )
    parser.add_argument(
        "--keeptemp",
        action="store_true",
        default=False,
        help="If set do not delete temp file directory.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Set syscall reporting to verbose.",
    )
    # HMMER options
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="Set number of cores available to hmmer software.",
    )
    parser.add_argument(
        "--maxeval",
        type=float,
        default=0.001,
        help="Maximum e-value allowed for valid hit. Default = 0.001",
    )
    parser.add_argument(
        "--maxdist",
        type=int,
        default=None,
        help="Maximum distance allowed between TIR candidates to consider valid pairing.",
    )
    parser.add_argument(
        "--nobias",
        action="store_true",
        default=False,
        help="Turn OFF bias correction of scores in nhmmer.",
    )
    parser.add_argument(
        "--matrix",
        type=str,
        default=None,
        help="Use custom DNA substitution matrix with nhmmer.",
    )
    parser.add_argument(
        "--mincov",
        type=float,
        default=0.5,
        help="Minimum valid hit length as prop of model length. Defaults to 0.5",
    )
    # Non-standard HMMER paths
    parser.add_argument(
        "--hmmpress",
        type=str,
        default="hmmpress",
        help="Set location of hmmpress if not in PATH.",
    )
    parser.add_argument(
        "--nhmmer",
        type=str,
        default="nhmmer",
        help="Set location of nhmmer if not in PATH.",
    )
    parser.add_argument(
        "--hmmbuild",
        type=str,
        default="hmmbuild",
        help="Set location of hmmbuild if not in PATH.",
    )
    args = parser.parse_args()
    return args


def missing_tool(tool_name):
    path = shutil.which(tool_name)
    if path is None:
        return [tool_name]
    else:
        return []


def main():
    """Do the work."""
    # Get cmd line args
    args = mainArgs()

    # Set up logging
    init_logging(loglevel=args.loglevel)

    # Check for required programs.
    tools = [args.hmmpress, args.nhmmer, args.hmmbuild]

    missing_tools = []
    for tool in tools:
        missing_tools += missing_tool(tool)
    if missing_tools:
        logging.warning(
            "Some tools required by tirmite could not be found: "
            + ", ".join(missing_tools)
        )
        logging.warning("You may need to install them to use all features.")

    # Create output and temp paths as required
    outDir, tempDir = dochecks(args)

    # Load reference genome
    logging.info("Loading genome from: %s " % args.genome)
    genome = importFasta(args.genome)

    # Import custom TIR hits from BEDfile.
    if args.pairbed:
        # Die if no input file
        if not glob.glob(os.path.abspath(args.pairbed)):
            logging.warning("BED file %s not found. Quitting." % args.pairbed)
            # Remove temp directory
            if not args.keeptemp:
                shutil.rmtree(tempDir)
            sys.exit(1)

        logging.info(
            "Skipping HMM search. Using custom TIRs from file: %s" % args.pairbed
        )

        # Import hits from BED file
        # Format: Chrm, start, end, name, evalue, strand
        hitTable = None
        logging.info("Loading custom TIR hits from: %s" % str(args.pairbed))
        hitTable = tirmite.import_BED(
            infile=args.pairbed, hitTable=hitTable, prefix=args.prefix
        )

        # Apply hit e-value filters
        logging.info("Filtering hits with e-value > %s" % str(args.maxeval))
        hitCount = len(hitTable.index)
        hitTable = tirmite.filterHitsEval(maxeval=args.maxeval, hitTable=hitTable)
        logging.info(
            "Excluded %s hits on e-value criteria."
            % str(hitCount - len(hitTable.index))
        )
        logging.info("Remaining hits: %s " % str(len(hitTable.index)))

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
            logging.info("Pairing is off. Reporting hits only.")
            # Remove temp directory
            if not args.keeptemp:
                shutil.rmtree(tempDir)
            sys.exit(1)

    # Else run nhmmer and load TIR hits.
    else:
        # If raw alignments provided, convert to stockholm format.
        if args.alnDir or args.alnFile:
            stockholmDir = tirmite.convertAlign(
                alnDir=args.alnDir,
                alnFile=args.alnFile,
                inFormat=args.alnFormat,
                tempDir=tempDir,
            )
        else:
            stockholmDir = None

        # If pre-built HMM provided, check correct format.
        if args.hmmFile:
            if os.path.splitext(os.path.basename(args.hmmFile))[1].lstrip(".") != "hmm":
                logging.warning("--hmmFile has non-hmm extension. Exiting.")
                # Remove temp directory
                if not args.keeptemp:
                    shutil.rmtree(tempDir)
                sys.exit(1)

        # Compose and run HMMER commands
        cmds, resultDir, hmmDB = cmdScriptHMMER(
            hmmDir=args.hmmDir,
            hmmFile=args.hmmFile,
            alnDir=stockholmDir,
            tempDir=tempDir,
            args=args,
        )
        run_cmd(cmds, verbose=args.verbose, tempDir=tempDir, keeptemp=args.keeptemp)

        # Die if no hits found
        if not glob.glob(os.path.join(os.path.abspath(resultDir), "*.tab")):
            logging.warning("No hits found in %s . Quitting." % resultDir)
            # Remove temp directory
            if not args.keeptemp:
                shutil.rmtree(tempDir)
            sys.exit(1)

        # Import hits from nhmmer result files
        hitTable = None
        modelCount = 0
        for resultfile in glob.glob(os.path.join(os.path.abspath(resultDir), "*.tab")):
            logging.info("Loading nhmmer hits from: %s " % resultfile)
            hitTable = tirmite.import_nhmmer(
                infile=resultfile, hitTable=hitTable, prefix=args.prefix
            )
            modelCount += 1

        logging.info(
            "Imported %s hits from %s models. "
            % (str(len(hitTable.index)), str(modelCount))
        )

        # Apply hit length filters
        logging.info("Filtering hits with < %s model coverage. " % str(args.mincov))
        hitCount = len(hitTable.index)
        hitTable = tirmite.filterHitsLen(
            hmmDB=hmmDB, mincov=args.mincov, hitTable=hitTable
        )
        logging.info(
            "Excluded %s hits on coverage criteria. "
            % str(hitCount - len(hitTable.index))
        )
        logging.info("Remaining hits: %s " % str(len(hitTable.index)))

        # Apply hit e-value filters
        logging.info("Filtering hits with e-value > %s" % str(args.maxeval))
        hitCount = len(hitTable.index)
        hitTable = tirmite.filterHitsEval(maxeval=args.maxeval, hitTable=hitTable)
        logging.info(
            "Excluded %s hits on e-value criteria."
            % str(hitCount - len(hitTable.index))
        )
        logging.info("Remaining hits: %s " % str(len(hitTable.index)))

        # Group hits by model and chromosome (hitsDict), and initiate hit tracker hitIndex to manage pair-searching
        hitsDict, hitIndex = tirmite.table2dict(hitTable)

        # If pairing is off, just report the hits
        if args.nopairing:
            tirmite.writeTIRs(
                outDir=outDir,
                hitTable=hitTable,
                maxeval=args.maxeval,
                genome=genome,
                padlen=args.padlen,
            )
            logging.info("Pairing is off. Reporting hits only.")
            # Remove temp directory
            if not args.keeptemp:
                shutil.rmtree(tempDir)
            sys.exit(1)

    # Run pairing on filtered TIR set

    # Populate hitIndex with acceptible candidate partners (compatible strand and distance.)
    hitIndex = tirmite.parseHits(
        hitsDict=hitsDict, hitIndex=hitIndex, maxDist=args.maxdist
    )

    # Run iterative pairing procedure
    hitIndex, paired, unpaired = tirmite.iterateGetPairs(
        hitIndex, stableReps=args.stableReps
    )

    # Write TIR hits to fasta for each pHMM
    logging.info("Writing all valid TIR hits to fasta.")
    tirmite.writeTIRs(
        outDir=outDir,
        hitTable=hitTable,
        maxeval=args.maxeval,
        genome=genome,
        prefix=args.prefix,
        padlen=args.padlen,
    )

    # Write paired TIR hits to fasta. Pairs named as element ID + L/R tag.
    if args.reportTIR in ["all", "paired"]:
        logging.info("Writing successfully paired TIRs to fasta.")
        tirmite.writePairedTIRs(
            outDir=outDir,
            paired=paired,
            hitIndex=hitIndex,
            genome=genome,
            prefix=args.prefix,
            padlen=args.padlen,
        )

    # Extract paired hit regions (candidate TEs / MITEs) elements are stored
    # as list of gffTup objects
    pairedEles = tirmite.fetchElements(paired=paired, hitIndex=hitIndex, genome=genome)

    # Write paired-TIR features to fasta
    logging.info("Writing TIR-elements to fasta.")
    tirmite.writeElements(outDir, eleDict=pairedEles, prefix=args.prefix)

    # Write paired features to gff3, optionally also report paired/unpaired TIRs
    if args.gffOut:
        # Get unpaired TIRs
        if args.reportTIR in ["all", "unpaired"]:
            # Unpaired TIR features are stored as list of gffTup objects
            unpairedTIRs = tirmite.fetchUnpaired(hitIndex=hitIndex)
        else:
            unpairedTIRs = None

        # Set gff file path
        if args.prefix:
            gffOutPath = os.path.join(outDir, args.prefix + ".gff3")
        else:
            gffOutPath = os.path.join(outDir, "tirmite_report.gff3")
        # Write gff3
        logging.info("Writing features to gff: %s " % gffOutPath)
        tirmite.gffWrite(
            outpath=gffOutPath,
            featureList=pairedEles,
            writeTIRs=args.reportTIR,
            unpaired=unpairedTIRs,
            prefix=args.prefix,
        )

    # Remove temp directory
    if not args.keeptemp:
        shutil.rmtree(tempDir)
