import argparse
import logging
import os
import shutil
import sys

from Bio import SeqIO
from pymummer import coords_file, nucmer  # alignment

from tirmite.logs import init_logging
from tirmite.runBlastn import makeBlast
from tirmite.utils import checkUniqueID, cleanID, getTimestring, manageTemp
from tirmite.wrapping import run_cmd


def mainArgs():
    parser = argparse.ArgumentParser(
        description="Extract terminal repeats from DNA transposons (TIRs). \
                    Optionally, compose synthetic MITES from complete DNA transposons.",
        prog="tsplit",
    )
    parser.add_argument(
        "--loglevel",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging level. Default: 'DEBUG'",
    )
    parser.add_argument(
        "-i",
        "--infile",
        type=str,
        required=True,
        default=None,
        help="Multifasta containing complete elements.",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        default=None,
        help="All output files begin with this string. \
                                        (Default: [infile name])",
    )
    parser.add_argument(
        "-d",
        "--outdir",
        type=str,
        default=None,
        help="Write output files to this directory. \
                                        (Default: cwd)",
    )
    parser.add_argument(
        "--splitmode",
        default="split",
        choices=["all", "split", "internal", "external", None],
        help='all= Report input sequence as well as internal and external segments. \
                                        split= Report internal and external segments after splitting. \
                                        internal = Report only internal segments \
                                        external = Report only terminal repeat segments. \
                                        If set to "None" then only synthetic MITES will be reported if --makemites is also set. \
                                        (Default: split)',
    )
    parser.add_argument(
        "--makemites",
        action="store_true",
        default=False,
        help="Attempt to construct synthetic MITE sequences from TIRs.",
    )
    parser.add_argument(
        "--keeptemp",
        action="store_true",
        default=False,
        help="If set do not remove temp directory on completion.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="If set, report progress.",
    )
    parser.add_argument(
        "-m",
        "--maxdist",
        type=int,
        default=2,
        help="Terminal repeat candidates must be no more than this many bases from end of input element. (Default: 2)\
                                        Note: Increase this value if you suspect that your element is nested within some flanking sequence.",
    )
    parser.add_argument(
        "--minid",
        type=float,
        default=80.0,
        help="Minimum identity between terminal repeat pairs. As float. \
                                        (Default: 80.0)",
    )
    parser.add_argument(
        "--minterm",
        type=int,
        default=10,
        help='Minimum length for a terminal repeat to be considered. \
                                        Equivalent to nucmer "--mincluster" \
                                        (Default: 10)',
    )
    parser.add_argument(
        "--minseed",
        type=int,
        default=5,
        help='Minimum length of a maximal exact match to be included in final match cluster. \
                                        Equivalent to nucmer "--minmatch". \
                                        (Default: 5)',
    )
    parser.add_argument(
        "--diagfactor",
        type=float,
        default=0.2,
        help="Maximum diagonal difference factor for clustering of matches within nucmer, \
                                        i.e. diagonal difference / match separation (default 0.20) \
                                        Note: Increase value for greater tolerance of indels between terminal repeats.",
    )
    parser.add_argument(
        "--method",
        default="nucmer",
        choices=["blastn", "nucmer"],
        help='Select alignment method: "blastn" or "nucmer".(Default: nucmer)',
    )
    args = parser.parse_args()
    return args


def missing_tool(tool_name):
    path = shutil.which(tool_name)
    if path is None:
        return [tool_name]
    else:
        return []


def tSplitchecks(args):
    """Housekeeping tasks: Create output files/dirs
    and temp dirs as required."""
    if not os.path.isfile(args.infile):
        print("Input sequence file does not exist. Quitting.")
        sys.exit(1)
    # Make outDir if does not exist else set to current dir.
    if args.outdir:
        absOutDir = os.path.abspath(args.outdir)
        if not os.path.isdir(args.outdir):
            os.makedirs(absOutDir)
        outDir = absOutDir
    else:
        outDir = os.getcwd()
    # Make temp directory
    tempDir = os.path.join(os.getcwd(), "temp_" + getTimestring())
    os.makedirs(tempDir)
    # Set prefix to infile basename if none
    if not args.prefix:
        prefix = os.path.splitext(os.path.basename(args.infile))[0]
    else:
        prefix = args.prefix
    # Create outfile paths
    outfile = prefix + "_tsplit_output.fasta"
    outpath = os.path.join(outDir, outfile)
    # Return full path to output file and temp directory
    return outpath, tempDir


## Fix: Do not load fasta into memory!
def importFasta2List(file):
    """Load elements from multifasta file. Check that seq IDs are unique."""
    # Read in elements from multifasta file, convert seqrecord iterator to list
    records = list(SeqIO.parse(file, "fasta"))
    # Check names are unique
    checkUniqueID(records)
    # If unique, return record list.
    return records


def getTIRs(
    elements=None,
    flankdist=2,
    minid=80,
    minterm=10,
    minseed=5,
    diagfactor=0.3,
    mites=False,
    report="split",
    temp=None,
    keeptemp=False,
    alignTool="nucmer",
    verbose=False,
):
    """
    Align elements to self and attempt to identify TIRs.
    Optionally attempt to construct synthetic MITEs from TIRs.
    """
    # Set temp directory to cwd if none.
    if not temp:
        temp = os.getcwd()
    # For each candidate TIR element
    for rec in elements:
        # Create temp paths for single element fasta and alignment coords
        tempFasta = os.path.join(temp, cleanID(rec.id) + ".fasta")
        tempCoords = tempFasta + ".coords"
        # Write current element to single fasta
        manageTemp(record=rec, tempPath=tempFasta, scrub=False)
        # Align to self with nucmer
        if alignTool == "nucmer":
            # Compose Nucmer script for current element vs self
            runner = nucmer.Runner(
                tempFasta,
                tempFasta,
                tempCoords,
                min_id=minid,
                min_length=minseed,
                diagfactor=diagfactor,
                mincluster=minterm,
                breaklen=200,
                maxmatch=True,
                simplify=False,
            )
            # Execute nucmer
            runner.run()
        elif alignTool == "blastn":
            # Alternatively, use blastn as search tool and write nucmer.coords-like output.
            cmds = makeBlast(seq=tempFasta, outfile=tempCoords, pid=minid)
            run_cmd(cmds, verbose=verbose)
        # Import coords file to iterator object
        file_reader = coords_file.reader(tempCoords)
        # Exclude hits to self. Also converts iterator output to stable list
        alignments = [hit for hit in file_reader if not hit.is_self_hit()]
        # Filter hits less than min length (Done internally for nucmer, not blastn.)
        alignments = [
            hit for hit in alignments if hit.ref_end - hit.ref_start >= minterm
        ]
        # Filter for hits on same strand i.e. tandem repeats / LTRs
        alignments = [hit for hit in alignments if not hit.on_same_strand()]
        # Filter for 5' repeats which begin within x bases of element start
        alignments = [hit for hit in alignments if hit.ref_start <= flankdist]
        # Scrub overlappying ref / query segments, and also complementary
        # 3' to 5' flank hits
        alignments = [hit for hit in alignments if hit.ref_end < hit.qry_end]
        # Sort largest to smallest dist between end of ref (subject) and start
        # of query (hit)
        # x.qry_end - x.ref_end = 5'end of right TIR - 3' end of left
        # TIR = length of internal segment
        # TIR pair with smallest internal segment (longest TIRs) is first in list.
        alignments = sorted(
            alignments, key=lambda x: (x.qry_end - x.ref_end), reverse=False
        )
        # If alignments exist after filtering report features using alignment
        # pair with largest internal segment i.e. first element in sorted list.
        if alignments:
            if verbose:
                [print(x) for x in alignments]
            if report in ["split", "external", "all"]:
                # yield TIR slice - append "_TIR"
                extSeg = rec[alignments[0].ref_start : alignments[0].ref_end + 1]
                extSeg.id = extSeg.id + "_TIR"
                extSeg.name = extSeg.id
                extSeg.description = "[" + rec.id + " TIR segment]"
                yield extSeg
            if report in ["split", "internal", "all"]:
                # yield internal slice - append "_I"
                intSeg = rec[alignments[0].ref_end : alignments[0].qry_end + 1]
                intSeg.id = intSeg.id + "_I"
                intSeg.name = intSeg.id
                intSeg.description = "[" + rec.id + " internal segment]"
                yield intSeg
            if report == "all":
                yield rec
            if mites:
                # Assemble TIRs into hypothetical MITEs
                synMITE = (
                    rec[alignments[0].ref_start : alignments[0].ref_end + 1]
                    + rec[alignments[0].qry_end : alignments[0].qry_start + 1]
                )
                synMITE.id = synMITE.id + "_synMITE"
                synMITE.name = synMITE.id
                synMITE.description = (
                    "[Synthetic MITE constructed from " + rec.id + " TIRs]"
                )
                yield synMITE
        else:
            # If alignment list empty after filtering print alert and continue
            print("No TIRs found for candidate element: %s" % rec.id)
        # Scrub single fasta and coords file for current element.
        if not keeptemp:
            manageTemp(tempPath=tempFasta, scrub=True)
            manageTemp(tempPath=tempCoords, scrub=True)


## Fix: Do not load fasta into memory!
def segWrite(outfile, segs=None):
    """
    Take a generator object yielding seqrecords and
    write each to outfile in fasta format.
    """
    seqcount = 0
    if segs:
        with open(outfile, "w") as handle:
            for seq in segs:
                seqcount += 1
                SeqIO.write(seq, handle, "fasta")
        if seqcount == 0:
            os.remove(outfile)


def main():
    """Do the work."""

    # Get cmd line args
    args = mainArgs()

    # Set up logging
    init_logging(loglevel=args.loglevel)

    # Check for required programs.
    tools = ["delta-filter", "nucmer", "show-coords"]
    if args.method == "blastn":
        tools.append("blastn")
    missing_tools = []
    for tool in tools:
        missing_tools += missing_tool(tool)
    if missing_tools:
        logging.warning(
            (
                "WARNING: Some tools required by Exterminate could not be found: "
                + ", ".join(missing_tools)
            )
        )
        logging.warning("You may need to install them to use all features.")

    # Create output paths as required
    outpath, tempdir = tSplitchecks(args)

    # Load elements to be screened
    ## Fix: Do not load fasta into memory!
    elements = importFasta2List(args.infile)

    # If TIR mode or if makeMITE mode enabled, search for inverted terminal repeats
    # Optionally construct synthetic MITE from TIRs if discovered
    segments = getTIRs(
        elements=elements,
        flankdist=args.maxdist,
        minterm=args.minterm,
        minseed=args.minseed,
        minid=args.minid,
        diagfactor=args.diagfactor,
        mites=args.makemites,
        report=args.splitmode,
        temp=tempdir,
        alignTool=args.method,
        keeptemp=args.keeptemp,
        verbose=args.verbose,
    )
    segWrite(outpath, segs=segments)

    # Remove temp directory
    if not args.keeptemp:
        shutil.rmtree(tempdir)


if __name__ == "__main__":
    main()
