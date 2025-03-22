import glob
import os
import re
from shlex import quote
import shutil

from tirmite.wrapping import run_cmd


def cleanID(s):
    """
    Remove non alphanumeric characters from string.
    Replace whitespace with underscores.
    """
    s = re.sub(r"[^\w\s]", "", s)
    s = re.sub(r"\s+", "_", s)
    return s


def _hmmbuild_command(
    exePath="hmmbuild", modelname=None, cores=None, inAlign=None, outdir=None
):
    """
    Construct the hmmbuild command.
    """
    # Make model name compliant
    modelname = cleanID(modelname)
    # Check for outdir
    if outdir:
        outdir = os.path.abspath(outdir)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
    else:
        outdir = os.getcwd()
    # Create hmm database dir
    outdir = os.path.join(outdir, "hmmDB")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # Base command
    command = quote(exePath) + " --dna -n " + modelname
    # Optional set cores
    if cores:
        command += " --cpu " + str(cores)
    if outdir:
        modelout = os.path.join(outdir, modelname + ".hmm")
    else:
        modelout = modelname + ".hmm"
    # Append output file name and source alignment, return command string
    command = " ".join(
        [command, quote(os.path.abspath(modelout)), quote(os.path.abspath(inAlign))]
    )
    return command, os.path.abspath(modelout)


def _hmmpress_command(exePath="hmmpress", hmmfile=None):
    """
    Construct the hmmbuild command.
    """
    # Base command
    command = quote(exePath) + " -f " + quote(os.path.abspath(hmmfile))
    return command


def _nhmmer_command(
    exePath="nhmmer",
    modelPath=None,
    genome=None,
    evalue=None,
    nobias=False,
    matrix=None,
    cores=None,
    outdir=None,
):
    """
    Construct the nhmmer command
    """
    # Get model hmm basename
    model_base = os.path.splitext(os.path.basename(modelPath))[0]
    # Check for outdir
    if outdir:
        outdir = os.path.abspath(outdir)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
    else:
        outdir = os.getcwd()
    # Make subdir for nhmmer tab results
    outdir = os.path.join(outdir, "nhmmer_results")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # Create resultfile name
    outfile = os.path.join(os.path.abspath(outdir), model_base + ".tab")
    command = quote(exePath) + " --tblout " + quote(outfile)
    if cores:
        command += " --cpu " + str(cores)
    if evalue:
        command += " -E " + str(evalue)
    if nobias:
        command += " --nobias"
    if matrix:
        command += " --mxfile " + quote(os.path.abspath(matrix))
    command += (
        " --noali --notextw --dna --max "
        + quote(os.path.abspath(modelPath))
        + " "
        + quote(os.path.abspath(genome))
    )
    return command, outdir


def cmdScriptHMMER(hmmDir=None, hmmFile=None, alnDir=None, tempDir=None, args=None):
    if tempDir:
        tempDir = os.path.abspath(tempDir)
        if not os.path.isdir(tempDir):
            os.makedirs(tempDir)
    else:
        tempDir = os.getcwd()
    # Create hmm database dir
    hmmDB = os.path.join(tempDir, "hmmDB")
    if not os.path.isdir(hmmDB):
        os.makedirs(hmmDB)
    # Copy all existing hmms to hmmDB
    if hmmDir:
        for hmm in glob.glob(os.path.join(hmmDir, "*.hmm")):
            shutil.copy2(hmm, hmmDB)
    if hmmFile:
        shutil.copy2(hmmFile, hmmDB)
    # Write cmds to list
    cmds = list()
    # Process alignments
    if alnDir:
        build_cmds = list()
        for alignment in glob.glob(os.path.join(alnDir, "*")):
            modelName = os.path.splitext(os.path.basename(alignment))[0]
            hmmbuildCmd, modelPath = _hmmbuild_command(
                exePath=args.hmmbuild,
                modelname=modelName,
                cores=args.cores,
                inAlign=alignment,
                outdir=tempDir,
            )
            build_cmds.append(hmmbuildCmd)
        run_cmd(build_cmds, verbose=args.verbose, keeptemp=args.keeptemp)
    # Press and write nhmmer cmd for all models in hmmDB directory
    for hmm in glob.glob(os.path.join(hmmDB, "*.hmm")):
        hmmPressCmd = _hmmpress_command(exePath=args.hmmpress, hmmfile=hmm)
        nhmmerCmd, resultDir = _nhmmer_command(
            exePath=args.nhmmer,
            nobias=args.nobias,
            matrix=args.matrix,
            modelPath=hmm,
            genome=args.genome,
            evalue=args.maxeval,
            cores=args.cores,
            outdir=tempDir,
        )
        cmds.append(hmmPressCmd)
        cmds.append(nhmmerCmd)
    # Return list of cmds and location of file result files
    return cmds, resultDir, hmmDB
