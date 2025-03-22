from collections import Counter
from datetime import datetime
import os
import re
import sys

from Bio import SeqIO


def dochecks(args):
    """
    Housekeeping tasks:
    Create output files/dirs and temp dirs as required.
    """
    # Check that specified files exist
    isfile(args.genome)
    if args.hmmFile:
        isfile(args.hmmFile)
    if args.alnFile:
        isfile(args.alnFile)
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
    # Return full path to output and temp directories
    return outDir, tempDir


def getTimestring():
    """
    Return int only string of current datetime with milliseconds.
    """
    (dt, micro) = datetime.utcnow().strftime("%Y%m%d%H%M%S.%f").split(".")
    dt = "%s%03d" % (dt, int(micro) / 1000)
    return dt


def isfile(path):
    if not os.path.isfile(path):
        print("Input file not found: %s" % path)
        sys.exit(1)


def cleanID(s):
    """
    Remove non alphanumeric characters from string.
    Replace whitespace with underscores.
    """
    s = re.sub(r"[^\w\s]", "", s)
    s = re.sub(r"\s+", "_", s)
    return s


# Fix: Do not load fasta into genome!
def manageTemp(record=None, tempPath=None, scrub=False):
    """
    Create single sequence fasta files or scrub temp files.
    """
    if scrub and tempPath:
        try:
            os.remove(tempPath)
        except OSError:
            pass
    else:
        with open(tempPath, "w") as f:
            SeqIO.write(record, f, "fasta")


# Fix: Do not load fasta into genome!
def checkUniqueID(records):
    """
    Check that IDs for input elements are unique.
    """
    seqIDs = [records[x].id for x in range(len(records))]
    IDcounts = Counter(seqIDs)
    duplicates = [k for k, v in IDcounts.items() if v > 1]
    if duplicates:
        print("Input sequence IDs not unique. Quiting.")
        print(duplicates)
        sys.exit(1)
    else:
        pass


# Fix: Do not load fasta into genome!
def importFasta(file):
    """
    Load elements from multifasta file. Check that seq IDs are unique.
    """
    # Read in elements from multifasta file, convert seqrecord iterator to list
    records = list(SeqIO.parse(file, "fasta"))
    # Check names are unique
    checkUniqueID(records)
    # If unique, return records as dict keyed by seq id
    recordsDict = dict()
    for rec in records:
        recordsDict[rec.id] = rec
    return recordsDict
