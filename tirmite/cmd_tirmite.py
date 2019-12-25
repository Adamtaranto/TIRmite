#!/usr/bin/env python

from __future__ import print_function
from tirmite import __version__
import os
import sys
import glob
import shutil
import tirmite
import argparse


def log(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

def mainArgs():
	'''Parse command line arguments.'''
	parser = argparse.ArgumentParser(
							description	=	'Map TIR-pHMM models to genomic sequences for annotation of MITES and complete DNA-Transposons.',
							prog		=	'tirmite'
							)
	parser.add_argument('--version', action='version',version='%(prog)s {version}'.format(version=__version__))
	# Input
	parser.add_argument('--genome',type=str,required=True,help='Path to target genome that will be queried with HMMs.')
	parser.add_argument('--hmmDir',type=str,default=None,help='Directory containing pre-prepared TIR-pHMMs.')
	parser.add_argument('--hmmFile',type=str,default=None,help='Path to single TIR-pHMM file. Incompatible with "--hmmDir".')
	parser.add_argument('--alnDir',type=str,default=None,help='Path to directory containing only TIR alignments to be converted to HMM.')
	parser.add_argument('--alnFile',type=str,default=None,help='Provide a single TIR alignment to be converted to HMM. Incompatible with "--alnDir".')
	parser.add_argument('--alnFormat',default='fasta',choices=["clustal","fasta","nexus","phylip","stockholm"],
						help='Alignments provided with "--alnDir" or "--alnFile" are all in this format.') 
	parser.add_argument('--pairbed',type=str,default=None,help='If set TIRmite will preform pairing on TIRs from custom bedfile only.')
	# Alternative search method
	#parser.add_argument('--useBowtie2',action='store_true',default=False,help='If set, map short TIR to genome with bowtie2. Potentially useful for very short though highly conserved TIRs where TIR-pHMM hits return high e-values.')
	#parser.add_argument('--btTIR',type=str,default=None,help='Fasta file containing a single TIR to be mapped with bowtie2.')
	#parser.add_argument('--bowtie2',type=str,default='bowtie2',help='Set location of bowtie2 if not in PATH.')
	#parser.add_argument('--bt2build',type=str,default='bowtie2-build',help='Set location of bowtie2-build if not in PATH.')
	#parser.add_argument('--samtools',type=str,default='samtools',help='Set location of samtools if not in PATH.')
	#parser.add_argument('--bedtools',type=str,default='bedtools',help='Set location of bedtools if not in PATH.')
	# Pairing heuristics
	parser.add_argument('--stableReps',type=int,default=0,help='Number of times to iterate pairing procedure when no additional pairs are found AND remaining unpaired hits > 0.')
	# Output and housekeeping
	parser.add_argument('--outdir',type=str,default=None,help='All output files will be written to this directory.')
	parser.add_argument('--prefix',type=str,default=None,help='Add prefix to all TIRs and Paired elements detected in this run. Useful when running same TIR-pHMM against many genomes.(Default = None)')
	parser.add_argument('--nopairing',action='store_true',default=False,help='If set, only report TIR-pHMM hits. Do not attempt pairing.')
	parser.add_argument('--gffOut',action='store_true',default=False,help='If set report features as prefix.gff3. File saved to outdir. Default: False')
	parser.add_argument('--reportTIR',default='all',choices=[None,'all','paired','unpaired'],help='Options for reporting TIRs in GFF annotation file.') 
	parser.add_argument('--padlen',type=int,default=None,help='Extract x bases either side of TIR when writing TIRs to fasta.')
	parser.add_argument('--keeptemp',action='store_true',default=False,help='If set do not delete temp file directory.')
	parser.add_argument('-v','--verbose',action='store_true',default=False,help='Set syscall reporting to verbose.')
	# HMMER options
	parser.add_argument('--cores',type=int,default=1,help='Set number of cores available to hmmer software.')
	parser.add_argument('--maxeval',type=float,default=0.001,help='Maximum e-value allowed for valid hit. Default = 0.001')
	parser.add_argument('--maxdist',type=int,default=None,help='Maximum distance allowed between TIR candidates to consider valid pairing.')
	parser.add_argument('--nobias',action='store_true',default=False,help='Turn OFF bias correction of scores in nhmmer.')
	parser.add_argument('--matrix',type=str,default=None,help='Use custom DNA substitution matrix with nhmmer.')
	parser.add_argument('--mincov',type=float,default=0.5,help='Minimum valid hit length as prop of model length. Defaults to 0.5')
	# Non-standard HMMER paths
	parser.add_argument('--hmmpress',type=str,default='hmmpress',help='Set location of hmmpress if not in PATH.')
	parser.add_argument('--nhmmer',type=str,default='nhmmer',help='Set location of nhmmer if not in PATH.')
	parser.add_argument('--hmmbuild',type=str,default='hmmbuild',help='Set location of hmmbuild if not in PATH.')
	args = parser.parse_args()
	return args

def missing_tool(tool_name):
    path = shutil.which(tool_name)
    if path is None:
        return [tool_name]
    else:
        return []

def main():
	'''Do the work.'''
	# Get cmd line args
	args = mainArgs()

	# Check for required programs.
	#tools = [args.hmmpress,args.nhmmer,args.hmmbuild,args.bowtie2,args.bt2build,args.samtools,args.bedtools]
	tools = [args.hmmpress,args.nhmmer,args.hmmbuild]
	
	missing_tools = []
	for tool in tools:
	    missing_tools += missing_tool(tool)
	if missing_tools:
	    log('WARNING: Some tools required by tirmite could not be found: ' +
	          ', '.join(missing_tools))
	    log('You may need to install them to use all features.')

	# Create output and temp paths as required
	outDir,tempDir = tirmite.dochecks(args)

	# Load reference genome
	log("Log: Loading genome from: %s " % args.genome)
	genome = tirmite.importFasta(args.genome)

	#if args.useBowtie2:
	#	# Check that input fasta exists
	#	tirmite.isfile(args.btTIR)
	#	btTIRname = tirmite.getbtName(args.btTIR)
	#	# Compose bowtie map and filter commands
	#	cmds = list()
	#	cmds.append(tirmite._bowtie2build_cmd(bt2Path=args.bt2build,genome=args.genome))
	#	cmds.append(tirmite._bowtie2_cmd(bt2Path=args.bowtie2,tirFasta=args.btTIR,cores=args.cores))
	#	bam2bed_cmds,mappedPath = tirmite._bam2bed_cmd(samPath=args.samtools,bedPath=args.bedtools,tempDir=tempDir)
	#	cmds += bam2bed_cmds
	#	# Run mapp and filter
	#	tirmite.run_cmd(cmds,verbose=args.verbose,keeptemp=args.keeptemp)
	#	# Import mapping locations
	#	hitTable = tirmite.import_mapped(infile=mappedPath,tirName=btTIRname,prefix=args.prefix)
	#else:

	# Import custom TIR hits from BEDfile.
	if args.pairbed:
		# Die if no input file
		if not glob.glob(os.path.abspath(args.pairbed)):
			log("WARNING: BED file %s not found. Quitting." % args.pairbed)
			# Remove temp directory
			if not args.keeptemp:
				shutil.rmtree(tempDir)
			sys.exit(1)

		log("Log: Skipping HMM search. Using custom TIRs from file: %s" %args.pairbed)
		
		# Import hits from BED file
		# Format: Chrm, start, end, name, evalue, strand
		hitTable = None
		log("Log: Loading custom TIR hits from: %s" % str(args.pairbed))
		hitTable = tirmite.import_BED(infile=args.pairbed,hitTable=hitTable,prefix=args.prefix)

		# Apply hit e-value filters
		log("Log: Filtering hits with e-value > %s" % str(args.maxeval))
		hitCount = len(hitTable.index)
		hitTable = tirmite.filterHitsEval(maxeval=args.maxeval, hitTable=hitTable)
		log("Log: Excluded %s hits on e-value criteria." % str(hitCount - len(hitTable.index)))
		log("Log: Remaining hits: %s " % str(len(hitTable.index)))

		# Group hits by model and chromosome (hitsDict), and initiate hit tracker hitIndex to manage pair-searching
		hitsDict,hitIndex = tirmite.table2dict(hitTable)

		# If pairing is off, just report the hits.
		if args.nopairing:
			tirmite.writeTIRs(outDir=outDir, hitTable=hitTable, maxeval=args.maxeval, genome=genome, padlen=args.padlen)
			log("Log: Pairing is off. Reporting hits only.")
			# Remove temp directory
			if not args.keeptemp:
				shutil.rmtree(tempDir)
			sys.exit(1)

	# Else run nhmmer and load TIR hits.
	else:
		# If raw alignments provided, convert to stockholm format.
		if args.alnDir or args.alnFile:
			stockholmDir = tirmite.convertAlign(alnDir=args.alnDir,alnFile=args.alnFile,inFormat=args.alnFormat,tempDir=tempDir)
		else:
			stockholmDir = None
		
		# If pre-built HMM provided, check correct format.
		if args.hmmFile:
			if os.path.splitext(os.path.basename(args.hmmFile))[1].lstrip('.') != 'hmm':
				log("WARNING: --hmmFile has non-hmm extension. Exiting.")
				# Remove temp directory
				if not args.keeptemp:
					shutil.rmtree(tempDir)
				sys.exit(1)

		# Compose and run HMMER commands
		cmds,resultDir,hmmDB = tirmite.cmdScript(hmmDir=args.hmmDir, hmmFile=args.hmmFile, alnDir=stockholmDir, tempDir=tempDir, args=args)
		tirmite.run_cmd(cmds, verbose=args.verbose, tempDir=tempDir, keeptemp=args.keeptemp)

		# Die if no hits found
		if not glob.glob(os.path.join(os.path.abspath(resultDir),'*.tab')):
			log("Log: No hits found in %s . Quitting." % resultDir)
			# Remove temp directory
			if not args.keeptemp:
				shutil.rmtree(tempDir)
			sys.exit(1)

		# Import hits from nhmmer result files
		hitTable = None
		modelCount = 0
		for resultfile in glob.glob(os.path.join(os.path.abspath(resultDir),'*.tab')):
			log("Log: Loading nhmmer hits from: %s " % resultfile)
			hitTable = tirmite.import_nhmmer(infile=resultfile,hitTable=hitTable,prefix=args.prefix)
			modelCount += 1

		log("Log: Imported %s hits from %s models. " % (str(len(hitTable.index)),str(modelCount)))
		
		# Apply hit length filters
		log("Log: Filtering hits with < %s model coverage. " % str(args.mincov))
		hitCount = len(hitTable.index)
		hitTable = tirmite.filterHitsLen(hmmDB=hmmDB, mincov=args.mincov, hitTable=hitTable)
		log("Log: Excluded %s hits on coverage criteria. " % str(hitCount - len(hitTable.index)))
		log("Log: Remaining hits: %s " % str(len(hitTable.index)))

		# Apply hit e-value filters
		log("Log: Filtering hits with e-value > %s" % str(args.maxeval))
		hitCount = len(hitTable.index)
		hitTable = tirmite.filterHitsEval(maxeval=args.maxeval, hitTable=hitTable)
		log("Log: Excluded %s hits on e-value criteria." % str(hitCount - len(hitTable.index)))
		log("Log: Remaining hits: %s " % str(len(hitTable.index)))

		# Group hits by model and chromosome (hitsDict), and initiate hit tracker hitIndex to manage pair-searching
		hitsDict,hitIndex = tirmite.table2dict(hitTable)

		# If pairing is off, just report the hits
		if args.nopairing:
			tirmite.writeTIRs(outDir=outDir, hitTable=hitTable, maxeval=args.maxeval, genome=genome, padlen=args.padlen)
			log("Log: Pairing is off. Reporting hits only.")
			# Remove temp directory
			if not args.keeptemp:
				shutil.rmtree(tempDir)
			sys.exit(1)

	# Run pairing on filtered TIR set

	# Populate hitIndex with acceptible candidate partners (compatible strand and distance.)
	hitIndex = tirmite.parseHits(hitsDict=hitsDict, hitIndex=hitIndex, maxDist=args.maxdist)

	# Run iterative pairing procedure
	hitIndex,paired,unpaired = tirmite.iterateGetPairs(hitIndex, stableReps=args.stableReps)

	# Write TIR hits to fasta for each pHMM
	log("Log: Writing all valid TIR hits to fasta.")
	tirmite.writeTIRs(outDir=outDir, hitTable=hitTable, maxeval=args.maxeval, genome=genome, prefix=args.prefix, padlen=args.padlen)

	# Write paired TIR hits to fasta. Pairs named as element ID + L/R tag.
	if args.reportTIR in ['all','paired']:
		log("Log: Writing successfully paired TIRs to fasta.")
		tirmite.writePairedTIRs(outDir=outDir, paired=paired, hitIndex=hitIndex, genome=genome, prefix=args.prefix, padlen=args.padlen)

	# Extract paired hit regions (candidate TEs / MITEs) elements are stored as list of gffTup objects
	pairedEles = tirmite.fetchElements(paired=paired, hitIndex=hitIndex, genome=genome)

	# Write paired-TIR features to fasta
	log("Log: Writing TIR-elements to fasta." )
	tirmite.writeElements(outDir, eleDict=pairedEles, prefix=args.prefix)

	# Write paired features to gff3, optionally also report paired/unpaired TIRs
	if args.gffOut:
		# Get unpaired TIRs
		if args.reportTIR in ['all','unpaired']:
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
		log("Log: Writing features to gff: %s " % gffOutPath)
		tirmite.gffWrite(outpath=gffOutPath, featureList=pairedEles, writeTIRs=args.reportTIR, unpaired=unpairedTIRs, prefix=args.prefix)

	# Remove temp directory
	if not args.keeptemp:
		shutil.rmtree(tempDir)
