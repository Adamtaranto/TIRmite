import tirmite
import argparse
import shutil
import glob
import sys
import os

def mainArgs():
	'''Parse command line arguments.'''
	parser = argparse.ArgumentParser(
							description	=	'Map TIR-pHMM models to genomic sequences for annotation of MITES and complete DNA-Transposons.',
							prog		=	'tirmite'
							)
	# Input
	parser.add_argument('--genome',type=str,required=True,help='Path to target genome that will be queried with HMMs.')
	parser.add_argument('--hmmDir',type=str,default=None,help='Directory containing pre-prepared TIR-pHMMs.')
	parser.add_argument('--hmmFile',type=str,default=None,help='Path to single TIR-pHMM file. Incompatible with "--hmmDir".')
	parser.add_argument('--alnDir',type=str,default=None,help='Path to directory containing only TIR alignments to be converted to HMM.')
	parser.add_argument('--alnFile',type=str,default=None,help='Provide a single TIR alignment to be converted to HMM. Incompatible with "--alnDir".')
	parser.add_argument('--alnFormat',default='fasta',choices=["clustal","emboss","fasta","fasta-m10","ig","maf","mauve","nexus","phylip","phylip-sequential","phylip-relaxed","stockholm"],
						help='Alignments provided with "--alnDir" or "--alnFile" are all in this format.') 
	# Pairing heuristics
	parser.add_argument('--stableReps',type=int,default=0,help='Number of times to iterate pairing procedure when no additional pairs are found AND remaining unpaired hits > 0.')
	# Output and housekeeping
	parser.add_argument('--outdir',type=str,default=None,help='All output files will be written to this directory.')
	parser.add_argument('--prefix',type=str,default=None,help='Add prefix to all TIRs and Paired elements detected in this run. Useful when running same TIR-pHMM against many genomes.(Default = None)')
	parser.add_argument('--nopairing',action='store_true',default=False,help='If set, only report TIR-pHMM hits. Do not attempt pairing.')
	parser.add_argument('--gffOut',type=str,default=None,help='GFF3 annotation filename. Do not write annotations if not set.')
	parser.add_argument('--reportTIR',default='all',choices=[None,'all','paired','unpaired'],help='Options for reporting TIRs in GFF annotation file.') 
	parser.add_argument('--keeptemp',action='store_true',default=False,help='If set do not delete temp file directory.')
	parser.add_argument('-v','--verbose',action='store_true',default=False,help='Set syscall reporting to verbose.')
	# HMMER options
	parser.add_argument('--cores',type=int,default=1,help='Set number of cores available to hmmer software.')
	parser.add_argument('--maxeval',type=float,default=0.001,help='Maximum e-value allowed for valid hit. Default = 0.001')
	parser.add_argument('--maxdist',type=int,default=None,help='Maximum distance allowed between TIR candidates to consider valid pairing.')
	parser.add_argument('--nobias',action='store_true',default=False,help='Turn OFF bias correction of scores in nhmmer.')
	parser.add_argument('--matrix',type=str,default=None,help='Use custom DNA substitution matrix with nhmmer.')
	# Non-standard HMMER paths
	parser.add_argument('--hmmpress',type=str,default='hmmpress',help='Set location of hmmpress if not in path.')
	parser.add_argument('--nhmmer',type=str,default='nhmmer',help='Set location of nhmmer if not in path.')
	parser.add_argument('--hmmbuild',type=str,default='hmmbuild',help='Set location of hmmbuild if not in path.')
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
	tools = [args.hmmpress,args.nhmmer,args.hmmbuild]
	missing_tools = []
	for tool in tools:
	    missing_tools += missing_tool(tool)
	if missing_tools:
	    print('WARNING: Some tools required by tirmite could not be found: ' +
	          ', '.join(missing_tools))
	    print('You may need to install them to use all features.')

	# Create output and temp paths as required
	outDir,tempDir = tirmite.dochecks(args)

	# Load reference genome
	genome = tirmite.importFasta(args.genome)

	# If raw alignments provided, convert to stockholm format.
	if args.alnDir or args.alnFile:
		stockholmDir = tirmite.convertAlign(alnDir=args.alnDir,alnFile=args.alnFile,inFormat=args.alnFormat,tempDir=tempDir)
	else:
		stockholmDir = None

	# Compose and run HMMER commands
	cmds,resultDir = tirmite.cmdScript(hmmDir=args.hmmDir, hmmFile=args.hmmFile, alnDir=stockholmDir, tempDir=tempDir, args=args)
	tirmite.run_cmd(cmds,verbose=args.verbose)

	# Die if no hits found
	if not glob.glob(os.path.join(os.path.abspath(resultDir),'*.tab')):
		print("No hits found in %s . Quitting." % resultDir)
		sys.exit(1)

	# Import hits from nhmmer result files
	hitTable = None
	for resultfile in glob.glob(os.path.join(os.path.abspath(resultDir),'*.tab')):
		hitTable = tirmite.import_nhmmer(infile=resultfile,hitTable=hitTable,prefix=args.prefix)

	# Group hits by model and chromosome (hitsDict), and initiate hit tracker hitIndex to manage pair-searching
	hitsDict,hitIndex = tirmite.table2dict(hitTable)

	# If pairing is off, just report the hits
	if args.nopairing:
		tirmite.writeTIRs(outDir=outDir, hitTable=hitTable, maxeval=args.maxeval, genome=genome)
		# Remove temp directory
		if not args.keeptemp:
			shutil.rmtree(tempDir)
		sys.exit(1)

	# Populate hitIndex with acceptible candidate partners (compatible strand and distance.)
	hitIndex = tirmite.parseHits(hitsDict=hitsDict, hitIndex=hitIndex, maxDist=args.maxdist)

	# Run iterative pairing procedure
	hitIndex,paired,unpaired = tirmite.iterateGetPairs(hitIndex, stableReps=args.stableReps)

	# Write TIR hits to fasta for each pHMM
	tirmite.writeTIRs(outDir=outDir, hitTable=hitTable, maxeval=args.maxeval, genome=genome)

	# Extract paired hit regions (candidate TEs / MITEs)
	pairedEles = tirmite.fetchElements(paired=paired, hitIndex=hitIndex, genome=genome)

	# Write paired TIR features to fasta
	tirmite.writeElements(outDir, eleDict=pairedEles)

	# Write paired features to gff3, optionally also report paired/unpaired TIRs
	if args.gffOut:
		# Get unpaired TIRs
		if args.reportTIR in ['all','unpaired']:
			unpairedTIRs = tirmite.fetchUnpaired(hitIndex=hitIndex)
		else:
			unpairedTIRs = None
	# Write gff3
		tirmite.gffWrite(outpath=args.gffOut, featureList=pairedEles, writeTIRs=args.reportTIR, unpaired=unpairedTIRs)

	# Remove temp directory
	if not args.keeptemp:
		shutil.rmtree(tempDir)

