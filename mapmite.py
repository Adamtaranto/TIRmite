#!/usr/bin/env python
#python 3
#mapmite.py
#Version 1.0 Adam Taranto, August 2017
#Contact, Adam Taranto, adam.taranto@anu.edu.au

########################################################################
# Map TIR-pHMM models to genomic sequences for annotation of MITES and #
# complete DNA-Transposons.                                            #
########################################################################

import os
import re
import sys
import glob
import shutil
import argparse
import tempfile
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from datetime import datetime
from collections import Counter
from collections import namedtuple
from operator import attrgetter

class Error (Exception): pass

def decode(x):
	try:
		s = x.decode()
	except:
		return x
	return s

def dochecks(args):
	"""Housekeeping tasks: Create output files/dirs and temp dirs as required."""
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
	tempDir = os.path.join(os.getcwd(),"temp_" + getTimestring())
	os.makedirs(tempDir)
	# Return full path to output and temp directories
	return outDir,tempDir

def isfile(path):
	if not os.path.isfile(path):
		print("Input file not found: %s" % path)
		sys.exit(1)

def cleanID(s):
	"""Remove non alphanumeric characters from string. Replace whitespace with underscores."""
	s = re.sub(r"[^\w\s]", '', s)
	s = re.sub(r"\s+", '_', s)
	return s

def getTimestring():
	"""Return int only string of current datetime with milliseconds."""
	(dt, micro) = datetime.utcnow().strftime('%Y%m%d%H%M%S.%f').split('.')
	dt = "%s%03d" % (dt, int(micro) / 1000)
	return dt

def checkUniqueID(records):
	"""Check that IDs for input elements are unique."""
	seqIDs = [records[x].id for x in range(len(records))]
	IDcounts = Counter(seqIDs)
	duplicates = [k for k, v in IDcounts.items() if v > 1]
	if duplicates:
		print("Input sequence IDs not unique. Quiting.")
		print(duplicates)
		sys.exit(1)
	else:
		pass

def manageTemp(record=None, tempPath=None, scrub=False):
	"""Create single sequence fasta files or scrub temp files."""
	if scrub and tempPath:
		try:
			os.remove(tempPath)
		except OSError:
			pass
	else:
		with open(tempPath, "w") as f:
			SeqIO.write(record, f, "fasta")

def importFasta(file):
	"""Load elements from multifasta file. Check that seq IDs are unique."""
	# Read in elements from multifasta file, convert seqrecord iterator to list.
	records = list(SeqIO.parse(file, "fasta"))
	# Check names are unique
	checkUniqueID(records)
	# If unique, return records as dict keyed by seq id
	recordsDict = dict()
	for rec in records:
		recordsDict[rec.id] = rec
	return recordsDict

def convertAlign(alnDir=None,alnFile=None,inFormat='fasta',tempDir=None):
	''' Convert input alignments into Stockholm format.'''
	# Construct out model path
	if tempDir:
		alnOutDir = os.path.join(os.path.abspath(tempDir), 'temp_aln')
	else:
		alnOutDir = os.path.join(os.getcwd(), 'temp_aln')
	# Create if does not exist
	if not os.path.isdir(alnOutDir):
		os.makedirs(alnOutDir)
	# Get list of alignment files to process
	if alnFile:
		alignments = [alnFile]
	else:
		alignments = glob.glob(alnDir)
	# Do conversion on each
	for infile in alignments:
		# Get basename
		inBase = os.path.splitext(os.path.basename(infile))[0]
		# Make outpath
		outAln = os.path.join(alnOutDir, inBase + '.stockholm')
		# Open files
		input_handle = open(infile, "rU")
		output_handle = open(outAln, "w")
		# Read alignment
		alignments = AlignIO.parse(input_handle, inFormat)
		# Write as stockholm
		AlignIO.write(alignments, output_handle, "stockholm")
		# Close handles
		output_handle.close()
		input_handle.close()
	return alnOutDir

def _hmmbuild_command(exePath="hmmbuild",modelname=None,cores=None,inAlign=None,outdir=None):
	'''Construct the hmmbuild command'''
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
	outdir = os.path.join(outdir,"hmmDB")
	if not os.path.isdir(outdir):
		os.makedirs(outdir)
	# Base command
	command = exePath + " --dna -n " + modelname
	# Optional set cores
	if cores:
		command += ' --cpu ' + str(cores)
	if outdir:
		modelout = os.path.join(outdir,modelname + ".hmm ")
	else:
		modelout = modelname + ".hmm "
	# Append output file name and source alignment, return command string
	command += ' ' + os.path.abspath(modelout) + os.path.abspath(inAlign)
	return command,os.path.abspath(modelout)

def _hmmpress_command(exePath="hmmpress", hmmfile=None):
	'''Construct the hmmbuild command'''
	# Base command
	command = exePath + " -f " + os.path.abspath(hmmfile)
	return command

def _nhmmer_command(exePath="nhmmer",modelPath=None,genome=None,evalue=None,nobias=False,matrix=None,cores=None,outdir=None):
	'''Construct the nhmmer command'''
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
	outdir = os.path.join(outdir,"nhmmer_results")
	if not os.path.isdir(outdir):
		os.makedirs(outdir)
	# Create resultfile name
	outfile = os.path.join(os.path.abspath(outdir),model_base + ".tab")
	command = exePath + " --tblout " + outfile
	if cores:
			command += ' --cpu ' + str(cores)
	if evalue:
			command += ' -E ' + str(evalue)
	if nobias:
		command += ' --nobias'
	if matrix:
		command += ' --mxfile ' + os.path.abspath(matrix)
	command += " --noali --notextw --dna --max " + os.path.abspath(modelPath) + " " + os.path.abspath(genome)
	return command,outdir

def _write_script(cmds,script):
	'''Write commands into a bash script'''
	f = open(script, 'w+')
	for cmd in cmds:
		print(cmd, file=f)
	f.close()

def cmdScript(hmmDir=None, hmmFile=None, alnDir=None, tempDir=None, args=None):
	if tempDir:
		tempDir = os.path.abspath(tempDir)
		if not os.path.isdir(tempDir):
			os.makedirs(tempDir)
	else:
		tempDir = os.getcwd()
	# Create hmm database dir
	hmmDB = os.path.join(tempDir,"hmmDB")
	if not os.path.isdir(hmmDB):
		os.makedirs(hmmDB)
	# Copy all existing hmms to hmmDB
	if hmmDir:
		for hmm in glob.glob(os.path.join(hmmDir,'*.hmm')):
			shutil.copy2(hmm,hmmDB)
	if hmmFile:
		shutil.copy2(hmmFile,hmmDB)
	# Write cmds to list
	cmds = list()
	# Process alignments 
	if alnDir:
		build_cmds = list()
		for alignment in glob.glob(os.path.join(alnDir,'*')):
			modelName = os.path.splitext(os.path.basename(alignment))[0]
			hmmbuildCmd,modelPath = _hmmbuild_command(exePath=args.hmmbuild,modelname=modelName,cores=args.cores,inAlign=alignment,outdir=tempDir)
			build_cmds.append(hmmbuildCmd)
		run_cmd(build_cmds,verbose=args.verbose)
	# Press and write nhmmer cmd for all models in hmmDB directory
	for hmm in glob.glob(os.path.join(hmmDB,'*.hmm')):
		hmmPressCmd = _hmmpress_command(exePath=args.hmmpress, hmmfile=hmm)
		nhmmerCmd,resultDir = _nhmmer_command(exePath=args.nhmmer,nobias=args.nobias,matrix=args.matrix,modelPath=hmm,genome=args.genome,evalue=args.maxeval,cores=args.cores,outdir=tempDir)
		cmds.append(hmmPressCmd)
		cmds.append(nhmmerCmd)
	# Return list of cmds and location of file result files
	return cmds,resultDir

def syscall(cmd, verbose=False):
	'''Manage error handling when making syscalls'''
	if verbose:
		print('Running command:', cmd, flush=True)
	try:
		output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
	except subprocess.CalledProcessError as error:
		print('The following command failed with exit code', error.returncode, file=sys.stderr)
		print(cmd, file=sys.stderr)
		print('\nThe output was:\n', file=sys.stderr)
		print(error.output.decode(), file=sys.stderr)
		raise Error('Error running command:', cmd)
	if verbose:
		print(decode(output))

def run_cmd(cmds,verbose=False):
	'''Write and excute HMMER script'''
	tmpdir = tempfile.mkdtemp(prefix='tmp.', dir=os.getcwd())
	original_dir = os.getcwd()
	os.chdir(tmpdir)
	script = 'run_jobs.sh'
	_write_script(cmds,script)
	syscall('bash ' + script, verbose=verbose)
	os.chdir(original_dir)
	shutil.rmtree(tmpdir)

def import_nhmmer(infile=None,hitTable=None,prefix=None):
	''' Read nhmmer tab files to pandas dataframe.'''
	hitRecords = list()
	with open(infile, "rU") as f:
		for line in f.readlines():
			li = line.strip()
			if not li.startswith("#"):
				li = li.split()
				if li[11] == '+':
					hitRecords.append({
					'target':li[0],
					'model':li[2],
					'hmmStart':li[4],
					'hmmEnd':li[5],
					'hitStart':li[6],
					'hitEnd':li[7],
					'strand':li[11],
					'evalue':li[12],
					'score':li[13],
					'bias':li[14]})
				elif li[11] == '-':
					hitRecords.append({
					'target':li[0],
					'model':li[2],
					'hmmStart':li[4],
					'hmmEnd':li[5],
					'hitStart':li[7],
					'hitEnd':li[6],
					'strand':li[11],
					'evalue':li[12],
					'score':li[13],
					'bias':li[14]})
	# Convert list of dicts into dataframe
	df = pd.DataFrame(hitRecords)
	# Reorder columns
	cols = ['model','target','hitStart','hitEnd','strand','evalue','score','bias','hmmStart','hmmEnd']
	df = df.ix[:, cols]
	if hitTable is not None:
		# If an existing table was passed, concatenate
		df = pd.concat([df,hitTable], ignore_index=True)
	# Sort hits by HMM, Chromosome, location, and strand
	df = df.sort_values(['model','target','hitStart','hitEnd','strand'], ascending=[True,True,True,True,True])
	# Reindex
	df = df.reset_index(drop=True)
	if prefix:
		df['model'] = str(prefix) + df['model'].astype(str)
	return df

def table2dict(hitTable):
	''' Convert pandas dataframe of nhmmer hits into dict[model][chrom] 
		and index[model] = [hitlist].withCandidates and pairing status'''
	# Set up empty dict
	hitsDict = dict()
	hitIndex = dict()
	# Populate keys from dataframe
	for hmm in hitTable.model.unique():
		hitsDict[hmm] = dict()
		hitIndex[hmm] = dict()
		for chr in hitTable[hitTable['model'] == hmm].target.unique():
			hitsDict[hmm][chr] = list()
	# Set up named tuple
	hitTup = namedtuple('Elem', ['model','target','hitStart','hitEnd','strand','idx'])
	# Add each record to dicts
	for row in hitTable.iterrows():
		record = hitTup(row[1].model,row[1].target,int(row[1].hitStart),int(row[1].hitEnd),row[1].strand,row[0])
		# Log hit for model on chromosome
		hitsDict[row[1].model][row[1].target].append(record)
		# Populate tracker
		hitIndex[hmm][row[0]] = {'rec':record,'partner':None,'candidates':list()}
	# Return master rec object and pairing tracker
	return hitsDict,hitIndex

def parseHits(hitsDict=None, hitIndex=None, maxDist=None):
	'''Populate hitIndex with pairing candidates'''
	if not maxDist:
		maxDist = float('inf')
	for hmm in hitIndex.keys():
		for UID in hitIndex[hmm].keys():
			ref = hitIndex[hmm][UID]['rec']
			if ref.strand == '+':
				for localhit in hitsDict[ref.model][ref.target]:
					if localhit.strand == '-' and localhit.hitStart >= ref.hitEnd and localhit.hitStart - ref.hitEnd <= maxDist:
						hitIndex[hmm][UID]['candidates'].append(localhit)
				# Sort candidate hit records from low to high on hitStart values
				hitIndex[hmm][UID]['candidates'] = sorted(hitIndex[hmm][UID]['candidates'], key=attrgetter('hitStart', 'hitEnd'))
			if ref.strand == '-':
				for localhit in hitsDict[ref.model][ref.target]:
					if localhit.strand == '+' and localhit.hitEnd <= ref.hitStart and ref.hitStart - localhit.hitEnd <= maxDist:
						hitIndex[hmm][UID]['candidates'].append(localhit)
				# Sort candidate hit records from high to low on hitEnd values
				hitIndex[hmm][UID]['candidates'] = sorted(hitIndex[hmm][UID]['candidates'], key=attrgetter('hitEnd','hitStart'), reverse=True)
	# hitIndex[model][idx].keys() == [rec,candidates,partner]
	return hitIndex

def isfirstUnpaired(ref=None, mate=None, model=None, index=None):
	'''Provided with a hitID (ref) and the ID of its nearest unpaired candidate partner 
	(mate), check if ref is also the nearest unpaired partner to mate.'''
	# Init result trackers
	found = None
	mateFUP = None
	# Scan candidate partners of 'mate' looking for ref
	for matePartner in index[model][mate]['candidates']:
		# If unpaired candidate is ref
		if not index[model][matePartner.idx]['partner'] and matePartner.idx == ref:
			found = set([matePartner.idx,mate])
			index[model][ref]['partner'] = mate
			index[model][mate]['partner'] = ref
			# return
			return found,index,mateFUP
		# If first unpaired candidate partner for mate is not ref
		elif not index[model][matePartner.idx]['partner']:
			# Return none and unchanged index
			mateFUP = matePartner.idx
			return found,index,mateFUP
		else:
			continue
	# If mate candidates include no unpaired reps, return unchanged index 
	# and blank 'found' and 'mateFUP' variables.
	return found,index,mateFUP

def getPairs(hitIndex=None, paired=None):
	'''Loop over all hit for all models and search for reciprocity within 
		2 degrees of the top unpaired candidate for each unpaired hit.'''
	# If pair tracker not given
	if not paired:
		# Create dict of empty lists, keyed by model name
		paired = dict()
		for model in hitIndex.keys():
			paired[model] = list()
	# For each HMM model
	for model in hitIndex.keys():
		# Ask each hit in genome
		for refID in hitIndex[model].keys():
			# If it has been asigned a partner
			if not hitIndex[model][refID]['partner']:
				# If not partnered, start checking candidate partners
				for Can1 in hitIndex[model][refID]['candidates']:
					# For a candidate that is also unpartnered
					if not hitIndex[model][Can1.idx]['partner']:
						# Check if unpartnered candidate is a reciprocal match for our hit
						found,hitIndex,mateFUP = isfirstUnpaired(ref=refID, mate=Can1.idx, model=model, index=hitIndex)
						if found:
							# If current hit is also the best return match of our candidate, store as pair.
							paired[model].append(found)
						elif mateFUP:
							# Else if not a return match, check candidate's first outbound match for reciprocity.
							found,hitIndex,mateFUP = isfirstUnpaired(ref=Can1.idx, mate=mateFUP, model=model, index=hitIndex)
							if found:
								# Store if found.
								paired[model].append(found)
	return hitIndex,paired

def countUnpaired(hitIndex):
	'''How many hits are still unpaired across all models.'''
	count = 0
	for model in hitIndex.keys():
		for hitID in hitIndex[model].keys():
			if not hitIndex[model][hitID]['partner']:
				count += 1
	return count

def listunpaired(hitIndex):
	'''Return list of all unpaired hit IDs'''
	unpaired = list()
	for model in hitIndex.keys():
		for hitID in hitIndex[model].keys():
			if not hitIndex[model][hitID]['partner']:
				unpaired.append(hitID)
	return unpaired

def iterateGetPairs(hitIndex, stableReps=0):
	''' Iterate pairing procedure for all models until no unpaired hits remain or
		number of reps without change is exceeded.'''
	# Init stable repeat counter
	reps = 0
	# Run initial pairing
	hitIndex,paired = getPairs(hitIndex=hitIndex)
	# Count remaining unpaired hits
	countUP = countUnpaired(hitIndex)
	# Iterate pairing procedure until either no unpaired remain 
	# OR max number of interations without new pairing is reached
	while countUP > 0 and reps < stableReps:
		# Re-run pairing procedure
		hitIndex,paired = getPairs(hitIndex=hitIndex, paired=paired)
		# Store previous unpaired hit count
		lastCountUP = countUP
		# Update unpaired hit count
		countUP = countUnpaired(hitIndex)
		# If no change in upaired hit count, iterate stable rep counter
		if lastCountUP == countUP:
			reps += 1
	# Get IDs of remaining unpaired hits
	unpaired = listunpaired(hitIndex)
	# Return results
	return hitIndex,paired,unpaired

def extractTIRs(model=None, hitTable=None, maxeval=0.001, genome=None):
	''' For significant hits in model, compose seqrecords.'''
	hitcount = 0
	seqList = list()
	for index,row in hitTable[hitTable['model'] == model].iterrows():
		if float(row['evalue']) <= maxeval:
			hitcount += 1
			hitrecord = genome[row['target']][int(row['hitStart'])-1:int(row['hitEnd'])]
			hitrecord.id = model + '_' + str(index)
			if row['strand'] == '-':
				hitrecord = hitrecord.reverse_complement(id=hitrecord.id + "_rc")
			hitrecord.name = hitrecord.id 
			hitrecord.description = '_'.join([row['target'],row['strand'],row['hitStart'],row['hitEnd'],'modelAlignment:' + row['hmmStart'],row['hmmEnd']])
			# Append record to list
			seqList.append(hitrecord)
		else:
			continue
	# Return seqrecord list and total hit count for model
	return seqList,hitcount 

def writeTIRs(outDir=None, hitTable=None, maxeval=0.001, genome=None):
	''' Write all hits per Model to a multifasta in the outdir'''
	if outDir:
		outDir = os.path.abspath(outDir)
		if not os.path.isdir(outDir):
				os.makedirs(outDir)
	else:
		outDir = os.getcwd()
	for model in hitTable['model'].unique():
		# List of TIR seqrecords, and count of hits
		seqList,hitcount = extractTIRs(model=model, hitTable=hitTable, maxeval=maxeval, genome=genome)
		outfile = os.path.join(outDir,model + "_hits_" + str(hitcount) + ".fasta")
		# Write extracted hits to model outfile
		with open(outfile, "w") as handle:
			for seq in seqList:
				SeqIO.write(seq, handle, "fasta")

def flipTIRs(x,y):
	''' Sort hits into left and right TIRs.'''
	left2right = sorted([x,y], key=attrgetter('hitStart', 'hitEnd'))
	return (left2right[0],left2right[1])

def fetchElements(paired=None, hitIndex=None, genome=None):
	''' Extract complete sequence of paired elements, 
		asign names and child TIRs for use in seq and GFF reporting.'''
	TIRelements = dict()
	gffTup = namedtuple('gffElem', ['model', 'chromosome', 'start', 'end', 'strand', 'type', 'id', 'leftHit' , 'rightHit','seq'])
	for model in paired.keys():
		TIRelements[model] = list()
		model_counter = 0
		for x,y in paired[model]:
			model_counter += 1
			x = hitIndex[model][x]['rec']
			y = hitIndex[model][y]['rec']
			leftHit,rightHit = flipTIRs(x,y)
			eleID = model + "_Element_" + str(model_counter)
			eleSeq = genome[leftHit.target][int(leftHit.hitStart)-1:int(rightHit.hitEnd)]
			eleSeq.id = eleID
			eleSeq.name = eleID
			eleSeq.description = '_'.join([leftHit.target,str(leftHit.hitStart),str(rightHit.hitEnd)]) + " len=" + str(rightHit.hitEnd-leftHit.hitStart)
			TIRelement = gffTup(model, leftHit.target , leftHit.hitStart , rightHit.hitEnd , leftHit.strand , "TIR_Element", eleID, leftHit,rightHit, eleSeq)
			TIRelements[model].append(TIRelement)
	# Return list of element info tuples
	return TIRelements

def writeElements(outDir, eleDict=None):
	''' Takes dict of extracted sequences keyed by model. Writes to fasta by model.'''
	for model in eleDict.keys():
		outfile = os.path.join(outDir,model + '_elements.fasta')
		with open(outfile, "w") as handle:
			for seq in eleDict[model]:
				SeqIO.write(seq.seq, handle, "fasta")

def fetchUnpaired(hitIndex=None):
	'''	Take list of unpaired hit IDs from listunpaired(),
		Compose TIR gff3 record. '''
	orphans = list()
	gffTup = namedtuple('gffElem', ['model', 'chromosome', 'start', 'end', 'strand', 'type', 'id', 'leftHit' , 'rightHit', 'seq'])
	for model in hitIndex.keys():
		for recID in hitIndex[model].keys():
			if not hitIndex[model][recID]['partner']:
				x = hitIndex[model][recID]['rec']
				orphan = gffTup(x.model, x.target , x.hitStart , x.hitEnd , x.strand , "orphan_TIR", x.idx, None , None, None)
				orphans.append(orphan)
	return orphans

def	gffWrite(outpath=None, featureList=list(), writeTIRs=True, unpaired=None):
	'''	Write predicted paired-TIR features (i.e. MITEs) from fetchElements() as GFF3.
		Optionally, write child TIRS and orphan TIRs to GFF3 also.'''
	# If path to output gff3 file not provided, set default location.
	if not outpath:
		outpath = os.path.join(os.getcwd(),"mapmite_features.gff3")
	# Unpack element dict to list
	allfeatures = list()
	for model in featureList.keys():
		for record in featureList[model]:
			allfeatures.append(record)
	# Add list of unpaired TIRs to main featureList if provided.
	if unpaired:
		allfeatures = allfeatures + unpaired
	# Sort features
	sortedFeatures = sorted(allfeatures, key=attrgetter('model', 'chromosome', 'start', 'end'))
	# Open GFF handle
	with open(outpath, 'w') as file:
		# Write headers
		file.write('##gff-version 3' + '\n')
		file.write('\t'.join(['#seqid','source','type','start','end','score','strand','phase','attributes']) + '\n')
		# Format features for GFF3
		for x in sortedFeatures:
			if x.type == "orphan_TIR" and writeTIRs in ['all','unpaired']:
				file.write('\t'.join(	[str(x.chromosome),"mapmite",x.type,str(x.start),str(x.end),".",x.strand,".",
										'ID=' + str(x.model) + "_" + str(x.id) + ';model=' + str(x.model) + ';']
										) + '\n'
										)
			if x.type == "TIR_Element":
				# Write Element line
				file.write('\t'.join(	[str(x.chromosome),"mapmite",x.type,str(x.start),str(x.end),".",x.strand,".",
										'ID=' + str(x.id) + ';model=' + str(x.model) + ';']
										) + '\n'
										)
				if writeTIRs in ['all','paired']:
					# Write left TIR line as child
					l = x.leftHit
					file.write('\t'.join(	[str(l.target),"mapmite","paired_TIR",str(l.hitStart),str(l.hitEnd),".",l.strand,".",
											'ID=' + str(x.model) + "_" + str(l.idx) + ';model=' + str(x.model) + ';Parent=' + str(x.id) + ';']
											) + '\n'
											)
					# Write right TIR line as child on neg strand
					r = x.rightHit
					file.write('\t'.join(	[str(r.target),"mapmite","paired_TIR",str(r.hitStart),str(r.hitEnd),".",r.strand,".",
											'ID=' + str(x.model) + "_" + str(r.idx) + ';model=' + str(x.model) + ';Parent=' + str(x.id) + ';']
											) + '\n'
											)


#gffTup fields: 'model', 'chromosome', 'start', 'end', 'strand', 'type', 'id', 'score','bias', 'evalue', 'leftHit' , 'rightHit', 'eleSeq'
#Types: "TIR_Element", "orphan_TIR"


def main(args):
	'''Do the work.'''
	# Create output and temp paths as required
	outDir,tempDir = dochecks(args)

	# Load reference genome
	genome = importFasta(args.genome)

	# If raw alignments provided, convert to stockholm format.
	if args.alnDir or args.alnFile:
		stockholmDir = convertAlign(alnDir=args.alnDir,alnFile=args.alnFile,inFormat=args.alnFormat,tempDir=tempDir)
	else:
		stockholmDir = None

	# Compose and run HMMER commands
	cmds,resultDir = cmdScript(hmmDir=args.hmmDir, hmmFile=args.hmmFile, alnDir=stockholmDir, tempDir=tempDir, args=args)
	run_cmd(cmds,verbose=args.verbose)

	# Die if no hits found
	if not glob.glob(os.path.join(os.path.abspath(resultDir),'*.tab')):
		print("No hits found in %s . Quitting." % resultDir)
		sys.exit(1)

	# Import hits from nhmmer result files
	hitTable = None
	for resultfile in glob.glob(os.path.join(os.path.abspath(resultDir),'*.tab')):
		hitTable = import_nhmmer(infile=resultfile,hitTable=hitTable,prefix=args.prefix)

	# Group hits by model and chromosome (hitsDict), and initiate hit tracker hitIndex to manage pair-searching
	hitsDict,hitIndex = table2dict(hitTable)

	# If pairing is off, just report the hits
	if args.nopairing:
		writeTIRs(outDir=outDir, hitTable=hitTable, maxeval=args.maxeval, genome=genome)
		# Remove temp directory
		if not args.keeptemp:
			shutil.rmtree(tempDir)
		sys.exit(1)

	# Populate hitIndex with acceptible candidate partners (compatible strand and distance.)
	hitIndex = parseHits(hitsDict=hitsDict, hitIndex=hitIndex, maxDist=args.maxdist)

	# Run iterative pairing procedure
	hitIndex,paired,unpaired = iterateGetPairs(hitIndex, stableReps=args.stableReps)

	# Write TIR hits to fasta for each pHMM
	writeTIRs(outDir=outDir, hitTable=hitTable, maxeval=args.maxeval, genome=genome)

	# Extract paired hit regions (candidate TEs / MITEs)
	pairedEles = fetchElements(paired=paired, hitIndex=hitIndex, genome=genome)

	# Write paired TIR features to fasta
	writeElements(outDir, eleDict=pairedEles)

	# Write paired features to gff3, optionally also report paired/unpaired TIRs
	if args.gffOut:
		# Get unpaired TIRs
		if args.reportTIR in ['all','unpaired']:
			unpairedTIRs = fetchUnpaired(hitIndex=hitIndex)
		else:
			unpairedTIRs = None
	# Write gff3
		gffWrite(outpath=args.gffOut, featureList=pairedEles, writeTIRs=args.reportTIR, unpaired=unpairedTIRs)

	# Remove temp directory
	if not args.keeptemp:
		shutil.rmtree(tempDir)

def mainArgs():
	'''Parse command line arguments.'''
	parser = argparse.ArgumentParser(
							description	=	'Map TIR-pHMM models to genomic sequences for annotation of MITES and complete DNA-Transposons.',
							prog		=	'mapmite'
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

if __name__== '__main__':
	args = mainArgs()
	main(args)


"""
# Notes:
## nhmmer output format (Warning: Fields delimited by random number of spaces. Convert to tabs.)
(1) target name: The name of the target sequence or profile.
(2) accession: The accession of the target sequence or profile, or ’-’ if none.

(3) query name: The name of the query sequence or profile.
(4) accession: The accession of the query sequence or profile, or ’-’ if none.

(5) hmmfrom: The position in the hmm at which the hit starts.
(6) hmm to: The position in the hmm at which the hit ends.

(7) alifrom: The position in the target sequence at which the hit starts.
(8) ali to: The position in the target sequence at which the hit ends.

(9) envfrom: The position in the target sequence at which the surrounding envelope starts. 
(10) env to: The position in the target sequence at which the surrounding envelope ends. 

(11) sq len: The length of the target sequence..

(12) strand: The strand on which the hit was found (“-” when alifrom > ali to).

(13) E-value: The expectation value (statistical significance) of the target, as above.

(14) score (full sequence): The score (in bits) for this hit. It includes the biased-composition correction.

(15) Bias (full sequence): The biased-composition correction, as above

(16) description of target: The remainder of the line is the target’s description line, as free text.

# Future
## Align signigicant TIR hits using query model as guide.
hmmalign -o outfile.aln \
--mapali hmmsourcealign.stckhlm \
--trim \
--dna \
--informat FASTA \
--outformat Stockholm \ #[Stockholm, Pfam, A2M, PSIBLAST]
<hmmfile> <seqfile>

"""