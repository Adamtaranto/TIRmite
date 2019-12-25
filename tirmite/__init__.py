#!/usr/bin/env python
#python 3
#tirmite.py
#Contact, Adam Taranto, adam.taranto@anu.edu.au

########################################################################
# Map TIR-pHMM models to genomic sequences for annotation of MITES and #
# complete DNA-Transposons.                                            #
########################################################################

import re
import os
import sys
import glob
import shutil
import tempfile
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from datetime import datetime
from collections import Counter
from collections import namedtuple
from operator import attrgetter
from .hmmer_wrappers import _hmmbuild_command,_hmmpress_command,_nhmmer_command
from .bowtie2_wrappers import _bowtie2build_cmd,_bowtie2_cmd,_bam2bed_cmd
from .runBlastn import makeBlast, run_blast
from pymummer import coords_file, alignment, nucmer

__version__ = "1.1.4"

class Error (Exception): pass

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

def tSplitchecks(args):
	"""Housekeeping tasks: Create output files/dirs and temp dirs as required."""
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
	tempDir = os.path.join(os.getcwd(),"temp_" + getTimestring())
	os.makedirs(tempDir)
	# Set prefix to infile basename if none
	if not args.prefix:
		prefix = os.path.splitext(os.path.basename(args.infile))[0]
	else:
		prefix = args.prefix
	# Create outfile paths
	outfile = prefix + "_tsplit_output.fasta"
	outpath = os.path.join(outDir,outfile)
	# Return full path to output file and temp directory
	return outpath,tempDir

def importFasta2List(file):
	"""Load elements from multifasta file. Check that seq IDs are unique."""
	# Read in elements from multifasta file, convert seqrecord iterator to list.
	records = list(SeqIO.parse(file, "fasta"))
	# Check names are unique
	checkUniqueID(records)
	# If unique, return record list.
	return records

class Error (Exception): pass

def decode(x):
	try:
		s = x.decode()
	except:
		return x
	return s

def cleanID(s):
	"""Remove non alphanumeric characters from string. Replace whitespace with underscores."""
	s = re.sub(r"[^\w\s]", '', s)
	s = re.sub(r"\s+", '_', s)
	return s

def _write_script(cmds,script):
	'''Write commands into a bash script'''
	f = open(script, 'w+')
	for cmd in cmds:
		print(cmd, file=f)
	f.close()

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

def run_cmd(cmds,verbose=False,tempDir=None,keeptemp=False):
	'''Write and excute HMMER script'''
	if not tempDir:
		tempDir = os.getcwd()
	tmpdir = tempfile.mkdtemp(prefix='tmp.', dir=tempDir)
	original_dir = os.getcwd()
	os.chdir(tmpdir)
	script = 'run_jobs.sh'
	_write_script(cmds,script)
	syscall('bash ' + script, verbose=verbose)
	os.chdir(original_dir)
	if not keeptemp:
		shutil.rmtree(tmpdir)

def isfile(path):
	if not os.path.isfile(path):
		print("Input file not found: %s" % path)
		sys.exit(1)

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

def getbtName(file):
	"""Load seqs from file. Check if more than one record. report first id lable."""
	# Read in elements from fasta file, convert seqrecord iterator to list.
	records = list(SeqIO.parse(file, "fasta"))
	# Check only one seq
	if len(records) > 1:
		print('Warning: %s contains multiple sequences! \n May result is incorrect TIR pairing.' % file)
	return cleanID(str(records[0].id))

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
		run_cmd(build_cmds,verbose=args.verbose,keeptemp=args.keeptemp)
	# Press and write nhmmer cmd for all models in hmmDB directory
	for hmm in glob.glob(os.path.join(hmmDB,'*.hmm')):
		hmmPressCmd = _hmmpress_command(exePath=args.hmmpress, hmmfile=hmm)
		nhmmerCmd,resultDir = _nhmmer_command(exePath=args.nhmmer,nobias=args.nobias,matrix=args.matrix,modelPath=hmm,genome=args.genome,evalue=args.maxeval,cores=args.cores,outdir=tempDir)
		cmds.append(hmmPressCmd)
		cmds.append(nhmmerCmd)
	# Return list of cmds and location of file result files
	return cmds,resultDir,hmmDB

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
	#if prefix:
	#	df['model'] = str(prefix) + '_' + df['model'].astype(str)
	return df


def import_BED(infile=None,hitTable=None,prefix=None):
	''' Read TIR bedfile to pandas dataframe.'''
	# Format: Chrm, start, end, name, evalue, strand
	hitRecords = list()
	with open(infile, "rU") as f:
		for line in f.readlines():
			li = line.strip()
			if not li.startswith("#"):
				li = li.split()
				hitRecords.append({
				'target':li[0],
				'model':li[3],
				'hmmStart':'NA',
				'hmmEnd':'NA',
				'hitStart':li[1],
				'hitEnd':li[2],
				'strand':li[5],
				'evalue':li[4],
				'score':'NA',
				'bias':'NA'})
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
	#if prefix:
	#	df['model'] = str(prefix) + '_' + df['model'].astype(str)
	return df

def import_mapped(infile=None,tirName="refTIR",hitTable=None,prefix=None):
	''' Read bowtie2 mapped TIR locations from bedfile to pandas dataframe.'''
	print('Bowtie2 mapped reads will not be filtered on e-value.')
	hitRecords = list()
	with open(infile, "rU") as f:
		for line in f.readlines():
			li = line.strip()
			if not li.startswith("#"):
				li = li.split()
				if li[3] == '+':
					hitRecords.append({
					'target':str(li[0]),
					'model':tirName,
					'hmmStart':'NA',
					'hmmEnd':'NA',
					'hitStart':int(li[1]),
					'hitEnd':int(li[2]),
					'strand':str(li[3]),
					'evalue':0,
					'score':'NA',
					'bias':'NA'})
				elif li[3] == '-':
					hitRecords.append({
					'target':str(li[0]),
					'model':tirName,
					'hmmStart':'NA',
					'hmmEnd':'NA',
					'hitStart':int(li[1]),
					'hitEnd':int(li[2]),
					'strand':str(li[3]),
					'evalue':0,
					'score':'NA',
					'bias':'NA'})
	# Convert list of dicts into dataframe
	df = pd.DataFrame(hitRecords)
	# Reorder columns
	cols = ['model','target','hitStart','hitEnd','strand','evalue','score','bias','hmmStart','hmmEnd']
	df = df.ix[:, cols]
	if hitTable is not None:
		# If an existing table was passed, concatenate
		df = pd.concat([df,hitTable], ignore_index=True)
	# Sort hits by model, Chromosome, location, and strand
	df = df.sort_values(['model','target','hitStart','hitEnd','strand'], ascending=[True,True,True,True,True])
	# Reindex
	df = df.reset_index(drop=True)
	#if prefix:
	#	df['model'] = str(prefix) + '_' + df['model'].astype(str)
	return df

def filterHitsLen(hmmDB=None, mincov=None, hitTable=None):
		modelLens = dict()
		for hmm in glob.glob(os.path.join(hmmDB,'*.hmm')):
			hmmLen = None
			hmmName = None
			with open(hmm, "rU") as f:
				for line in f.readlines():
					li = line.strip()
					if li.startswith("LENG"):
						hmmLen = int(li.split()[1])
					if li.startswith("NAME"):
						hmmName = str(li.split()[1])
				if hmmLen and hmmName:
					modelLens[hmmName] = hmmLen
		for model in modelLens.keys():
			minlen = modelLens[model] * mincov
			hitTable = hitTable.ix[~((hitTable['model'] == model) & ((hitTable['hitEnd'].astype(int) - hitTable['hitStart'].astype(int)) + 1 < minlen))]
		return hitTable

def filterHitsEval(maxeval=None, hitTable=None):
	''' Filter hitTable df to remove hits with e-value in excess of --maxeval.
	'''
	hitTable = hitTable.ix[((hitTable['evalue'].astype(float)) < float(maxeval))]
	return hitTable

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
	hitTup = namedtuple('Elem', ['model','target','hitStart','hitEnd','strand','idx','evalue'])
	# Add each record to dicts
	for row in hitTable.iterrows():
		record = hitTup(row[1].model,row[1].target,int(row[1].hitStart),int(row[1].hitEnd),row[1].strand,row[0],row[1].evalue)
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

def extractTIRs(model=None, hitTable=None, maxeval=0.001, genome=None, padlen=None):
	''' For significant hits in model, compose seqrecords.'''
	# Note: Padding not yet enabled for TIR extraction.
	hitcount = 0
	seqList = list()
	for index,row in hitTable[hitTable['model'] == model].iterrows():
		if float(row['evalue']) <= maxeval:
			hitcount += 1
			if padlen:
				hitrecord = genome[row['target']][int(row['hitStart'])-1-padlen:int(row['hitStart'])-1].lower() + \
				genome[row['target']][int(row['hitStart'])-1:int(row['hitEnd'])] + \
				genome[row['target']][int(row['hitEnd']):int(row['hitEnd'])+padlen].lower()
			else:
				hitrecord = genome[row['target']][int(row['hitStart'])-1:int(row['hitEnd'])]
			hitrecord.id = model + '_' + str(index)
			if row['strand'] == '-':
				hitrecord = hitrecord.reverse_complement(id=hitrecord.id + "_rc")
			hitrecord.name = hitrecord.id 
			hitrecord.description = '_'.join(['[' + str(row['target']) + ':' + str(row['strand']),str(row['hitStart']),str(row['hitEnd']) + ' modelAlignment:' + row['hmmStart'],row['hmmEnd'] + ' E-value:' + str(row['evalue']) + ']'])
			# Append record to list
			seqList.append(hitrecord)
		else:
			continue
	# Return seqrecord list and total hit count for model
	return seqList,hitcount 

def writeTIRs(outDir=None, hitTable=None, maxeval=0.001, genome=None, prefix=None, padlen=None):
	''' Write all hits per Model to a multifasta in the outdir'''
	# Note: Padding not yet enabled for TIR extraction.
	if prefix:
		prefix = cleanID(prefix) + '_'
	else:
		prefix = ''
	if outDir:
		outDir = os.path.abspath(outDir)
		if not os.path.isdir(outDir):
				os.makedirs(outDir)
	else:
		outDir = os.getcwd()
	for model in hitTable['model'].unique():
		# List of TIR seqrecords, and count of hits
		seqList,hitcount = extractTIRs(model=model, hitTable=hitTable, maxeval=maxeval, genome=genome, padlen=padlen)
		outfile = os.path.join(outDir, prefix + model + "_hits_" + str(hitcount) + ".fasta")
		# Write extracted hits to model outfile
		with open(outfile, "w") as handle:
			for seq in seqList:
				seq.id = prefix + str(seq.id)
				SeqIO.write(seq, handle, "fasta")

#CS10_Chromosome_02_+_88294_88353_modelAlignment:1_60

def flipTIRs(x,y):
	''' Sort hits into left and right TIRs.'''
	left2right = sorted([x,y], key=attrgetter('hitStart', 'hitEnd'))
	return (left2right[0],left2right[1])

def fetchElements(paired=None, hitIndex=None, genome=None):
	''' Extract complete sequence of paired elements, 
		asign names and child TIRs for use in seq and GFF reporting.'''
	TIRelements = dict()
	gffTup = namedtuple('gffElem', ['model', 'chromosome', 'start', 'end', 'strand', 'type', 'id', 'leftHit' , 'rightHit','seq','evalue'])
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
			eleSeq.description = '_'.join([ '[' + leftHit.target + ':' + str(leftHit.hitStart),str(rightHit.hitEnd)]) + " len=" + str(rightHit.hitEnd-leftHit.hitStart) + ']'
			TIRelement = gffTup(model, leftHit.target , leftHit.hitStart , rightHit.hitEnd , leftHit.strand , "TIR_Element", eleID, leftHit,rightHit, eleSeq,'NA')
			TIRelements[model].append(TIRelement)
	# Return list of element info tuples
	return TIRelements

def writeElements(outDir, eleDict=None, prefix=None):
	''' Takes dict of extracted sequences keyed by model. Writes to fasta by model.'''
	if prefix:
		prefix = cleanID(prefix) + '_'
	else:
		prefix = ''
	for model in eleDict.keys():
		outfile = os.path.join(outDir,prefix + model + '_elements.fasta')
		with open(outfile, "w") as handle:
			for element in eleDict[model]:
				element.seq.id = prefix + str(element.seq.id)
				SeqIO.write(element.seq, handle, "fasta")

def writePairedTIRs(outDir=None,paired=None, hitIndex=None, genome=None, prefix=None, padlen=None):
	''' Extract TIR sequence of paired hits, write to fasta.'''
	# Note: Sequence padding not yet enabled for paired TIRs
	TIRpairs = dict()
	gffTup = namedtuple('gffElem', ['model', 'chromosome', 'start', 'end', 'strand', 'type', 'id', 'leftHit' , 'rightHit','seq','evalue'])
	for model in paired.keys():
		TIRpairs[model] = list()
		model_counter = 0
		for x,y in paired[model]:
			model_counter += 1
			x = hitIndex[model][x]['rec']
			y = hitIndex[model][y]['rec']
			leftHit,rightHit = flipTIRs(x,y)
			eleID = model + "_TIRpair_" + str(model_counter)
			# If padlen set, extract hit with x bases either side
			if padlen:
				eleSeqLeft = genome[leftHit.target][int(leftHit.hitStart)-1-padlen:int(leftHit.hitStart)-1].lower() + \
				genome[leftHit.target][int(leftHit.hitStart)-1:int(leftHit.hitEnd)] + \
				genome[leftHit.target][int(leftHit.hitEnd):int(leftHit.hitEnd)+padlen].lower()
				eleSeqRight = genome[leftHit.target][int(rightHit.hitStart)-1-padlen:int(rightHit.hitStart)-1].lower() + \
				genome[leftHit.target][int(rightHit.hitStart)-1:int(rightHit.hitEnd)] + \
				genome[leftHit.target][int(rightHit.hitEnd):int(rightHit.hitEnd)+padlen].lower()
			else:
				eleSeqLeft = genome[leftHit.target][int(leftHit.hitStart)-1:int(leftHit.hitEnd)]
				eleSeqRight = genome[leftHit.target][int(rightHit.hitStart)-1:int(rightHit.hitEnd)]	
			eleSeqRight = eleSeqRight.reverse_complement()
			eleSeqLeft.id = eleID + "_L"
			eleSeqLeft.name = eleID + "_L"
			eleSeqLeft.description = '_'.join([ '[' + leftHit.target + ':' + str(leftHit.hitStart),str(leftHit.hitEnd)]) + ']'
			eleSeqRight.id = eleID + "_R"
			eleSeqRight.name = eleID + "_R"
			eleSeqRight.description = '_'.join([ '[' + leftHit.target + ':' + str(rightHit.hitEnd),str(rightHit.hitStart)]) + ']'
			TIRleft = gffTup(model, leftHit.target , leftHit.hitStart , leftHit.hitEnd , leftHit.strand , "TIR", eleSeqLeft.id, leftHit,rightHit, eleSeqLeft,'NA')
			TIRright = gffTup(model, leftHit.target , rightHit.hitStart , rightHit.hitEnd , leftHit.strand , "TIR", eleSeqRight.id, leftHit,rightHit, eleSeqRight,'NA')
			TIRpairs[model].append(TIRleft)
			TIRpairs[model].append(TIRright)
	if prefix:
		prefix = cleanID(prefix) + '_'
	else:
		prefix = ''
	for model in TIRpairs.keys():
		outfile = os.path.join(outDir,prefix + model + '_paired_TIR_hits_' + str(model_counter * 2) + '.fasta')
		with open(outfile, "w") as handle:
			for element in TIRpairs[model]:
				element.seq.id = prefix + str(element.seq.id)
				SeqIO.write(element.seq, handle, "fasta")

def fetchUnpaired(hitIndex=None):
	'''	Take list of unpaired hit IDs from listunpaired(),
		Compose TIR gff3 record. '''
	orphans = list()
	gffTup = namedtuple('gffElem', ['model', 'chromosome', 'start', 'end', 'strand', 'type', 'id', 'leftHit' , 'rightHit', 'seq','evalue'])
	for model in hitIndex.keys():
		for recID in hitIndex[model].keys():
			if not hitIndex[model][recID]['partner']:
				x = hitIndex[model][recID]['rec']
				orphan = gffTup(x.model, x.target , x.hitStart , x.hitEnd , x.strand , "orphan_TIR", x.idx, None , None, None, x.evalue)
				orphans.append(orphan)
	return orphans

def	gffWrite(outpath=None, featureList=list(), writeTIRs=True, unpaired=None, suppressMeta=False, prefix=None):
	'''	Write predicted paired-TIR features (i.e. MITEs) from fetchElements() as GFF3.
		Optionally, write child TIRS and orphan TIRs to GFF3 also.'''
	if prefix:
		prefix = cleanID(prefix) + '_'
	else:
		prefix = ''
	# If path to output gff3 file not provided, set default location.
	if not outpath:
		outpath = os.path.join(os.getcwd(),"tirmite_features.gff3")
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
				file.write('\t'.join(	[str(x.chromosome),"tirmite",x.type,str(x.start),str(x.end),".",x.strand,".",
										'ID=' + prefix + str(x.model) + "_" + str(x.id) + ';model=' + str(x.model) + ';evalue=' + str(x.evalue) + ';']
										) + '\n'
										)
			if x.type == "TIR_Element":
				# Write Element line
				file.write('\t'.join(	[str(x.chromosome),"tirmite",x.type,str(x.start),str(x.end),".",x.strand,".",
										'ID=' + prefix + str(x.id) + ';model=' + str(x.model) + ';']
										) + '\n'
										)
				if writeTIRs in ['all','paired']:
					# Write left TIR line as child
					l = x.leftHit
					file.write('\t'.join(	[str(l.target),"tirmite","paired_TIR",str(l.hitStart),str(l.hitEnd),".",l.strand,".",
											'ID=' + prefix + str(x.model) + "_" + str(l.idx) + ';model=' + str(x.model) + ';Parent=' + str(x.id) + ';evalue=' + str(l.evalue) + ';']
											) + '\n'
											)
					# Write right TIR line as child on neg strand
					r = x.rightHit
					file.write('\t'.join(	[str(r.target),"tirmite","paired_TIR",str(r.hitStart),str(r.hitEnd),".",r.strand,".",
											'ID=' + prefix + str(x.model) + "_" + str(r.idx) + ';model=' + str(x.model) + ';Parent=' + str(x.id) + ';evalue=' + str(r.evalue) + ';']
											) + '\n'
											)

def getLTRs(elements=None, flankdist=10, minid=80, minterm=10, minseed=5, diagfactor=0.3, report='split', temp=None,keeptemp=False,alignTool='nucmer',verbose=False):
	"""Align elements to self and attempt to identify LTRs."""
	# Set temp directory to cwd if none.
	if not temp:
		temp = os.getcwd()
	# For each candidate LTR element
	for rec in elements:
		# Create temp paths for single element fasta and alignment coords
		tempFasta = os.path.join(temp, cleanID(rec.id) + '.fasta')
		tempCoords = tempFasta + '.coords'
		# Write current element to single fasta
		manageTemp(record=rec, tempPath=tempFasta, scrub=False)
		# Align to self with nucmer
		if alignTool == 'nucmer':
			# Compose Nucmer script for current element vs self
			runner = nucmer.Runner(	tempFasta, tempFasta, tempCoords,
									min_id		=	minid, 
									min_length	=	minseed,
									diagfactor	=	diagfactor,
									mincluster	=	minterm,
									breaklen	=	200,
									maxmatch	=	True,
									simplify	=	False
									)
			# Execute nucmer
			runner.run()
		elif alignTool == 'blastn':
			# Alternatively, use blastn as search tool and write nucmer.coords-like output.
			cmds = makeBlast(seq=tempFasta, outfile=tempCoords, pid=minid)
			run_blast(cmds, verbose=verbose)
		# Import coords file to iterator object
		file_reader = coords_file.reader(tempCoords)
		# Exclude hits to self. Also converts iterator output to stable list
		alignments = [hit for hit in file_reader if not hit.is_self_hit()]
		# Filter hits less than min length (Done internally for nucmer, not blastn.)
		alignments = [hit for hit in alignments if hit.ref_end - hit.ref_start >= minterm]
		# Filter for hits on same strand i.e. tandem repeats / LTRs
		alignments = [hit for hit in alignments if hit.on_same_strand()]
		# Filter for 5' repeats which begin within x bases of element start
		alignments = [hit for hit in alignments if hit.ref_start <= flankdist]
		# Filter for 5' repeats whose 3' match ends within x bases of element end
		alignments = [hit for hit in alignments if len(rec) - hit.qry_end <= flankdist]
		# Scrub overlappying ref / query segments, and also complementary 3' to 5' flank hits
		alignments = [hit for hit in alignments if hit.ref_end < hit.qry_start]
		# Sort largest to smallest dist between end of ref (subject) and start of query (hit)
		# x.qry_start (3') - x.ref_end (5') = Length of internal segment
		#Note: Need to check that this sorts correctlly for alignments is same orientation.
		print("Warning: Check candidate pairs are correctly sorted.")
		alignments = sorted(alignments, key=lambda x: (x.qry_start - x.ref_end), reverse=True)
		# If alignments exist after filtering report features using alignment pair with largest 
		# internal segment i.e. first element in sorted list.
		if alignments:
			if verbose:
				[print(x) for x in alignments]
			if report == 'all':
				# yield original element
				yield rec
			if report in ['split','external']:
				# yield LTR slice - append "_LTR"
				extSeg = rec[alignments[0].ref_start:alignments[0].ref_end + 1]
				extSeg.id = extSeg.id + "_LTR"
				extSeg.name = extSeg.id
				extSeg.description = "[" + rec.id + " LTR segment]"
				yield extSeg
			if report in ['split','internal']:
				# yield internal slice - append "_I"
				intSeg = rec[alignments[0].ref_end:alignments[0].qry_start + 1] 
				intSeg.id = intSeg.id + "_I"
				intSeg.name = intSeg.id
				intSeg.description = "[" + rec.id + " internal segment]"
				yield intSeg
		else:
			# If alignment list empty after filtering print alert and continue
			print('No LTRs found for candidate element: %s' % rec.id)
		# Scrub single fasta and coords file for current element.
		if not keeptemp:
			manageTemp(tempPath=tempFasta, scrub=True)
			manageTemp(tempPath=tempCoords, scrub=True)

def getTIRs(elements=None, flankdist=2, minid=80, minterm=10, minseed=5, diagfactor=0.3, mites=False, report='split', temp=None,keeptemp=False,alignTool='nucmer',verbose=False):
	""" 
	Align elements to self and attempt to identify TIRs. 
	Optionally attempt to construct synthetic MITEs from TIRs.
	"""
	# Set temp directory to cwd if none.
	if not temp:
		temp = os.getcwd()
	# For each candidate LTR element
	for rec in elements:
		# Create temp paths for single element fasta and alignment coords
		tempFasta = os.path.join(temp, cleanID(rec.id) + '.fasta')
		tempCoords = tempFasta + '.coords'
		# Write current element to single fasta
		manageTemp(record=rec, tempPath=tempFasta, scrub=False)
		# Align to self with nucmer
		if alignTool == 'nucmer':
			# Compose Nucmer script for current element vs self
			runner = nucmer.Runner(	tempFasta, tempFasta, tempCoords,
									min_id		=	minid,
									min_length	=	minseed,
									diagfactor	=	diagfactor,
									mincluster	=	minterm,
									breaklen	=	200,
									maxmatch	=	True,
									simplify	=	False
									)
			# Execute nucmer
			runner.run()
		elif alignTool == 'blastn':
			# Alternatively, use blastn as search tool and write nucmer.coords-like output.
			cmds = makeBlast(seq=tempFasta, outfile=tempCoords, pid=minid)
			run_blast(cmds,verbose=verbose)
		# Import coords file to iterator object
		file_reader = coords_file.reader(tempCoords)
		# Exclude hits to self. Also converts iterator output to stable list
		alignments = [hit for hit in file_reader if not hit.is_self_hit()]
		# Filter hits less than min length (Done internally for nucmer, not blastn.)
		alignments = [hit for hit in alignments if hit.ref_end - hit.ref_start >= minterm]
		# Filter for hits on same strand i.e. tandem repeats / LTRs
		alignments = [hit for hit in alignments if not hit.on_same_strand()]
		# Filter for 5' repeats which begin within x bases of element start
		alignments = [hit for hit in alignments if hit.ref_start <= flankdist]
		# Scrub overlappying ref / query segments, and also complementary 3' to 5' flank hits
		alignments = [hit for hit in alignments if hit.ref_end < hit.qry_end]
		# Sort largest to smallest dist between end of ref (subject) and start of query (hit)
		# x.qry_end - x.ref_end = 5'end of right TIR - 3' end of left TIR = length of internal segment
		# TIR pair with smallest internal segment (longest TIRs) is first in list.
		alignments = sorted(alignments, key=lambda x: (x.qry_end - x.ref_end), reverse=False)
		# If alignments exist after filtering report features using alignment pair with largest 
		# internal segment i.e. first element in sorted list.
		if alignments:
			if verbose:
				[print(x) for x in alignments]
			if report in ['split','external','all']:
				# yield TIR slice - append "_TIR"
				extSeg = rec[alignments[0].ref_start:alignments[0].ref_end + 1]
				extSeg.id = extSeg.id + "_TIR"
				extSeg.name = extSeg.id
				extSeg.description = "[" + rec.id + " TIR segment]"
				yield extSeg
			if report in ['split','internal','all']:
				# yield internal slice - append "_I"
				intSeg = rec[alignments[0].ref_end:alignments[0].qry_end + 1]
				intSeg.id = intSeg.id + "_I"
				intSeg.name = intSeg.id
				intSeg.description = "[" + rec.id + " internal segment]"
				yield intSeg
			if report == 'all':
				yield rec
			if mites:
				# Assemble TIRs into hypothetical MITEs
				synMITE = rec[alignments[0].ref_start:alignments[0].ref_end + 1] + rec[alignments[0].qry_end:alignments[0].qry_start + 1]
				synMITE.id = synMITE.id + "_synMITE"
				synMITE.name = synMITE.id
				synMITE.description = "[Synthetic MITE constructed from " + rec.id + " TIRs]"
				yield synMITE
		else:
			# If alignment list empty after filtering print alert and continue
			print('No TIRs found for candidate element: %s' % rec.id)
		# Scrub single fasta and coords file for current element.
		if not keeptemp:
			manageTemp(tempPath=tempFasta, scrub=True)
			manageTemp(tempPath=tempCoords, scrub=True)

def segWrite(outfile,segs=None):
	"""Take a generator object yielding seqrecords and write each to outfile in fasta format."""
	seqcount = 0
	if segs:
		with open(outfile, "w") as handle:
			for seq in segs:
				seqcount += 1
				SeqIO.write(seq, handle, "fasta")
		if seqcount == 0:
			os.remove(outfile)

#gffTup fields: 'model', 'chromosome', 'start', 'end', 'strand', 'type', 'id', 'score','bias', 'evalue', 'leftHit' , 'rightHit', 'eleSeq'
#Types: "TIR_Element", "orphan_TIR"

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

"""
# Useful attributes of pymummer objects:
[x.ref_start for x in alignments]
[x.ref_end for x in alignments]
[x.qry_start for x in alignments]
[x.qry_end for x in alignments]
[x.hit_length_ref for x in alignments]
[x.hit_length_qry for x in alignments]
[x.percent_identity for x in alignments]
[x.ref_length for x in alignments]
[x.qry_length for x in alignments]
[x.frame for x in alignments]
[x.ref_name for x in alignments]
[x.qry_name for x in alignments]

## Don't use these, bizzaro format. Not indexed to 0. Cannot sort as ints.
#coord.reverse_query()
#coord.reverse_reference()
#coord.qry_coords()
#coord.ref_coords()
"""
