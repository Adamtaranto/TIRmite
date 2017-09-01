#!/usr/bin/env python
#python 3
#tirmite.py
#Version 1.0 Adam Taranto, August 2017
#Contact, Adam Taranto, adam.taranto@anu.edu.au

########################################################################
# Map TIR-pHMM models to genomic sequences for annotation of MITES and #
# complete DNA-Transposons.                                            #
########################################################################

import os
import sys
import glob
import shutil
from .hmmer_wrappers import *
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from datetime import datetime
from collections import Counter
from collections import namedtuple
from operator import attrgetter

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
										'ID=' + str(x.model) + "_" + str(x.id) + ';model=' + str(x.model) + ';']
										) + '\n'
										)
			if x.type == "TIR_Element":
				# Write Element line
				file.write('\t'.join(	[str(x.chromosome),"tirmite",x.type,str(x.start),str(x.end),".",x.strand,".",
										'ID=' + str(x.id) + ';model=' + str(x.model) + ';']
										) + '\n'
										)
				if writeTIRs in ['all','paired']:
					# Write left TIR line as child
					l = x.leftHit
					file.write('\t'.join(	[str(l.target),"tirmite","paired_TIR",str(l.hitStart),str(l.hitEnd),".",l.strand,".",
											'ID=' + str(x.model) + "_" + str(l.idx) + ';model=' + str(x.model) + ';Parent=' + str(x.id) + ';']
											) + '\n'
											)
					# Write right TIR line as child on neg strand
					r = x.rightHit
					file.write('\t'.join(	[str(r.target),"tirmite","paired_TIR",str(r.hitStart),str(r.hitEnd),".",r.strand,".",
											'ID=' + str(x.model) + "_" + str(r.idx) + ';model=' + str(x.model) + ';Parent=' + str(x.id) + ';']
											) + '\n'
											)


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