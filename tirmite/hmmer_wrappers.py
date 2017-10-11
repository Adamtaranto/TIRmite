import re 
import os
from shlex import quote

def cleanID(s):
	"""Remove non alphanumeric characters from string. Replace whitespace with underscores."""
	s = re.sub(r"[^\w\s]", '', s)
	s = re.sub(r"\s+", '_', s)
	return s

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
	command = quote(exePath) + " --dna -n " + modelname
	# Optional set cores
	if cores:
		command += ' --cpu ' + str(cores)
	if outdir:
		modelout = os.path.join(outdir,modelname + ".hmm")
	else:
		modelout = modelname + ".hmm"
	# Append output file name and source alignment, return command string
	command = ' '.join([command,quote(os.path.abspath(modelout)),quote(os.path.abspath(inAlign))])
	return command,os.path.abspath(modelout)

def _hmmpress_command(exePath="hmmpress", hmmfile=None):
	'''Construct the hmmbuild command'''
	# Base command
	command = quote(exePath) + " -f " + quote(os.path.abspath(hmmfile))
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
	command = quote(exePath) + " --tblout " + quote(outfile)
	if cores:
			command += ' --cpu ' + str(cores)
	if evalue:
			command += ' -E ' + str(evalue)
	if nobias:
		command += ' --nobias'
	if matrix:
		command += ' --mxfile ' + quote(os.path.abspath(matrix))
	command += " --noali --notextw --dna --max " + quote(os.path.abspath(modelPath)) + " " + quote(os.path.abspath(genome))
	return command,outdir