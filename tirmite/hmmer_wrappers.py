import re
import os
import sys
import glob
import shutil
import tempfile
import subprocess

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