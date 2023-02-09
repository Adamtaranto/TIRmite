import os
import sys
import shutil
import tempfile
import subprocess
from shlex import quote

class Error (Exception): pass

def decode(x):
	try:
		s = x.decode()
	except:
		return x
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

def makeBlast(seq=None, outfile=None, pid=60):
	cmd = 'blastn -word_size 4 -outfmt "6 qstart qend sstart send length positive pident qlen slen qframe sframe qseqid sseqid" -query ' + quote(str(seq)) + ' -subject ' + quote(str(seq)) + ' -out ' + quote(str(outfile)) + ' -perc_identity ' + str(pid)
	return [cmd]

def run_blast(cmds,verbose=False):
	'''Write and excute HMMER script'''
	tmpdir = tempfile.mkdtemp(prefix='tmp.', dir=os.getcwd())
	original_dir = os.getcwd()
	os.chdir(tmpdir)
	script = 'run_jobs.sh'
	_write_script(cmds,script)
	syscall('bash ' + script, verbose=verbose)
	os.chdir(original_dir)
	shutil.rmtree(tmpdir)


'''
Recreate pymummer-like coords file output from blast.

#Blast field, coords field, Description

qstart			[S1] 	Start of the alignment region in the reference sequence 
qend			[E1] 	End of the alignment region in the reference sequence 
sstart			[S2] 	Start of the alignment region in the query sequence 
send 			[E2] 	End of the alignment region in the query sequence 
length			[LEN 1] Length of the alignment region in the reference sequence
positive		[LEN 2] Length of the alignment region in the query sequence (#"positive" is just a filler as blast won't all repeated fields)
pident			[% IDY] Percent identity of the alignment 
qlen			[LEN R] Length of the reference sequence 
slen			[LEN Q] Length of the query sequence 
qframe sframe	[FRM] 	Reading frame for the reference AND query sequence alignments respectively 
qseqid sseqid	[TAGS] 	The reference AND query FastA IDs respectively. All output coordinates and lengths are relative to the forward strand of the reference DNA sequence.
'''