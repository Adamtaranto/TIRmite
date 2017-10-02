import os
from shlex import quote

def _bowtie2build_cmd(bt2Path="bowtie2-build",IdxPath="db/GenIdx",genome=None):
	'''Construct the bowtie2-build command'''
	# Base command
	cmd = ' '.join(quote(bt2Path),quote(os.path.abspath(genome)),IdxPath)
	return cmd

def _bowtie2_cmd(bt2Path="bowtie2",tirFasta=None,IdxPath="db/GenIdx",cores=None):
	'''Construct commands for bowtie2 mapping.'''
	# bowtie2 -x genidx -f -a --very-sensitive-local -U TIR.fa --al alignments.bam 
	# Base command
	cmd = ''.join(quote(bt2Path),'--al alignments.bam -f -a --very-sensitive-local -x',IdxPath,'-U',quote(os.path.abspath(tirFasta)))
	# Optional set cores
	if cores:
		cmd += ' --threads ' + str(cores)
	return cmd

def _bam2bed_cmd(samPath="samtools",bedPath="bedtools",tempDir=None):
	''' Filtering mapped reads with bedtools and samtools.
	# Fwd hits
	samtools view -b -F 0x10 alignments.bam | bedtools bamtobed -i stdin | awk -v OFS='\t' '{print $1,$2,$3,"+"}' > mapped.bed 
	# Rev hits
	samtools view -b -f 0x10 alignments.bam | bedtools bamtobed -i stdin | awk -v OFS='\t' '{print $1,$2,$3,"-"}' >> mapped.bed 
	'''
	# Base command
	mappedPath = os.path.join(tempDir,'bowtie2mappedTIR.bed')
	cmds = list() 
	cmds.append(' '.join(quote(samPath),"view -b -F 0x10 alignments.bam |",quote(bedPath),"bamtobed -i stdin | awk -v OFS='\\t' '{print $1,$2,$3,\"+\"}' >",quote(mappedPath)))
	cmds.append(' '.join(quote(samPath),"view -b -f 0x10 alignments.bam |",quote(bedPath),"bamtobed -i stdin | awk -v OFS='\\t' '{print $1,$2,$3,\"+\"}' >>",quote(mappedPath)))
	return cmds,mappedPath