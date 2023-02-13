from shlex import quote


def makeBlast(seq=None, outfile=None, pid=60):
    cmd = (
        'blastn -word_size 4 -outfmt "6 qstart qend sstart send length positive pident qlen slen qframe sframe qseqid sseqid" -query '
        + quote(str(seq))
        + " -subject "
        + quote(str(seq))
        + " -out "
        + quote(str(outfile))
        + " -perc_identity "
        + str(pid)
    )
    return [cmd]


"""
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
"""
