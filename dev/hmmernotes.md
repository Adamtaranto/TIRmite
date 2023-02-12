
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
hmmalign -o outfile.aln
--mapali hmmsourcealign.stckhlm
--trim
--dna
--informat FASTA
--outformat Stockholm #[Stockholm, Pfam, A2M, PSIBLAST]
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
