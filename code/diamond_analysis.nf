
params.resdir="/pasteur/homes/mamorel/GEVA/HCV/Data/assembly/"


resdir=file(params.resdir)
resdir.with {mkdirs()}


//inputChannel = Channel.fromPath( '/pasteur/homes/mamorel/GEVA/HCV/Data/assembly/*scaffolds.fasta' )
inputChannel = Channel.fromPath( '/pasteur/homes/mamorel/GEVA/HCV/Data/raw/hybridemapping/HCV/scaffolds_HCV_11A*.fasta' )

process interesting_contigs{

	input: 
	file(scaffolds) from inputChannel
	output:
	file("big*.fasta") into selected_contigs

	script:
	"""
#!~/miniconda3/bin/python
from Bio import SeqIO
NAME = "${scaffolds}".split('assembly')[0]
to_keep=[]
for record in SeqIO.parse("${scaffolds}", "fasta"):
    if len(record) > 250:
        to_keep.append(record)
	
SeqIO.write(to_keep, "".join(["big", NAME, 'scaffolds.fasta']), "fasta")
	"""

}

process diamond{
    publishDir "${resdir}", pattern: '*matches.m8*', mode: 'copy'
	module 'diamond/0.9.24'
	input:
	file(contigs) from selected_contigs
	output:
	file("*matches.m8") into blast_results
	shell:
	'''
	DB='/pasteur/services/policy01/banques/prod/rel/nrprot/current/diamond/0.9/nr.dmnd'
	diamond blastx -d $DB -q !{contigs} -p 6 -k 10 -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle -o !{contigs}matches.m8 
	'''
}
