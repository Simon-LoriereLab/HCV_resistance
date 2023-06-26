
// where the outputs will be stored
params.resdir="/pasteur/homes/mamorel/GEVA/HCV/Data/assembly/"


resdir=file(params.resdir)
resdir.with {mkdirs()}

//retrieve the samples per pair (paired end)
Channel
    .fromFilePairs("/pasteur/homes/mamorel/GEVA/HCV/Data/raw/hybridemapping/secondpart/HCV_6B*_R{1,2}.fastq.gz", checkIfExists:true)
    .set { samples_ch }

//meta spades
process denovo_assembly {
	publishDir "${resdir}", pattern: '*_assembly*', mode: 'copy'
	module 'SPAdes/3.12.0:bwa/0.7.7'
	input:
	set sampleId, file(reads) from samples_ch
	output:
	set sampleId, file("*_assembly/scaffolds.fasta") into metaspades_output

	shell:
	'''
	spades.py --meta -t 4 --pe1-1 !{sampleId}_R1.fastq.gz --pe1-2 !{sampleId}_R2.fastq.gz -o !{sampleId}_assembly
	'''
}

//python script to retrieve only the contigs that are more than 250 bp
process interesting_contigs{

	input: 
	set sampleId, file(scaffolds) from metaspades_output
	output:
	set sampleId, file("big_scaffolds.fasta") into selected_contigs

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

//diamond blast to analyse interesting contigs 
process diamond{
	module 'diamond/0.9.24'
	input:
	set sampleId, file(contigs) from selected_contigs
	output:
	set sampleId, file("*matches.m8") into blast_results
	shell:
	'''
	DB='/pasteur/services/policy01/banques/prod/rel/nrprot/current/diamond/0.9/nr.dmnd'
	diamond blastx -d $DB -q !{contigs} -p 6 -k 10 -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle -o !{sampleId}matches.m8 
	'''
}

blast_results.subscribe{ID, table ->  table.copyTo(file("${resdir}"));}

