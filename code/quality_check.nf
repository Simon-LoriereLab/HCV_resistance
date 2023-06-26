

// where the outputs will be stored
params.resdir="/pasteur/homes/mamorel/GEVA/Marie/Data/assembly/denovo/"
params.adapter = "/pasteur/homes/mamorel/GEVA/HCV/Data/raw/illumina_adapter.fasta"
params.pairs = "/pasteur/homes/mamorel/GEVA/Marie/Data/raw/NGS/HCV_11B*_R{1,2}.fastq.gz"

resdir=file(params.resdir)
resdir.with {mkdirs()}
adapter = file(params.adapter)

//retrieve the samples per pair (paired end)

//"/pasteur/homes/mamorel/GEVA/Marie/Data/raw/NGS/HCV_130*_R{1,2}.fastq.gz"
//"/pasteur/homes/mamorel/GEVA/HCV/Data/raw/hybridemapping/thirdpart/HCV*_R{1,2}.fastq.gz"
Channel
    .fromFilePairs(params.pairs, checkIfExists:true)
    .set { samples_ch}


//trim the fastq files and take as output the trimmed paired fq files 

process quality_control {
    publishDir "${resdir}", pattern: '*fastqc.zip', mode: 'copy'
    module 'fastqc/0.10.1'
    input:
    set sampleId, file(reads) from samples_ch
    output:
    set sampleId, file(reads) into quality_ch
    file "*fastqc.zip" into results
    shell:
    '''
    fastqc !{sampleId}_R1.fastq.gz !{sampleId}_R2.fastq.gz
    '''

}


process trimming {
	module 'Trimmomatic/0.36'
	input:
	set sampleId, file(reads) from quality_ch
	file (adapt) from adapter
	output:
	set sampleId, file("*P.fq.gz") into Cleanreads
	shell:
	'''
	 Trimmomatic PE -threads 4 !{sampleId}_R1.fastq.gz\
	 !{sampleId}_R2.fastq.gz !{sampleId}_1P.fq.gz !{sampleId}_1U.fq.gz !{sampleId}_2P.fq.gz \
	 !{sampleId}_2U.fq.gz ILLUMINACLIP:!{adapt}:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
	'''
	
}



//meta spades
process denovo_assembly {
	publishDir "${resdir}", pattern: '*assembly/scaffolds.fasta', mode: 'copy'
	module 'SPAdes/3.12.0:bwa/0.7.7'
	input:
	set sampleId, file(reads) from Cleanreads
	output:
	set sampleId, file("*_assembly/scaffolds.fasta") into metaspades_output

	shell:
	'''
	spades.py --meta -t 4 --pe1-1 !{sampleId}_1P.fq.gz --pe1-2 !{sampleId}_2P.fq.gz -o !{sampleId}_assembly
	'''
}

//python script to retrieve only the contigs that are more than 250 bp
process interesting_contigs{
    conda '/pasteur/homes/mamorel/miniconda3/envs/jupyter-notebook'
	input: 
	set sampleId, file(scaffolds) from metaspades_output
	output:
	set sampleId, file("big_scaffolds.fasta") into selected_contigs

	script:
	"""
#!/usr/bin/env python
from Bio import SeqIO
to_keep=[]
for record in SeqIO.parse("${scaffolds}", "fasta"):
    if len(record) > 250:
        to_keep.append(record)
	
SeqIO.write(to_keep, "big_scaffolds.fasta", "fasta")
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

'''
#for i in `find . -wholename "*assembly/scaffolds.fasta"` ; do a=`echo ${i#*/*/*/} | cut -d "/" -f 1` ; cp ${i} ~/GEVA/HCV/Data/assembly/denovo_analysis/thirdpart/${a}.scaffolds.fasta ; done
'''