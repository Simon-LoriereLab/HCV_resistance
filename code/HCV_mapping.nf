params.refsequence = "/pasteur/homes/mamorel/GEVA/Marie/Data/assembly/mapping/second_mapping"
params.pairs = "/pasteur/homes/mamorel/GEVA/Marie/Data/raw/HCV_130*_R{1,2}.fastq.gz"
//"/pasteur/homes/mamorel/GEVA/HCV/Data/raw/hybridemapping/thirdpart/HCV128*_R{1,2}.fastq.gz"


//"/pasteur/homes/mamorel/GEVA/Marie/Data/raw/NGS/HCV_130*_R{1,2}.fastq.gz"

params.resdir="/pasteur/homes/mamorel/GEVA/Marie/Data/assembly/mapping/second_mapping/"
//"/pasteur/homes/mamorel/GEVA/HCV/Data/assembly/mapping/first_consensus_building/hybride/"
// "/pasteur/homes/mamorel/GEVA/Marie/Data/mapping/hybride/"

params.adapter = "/pasteur/homes/mamorel/GEVA/HCV/Data/raw/illumina_adapter.fasta"

refsequence = file(params.refsequence)
adapter = file(params.adapter)
//pairs = file(params.pairs)
resdir=file(params.resdir)
resdir.with {mkdirs()}

//import files two by two according to their pair ID
//In samples_ch we get one single identifier for the pair and the two files of the pair
Channel
    .fromFilePairs(params.pairs, checkIfExists:true)
    .set { samples_ch}


//trim the fastq files and take as output the trimmed paired fq files 
process trimming {
	module 'Trimmomatic/0.36'
	input:
	set sampleId, file(reads) from samples_ch
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


//mapping not too stringeant on the hybride 
process clc_assembly {
	publishDir "${resdir}", pattern: '*.cas', mode: 'copy'
	module 'clc-assembly-cell/5.1.0'
	input:
	set sampleId, file(reads) from Cleanreads
	file(ref) from refsequence
	output:
	set sampleId, file(ref), file(reads), file("*.cas") into First_mapping

	shell:
	'''
	clc_mapper -o !{sampleId}.cas -a local \
	-t 1 -r random -l 0.7 -s 0.7 -x 2 -g 3 -e 3 -d !{ref} \
	 -q -p fb ss 70 1200 -i !{sampleId}_1P.fq.gz !{sampleId}_2P.fq.gz

	'''
}

//first consensus extraction
process extract_consensus{
	publishDir "${resdir}", pattern: '*.consensus', mode: 'copy'
	module 'clc-assembly-cell/5.1.0'
	input : 
	set sampleId, file(ref), file(reads), file(casfile) from First_mapping 
	output:
	set sampleId, file(reads), file("*.consensus") into First_consensus
	shell:
	'''
	clc_extract_consensus -a !{casfile} -o !{sampleId}.consensus
	'''
}

// mapping on the consensus 
process clc_assembly_bis {
	//publishDir "${resdir}", pattern: '*.cas', mode: 'copy'
	module 'clc-assembly-cell/5.1.0'
	input:
	set sampleId, file(reads), file(consensus) from First_consensus
	output:
	set sampleId, file(reads), file(consensus), file("extract*.cas") into Second_mapping 

	shell:
	'''
	clc_mapper -o extract_!{sampleId}.cas -a local \
	-t 1 -r random -l 0.9 -s 0.9 -x 2 -g 3 -e 3 -d !{consensus} \
	 -q -p fb ss 70 1200 -i !{sampleId}_1P.fq.gz !{sampleId}_2P.fq.gz
	'''
}

process extract_consensus_bis{
	publishDir "${resdir}", pattern: '*.consensus', mode: 'copy'
	module 'clc-assembly-cell/5.1.0'
	input : 
	set sampleId, file(reads), file(ref), file(casfile) from Second_mapping
	output:
	set sampleId, file(reads), file("*.consensus") into Second_consensus
	shell:
	'''
	clc_extract_consensus -a !{casfile} -o extract!{sampleId}.consensus
	'''
}

// mapping on the second consensus 
process clc_assembly_final {
	//publishDir "${resdir}", pattern: '*.cas', mode: 'copy'
	module 'clc-assembly-cell/5.1.0'
	input:
	set sampleId, file(reads), file(consensus) from Second_consensus
	output:
	set sampleId, file(reads), file(consensus), file("final*.cas") into Final_mapping
	shell:
	'''
	clc_mapper -o final_!{sampleId}.cas -a local \
	-t 1 -r random -l 0.9 -s 0.9 -x 2 -g 3 -e 3 -d !{consensus} \
	 -q -p fb ss 70 1200 -i !{sampleId}_1P.fq.gz !{sampleId}_2P.fq.gz
	'''
}

process clc_info_final {
	publishDir "${resdir}", pattern: 'final*',  mode: 'copy'
	module 'clc-assembly-cell/5.1.0'
	input:
	set sampleID, file(reads), file(consensus),file (casfile) from Final_mapping
	output : 
	set sampleID, file(reads), file (consensus), file (casfile), file("final*.mapping_info"), file("final_*.bam") into Clc_info
	shell:
	'''
	clc_mapping_info -c -p fb ss 70 1200 -e 10 !{casfile} > final_!{sampleID}.mapping_info
	clc_cas_to_sam -a !{casfile} -o final_!{sampleID}.bam -f 33 -u
	'''
}

process bam_sorting_final{
	publishDir "${resdir}",pattern: 'final*.sorted.ba*',  mode: 'copy'
	module 'samtools/1.3'
	input:
	set sampleID, file(reads), file(consensus), file(casfile), file(info), file(bam) from Clc_info
	output: 
	set sampleID, file(consensus),  file(casfile), file("final*sorted.bam"), file("final*sorted.bam.bai") into LofreqChannel, VphaserChannel, IvarChannel
	shell:
	'''
	samtools sort -@ 4 !{bam} -o final!{sampleID}.sorted.bam
	samtools index final!{sampleID}.sorted.bam 
	'''
}

process consensus_final {
	publishDir "${resdir}",  mode: 'copy'
	module 'samtools/1.3'
	
	input:
	set sampleId, file(consensus), file(casfile), file(sortedbam), file(sortedbai) from IvarChannel
	output:
	file("*cns.fq") into Subscribeivar 
	shell:
	'''
	samtools flagstat !{sortedbam} > stats
	samtools mpileup -uf !{consensus} !{sortedbam} | bcftools call -c | vcfutils.pl vcf2fq > !{sampleId}cns.fq
	'''

}


"""

// low frequency variant extraction
process lofreqcalling{
	module 'samtools/1.3:lofreq/2.1.3.1'
	publishDir "${resdir}", mode: 'copy'
	input : 
	set sampleId, file(consensus), file(casfile), file(sortedbam), file(sortedbai) from LofreqChannel
	output:
	file("lofreq*") into Subscribelofreq
	shell:
	'''
	lofreq call -f !{consensus} -o lofreq!{sampleId}.vcf !{sortedbam}
	'''
}

process vphasercalling{
	module 'VPhaser-2'
	publishDir "${resdir}", mode: 'copy'
	input: 
	set sampleId, file(consensus), file(casfile), file(sortedbam), file(sortedbai) from VphaserChannel
	output:
	file("vphaser_calling*/*") into Subscribevphaser
	shell:
	'''
	mkdir vphaser_calling!{sampleId}
	variant_caller -i !{sortedbam} -o vphaser_calling!{sampleId}
	'''
}

process ivarcalling{
	module 'ivar/1.0:samtools/1.3'
	publishDir "${resdir}", mode: 'copy'
	input:
	set sampleId, file(consensus), file(casfile), file(sortedbam), file(sortedbai) from IvarChannel
	output:
	file("ivar*Qual*") into Subscribeivar 
	shell:
	'''
	#classique
	samtools mpileup -A -d 600000 --reference !{consensus} -B -Q 0 !{sortedbam} | ivar variants -p ivar_!{sampleId}Qual20 -q 20 -t 0.02 
	#all quality and freq 1% 
	samtools mpileup -A -d 600000 --reference !{consensus} -B -Q 0 !{sortedbam} | ivar variants -p ivar_!{sampleId}Qual0 -q 0 -t 0.01 
	#all quality and freq 0.1%
	samtools mpileup -A -d 600000 --reference  !{consensus} -B -Q 0 !{sortedbam} | ivar variants -p ivar_!{sampleId}Qual0max -q 0 -t 0.001
	'''
}


"""


