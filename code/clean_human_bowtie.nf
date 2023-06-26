nextflow.enable.dsl=1
//"/pasteur/zeus/projets/p02/GEVA/To_sort/GEVAmetagenomics/users/Marie/Data/raw/NGS/"
//"/pasteur/zeus/projets/p02/GEVA/To_sort/GEVAmetagenomics/users/Marie/Data/raw/NGS/2020-10-29_lastHCV/"
//"/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/GEVA/Data/raw/hybridemapping/control"


params.refsequence = "/pasteur/zeus/services/p01/banques-prod/prod/rel/hg38/current/fasta/3.6/hg38.fa"
params.pairs = "/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/GEVA/Data/raw/hybridemapping/control/HCV*_R{1,2}.fastq.gz"

params.resdir="/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/GEVA/Data/clean_human/"

params.adapter = "/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/GEVA/Data/raw/illumina_adapter.fasta"

refsequence = file(params.refsequence)
adapter = file(params.adapter)
resdir=file(params.resdir)
resdir.with {mkdirs()}

//import files two by two according to their pair ID
//In samples_ch we get one single identifier for the pair and the two files of the pair
Channel
    .fromFilePairs(params.pairs, checkIfExists:true)
    .set { samples_ch}


//trim the fastq files and take as output the trimmed paired fq files 
process trimming {
    module 'graalvm/ce-java8-20.0.0'
	module 'Trimmomatic/0.39'
	input:
	set val(sampleId), file(reads) from samples_ch
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

process build_reference {
    module 'bowtie2/2.3.5.1'
    input: 
    file refsequence
    output: 
    file "host_DB*" into Human_ref
    script:
    """
    bowtie2-build $refsequence host_DB
    """
}

process mapping {
    module 'bowtie2/2.3.5.1'
    input: 
    file host_DB from Human_ref
    set sampleId, file(reads) from Cleanreads
    output:
    set sampleId, file("*_mapped_and_unmapped.sam") into MappedReadds
    script:
    """
    bowtie2 -p 8 -x host_DB -1 ${sampleId}_1P.fq.gz -2 ${sampleId}_2P.fq.gz -S ${sampleId}_mapped_and_unmapped.sam
    """
}

process sam2bam {
    module 'samtools/1.4'
    input: 
    set sampleId, file(sam) from MappedReadds
    output: 
    set sampleId, file("*_mapped_and_unmapped.bam") into BamFile
    script:
    """
    samtools view -bS ${sam} > ${sampleId}_mapped_and_unmapped.bam
    """
}


process filter_unmapped {
    module 'samtools/1.4'
    input:
    set sampleId, file(bam) from  BamFile
    output:
    set sampleId, file("*_bothReadsUnmapped.bam") into FilterReads
    script:
    """
    samtools view -b -f 12 -F 256 ${bam} > ${sampleId}_bothReadsUnmapped.bam 
    """
}

process split_fastq{
    module 'samtools/1.4'
    publishDir "${resdir}", mode: 'copy'
    input:
    set sampleId, file(unmaped_bam) from FilterReads
    output:
    set file("*host_removed_R2.fastq.gz"), file("*host_removed_R1.fastq.gz")
    script:
    """
    samtools sort -n -m 5G -@ 2 ${unmaped_bam} -o ${sampleId}_bothReadsUnmapped_sorted.bam
    samtools fastq -@ 8 ${sampleId}_bothReadsUnmapped_sorted.bam -1 ${sampleId}_host_removed_R1.fastq.gz -2 ${sampleId}_host_removed_R2.fastq.gz -0 /dev/null -s /dev/null -n
    """
}