nextflow.enable.dsl=1
//"/pasteur/zeus/projets/p02/GEVA/To_sort/GEVAmetagenomics/users/Marie/Data/raw/NGS/"
//"/pasteur/zeus/projets/p02/GEVA/To_sort/GEVAmetagenomics/users/Marie/Data/raw/NGS/2020-10-29_lastHCV/"
//"/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/GEVA/Data/raw/hybridemapping/control"


params.refsequence = "/pasteur/zeus/services/p01/banques-prod/prod/rel/hg38/current/fasta/3.6/hg38.fa"
params.pairs = "/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/GEVA/Data/raw/hybridemapping/treatment_failure/NGS/HCV_01*_R{1,2}.fastq.gz"

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


//mapping not too stringeant on the human genome
process clc_assembly {
	publishDir "${resdir}", pattern: '*.cas', mode: 'copy'
	module 'clc-assembly-cell/5.1.0'
	input:
	set val(sampleId), file(reads) from Cleanreads
	file(ref) from refsequence
	output:
	set val(sampleId), file(reads), file("*.cas") into First_mapping

	shell:
	'''
	clc_mapper -o !{sampleId}.cas -a local \
	-t 1 -r random -l 0.9 -s 0.9 -x 2 -g 3 -e 3 -d !{ref} \
	 -q -p fb ss 70 1200 -i !{sampleId}_1P.fq.gz !{sampleId}_2P.fq.gz

	'''
}


process clc_info_final {
	publishDir "${resdir}", pattern: 'final*',  mode: 'copy'
	module 'clc-assembly-cell/5.1.0'
	input:
	set val(sampleID), file(reads), file (casfile) from First_mapping
    file(ref) from refsequence
	output : 
	set val(sampleID), file (casfile), file("final*.mapping_info"), file("final_*.bam") into Clc_info
	shell:
	'''
	clc_mapping_info -c -p fb ss 70 1200 -e 10 !{casfile} > final_!{sampleID}.mapping_info
	clc_cas_to_sam -a !{casfile} -o final_!{sampleID}.bam -f 33 -u
	'''
}

process bam_index{
	module 'samtools/1.4'
	input:
	set val(sampleID), file(casfile), file(info), file(bam) from Clc_info
	output: 
	set val(sampleID), file("*sorted.bam"), file("*sorted.bam.bai") into bam_files_paired_map_map, bam_files_paired_unmap_unmap, bam_files_paired_unmap_map, bam_files_paired_map_unmap
    //samtools sort -@ 4 !{bam} -o final!{sampleID}.sorted.bam
	shell:
	'''
    samtools sort -@ 4 !{bam} -o final!{sampleID}.sorted.bam
    samtools index final!{sampleID}.sorted.bam
	'''
}

/*
 * Step 2a: Handle paired-end bams
 */
process pairedEndMapMap{
    module 'samtools/1.4'
    input:
    set val(sampleID), file(bam), file(bai) from bam_files_paired_map_map

    output:
    set val(sampleID), file( '*.map_map.bam') into map_map_bam

    script:
    """
    samtools view -b -f1 -F12 $bam -@ 4 -o ${sampleID}.map_map.bam
    """
}

process pairedEndUnmapUnmap{
  module 'samtools/1.4'
  input:
  set val(sampleID), file(bam), file(bai) from bam_files_paired_unmap_unmap

  output:
  set val(sampleID), file('*.unmap_unmap.bam') into unmap_unmap_bam

  script:
  """
  samtools view -b -f12 -F256 $bam --@ 4 -o ${sampleID}.unmap_unmap.bam
  """
}

process pairedEndUnmapMap{
  module 'samtools/1.4'
  input:
  set val(sampleID), file(bam), file(bai) from bam_files_paired_unmap_map

  output:
  set val(sampleID), file( '*.unmap_map.bam') into unmap_map_bam


  script:
  """
  samtools view -b -f4 -F264 $bam -@ 4 -o ${sampleID}.unmap_map.bam
  """
}

process pairedEndMapUnmap{
  module 'samtools/1.4'
  input:
  set val(sampleID), file(bam), file(bai) from bam_files_paired_map_unmap

  output:
  set val(sampleID), file( '*.map_unmap.bam') into map_unmap_bam

  script:
  """
  samtools view -b -f8 -F260 $bam  -@ 4 -o ${sampleID}.map_unmap.bam
  """
}

unmap_unmap_bam.join(map_unmap_bam, remainder: true)
                .join(unmap_map_bam, remainder: true)
                .set{ all_unmapped_bam }

process mergeUnmapped{
  module 'samtools/1.4'
  input:
  set val(sampleID), file(unmap_unmap), file (map_unmap),  file(unmap_map) from all_unmapped_bam

  output:
  set val(sampleID), file('*.merged_unmapped.bam') into merged_unmapped

  script:
  """
  samtools merge ${sampleID}.merged_unmapped.bam $unmap_unmap $map_unmap $unmap_map  -@ 4
  """
}

process sortExtractMapped{
  module 'samtools/1.4'

  publishDir "${params.outdir}/reads", mode: 'copy',
        saveAs: { filename ->
            if (filename.indexOf(".fq.gz") > 0) filename
            else null
        }

  input:
  set val(sampleID), file(all_map_bam) from map_map_bam

  output:
  set val(sampleID), file('*_mapped.fq.gz') into reads_mapped

  script:  
  """
  samtools collate -Ou -@,4 $all_map_bam \
    | samtools fastq -@,4 -1 ${sampleID}_R1_mapped.fq.gz -2 ${sampleID}_R2_mapped.fq.gz -s ${sampleID}_mapped_singletons.fq.gz -n 
  """
}

process sortExtractUnmapped{
  publishDir "${params.outdir}/reads", mode: 'copy',
        saveAs: { filename ->
            if (filename.indexOf(".fq.gz") > 0) filename
            else null
        }
  module 'samtools/1.4'
  input:
  set val(sampleID), file(all_unmapped) from merged_unmapped

  output:
  set val(sampleID), file('*_unmapped.fq.gz') into reads_unmapped

  script:  
  """
  samtools collate -Ou $all_unmapped -@,8 \
      | samtools fastq -@,8 -1 ${sampleID}_R1_unmapped.fq.gz -2 ${sampleID}_R2_unmapped.fq.gz -s ${sampleID}_unmapped_singletons.fq.gz -n
  """
}

reads_mapped.join(reads_unmapped, remainder: true)
            .map{
              row -> tuple(row[0], row[1][0], row[1][1], row[2][0], row[2][1])
            }
            .set{ all_fastq }

process joinMappedAndUnmappedFastq{
  publishDir "${params.outdir}/reads", mode: 'copy',
        saveAs: { filename ->
            if (filename.indexOf(".fq.gz") > 0) filename
            else null
        }

  input:
  set val(sampleID), file(mapped_fq1), file(mapped_fq2), file(unmapped_fq1), file(unmapped_fq2) from all_fastq.filter{ it.size()>0 }

  output:
  set file('*1.fq.gz'), file('*2.fq.gz') into read_qc


  script:
  """
  cat $unmapped_fq1 >> $mapped_fq1
  mv $mapped_fq1 ${sampleID}.1.fq.gz
  cat $unmapped_fq2 >> $mapped_fq2
  mv $mapped_fq2 ${sampleID}.2.fq.gz
  """
}






