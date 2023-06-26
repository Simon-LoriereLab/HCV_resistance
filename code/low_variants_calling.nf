
params.path = "/pasteur/homes/mamorel/GEVA/Marie/Data/assembly/mapping/second_mapping/others/"
params.sampleId = "HCV128_S26"

params.consensus = params.path+"extract"+params.sampleId+".consensus"
params.casfile = params.path+"final_"+params.sampleId+".cas"
params.sortedbam = params.path+"final"+params.sampleId+".sorted.bam"
params.sortedbai = params.path+"final"+params.sampleId+".sorted.bam.bai"
params.resdir = "/pasteur/homes/mamorel/GEVA/Marie/Data/variant_calling/others/"


resdir=file(params.resdir)
resdir.with {mkdirs()}

sampleId = params.sampleId

consensus = file(params.consensus)
casfile = file(params.casfile)
sortedbam = file(params.sortedbam)
sortedbai = file(params.sortedbai)


// low frequency variant extraction
process lofreqcalling{
	module 'samtools/1.3:lofreq/2.1.3.1'
	publishDir "${resdir}", mode: 'copy'
	input : 
    val sampleId
    file consensus 
    file casfile
    file sortedbam
    file sortedbai
	output:
	file("lofreq*") into Subscribelofreq
	shell:
	'''
	lofreq call -f !{consensus} -o lofreq!{sampleId}.vcf !{sortedbam}
	'''
}


process ivarcalling{
	module 'ivar/1.0:samtools/1.3'
	publishDir "${resdir}", mode: 'copy'
	input:
    val sampleId
    file consensus 
    file casfile
    file sortedbam
    file sortedbai
    //val sampleId
	//set file(consensus), file(casfile), file(sortedbam), file(sortedbai) from IvarChannel
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