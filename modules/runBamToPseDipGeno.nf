/*
*  samtools modules and workflows
*  if run with the option "with-singularity" then the tool will be run in container, otherwise it will create conda environment. 
*/

params.CONTAINER = "biocontainers/samtools:v1.9-4-deb_cv1"

/*
 * sort bam files by read name
 */




process gatkPileup{
    container params.CONTAINER
    conda "$baseDir/conda/gatk4.yaml"
    tag { "pileup_gatk" }
    label "oneCpu"
    							
    input:
	    tuple val(sample), path(trimFiltBamIdx), path(trimFiltBam)
	    path(faIdx)
        path(faDict)
        path(reference)

    output:									
   	    tuple val(sample), path("*.pileup.txt")

    script:
										
	    def tmpSamtoolsFolder = params.tmpSamtoolsFolder

	"""

	if [ ! -d ${tmpSamtoolsFolder}/${sample} ];then mkdir ${tmpSamtoolsFolder}/${sample};fi

	gatk --java-options "-Xmx${task.memory.toGiga()}G -Xms${task.memory.toGiga()}G" Pileup -R ${reference} -I ${trimFiltBam} -O ${sample}.pileup.txt
	
	if [ -d ${tmpSamtoolsFolder}/${sample} ];then rm -r ${tmpSamtoolsFolder}/${sample};fi

    	"""
}


process pileupToPed{
    tag { "pileup_gatk" }
    label "oneCpu"    							
    publishDir(params.pedResultsOut, pattern:"*.{ped,map}",mode:"move")

    input:
	tuple val(sample), path(pileupIn)

    output:									
   	tuple val(sample), path("*.{ped,map}")

    script:
										
	def minBaseQual = params.minBaseQual

	"""

	$baseDir/bin/pileup2plink ${pileupIn} ${minBaseQual} ${sample}
	

    	"""
}

workflow RUNBAMTOPSEDIPGENO {
 
    take: 
        trimFiltBam
        faIdx
        faDict
        reference
    
    main:

	    pileupOut  = gatkPileup(trimFiltBam, faIdx, faDict, reference)
	    pileupToPed(pileupOut)
}



