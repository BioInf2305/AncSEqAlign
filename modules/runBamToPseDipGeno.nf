/*
*  samtools modules and workflows
*  if run with the option "with-singularity" then the tool will be run in container, otherwise it will create conda environment. 
*/

params.CONTAINER = "biocontainers/samtools:v1.9-4-deb_cv1"

/*
 * sort bam files by read name
 */





process samtoolsFaidx{
    container params.CONTAINER
    conda "$baseDir/conda/samtools.yaml"
    tag { "faidx_samtools" }
    label "OneCpu"
    							
    input:
	tuple val(refPrefix), path(refFile)

    output:									
   	path("*.fai")

    script:
										

	"""

	samtools faidx ${refFile}


    	"""
}

process samtoolsBamIndex{
    container params.CONTAINER
    conda "$baseDir/conda/samtools.yaml"
    tag { "index_samtools" }
    label "OneCpu"
    							
    input:
	path(trimFiltBam)

    output:									
   	path("*.bai")

    script:
										

	"""

	samtools index ${trimFiltBam}


    	"""
}

process picardIndex{
    container params.CONTAINER
    conda "$baseDir/conda/picard.yaml"
    tag { "createSeqDict_picard" }
    label "oneCpu"
    							
    input:
	tuple val(refPrefix), path(refFile)

    output:									
   	path("*.dict")

    script:
										
 	def sampleName = params.sampleName
	def tmpSamtoolsFolder = params.tmpSamtoolsFolder

	"""

	if [ ! -d ${tmpSamtoolsFolder}/${sampleName} ];then mkdir ${tmpSamtoolsFolder}/${sampleName};fi

	picard CreateSequenceDictionary TMP_DIR=${tmpSamtoolsFolder}/${sampleName} R=${refFile} O=${refPrefix}.dict
	
	if [ -d ${tmpSamtoolsFolder}/${sampleName} ];then rm -r ${tmpSamtoolsFolder}/${sampleName};fi

    	"""
}


process gatkPileup{
    container params.CONTAINER
    conda "$baseDir/conda/gatk4.yaml"
    tag { "pileup_gatk" }
    label "oneCpu"
    							
    input:
	path(trimFiltBam)
	path(bamIndex)
	tuple val(refPrefix), path(refFile)
	path(picardDict)
	path(faIdx)

    output:									
   	path("*.pileup.txt")

    script:
										
 	def sampleName = params.sampleName
	def tmpSamtoolsFolder = params.tmpSamtoolsFolder

	"""

	if [ ! -d ${tmpSamtoolsFolder}/${sampleName} ];then mkdir ${tmpSamtoolsFolder}/${sampleName};fi

	gatk --java-options "-Xmx${task.memory.toGiga()}G -Xms${task.memory.toGiga()}G" Pileup -R ${refFile} -I ${trimFiltBam} -O ${sampleName}.pileup.txt
	
	if [ -d ${tmpSamtoolsFolder}/${sampleName} ];then rm -r ${tmpSamtoolsFolder}/${sampleName};fi

    	"""
}


process pileupToPed{
    tag { "pileup_gatk" }
    label "oneCpu"    							
    publishDir(params.pedResultsOut, pattern:"*.{ped,map}",mode:"move")

    input:
	path(pileupIn)

    output:									
   	path("*.{ped,map}")

    script:
										
 	def sampleName  = params.sampleName
	def minBaseQual = params.minBaseQual

	"""

	$baseDir/bin/pileup2plink ${pileupIn} ${minBaseQual} ${sampleName}
	

    	"""
}

workflow RUNBAMTOPSEDIPGENO {
 
    take: 
	trimFiltBam
	refTupleGatk
    
    main:
	bamIndex = samtoolsBamIndex(trimFiltBam)
	if (params.createGatkIndex == "Yes"){
		picardDict = picardIndex(refTupleGatk)
		faIdx      = samtoolsFaidx(refTupleGatk)
		pileupOut  = gatkPileup(trimFiltBam, bamIndex, refTupleGatk, picardDict, faIdx)
		pileupToPed(pileupOut)
		
	}
}



