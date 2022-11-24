/*
*  samtools modules and workflows
*  if run with the option "with-singularity" then the tool will be run in container, otherwise it will create conda environment. 
*/

params.CONTAINER = "biocontainers/samtools:v1.9-4-deb_cv1"

/*
 * sort bam files by read name
 */

process runBwaAlnAlignment {
    container params.CONTAINER
    conda "$baseDir/conda/bwa.yaml"
    tag { "runBwa_${sampleName}" }
    label "sixteenCpus"
    							
    input:
    	path(rawFastqFiles)
        path(refBwaIdx)
        path(reference)

    output:									
   	    path("*.sam")

    script:									
	
 	    def seedLengthValue = params.seedLengthValue
 	    def fastqBaseName              = rawFastqFiles.baseName
            sampleName                 = fastqBaseName.split("\\.")[0].split("_")[0]+"_"+fastqBaseName.split("\\.")[0].split("_")[1]
            def finalSampleName = fastqBaseName.split("\\.")[0].split("_")[0] 


	"""
	
	bwa aln -t ${task.cpus} -l ${seedLengthValue} ${reference} ${rawFastqFiles} > ${sampleName}.sai

	bwa samse -r '@RG\\tID:${sampleName}\\tSM:${finalSampleName}\\tPL:ILLUMINA' ${reference} ${sampleName}.sai ${rawFastqFiles} > ${sampleName}.sam

    	"""
}


process convertSamToBam{
    container params.CONTAINER
    conda "$baseDir/conda/samtools.yaml"
    tag { "runSamtools_${sampleName}" }
    label "twentyCpus"
    							
    input:
        path(samFile)

    output:									
   	    tuple val(sampleName), path("*.sorted.bam"),path("*.sorted*.bai")
    script:
						
 	    def samBaseName              = samFile.baseName
            sampleName                 = samBaseName.split("\\.")[0].split("_")[0]+"_"+samBaseName.split("\\.")[0].split("_")[1]
	    def tmpSamtoolsFolder = params.tmpSamtoolsFolder
            //finalSampleName = samBaseName.split("\\.")[0].split("_")[0] 
 	
	"""
	if [ ! -d ${tmpSamtoolsFolder}/${sampleName} ];then mkdir ${tmpSamtoolsFolder}/${sampleName};fi

	samtools view -@ ${task.cpus} -O BAM -o ${sampleName}.bam ${samFile} 
	
	samtools sort -T ${tmpSamtoolsFolder}/${sampleName} -O BAM -o ${sampleName}.sorted.bam -m 1G -@ ${task.cpus} ${sampleName}.bam

	samtools index ${sampleName}.sorted.bam


	if [ -d ${tmpSamtoolsFolder}/${sampleName} ];then rm -r ${tmpSamtoolsFolder}/${sampleName};fi

    	"""
}


process removeDuplicates{
    container params.CONTAINER
    conda "$baseDir/conda/picard.yaml"
    tag { "removeDuplicate_picard" }
    label "oneCpu"
    							
    input:
        tuple val(sampleName), path(sortedMergedBamFile), path(soredMergedBamIndex)
    output:									
        tuple val(finalSampleName), path("*rmDup.bam"), path("*rmDup*.bai")

    script:
										
	    def tmpSamtoolsFolder = params.tmpSamtoolsFolder
            finalSampleName = sampleName.split("\\_")[0]

	"""

	if [ ! -d ${tmpSamtoolsFolder}/${sampleName} ];then mkdir ${tmpSamtoolsFolder}/${sampleName};fi

	picard "-Xmx${task.memory.toGiga()}G " MarkDuplicates TMP_DIR=${tmpSamtoolsFolder}/${sampleName} I=${sortedMergedBamFile} O=${sampleName}.sorted.merged.rmDup.bam AS=true REMOVE_DUPLICATES=true METRICS_FILE=${sampleName}.rmDupMetrics.txt VALIDATION_STRINGENCY=LENIENT

	picard BuildBamIndex I=${sampleName}.sorted.merged.rmDup.bam
	
	if [ -d ${tmpSamtoolsFolder}/${sampleName} ];then rm -r ${tmpSamtoolsFolder}/${sampleName};fi

    	"""
}

process mergeSortedBam{
    container params.CONTAINER
    conda "$baseDir/conda/samtools.yaml"
    tag { "runSamtools_merge" }
    label "twentyCpus"
    							
    input:
        tuple val(sampleName), path(sortedBamFiles), path(bamIndexFiles)

    output:									
   	    tuple val(sampleName), path("*.sorted.merged.bam"), path("*.sorted.merged*.bai")

    script:
										
	    def tmpSamtoolsFolder = params.tmpSamtoolsFolder

	"""

	if [ ! -d ${tmpSamtoolsFolder}/${sampleName} ];then mkdir ${tmpSamtoolsFolder}/${sampleName};fi

	samtools merge -@ ${task.cpus} -o ${sampleName}.bam ${sortedBamFiles}

	samtools sort -T ${tmpSamtoolsFolder}/${sampleName} -m 1G -O BAM -o ${sampleName}.sorted.merged.bam -@ ${task.cpus} ${sampleName}.bam

	samtools index ${sampleName}.sorted.merged.bam
	
	if [ -d ${tmpSamtoolsFolder}/${sampleName} ];then rm -r ${tmpSamtoolsFolder}/${sampleName};fi

    	"""
}


process filterBam{
    container params.CONTAINER
    conda "$baseDir/conda/samtools.yaml"
    publishDir(params.filterBamOut, pattern:"*.sorted.merged.rmDup.filt.{bam,bam.bai}",mode:"copy")
    tag { "filterBam" }
    label "sixteenCpus"
    							
    input:
    	tuple val(sampleName), path(sortedMergedRmDupBamFile), path(sortedMergedRmDupBamIdx)

    output:									
   	    tuple val(sampleName), path("*.sorted.merged.rmDup.filt.bam"), path("*.sorted.merged.rmDup.filt.*.bai")
	

    script:
										
	    def tmpSamtoolsFolder = params.tmpSamtoolsFolder
	    def minMapQ           = params.minMapQ

	"""

	if [ ! -d ${tmpSamtoolsFolder}/${sampleName} ];then mkdir ${tmpSamtoolsFolder}/${sampleName};fi

	samtools index ${sortedMergedRmDupBamFile}

	samtools view -@ ${task.cpus} -F4 -q ${minMapQ} -O BAM -o ${sampleName}.sorted.merged.rmDup.filt.bam ${sortedMergedRmDupBamFile}

	samtools index ${sampleName}.sorted.merged.rmDup.filt.bam

	if [ -d ${tmpSamtoolsFolder}/${sampleName} ];then rm -r ${tmpSamtoolsFolder}/${sampleName};fi

    	"""
}

process trimBam{
    container params.CONTAINER
    conda "$baseDir/conda/bamutil.yaml"
    tag { "trimmingBam" }
    label "oneCpu"
    publishDir(params.filterBamOut, pattern:"*.sorted.merged.rmDup.filt.trimmed.bam",mode:"copy")
    							
    input:
        tuple val(sampleName), path(filteredBamFile), path(filteredBamIndex)

    output:									
   	    tuple val(sampleName), path("*.sorted.merged.rmDup.filt.trimmed.bam")
	

    script:
										
	    def numTrimBase       = params.numTrimBase

	"""
	
	    bam trimBam ${filteredBamFile} ${sampleName}.sorted.merged.rmDup.filt.trimmed.bam ${numTrimBase}

    	"""
}

process trimBamIndex {
    container params.CONTAINER
    conda "$baseDir/conda/samtools.yaml"
    tag { "runBwa_${sample}" }
    label "oneCpu"
    publishDir(params.filterBamOut, pattern:"*.sorted.merged.rmDup.filt.trimmed.bam.bai",mode:"copy")
    							
    input:
	    tuple val(sample), path(bamFile)

    output:									
   	    tuple val(sample), path("*.bai")

    script:									
	
	"""
	
	samtools index ${bamFile}

    """
}

workflow RUNALIGNMENT {
 
    take: 
	rawFastqFiles
        refBwaIdx
        reference
    
    main:
        rawSamFilesTuple = runBwaAlnAlignment( rawFastqFiles, refBwaIdx, reference )
        bamAndIndexFiles = convertSamToBam( rawSamFilesTuple )
        remDupBamFile = removeDuplicates( bamAndIndexFiles )
        remDupBamIdxGrouped = remDupBamFile.groupTuple()
        mergedBamIdx = mergeSortedBam( remDupBamIdxGrouped )
        filteredBamOut = filterBam( mergedBamIdx )
        filtTrimBam = trimBam( filteredBamOut )
        filtTrimIdxAndBam = trimBamIndex(filtTrimBam).combine(filtTrimBam, by:0)
   emit:
	filtTrimIdxAndBam = filtTrimIdxAndBam
	
}
