/*
*  samtools modules and workflows
*  if run with the option "with-singularity" then the tool will be run in container, otherwise it will create conda environment. 
*/

params.CONTAINER = "biocontainers/samtools:v1.9-4-deb_cv1"

/*
 * sort bam files by read name
 */

process runBwaIndex {
    container params.CONTAINER
    conda "$baseDir/conda/bwa.yaml"
    tag { "runBwa_${sample}" }
    label "oneCpu"
    							
    input:
	tuple val(refPrefix), path(refFiles)

    output:									
   	tuple val(refPrefix), path("*.{amb,ann,bwt,pac,sa}")

    script:									
	
	"""
	bwa index -p ${refPrefix} ${refFiles}

    	"""
}


process runBwaAlnAlignment {
    container params.CONTAINER
    conda "$baseDir/conda/bwa.yaml"
    tag { "runBwa_${idTag}" }
    label "sixteenCpus"
    							
    input:
    	tuple val(idTag), path(rawFastqFiles), val(refPrefix), path(refFiles)

    output:									
   	tuple val(idTag), path("*.sam")

    script:									
	
 	def seedLengthValue = params.seedLengthValue
	def sampleName      = params.sampleName

	"""
	
	bwa aln -t ${task.cpus} -l ${seedLengthValue} ${refPrefix} ${rawFastqFiles} > ${idTag}.sai

	bwa samse -r '@RG\\tID:${idTag}\\tSM:${sampleName}\\tPL:ILLUMINA' ${refPrefix} ${idTag}.sai ${rawFastqFiles} > ${idTag}.sam

    	"""
}


process convertSamToBam{
    container params.CONTAINER
    conda "$baseDir/conda/samtools.yaml"
    tag { "runSamtools_${sample}" }
    label "sixteenCpus"
    							
    input:
    	tuple val(sample), path(samFile)

    output:									
   	path("*.sorted.bam"), emit: sortedBamFile
	path("*.bai"), emit: bamIndexFile

    script:
										
	def tmpSamtoolsFolder = params.tmpSamtoolsFolder
 	
	"""
	if [ ! -d ${tmpSamtoolsFolder}/${sample} ];then mkdir ${tmpSamtoolsFolder}/${sample};fi

	samtools view -@ ${task.cpus} -O BAM -o ${sample}.bam ${samFile} 
	
	samtools sort -T ${tmpSamtoolsFolder}/${sample} -O BAM -o ${sample}.sorted.bam -m "${task.memory.toGiga()}G" -@ ${task.cpus} ${sample}.bam

	samtools index ${sample}.sorted.bam

	if [ -d ${tmpSamtoolsFolder}/${sample} ];then rm -r ${tmpSamtoolsFolder}/${sample};fi

    	"""
}


process mergeSortedBam{
    container params.CONTAINER
    conda "$baseDir/conda/samtools.yaml"
    tag { "runSamtools_merge" }
    label "sixteenCpus"
    							
    input:
    	path(sortedBamFiles)
	path(bamFileIndex)

    output:									
   	path("*.sorted.merged.bam"), emit: mergedBamFile
	path("*.sorted.merged*.bai"), emit: mergedBamIndex

    script:
										
 	def sampleName = params.sampleName
	def tmpSamtoolsFolder = params.tmpSamtoolsFolder

	"""

	if [ ! -d ${tmpSamtoolsFolder}/${sampleName} ];then mkdir ${tmpSamtoolsFolder}/${sampleName};fi

	samtools merge -@ ${task.cpus} -o ${sampleName}.bam ${sortedBamFiles}

	samtools sort -T ${tmpSamtoolsFolder}/${sampleName} -m ${task.memory.toGiga()}G -O BAM -o ${sampleName}.sorted.merged.bam -@ ${task.cpus} ${sampleName}.bam

	samtools index ${sampleName}.sorted.merged.bam
	
	if [ -d ${tmpSamtoolsFolder}/${sampleName} ];then rm -r ${tmpSamtoolsFolder}/${sampleName};fi

    	"""
}

process removeDuplicates{
    container params.CONTAINER
    conda "$baseDir/conda/picard.yaml"
    tag { "removeDuplicate_picard" }
    label "oneCpu"
    							
    input:
    	path(sortedMergedBamFiles)
	path(sortedMergedBamFileIndex)

    output:									
   	path("*.sorted.merged.rmDup.bam")

    script:
										
 	def sampleName = params.sampleName
	def tmpSamtoolsFolder = params.tmpSamtoolsFolder

	"""

	if [ ! -d ${tmpSamtoolsFolder}/${sampleName} ];then mkdir ${tmpSamtoolsFolder}/${sampleName};fi

	picard "-Xmx${task.memory.toGiga()}G " MarkDuplicates TMP_DIR=${tmpSamtoolsFolder}/${sampleName} I=${sortedMergedBamFiles} O=${sampleName}.sorted.merged.rmDup.bam AS=true REMOVE_DUPLICATES=true METRICS_FILE=${sampleName}.rmDupMetrics.txt VALIDATION_STRINGENCY=LENIENT
	
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
    	path(sortedMergedRmDupBamFile)

    output:									
   	path("*.sorted.merged.rmDup.filt.bam"), emit: mergedRmDupFiltBamFile
	path("*.sorted.merged.rmDup.filt.*.bai"), emit: mergedRmDupFiltBamIndex
	

    script:
										
 	def sampleName        = params.sampleName
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
    	path(filteredBamFile)
	path(filteredBamIndex)

    output:									
   	path("*.sorted.merged.rmDup.filt.trimmed.bam")
	

    script:
										
 	def sampleName        = params.sampleName
	def numTrimBase       = params.numTrimBase

	"""
	
	bam trimBam ${filteredBamFile} ${sampleName}.sorted.merged.rmDup.filt.trimmed.bam ${numTrimBase}

    	"""
}

workflow RUNALIGNMENT {
 
    take: 
	refTupleFastqTuple
    
    main:
	rawSamFilesTuple = runBwaAlnAlignment( refTupleFastqTuple )
	convertSamToBam( rawSamFilesTuple )
	mergeSortedBam( convertSamToBam.out.sortedBamFile.collect(), convertSamToBam.out.bamIndexFile.collect() )
	remDupBamFile = removeDuplicates( mergeSortedBam.out.mergedBamFile, mergeSortedBam.out.mergedBamIndex )
	filterBam( remDupBamFile )
	filtTrimBam = trimBam( filterBam.out.mergedRmDupFiltBamFile, filterBam.out.mergedRmDupFiltBamIndex )

   emit:
	filtTrimBam = filtTrimBam
	
}



