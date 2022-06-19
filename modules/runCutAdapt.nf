/*
*  samtools modules and workflows
*  if run with the option "with-singularity" then the tool will be run in container, otherwise it will create conda environment. 
*/

params.CONTAINER = "biocontainers/samtools:v1.9-4-deb_cv1"

/*
 * sort bam files by read name
 */


process filterCutadapt {
    container params.CONTAINER
    conda "$baseDir/conda/cutadapt.yaml"
    //publishDir(params.fqOutputDir, pattern:"*trimmed".gz",mode:"cp")
    tag { "cutadapt_${sample}" }
    label "oneCpu"
    							
    input:
    	tuple val(sample), path(rawFastqFile)

    output:									
   	tuple val(sample), path("*trimmed*")

    script:									
	
	def illuAdapters = params.illuAdapters
	def minAdapterOverlap = params.minAdapterOverlap
	def minReadLength = params.minReadLength
	def command = ""
	command = command + " " + "-a "+ params.illuAdapters +" -O "+ params.minAdapterOverlap + " -m "+ minReadLength + " "
	

 	"""
	
	cutadapt ${command} -o ${sample}_cutadapt_trimmed.O${params.minAdapterOverlap}m${minReadLength}.fastq.gz ${rawFastqFile} 


    	"""
}


workflow FILTERFASTQ {
 
    take: 
    	fastqFiles
    
    main:
	filteredFastqFiles = filterCutadapt(fastqFiles)

    emit:
	filteredFastqFiles = filteredFastqFiles
}



