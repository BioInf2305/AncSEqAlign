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
    tag { "cutadapt_${rawFastqFile}" }
    label "oneCpu"
    							
    input:
        path(rawFastqFile)

    output:									
   	    path("*trimmed*")

    script:									
	
 	    def fastqBaseName = rawFastqFile.baseName
        def sampleName = fastqBaseName.split("\\.")[0]
        println sampleName
        

	def illuAdapters = params.illuAdapters
	def minAdapterOverlap = params.minAdapterOverlap
	def minReadLength = params.minReadLength
	def command = ""
	command = command + " " + "-a "+ params.illuAdapters +" -O "+ params.minAdapterOverlap + " -m "+ minReadLength + " "
	

 	"""
	
	cutadapt ${command} -o ${sampleName}_cutadapt_trimmed_O${params.minAdapterOverlap}_m${minReadLength}.fastq.gz ${rawFastqFile} 


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



