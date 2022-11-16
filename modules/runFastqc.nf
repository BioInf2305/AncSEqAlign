/*
*  samtools modules and workflows
*  if run with the option "with-singularity" then the tool will be run in container, otherwise it will create conda environment. 
*/

params.CONTAINER = "biocontainers/samtools:v1.9-4-deb_cv1"

/*
 * sort bam files by read name
 */


process runFastqc {
    container params.CONTAINER
    conda "$baseDir/conda/fastqc.yaml"
    publishDir(params.fastqcResultsOut, pattern:"*.html",mode:"move")
    tag { "runFastqc_${sample}" }
    label "oneCpu"
    							
    input:
        path(rawFastqFile)

    output:									
        path("*.html")

    script:									
	
 	    def fastqBaseName              = rawFastqFile.baseName
        def sampleName                 = fastqBaseName.split("\\.")[0]
	    def minLength                  = params.minLength
	    def contaminantsFile           = params.contaminantsFile
	    def adaptersFile               = params.adaptersFile
	    def limitsFile                 = params.limitsFile
	    def kmers                      = params.kmers
	    def fastqcTmpDir               = params.fastqcTmpDir
	command = "" 
	
	
	if ( minLength != "0" ) {
	
	command = command + " --min_length "+minLength
	
	}
	
	if ( contaminantsFile != "No" ){
	
	command = command + " -c "+contaminantsFile
	
	}

	if ( adaptersFile != "No" ){
	
	command = command + " -a "+adaptersFile
	
	}

	if ( limitsFile != "No" ){
	
	command = command + " -l "+limitFile
	
	}

	
	command = command + " -k "+kmers+ " -d "+ fastqcTmpDir
	

 	"""
	
	fastqc ${rawFastqFile} ${command}

    	"""
}



workflow RUNFASTQC {
 
    take: 
    	fastqFiles
    
    main:
	runFastqc(fastqFiles)
}



