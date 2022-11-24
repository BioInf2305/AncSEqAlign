/*
*  samtools modules and workflows
*  if run with the option "with-singularity" then the tool will be run in container, otherwise it will create conda environment. 
*/

params.CONTAINER = "biocontainers/samtools:v1.9-4-deb_cv1"

/*
 * sort bam files by read name
 */


process runMultiBamqc {
    //container params.CONTAINER
    conda "$baseDir/conda/multiqc.yaml"
    publishDir(params.multiqcBam, pattern:"multiqc_*",mode:"move")

    tag { "runMultiqc_bam" }
    label "oneCpu"
    							
    input:
        path(bamQcFolders)

    output:									
	file "multiqc_report.html"
	file "multiqc_data"

    script:									
	
	"""	

	multiqc .

    	"""
}



workflow RUNMULTIBAMQC {
 
    take: 
    	bamQcFolders
    
    main:
	runMultiBamqc(bamQcFolders)
    
}
