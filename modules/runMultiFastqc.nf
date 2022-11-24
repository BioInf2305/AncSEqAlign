/*
*  samtools modules and workflows
*  if run with the option "with-singularity" then the tool will be run in container, otherwise it will create conda environment. 
*/

params.CONTAINER = "biocontainers/samtools:v1.9-4-deb_cv1"

/*
 * sort bam files by read name
 */


process runMultiFastqc {
    conda "$baseDir/conda/multiqc.yaml"
    publishDir(params.multiqcFastqc, pattern:"*.html",mode:"move")

    tag { "runMultiqc_fastq" }
    label "oneCpu"
    							
    input:
        path(fastqcZip)

    output:									
	file "multiqc_report.html"
	file "multiqc_data"

    script:									
	
	"""	

	multiqc .

    	"""
}



workflow RUNMULTIFASTQC {
 
    take: 
    	fastqcZip
    
    main:
	runMultiFastqc(fastqcZip)
    
}
