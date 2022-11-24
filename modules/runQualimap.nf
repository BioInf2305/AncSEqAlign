/*
*  samtools modules and workflows
*  if run with the option "with-singularity" then the tool will be run in container, otherwise it will create conda environment. 
*/

params.CONTAINER = "biocontainers/samtools:v1.9-4-deb_cv1"

/*
 * sort bam files by read name
 */


process runBamQc {
    container params.CONTAINER
    conda "$baseDir/conda/qualimap.yaml"
    publishDir(params.qualimapResultsOut, mode:"copy")
    tag { "runQualimap_${sample}" }
    label "oneCpu"
    							
    input:
        tuple val(sample), path(bamIdx), path(bamFile)

    output:									
        path(sample) 

    script:									
	
	

 	"""
	
	qualimap bamqc -bam ${bamFile} -outdir ${sample}

    	"""
}



workflow RUNQUALIMAP {
 
    take: 
    	trimFiltBamFile
    
    main:
	outBamQc = runBamQc(trimFiltBamFile)

    emit:
	outBamQc = outBamQc
    
}



