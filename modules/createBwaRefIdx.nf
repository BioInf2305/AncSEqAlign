
process runBwaIndex {
    container params.CONTAINER
    conda "$baseDir/conda/bwa.yaml"
    tag { "runBwa_refIndex" }
    label "oneCpu"
    							
    input:
	    path(reference)

    output:									
   	    path("*.{amb,ann,bwt,pac,sa}")

    script:									
	
	    """
	        bwa index ${reference}

    	"""
}

workflow CREATEBWAREFIDX {
 
    take: 
	    refFile
    main:
        bwaFaIdx = runBwaIndex(refFile)
    emit:
        bwaFaIdx 
}
