process samtoolsFaidx{
    container params.CONTAINER
    conda "$baseDir/conda/samtools.yaml"
    tag { "faidx_samtools" }
    label "OneCpu"
    							
    input:
	    path(refFile)

    output:									
   	    path("*.fai")

    script:
        refBaseName = refFile.baseName

	"""

	samtools faidx ${refFile}


    """
}


workflow CREATEREFINDICES {
 
    take: 
	    refFile
    main:
        faIdx = samtoolsFaidx(refFile)
    emit:
        faIdx
}
