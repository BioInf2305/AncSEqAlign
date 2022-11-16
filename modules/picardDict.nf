
process picardIndex{
    container params.CONTAINER
    conda "$baseDir/conda/picard.yaml"
    tag { "createSeqDict_picard" }
    label "oneCpu"
    							
    input:
        path(refFile)

    output:									
   	    path("*.dict")

    script:
        def refBaseName = refFile.baseName
	    def tmpSamtoolsFolder = params.tmpSamtoolsFolder

	"""

	if [ ! -d ${tmpSamtoolsFolder}/${refBaseName} ];then mkdir ${tmpSamtoolsFolder}/${refBaseName};fi

	picard CreateSequenceDictionary TMP_DIR=${tmpSamtoolsFolder}/${refBaseName} R=${refFile} O=${refBaseName}.dict
	
	if [ -d ${tmpSamtoolsFolder}/${refBaseName} ];then rm -r ${tmpSamtoolsFolder}/${refBaseName};fi

    	"""
}

workflow CREATEREFDICT {
 
    take: 
	    refFile
    main:
        picardIdx = picardIndex(refFile)
    emit:
        picardIdx
}
