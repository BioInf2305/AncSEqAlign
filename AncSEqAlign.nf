#!/usr/bin/env nextflow


nextflow.enable.dsl=2



version                 = "1.0"
// this prevents a warning of undefined parameter
params.help             = false

// this prints the input parameters

fastqFilesPath  = params.fastqFiles
refFilesPath    = params.refPrefixFiles
refFilesGatk    = params.refFilesGatk

Channel
    .fromFilePairs( fastqFilesPath, size : 1 )
    .set { fastqFileTuple }

Channel
    .fromFilePairs( refFilesPath, size: -1 )
    .set{ refFilesTuple }

Channel
    .fromFilePairs( refFilesGatk, size: -1 )
    .set{ refFilesGatkTuple }

include { FILTERFASTQ } from "${baseDir}/modules/runCutAdapt" addParams(

	illuAdapters              : params.illuAdapters,
	minAdapterOverlap         : params.minAdapterOverlap,
	minReadLength             : params.minReadLength
	
	)

include { RUNFASTQC } from "${baseDir}/modules/runFastqc" addParams(

	fastqcResultsOut          : params.fastqcResultsOut,
	minLength                 : params.minLength,
	contaminantsFile          : params.contaminantsFile,
	adaptersFile              : params.adaptersFile,
	limitsFile                : params.limitsFile,
	kmers                     : params.kmers,
	fastqcTmpDir              : params.fastqcTmpDir

	)

include { RUNALIGNMENT } from "${baseDir}/modules/runAlignment" addParams(

	createBwaIndex            : params.createBwaIndex,
	algorithmBwa              : params.algorithmBwa,
	seedLengthValue           : params.seedLengthValue,
	sampleName                : params.sampleName,
	minMapQ                   : params.minMapQ,
	filterBamOut              : params.filterBamOut
	
	)

include {RUNBAMTOPSEDIPGENO} from "${baseDir}/modules/runBamToPseDipGeno" addParams(


	createGatkIndex           : params.createGatkIndex,
	minBaseQual               : params.minBaseQual


	)

workflow {
	if ( params.filterFastq == "Yes" ){
		filteredFastqFiles = FILTERFASTQ( fastqFileTuple )
	}
	else{
		filteredFastqFiles = fastqFileTuple
	}
	RUNFASTQC(filteredFastqFiles)
	if (params.alignment == "Yes" ){
		refFilesFastqCombined = filteredFastqFiles.combine(refFilesTuple)
		trimmedFiltBamFile = RUNALIGNMENT( refFilesFastqCombined )
		if (params.pseudoDiploGeno == "Yes" ){
			RUNBAMTOPSEDIPGENO(trimmedFiltBamFile, refFilesGatkTuple)
		}
	}
}




workflow.onComplete { 
    println ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
    }
