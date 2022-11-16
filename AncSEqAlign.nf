#!/usr/bin/env nextflow


nextflow.enable.dsl=2



version                 = "1.0"
// this prevents a warning of undefined parameter
params.help             = false

// this prints the input parameters

fastqFilesPath  = params.fastqFiles
reference       = params.reference

Channel
    .fromFilePairs( fastqFilesPath, size : -1 )
    .map{it[1]}
    .flatten()
    .set{ fastqFile }



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
	minMapQ                   : params.minMapQ,
	filterBamOut              : params.filterBamOut
	
	)

include {RUNBAMTOPSEDIPGENO} from "${baseDir}/modules/runBamToPseDipGeno" addParams(


	createGatkIndex           : params.createGatkIndex,
	minBaseQual               : params.minBaseQual


	)

include {CREATEREFINDICES} from "${baseDir}/modules/samtoolsFaidx" 


include {CREATEREFDICT} from "${baseDir}/modules/picardDict" 


include {CREATEBWAREFIDX} from "${baseDir}/modules/createBwaRefIdx" 

workflow {
    
    faIdx = CREATEREFINDICES(reference)

    faDict = CREATEREFDICT(reference) 

    bwaRefIdx = CREATEBWAREFIDX(reference)

	if ( params.filterFastq == "Yes" ){
        fastqFile.view()   
		filteredFastqFiles = FILTERFASTQ( fastqFile )
	}
	else{
		filteredFastqFiles = fastqFile
	}
	RUNFASTQC(filteredFastqFiles)
	if (params.alignment == "Yes" ){
		trimmedFiltBamFile = RUNALIGNMENT( filteredFastqFiles, bwaRefIdx , reference)
		if (params.pseudoDiploGeno == "Yes" ){
		    RUNBAMTOPSEDIPGENO(trimmedFiltBamFile, faIdx, faDict, reference)
		}
	}
}



workflow.onComplete { 
    println ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
    }
