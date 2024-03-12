#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    phac-nml/viralassembly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/phac-nml/viralassembly
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NANOPORE } from './workflows/nanopore.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Inital workflow and param checks
//
WorkflowMain.initialise(workflow, params, log)

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow VIRALASSEMBLY {
    // Fastq pass directory input
    if ( params.fastq_pass ) {
        // Create fastqs channel based on if barcode dirs or .fastq files found in input
        nanoporeBarcodeDirs = file("${params.fastq_pass}/barcode*", type: 'dir', maxdepth: 1 )
        nanoporeFastqs = file("${params.fastq_pass}/*.fastq*", type: 'file', maxdepth: 1)
        if ( nanoporeBarcodeDirs ) {
            Channel.fromPath( nanoporeBarcodeDirs )
                .filter( ~/.*barcode[0-9]{1,4}$/ )
                .map{ dir -> [ [id: dir.baseName], file(dir) ] }
                .branch{
                    pass: it[1].listFiles().size() >= 1
                    empty: it[1].listFiles().size() == 0
                }.set{ ch_fastqs }
        } else if ( (nanoporeFastqs && params.variant_caller == 'medaka') || (nanoporeFastqs && params.variant_caller == 'clair3') ) {
            Channel.fromPath( nanoporeFastqs )
                .map{ fastq -> [ [id: fastq.baseName.replaceAll(~/\.fastq.*$/, '')], file(fastq) ] }
                .branch{
                    pass: it[1].countFastq() >= 1
                    empty: it[1].countFastq() == 0
                }.set{ ch_fastqs }
        // Failing conditions
        } else if ( nanoporeFastqs && params.variant_caller == 'nanopolish' ) {
            log.error("ERROR: Found fastq files but nanopolish requires barcoded directories as input")
            System.exit(1)
        } else {
            log.error("ERROR: Couldn't detect any fastq files in --fastq_pass ${params.fastq_pass}")
            System.exit(1)
        }
    }
    // Input CSV file with path to fastq folders or singular reads
    // -- May need to adjust logic for better error tracking later
    else {
        Channel.fromPath( params.input, type: 'file', checkIfExists: true )
            .splitCsv(header:true)
            .map{ row -> [ [id: row.sample], file(row.reads, type: 'dir', checkIfExists: true) ] }
            .branch{ 
                    pass: it[1].listFiles().size() >= 1
                    empty: it[1].listFiles().size() == 0
                }.set{ ch_fastqs }
    }

    main:

    //
    // WORKFLOW: Run pipeline
    //
    NANOPORE (
        ch_fastqs.pass,
        ch_fastqs.empty
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // WORKFLOW: Run main workflow
    //
    VIRALASSEMBLY()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
