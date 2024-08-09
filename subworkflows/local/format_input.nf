/*
    Subworkflow to format input fastq files/folders from either directories or samplesheet
        Want to concat samples to one fastq file when using either the csv or directory
            For the CSV, we can do a barcode dir and a fastq file and then combine
                Somehow?

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { samplesheetToList } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INITIALIZE CHANNELS FROM PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow FORMAT_INPUT {
    //
    // Fastq pass directory input
    //
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
        } else if ( nanoporeFastqs ) {
            Channel.fromPath( nanoporeFastqs )
                .map{ fastq -> [ [id: fastq.baseName.replaceAll(~/\.fastq.*$/, '')], file(fastq) ] }
                .branch{
                    pass: it[1].countFastq() >= 1
                    empty: it[1].countFastq() == 0
                }.set{ ch_fastqs }
        // Failing conditions
        } else {
            log.error("ERROR: Couldn't detect any fastq files in --fastq_pass ${params.fastq_pass}")
            System.exit(1)
        }
    }

    //
    // Input CSV file with path to fastq barcode directories or single reads
    //
    else {
        // This works for directories
        Channel.fromPath( params.input, type: 'file', checkIfExists: true )
            .splitCsv(header:true)
            .map{ row -> [ [id: row.sample], file(row.reads, type: 'dir', checkIfExists: true) ] }
            .branch{
                    pass: it[1].listFiles().size() >= 1
                    empty: it[1].listFiles().size() == 0
            }.set{ ch_fastqs }
    }

    emit:
    fastq_pass = ch_fastqs.pass     // 
    fastq_empty = ch_fastqs.empty   // 
}
