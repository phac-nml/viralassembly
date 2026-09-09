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
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow FORMAT_INPUT {
    main:
    //
    // Fastq pass directory input
    //
    if ( params.fastq_pass && params.platform == 'nanopore' ) {
        // Create fastqs channel based on if barcode dirs or .fastq/.fq files found in input
        nanoporeBarcodeDirs = file("${params.fastq_pass}/barcode*", type: 'dir', maxdepth: 1 )
        nanoporeFastqs = file("${params.fastq_pass}/*.{fastq,fq}{,.gz}", type: 'file', maxdepth: 1)
        // Barcode DIRs
        if ( nanoporeBarcodeDirs ) {
            channel.fromPath( nanoporeBarcodeDirs )
                .filter( ~/.*barcode[0-9]{1,4}$/ )
                .map{ dir -> [ [id: dir.baseName], file(dir) ] }
                .branch{ _meta, dir ->
                    pass: dir.listFiles().size() >= 1
                    empty: dir.listFiles().size() == 0
                }.set{ ch_fastqs }
        // FASTQS
        } else if ( nanoporeFastqs ) {
            channel.fromPath( nanoporeFastqs )
                .map{ fastq -> [ [id: fastq.baseName.replaceAll(~/\.(fastq|fq)(\.gz)?$/, '')], file(fastq) ] }
                .branch{ _meta, fastq ->
                    pass: fastq.countFastq() >= 1
                    empty: fastq.countFastq() == 0
                }.set{ ch_fastqs }
        // Failing to detect
        } else {
            log.error("Couldn't detect any barcode directories or fastq files in --fastq_pass ${params.fastq_pass}")
            System.exit(1)
        }
    } else if ( params.fastq_pass && params.platform == 'illumina' ) {
        // Create fastqs channel based on if .fastq/.fq files found in input
        illuminaFastqs = file("${params.fastq_pass}/*.{fastq,fq}{,.gz}", type: 'file', maxdepth: 1)
        if (illuminaFastqs){
            channel.fromFilePairs("${params.fastq_pass}/*_{R1,R2}*.{fastq,fq}{,.gz}", size: 2)
                .map{ id, fastqs -> [ [id: id], fastqs ] }
                .branch{ _meta, fastqs ->
                    pass: fastqs[0].countFastq() >= 1 && fastqs[1].countFastq() >= 1
                    empty: fastqs[0].countFastq() == 0 && fastqs[1].countFastq() == 0
                }.set{ ch_fastqs }
        // Failing to detect
        } else {
            log.error("Couldn't detect any fastq files in --fastq_pass ${params.fastq_pass}")
            System.exit(1)
        }
    }

    //
    // Input CSV file with path to fastq single reads to comply with IRIDA Next (eventually)
    //
    else {
        // Using the samplesheet, were adding fastq_2 just as a thingy for now
        def processedIDs = [] as Set
        Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))
            .map { meta, fastq_1, fastq_2 ->
                // Meta ID assignment
                if (!meta.id) {
                    meta.id = meta.irida_id
                } else {
                    meta.id = meta.id.replaceAll(/[^A-Za-z0-9_.\-]/, '_')
                }

                // Ensure ID is unique by appending meta.irida_id if needed
                //  Note that nextflow does not like while loops
                while (processedIDs.contains(meta.id)) {
                    meta.id = "${meta.id}_${meta.irida_id}"
                }
                // Add the ID to the set of processed IDs
                processedIDs << meta.id

                // File assignment
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
            }
            .groupTuple()
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map { meta, fastqs ->
                return [ meta, fastqs.flatten() ]
            }.set { ch_tmp_fastqs }

        ch_tmp_fastqs
            .map { meta, fastqs ->
                // Temp to match above filtering - 2025-12-12
                def count = 0
                for (fastq in fastqs) {
                    count += fastq.countFastq()
                }
                return [ meta, fastqs, count ]
            }
            .branch { meta, fastq, count ->
                empty: count == 0
                    return [ meta, fastq ]
                pass: count >= 1
                    return [ meta, fastq ]
            }
            .set{ ch_fastqs }
    }

    emit:
    pass  = ch_fastqs.pass    // channel: [ val(meta), file(fastq) ]
    empty = ch_fastqs.empty   // channel: [ val(meta), file(fastq) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
