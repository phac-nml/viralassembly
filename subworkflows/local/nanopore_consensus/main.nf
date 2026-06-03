/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Read QC
include { CHOPPER                   } from '../../../modules/local/chopper/main'
include { NANOSTAT                  } from '../../../modules/local/nanostat/main'
include { RENAME_FASTQ              } from '../../../modules/local/custom/utils.nf'

// Clair3 model
include { GET_MODEL                 } from '../../../modules/local/get_model/main'

// Artic related
include { ARTIC_GUPPYPLEX           } from '../../../modules/local/artic/guppyplex/main'
include { ARTIC_MINION              } from '../../../modules/local/artic/minion/main'

// Subworkflows
include { WF_NANOPORE_AMPLICON      } from '../../../subworkflows/local/nanopore_amplicon'
include { WF_NANOPORE_SHOTGUN       } from '../../../subworkflows/local/nanopore_shotgun'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow WF_NANOPORE_CONSENSUS {
    take:
    ch_fastqs       // channel: [ val(meta), file(fastqs) ]
    ch_reference    // channel: [ val(meta), file(reference) ]
    ch_fai          // channel: [ file(fai) ]
    ch_refstats     // channel: [ file(refstats) ]
    ch_primer_bed   // channel: [ file(primer.bed) ]
    ch_amplicon_bed // channel: [ file(amplicon.bed) ]

    main:

    // Nanopolish required channels, will be ignored when running clair3 or medaka but still passed to the workflow
    ch_fast5s = params.fast5_pass ? file(params.fast5_pass, type: 'dir', checkIfExists: true) : []
    ch_seqsum = params.sequencing_summary ? file(params.sequencing_summary, type: 'file', checkIfExists: true) : []

    // Tool version tracking
    ch_versions = channel.empty()

    // Clair3 model
    ch_clair3_model = channel.empty()
    if ( params.clair3_local_model ) {
        ch_clair3_model = file(params.clair3_local_model, checkIfExists: true)
    } else {
        GET_MODEL(params.clair3_model)
        ch_clair3_model = GET_MODEL.out.model
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Read QC and Statistics
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    ARTIC_GUPPYPLEX(
        ch_fastqs
    )
    ch_fastqs = ARTIC_GUPPYPLEX.out.fastq
    ch_versions = ch_versions.mix(ARTIC_GUPPYPLEX.out.versions)


    // Chopper may be useless as we already filter based on length earlier
    //  But it also does add quality filtering
    CHOPPER(
        ch_fastqs
    )
    ch_versions = ch_versions.mix(CHOPPER.out.versions)
    ch_fastqs = CHOPPER.out.fastq

    // Pass/fail reads based on count after length and quality filtering
    ch_fastqs
        .branch{ _meta, fastq ->
            pass: fastq.countFastq() >= params.min_reads
            empty: fastq.countFastq() < params.min_reads
        }.set{ ch_filtered_fastqs }

    // Stats on final reads
    NANOSTAT(
        ch_filtered_fastqs.pass
    )
    ch_versions = ch_versions.mix(NANOSTAT.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Chose which pipeline to run based on input params
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    if ( !params.primer_bed ) {
        WF_NANOPORE_SHOTGUN(
            ch_filtered_fastqs.pass,
            ch_fast5s,
            ch_seqsum,
            ch_reference,
            ch_fai,
            ch_refstats,
            ch_clair3_model
        )
        ch_consensus = WF_NANOPORE_SHOTGUN.out.consensus
        ch_bam = WF_NANOPORE_SHOTGUN.out.bam
        ch_vcf = WF_NANOPORE_SHOTGUN.out.vcf
        ch_versions = ch_versions.mix(WF_NANOPORE_SHOTGUN.out.versions)
    } else if ( !params.use_artic_tool ) {
        WF_NANOPORE_AMPLICON(
            ch_filtered_fastqs.pass,
            ch_fast5s,
            ch_seqsum,
            ch_reference,
            ch_fai,
            ch_refstats,
            ch_primer_bed,
            ch_amplicon_bed,
            ch_clair3_model
        )
        ch_consensus = WF_NANOPORE_AMPLICON.out.consensus
        ch_bam = WF_NANOPORE_AMPLICON.out.bam
        ch_vcf = WF_NANOPORE_AMPLICON.out.vcf
        ch_versions = ch_versions.mix(WF_NANOPORE_AMPLICON.out.versions)
    } else {
        ARTIC_MINION(
            ch_filtered_fastqs.pass,
            ch_reference.collect{ _meta, ref -> ref },
            ch_primer_bed,
            ch_clair3_model.ifEmpty([])
        )
        ch_consensus = ARTIC_MINION.out.consensus
        ch_bam = ARTIC_MINION.out.bam
        ch_vcf = ARTIC_MINION.out.vcf
        ch_versions = ch_versions.mix(ARTIC_MINION.out.versions)
    }

    emit:
    consensus             = ch_consensus                // channel: [ val(meta), file(consensus) ]
    bam                   = ch_bam                      // channel: [ val(meta), file(bam), file(bai) ]
    vcf                   = ch_vcf                      // channel: [ val(meta), file(vcf) ]
    empty_filtered_fastqs = ch_filtered_fastqs.empty    // channel: [ val(meta), file(fastqs) ]
    stats                 = NANOSTAT.out.stats          // channel: [ val(meta), file(nanostat.txt) ]
    versions              = ch_versions                 // channel: [ path(versions.yml) ]
}
