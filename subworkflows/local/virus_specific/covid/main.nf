/*
    Subworkflow to run covid specific tools
        1. Identifies covid lineages

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PANGOLIN_UPDATEDATA } from '../../../../modules/local/pangolin/updatedata/main'
include { PANGOLIN_RUN        } from '../../../../modules/local/pangolin/run/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow WF_VIRUS_COVID {
    take:
    ch_consensus        // channel: [ val(meta), path(consensus), path(bam) ]

    main:
    ch_versions = Channel.empty()

    //
    // Lineage analysis with Pangolin
    //
    ch_pangolin_report = Channel.empty()

    if (!params.skip_pangolin) {
        if (!params.pango_database) {
            PANGOLIN_UPDATEDATA('pangolin_db')
            pango_database = PANGOLIN_UPDATEDATA.out.db
            ch_versions   = ch_versions.mix(PANGOLIN_UPDATEDATA.out.versions.first())
        } else {
            pango_database = Channel.value(file(params.pango_database, type: 'dir'))
        }

        PANGOLIN_RUN (
            ch_consensus,
            pango_database
        )
        ch_pangolin_report = PANGOLIN_RUN.out.report
        ch_versions        = ch_versions.mix(PANGOLIN_RUN.out.versions.first())
    }

    emit:
    pangolin_report  = ch_pangolin_report            // channel: [ val(meta), [ csv ] ]
    versions         = ch_versions                   // channel: [ path(versions.yml) ]
}