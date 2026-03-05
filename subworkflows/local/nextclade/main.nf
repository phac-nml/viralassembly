/*
    Subworkflow to run nextclade on consensus files
        1. Gets the nextstrain or custom dataset to use
        2. Analyzes sequences to identify mutations and assess sequence quality

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { NEXTCLADE_DATASETGET  } from '../../../modules/local/nextclade/datasetget/main'
include { NEXTCLADE_RUN         } from '../../../modules/local/nextclade/run/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow WF_NEXTCLADE {
    take:
    ch_consensus  // channel: [ val(meta), file(consensus) ]

    main:
    // Version tracking
    ch_versions = channel.empty()

    if (! params.nextclade_dataset_dir ) {
        ch_dataset_name = params.nextclade_dataset_name ? params.nextclade_dataset_name : ''
        ch_dataset_tag = params.nextclade_dataset_tag ? params.nextclade_dataset_tag : ''
        ch_virus_name = params.virus_name ? params.virus_name : ''

        NEXTCLADE_DATASETGET(
            ch_dataset_name,
            ch_dataset_tag,
            ch_virus_name
        )
        ch_versions = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)
    }
    ch_nextclade_dataset = params.nextclade_dataset_dir ? 
        Channel.value(file(params.nextclade_dataset_dir, type: 'dir', checkIfExists: true)) : NEXTCLADE_DATASETGET.out.dataset

    NEXTCLADE_RUN(
        ch_consensus,
        ch_nextclade_dataset
    )
    ch_versions = ch_versions.mix(NEXTCLADE_RUN.out.versions)

    emit:
    csv = NEXTCLADE_RUN.out.csv
    versions = ch_versions
}
