/*
    Subworkflow to run nextclade on consensus files
        1. Gets the nextstrain or custom dataset to use
        2. Analyzes sequences to identify mutations and assess sequence quality

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { NEXTCLADE_SORT        } from '../../../modules/local/nextclade/sort/main'
include { NEXTCLADE_DATASETGET  } from '../../../modules/local/nextclade/datasetget/main'
include { NEXTCLADE_RUN         } from '../../../modules/local/nextclade/run/main'
include { COLLATE_CSVS          } from '../../../modules/local/custom/collate_csvs'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow WF_NEXTCLADE {
    take:
    ch_consensus  // channel: [ val(meta), file(consensus) ]
    segmented     // boolean: Boolean whether virus is segmented or not

    main:
    // Version tracking
    ch_versions = channel.empty()

    // Split multi-FASTA (segmented viruses)
    ch_consensus = ch_consensus.flatMap { meta, fasta ->
        if ( segmented ) {
            return fasta
                .splitFasta(by: 1, file: true)
                .collect { split_fasta ->
                    tuple(meta + [file_name: fasta.name], split_fasta)
                }
        } else {
            return [ tuple(meta + [file_name: fasta.name], fasta) ]
        }
    }

    // Figure out the best dataset to use with Nextclade Sort or use a specified dataset
    if ( segmented || !(params.nextclade_dataset_name || params.nextclade_dataset_dir) ) {
        // MODULE: Define the most appropriate dataset for a sample or segment
        NEXTCLADE_SORT(
            ch_consensus
        )
        ch_dataset_tag = ''
        ch_dataset_name = NEXTCLADE_SORT.out.dataset_name
            .map{ _meta, _fasta, stdout ->
                stdout.toString().trim()
            }
            .filter { it }
            .unique()
    } else if ( !(segmented) && !(params.nextclade_dataset_dir) ) {
        ch_dataset_name = params.nextclade_dataset_name ? params.nextclade_dataset_name : ''
        ch_dataset_tag = params.nextclade_dataset_tag ? params.nextclade_dataset_tag : ''
    }

    // MODULE: Download the specified dataset
    NEXTCLADE_DATASETGET(
        ch_dataset_name,
        ch_dataset_tag
    )
    ch_versions = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)
    ch_nextclade_dataset = params.nextclade_dataset_dir ?
        Channel.value(file(params.nextclade_dataset_dir, type: 'dir', checkIfExists: true)) : NEXTCLADE_DATASETGET.out.dataset

    ch_nextclade_dataset.ifEmpty {
        log.warn("There were no matching nextclade datasets for this virus. Skipping nextclade.")
    }
    // Define input for Nextclade run
    if ( segmented ) {
        ch_nextclade_run_input = NEXTCLADE_SORT.out.dataset_name
            .map{meta, fasta, stdout ->
                tuple(stdout.trim(), meta, fasta)
            }
            .combine(ch_nextclade_dataset, by: 0)
            .map{ _dataset_name, meta, fasta, dataset_path ->
                tuple(meta, fasta, dataset_path)
            }
    } else {
        ch_nextclade_run_input = ch_consensus
            .combine(ch_nextclade_dataset.map{ _dataset_name, dataset_path -> dataset_path })
    }

    // MODULE: Run nextclade to determine QC issues in consensus sequences
    NEXTCLADE_RUN(
        ch_nextclade_run_input
    )
    ch_versions = ch_versions.mix(NEXTCLADE_RUN.out.versions)

    // Collate all results in a single output if virus is segmented
    if ( segmented ) {
        ch_csvs = NEXTCLADE_RUN.out.csv
            .groupTuple()

        // MODULE: Collate each segments output into a final csv file
        COLLATE_CSVS(
            ch_csvs
        )
        ch_nextclade_csv = COLLATE_CSVS.out.final_csv
    } else {
        ch_nextclade_csv = NEXTCLADE_RUN.out.csv
    }

    emit:
    csv      = ch_nextclade_csv // channel: [ val(meta), file(csv) ]
    versions = ch_versions      // channel: [ path(versions.yml) ]
}
