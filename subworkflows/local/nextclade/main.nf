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

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Segmented Virus Handling
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Split multi-FASTA (segmented viruses)
    //  and add in nextclade_file_name to meta for the final nextcalde output name
    ch_consensus = ch_consensus.flatMap { meta, fasta ->
        if ( segmented ) {
            return fasta
                .splitFasta(by: 1, file: true)
                .collect { split_fasta ->
                    tuple(meta + [nextclade_file_name: fasta.name.replace('.consensus.fasta','')], split_fasta)
                }
        } else {
            return [ tuple(meta + [nextclade_file_name: fasta.name.replace('.consensus.fasta','')], fasta) ]
        }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Dataset Specification & Download
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Function to run Nextclade Sort
    def run_nextclade_sort = !params.nextclade_dataset_dir &&
        (segmented || !params.nextclade_dataset_name)

    // Create / Download / Check nextclade dataset directory(s)
    ch_dataset_name = params.nextclade_dataset_name ? params.nextclade_dataset_name : ''
    ch_dataset_tag = params.nextclade_dataset_tag ? params.nextclade_dataset_tag : ''
    if ( params.nextclade_dataset_dir ) {
        ch_nextclade_dataset = Channel.value(
            tuple('', file(params.nextclade_dataset_dir, type: 'dir', checkIfExists: true))
        )
    } else {
        // Figure out the best dataset to use with Nextclade Sort
        if ( run_nextclade_sort ) {
            NEXTCLADE_SORT(
                ch_consensus
            )

            // WARN that a sample didn't have hits
            //  therefore there will be no nextclade output
            NEXTCLADE_SORT.out.dataset_name.map { meta, fasta, dataset ->
                if (dataset.isEmpty() ) {
                    log.warn("${meta.id} has no matching nextclade datasets. No nextclade output will be available for this sample.")
                }
                tuple(meta, fasta, dataset)
            }
            ch_dataset_tag = ''
            ch_dataset_name = NEXTCLADE_SORT.out.dataset_name
                .map{ _meta, _fasta, dataset ->
                    dataset.toString().trim()
                }
                .filter { it }
                .unique()
        }
        // Actually download the dataset
        NEXTCLADE_DATASETGET(
            ch_dataset_name,
            ch_dataset_tag
        )
        ch_versions = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)
        ch_nextclade_dataset = NEXTCLADE_DATASETGET.out.dataset

    }

    // Mention if there were no datasets
    //  May need an easier way to see this for troubleshooting later but ok for now
    ch_nextclade_dataset.ifEmpty {
        log.warn("There were no matching nextclade datasets for this virus. Skipping nextclade.")
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Run Nextclade
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Define input for Nextclade run
    if ( run_nextclade_sort ) {
        ch_nextclade_run_input = NEXTCLADE_SORT.out.dataset_name
            .map{meta, fasta, dataset ->
                tuple(dataset.trim(), meta, fasta)
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

    // Collate all results in a single output for both segmented and non-segmented viruses
    COLLATE_CSVS(
        NEXTCLADE_RUN.out.csv
            .groupTuple()
    )

    emit:
    csv      = COLLATE_CSVS.out.final_csv // channel: [ val(meta), file(csv) ]
    versions = ch_versions                // channel: [ path(versions.yml) ]
}
