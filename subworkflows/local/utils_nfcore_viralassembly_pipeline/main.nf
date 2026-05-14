//
// Subworkflow with functionality specific to the phac-nml/viralassembly pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args // array: List of positional nextflow CLI args
    outdir            // string: The output directory where the results will be saved

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Check logic required for the pipeline to function properly
    //  Stuff like files for different inputs, models, etc
    //
    // Variant Callers - Clair3 is default but we only allow these 3 (and probably will remove them later for just C3)
    if ( ! ['medaka', 'nanopolish', 'clair3'].contains(params.variant_caller) ) {
        log.error("Please provide an input for --variant_caller with any of [ 'clair3', 'nanopolish', 'medaka' ]")
        System.exit(1)
    }

    //-- Data Inputs
    if ( !params.input && !params.fastq_pass ) {
        log.error("Please provide input data with either: '--input input.csv' or '--fastq_pass fastq_dir/'")
        System.exit(1)
    } else if ( params.input && params.fastq_pass ) {
        log.error("Please provide input data with either: '--input input.csv' or '--fastq_pass fastq_dir/' but not both")
        System.exit(1)
    } else if ( params.variant_caller == 'nanopolish' ) {
        if ( ! params.fast5_pass || ! params.sequencing_summary ) {
            log.error("Please pass both '--fast5_pass fast5_dir/' and '--sequencing_summary seqsum.txt' to run nanopolish")
            System.exit(1)
        }
    }

    // Multiple Nextclade Inputs
    if (! params.skip_nextclade && (
        (params.nextclade_dataset_dir && params.nextclade_dataset_name))) {
            log.error("Please pass only one of the following to run nextclade: `--nextclade_dataset_name', '--nextclade_dataset_dir', or '--virus_name'")
            System.exit(1)
    }

    // Skip nextlade while providing a nextclade dataset name or directory
    if ( params.skip_nextclade && params.nextclade_dataset_name||params.nextclade_dataset_dir ) {
        log.error("Remove `--skip_nextclade' to run nextclade or remove '--nextclade_dataset_name' to skip netxclade.")
        System.exit(1)
    } else if ( params.skip_nextclade && params.nextclade_dataset_dir ) {
        log.error("Remove `--skip_nextclade' to run nextclade or remove '--nextclade_dataset_dir' to skip netxclade.")
        System.exit(1)
    }

    //
    // Summarize and Validate Params
    //
    if (validate_params) {
        validateParameters()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    outdir          // path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output

    main:

    //
    // Completion email and summary
    //
    workflow.onComplete {

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}
