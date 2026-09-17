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
    //-- Data Inputs
    if ( !params.input && !params.fastq_pass ) {
        log.error("Please provide input data with either: '--input input.csv' or '--fastq_pass fastq_dir/'")
        System.exit(1)
    } else if ( params.input && params.fastq_pass ) {
        log.error("Please provide input data with either: '--input input.csv' or '--fastq_pass fastq_dir/' but not both")
        System.exit(1)
    }

    // Multiple Nextclade Inputs
    if (! params.skip_nextclade && (
        (params.nextclade_dataset_dir && params.nextclade_dataset_name))) {
            log.error("Please pass only one of the following to run nextclade: `--nextclade_dataset_name' or '--nextclade_dataset_dir'")
            System.exit(1)
    }

    //
    // Summarize and Validate Params
    //
    if (validate_params) {
        validateParameters()
    }

    // Check if reference is segmented
    fasta = file(params.reference, type: 'file', checkIfExists: true)
    segmented = isSegmented(fasta)

    // Nextclade input when virus is segmented
    if ( segmented && (params.nextclade_dataset_dir || params.nextclade_dataset_name) ) {
        log.error("The reference FASTA used is a segmented virus. Please remove the 'nextclade_dataset_dir' or 'nextclade_dataset_name' argument as the pipeline will assign the appropriate nextclade dataset to each segment.")
        System.exit(1)
    }

    emit:
    segmented = segmented       // boolean: If virus is segmented

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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check if the virus is segmented
def isSegmented(reference) {
    def input = reference.name.endsWith('.gz')
        ? new java.util.zip.GZIPInputStream(reference.newInputStream())
        : reference.newInputStream()

    def segmented

    input.withReader { reader ->
        segmented = reader
            .readLines()
            .findAll { it.startsWith('>') }
            .size() > 1
    }

    return segmented
}

// Get FASTA header
def fastaHeaderId(fasta) {
    def headers = fasta.readLines()
        .findAll { it.startsWith('>') }
        .collect { it.substring(1).tokenize()[0] }

    return headers
}
