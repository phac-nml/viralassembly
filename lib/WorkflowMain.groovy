//
// This file holds several functions specific to the main.nf workflow in the f/test pipeline
//

import nextflow.Nextflow

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            // TODO nf-core: Add Zenodo DOI for pipeline after first release
            //"* The pipeline\n" +
            //"  https://doi.org/10.5281/zenodo.XXXXXXX\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    //
    // Generate help string
    //
    public static String help(workflow, params, log) {
        def command = "nextflow run ${workflow.manifest.name} -profile mamba --variant_caller 'medaka' [ --input samplesheet.csv || --fastq_pass <FASTQ_DIR/> ]"
        def help_string = ''
        help_string += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        help_string += NfcoreSchema.paramsHelp(workflow, params, command)
        //help_string += '\n' + citation(workflow) + '\n'
        help_string += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return help_string
    }

    //
    // Generate parameter summary log string
    //
    public static String paramsSummaryLog(workflow, params, log) {
        def summary_log = ''
        summary_log += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        summary_log += NfcoreSchema.paramsSummaryLog(workflow, params)
        // summary_log += '\n' + citation(workflow) + '\n'
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return summary_log
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {

        // Print workflow version and exit on --version
        if (params.version) {
            String workflow_version = NfcoreTemplate.version(workflow)
            log.info "${workflow.manifest.name} ${workflow_version}"
            System.exit(0)
        }

        // Print workflow help statement and exit on --help
        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

        // Print parameter summary log before starting anything
        log.info paramsSummaryLog(workflow, params, log)

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)

        // Check that conda channels are set-up correctly
        if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Can't run singularity with the custom report at the moment
        if ( workflow.containerEngine == 'singularity' && params.custom_report ) {
            log.error("Custom report is not currently compatible with singularity")
            System.exit(1)
        }

        // Input logic checks
        //-- Variant caller check
        if ( ! ['medaka', 'nanopolish', 'clair3'].contains(params.variant_caller) ) {
            log.error("Please provide an input for --variant_caller with any of ['nanopolish', 'medaka', 'clair3']")
            System.exit(1)
        }
        //-- Data inputs
        if ( ! params.input && ! params.fastq_pass ) {
            log.error("Please provide input data with either: '--input samplesheet.csv' or '--fastq_pass fastq_dir/'")
            System.exit(1)
        } else if ( params.input && params.fastq_pass ) {
            log.error("Please provide input data with either: '--input samplesheet.csv' or '--fastq_pass fastq_dir/'")
            System.exit(1)
        } else if ( params.variant_caller == 'nanopolish' ) {
            if ( ! params.fast5_pass || ! params.sequencing_summary ) {
                log.error("Please pass both '--fast5_pass fast5_dir/' and '--sequencing_summary' to run nanopolish")
                System.exit(1)
            }
        }
    }

    //
    // Get attribute from genome config file e.g. fasta
    //
    public static Object getGenomeAttribute(params, attribute) {
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                return params.genomes[ params.genome ][ attribute ]
            }
        }
        return null
    }
}
