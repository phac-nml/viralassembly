/*
    Subworkflow to create the custom report visualizations
        Included Options:
            Custom HTML Interactive
            Custom PDF Static        - *inprogress eventually*

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Bedtools
include { BEDTOOLS_COVERAGE_GENOME_BED      } from '../../modules/local/bedtools/coverage/main'
include { BEDTOOLS_COVERAGE_AMPLICON_BED    } from '../../modules/local/bedtools/coverage/main'

// Visualization
include { CREATE_READ_VARIATION_CSV         } from '../../modules/local/visualization/main'
include { CREATE_VARIANT_TSV                } from '../../modules/local/visualization/main'
include { COMBINE_AMPLICON_COVERAGE         } from '../../modules/local/visualization/main'
include { CREATE_AMPLICON_COMPLETENESS      } from '../../modules/local/visualization/main'
include { CREATE_ALL_SAMPLE_SUMMARY_REPORT  } from '../../modules/local/visualization/main'

// Software Version Dump
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow WF_CREATE_CUSTOM_REPORT {
    take:
    ch_consensus        // channel: [ val(meta), path(consensus) ]
    ch_bam              // channel: [ val(meta), path(bam), path(bai) ]
    ch_vcf              // channel: [ val(meta), path(vcf) ]
    ch_reference        // channel: [ path(reference) ]
    ch_genome_bed       // channel: [ path(genome.bed) ]
    ch_amplicon_bed     // channel: [ path(amplicon_bed) ]
    ch_qc_csv           // channel: [ path(qc_csv) ]
    ch_versions         // channel: [ path(versions) ]

    main:
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Variant analysis
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    CREATE_READ_VARIATION_CSV(
        ch_bam,
        ch_reference
    )
    ch_versions = ch_versions.mix(CREATE_READ_VARIATION_CSV.out.versions)

    CREATE_VARIANT_TSV(
        ch_vcf
    )
    ch_versions = ch_versions.mix(CREATE_VARIANT_TSV.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Genome coverage
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    BEDTOOLS_COVERAGE_GENOME_BED(
        ch_bam,
        ch_genome_bed
    )
    ch_versions = ch_versions.mix(BEDTOOLS_COVERAGE_GENOME_BED.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Amplicon coverages
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    ch_amplicon_coverage = Channel.empty()
    ch_amplicon_completeness = Channel.empty()
    if ( ! params.reference ) {
        // Coverage
        BEDTOOLS_COVERAGE_AMPLICON_BED(
            ch_bam,
            ch_amplicon_bed
        )
        ch_versions = ch_versions.mix(BEDTOOLS_COVERAGE_AMPLICON_BED.out.versions)

        //  Custom Output - should only run with that later
        COMBINE_AMPLICON_COVERAGE(
            BEDTOOLS_COVERAGE_AMPLICON_BED.out.amplicon_coverage
                .collect{ it[1] }
        )
        ch_amplicon_coverage = COMBINE_AMPLICON_COVERAGE.out.amplicons

        // Completeness
        CREATE_AMPLICON_COMPLETENESS(
            ch_consensus,
            ch_amplicon_bed
        )
        CREATE_AMPLICON_COMPLETENESS.out.amplicon_completeness
            .collectFile(keepHeader: true, sort: { it.baseName }, skip: 1, name: 'merged_amplicon_completeness.csv')
            .set { ch_amplicon_completeness }

        ch_versions = ch_versions.mix(CREATE_AMPLICON_COMPLETENESS.out.versions)
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Tracking versions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions
            .unique()
            .collectFile(name: 'collated_versions.yml')
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Final Report
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    CREATE_ALL_SAMPLE_SUMMARY_REPORT(
        CREATE_READ_VARIATION_CSV.out.csv
            .collect{ it[1] },
        CREATE_VARIANT_TSV.out.tsv
            .collect{ it[1] },
        BEDTOOLS_COVERAGE_GENOME_BED.out.cov_bed
            .collect{ it[1] },
        ch_amplicon_coverage
            .ifEmpty([]),
        ch_amplicon_completeness
            .ifEmpty([]),
        ch_qc_csv,
        CUSTOM_DUMPSOFTWAREVERSIONS.out.yml
    )
}
