/*
    Subworkflow to create the custom report visualizations
        Included:
            MultiQC Overall HTML Interactive Report
            MultiQC Sample  HTML Interactive Report

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_STATS                 } from '../../modules/local/bcftools/stats/main'
include { BEDTOOLS_COVERAGE_AMPLICON_BED } from '../../modules/local/bedtools/coverage/main'
include { QUALIMAP_BAMQC                 } from '../../modules/local/qualimap/bamqc/main'
include { SAMTOOLS_FLAGSTAT              } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_REHEADER              } from '../../modules/local/samtools/reheader/main'

// Visualization
include { CREATE_READ_VARIATION_CSV      } from '../../modules/local/visualization/main'
include { CREATE_VARIANT_TSV             } from '../../modules/local/visualization/main'
include { COMBINE_AMPLICON_COVERAGE      } from '../../modules/local/visualization/main'
include { CSVTK_SAMPLE_AMPLICON_DEPTH    } from '../../modules/local/visualization/main'
include { CREATE_AMPLICON_COMPLETENESS   } from '../../modules/local/visualization/main'

// Software Version Dump
include { CUSTOM_DUMPSOFTWAREVERSIONS    } from '../../modules/nf-core/custom/dumpsoftwareversions/main'

// MultiQC
include { MULTIQC_SAMPLE    } from '../../modules/local/multiqc/main'
include { MULTIQC_OVERALL   } from '../../modules/local/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INITIALIZE CHANNELS FROM PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
ch_multiqc_overall_conf = file(params.multiqc_config_overall, checkIfExists: true)
ch_multiqc_sample_conf = file(params.multiqc_config_sample, checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow WF_CREATE_MULTIQC_REPORTS {
    take:
    ch_consensus        // channel: [ val(meta), path(consensus) ]
    ch_bam              // channel: [ val(meta), path(bam), path(bai) ]
    ch_vcf              // channel: [ val(meta), path(vcf) ]
    ch_sample_csv       // channel: [ val(meta), path(csv) ]
    ch_snpeff_csv       // channel: [ val(meta), path(csv) ] || empty
    ch_reference        // channel: [ path(reference) ]
    ch_amplicon_bed     // channel: [ path(amplicon_bed) ] || empty
    ch_overall_qc_csv   // channel: [ path(qc_csv) ]
    ch_versions         // channel: [ path(versions) ]

    main:
    // Sample variant analysis
    CREATE_READ_VARIATION_CSV(
        ch_bam,
        ch_reference
    )
    ch_versions = ch_versions.mix(CREATE_READ_VARIATION_CSV.out.versions)

    CREATE_VARIANT_TSV(
        ch_vcf
    )
    ch_versions = ch_versions.mix(CREATE_VARIANT_TSV.out.versions)

    // Amplicon analysis
    ch_amplicon_completeness = Channel.empty()
    if ( ! params.reference_no_scheme ) {
        // Coverage
        BEDTOOLS_COVERAGE_AMPLICON_BED(
            ch_bam,
            ch_amplicon_bed
        )
        ch_versions = ch_versions.mix(BEDTOOLS_COVERAGE_AMPLICON_BED.out.versions)

        // Per-sample amplicon depth
        CSVTK_SAMPLE_AMPLICON_DEPTH(
            BEDTOOLS_COVERAGE_AMPLICON_BED.out.amplicon_coverage
        )
        ch_sample_amplicon_depth = CSVTK_SAMPLE_AMPLICON_DEPTH.out.tsv

        ch_versions = ch_versions.mix(CSVTK_SAMPLE_AMPLICON_DEPTH.out.versions)

        // Completeness
        CREATE_AMPLICON_COMPLETENESS(
            ch_consensus,
            ch_amplicon_bed
        )
        CREATE_AMPLICON_COMPLETENESS.out.amplicon_completeness
            .collectFile(keepHeader: true, sort: { it.baseName }, skip: 1, name: 'merged_amplicon_completeness.csv')
            .set { ch_amplicon_completeness }

        ch_versions = ch_versions.mix(CREATE_AMPLICON_COMPLETENESS.out.versions)
    } else {
        // Create empty amplicon depth tuple with sample meta values to still get sample mqc reports
        ch_sample_amplicon_depth = ch_bam.map{ it -> tuple(it[0], []) }
    }

    // Stats from tools
    BCFTOOLS_STATS(
        ch_vcf,
        ch_reference
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)

    SAMTOOLS_FLAGSTAT(
        ch_bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    // If not using a scheme, need to correct the headers for qualimap by removing the empty RG
    // BAM channel also no longer needs bai file
    ch_bam = ch_bam.map { it -> tuple(it[0], it[1])}
    if ( ! params.reference_no_scheme ) {
        SAMTOOLS_REHEADER(
            ch_bam,
            "-c 'grep -v ^@RG'"
        )
        ch_bam = SAMTOOLS_REHEADER.out.bam
        ch_versions = ch_versions.mix(SAMTOOLS_REHEADER.out.versions)
    }
    QUALIMAP_BAMQC(
        ch_bam
    )
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions)

    // Tracking versions
    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // Final Reports
    MULTIQC_SAMPLE(
        ch_multiqc_sample_conf,
        ch_sample_csv
            .join(CREATE_READ_VARIATION_CSV.out.csv, by: [0])
            .join(CREATE_VARIANT_TSV.out.tsv, by: [0])
            .join(QUALIMAP_BAMQC.out.results, by: [0])
            .join(ch_sample_amplicon_depth, by: [0])
    )

    MULTIQC_OVERALL(
        ch_multiqc_overall_conf,
        ch_sample_amplicon_depth
            .collect{ it[1] }
            .ifEmpty([]),
        ch_amplicon_completeness
            .ifEmpty([]),
        BCFTOOLS_STATS.out.stats
            .collect{ it[1] },
        SAMTOOLS_FLAGSTAT.out.flagstat
            .collect{ it[1] },
        QUALIMAP_BAMQC.out.results
            .collect{ it[1] },
        ch_snpeff_csv
            .collect{ it[1] }
            .ifEmpty([]),
        ch_overall_qc_csv,
        CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml
    )
}
