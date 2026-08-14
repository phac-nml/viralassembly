/*
    Subworkflow to annotate VCF file using SnpEFF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SNPEFF_ANNOTATE   } from '../../../modules/local/snpeff/annotate/main'
include { ZIP_AND_INDEX_VCF } from '../../../modules/local/artic/zip_and_index/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow WF_SNPEFF_ANNOTATE {
    take:
    ch_vcf              // channel: [ val(meta), file(vcf) ]
    ch_snpeff_db        // channel: [ val(genome), path(db) ]
    ch_snpeff_config    // channel: [ path(config) ]
    level               // Flag used to adjust name for minor variant vcf

    main:
    // Version tracking
    ch_versions = channel.empty()

    SNPEFF_ANNOTATE(
        ch_vcf,
        ch_snpeff_db,
        ch_snpeff_config.ifEmpty([]),
        level
    )
    ch_versions = ch_versions.mix(SNPEFF_ANNOTATE.out.versions)

    // Zip and index vcf to match pass vcf
    ZIP_AND_INDEX_VCF(
        SNPEFF_ANNOTATE.out.vcf
    )
    ch_versions = ch_versions.mix(ZIP_AND_INDEX_VCF.out.versions)

    // Remove tabix index from vcf as it is not needed and won't match the normal artic steps as output
    ZIP_AND_INDEX_VCF.out.vcf
        .map { meta, vcf, _tbi -> [ meta, vcf ] }
        .set { ch_ann_vcf }

    emit:
    vcf         = ch_ann_vcf                // channel: [ val(meta), file(vcf) ]
    csv         = SNPEFF_ANNOTATE.out.csv   // channel: [ val(meta), file(csv) ]
    versions    = ch_versions               // channel: [ path(versions.yml) ]
}
