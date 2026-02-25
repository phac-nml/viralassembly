/*
    Subworkflow to annotate VCF file using SnpEFF
        1. Checks if a database is available
        2. If it is downloads it, otherwise attempts to make it from NCBI refseq genbank file
        3. Annotates VCF file using downloaded or made DB

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SNPEFF_DATABASE   } from '../../../modules/local/snpeff/database/main'
include { SNPEFF_ANNOTATE   } from '../../../modules/local/snpeff/annotate/main'
include { ZIP_AND_INDEX_VCF } from '../../../modules/local/artic_subcommands/zip_and_index/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow WF_SNPEFF_ANNOTATE {
    take:
    ch_vcf          // channel: [ val(meta), file(vcf) ]
    ch_reference    // channel: [ path(reference) ]

    main:
    // Version tracking
    ch_gff = params.gff ? file(params.gff, type: 'file', checkIfExists: true) : []
    ch_versions = channel.empty()

    // Get reference id
    ch_reference.splitFasta( record: [ id: true ] )
        .map{ record -> record.id.toString() }
        .collect() // To collect segmented and turn to a value channel
        .set{ ch_ref_ids }

    SNPEFF_DATABASE(
        ch_ref_ids,
        ch_reference,
        ch_gff
    )
    ch_versions = ch_versions.mix(SNPEFF_DATABASE.out.versions)

    SNPEFF_ANNOTATE(
        ch_vcf,
        SNPEFF_DATABASE.out.db,
        SNPEFF_DATABASE.out.config
            .ifEmpty([])
    )
    ch_versions = ch_versions.mix(SNPEFF_DATABASE.out.versions)

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
    vcf = ch_ann_vcf
    csv = SNPEFF_ANNOTATE.out.csv
    versions = ch_versions
}
