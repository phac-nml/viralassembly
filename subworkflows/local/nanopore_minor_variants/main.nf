/*
    Subworkflow to run pipeline steps for calling minor variants from nanopore data

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLAIRSTO_VARIANTS            } from '../../../modules/local/nanopore_minor_variants/main'
include { CAT_VCF                      } from '../../../modules/local/nanopore_minor_variants/main'
include { CAT_VCF as CAT_FULL_VCF      } from '../../../modules/local/nanopore_minor_variants/main'
include { DEDUP_VCFS                   } from '../../../modules/local/nanopore_minor_variants/main'
include { FIX_VCF                      } from '../../../modules/local/nanopore_minor_variants/main'
include { PASS_VCF                     } from '../../../modules/local/nanopore_minor_variants/main'
include { ZIP_AND_INDEX_VCF            } from '../../../modules/local/artic_subcommands/zip_and_index/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow WF_NANOPORE_MINOR_VARIANTS {
    take:
    ch_bam       // channel: [ val(meta), file(bam), file(bai) ]
    ch_reference    // channel: [ file(reference) ]
    ch_ref_fai      // channel: [ file(reference.fai) ]
    ch_con_vcf          // channel: [  val(meta), file(vcf) ]
    ch_clairSTO_model // channel:  val(model_name)

    main:

    // Version tracking
    ch_versions = channel.empty()

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Variant Calling
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    CLAIRSTO_VARIANTS(
        ch_bam,
        ch_reference,
        ch_ref_fai,
        ch_clairSTO_model
    )
    ch_versions = ch_versions.mix(CLAIRSTO_VARIANTS.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Vcf reformatting
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Combine indels and snvs (printed to different vcfs by clairsto)
    CAT_VCF(
        CLAIRSTO_VARIANTS.out.vcf
    )
    ch_versions = ch_versions.mix(CAT_VCF.out.versions)

    ZIP_AND_INDEX_VCF(
        CAT_VCF.out.vcf
    )
    ch_versions = ch_versions.mix(ZIP_AND_INDEX_VCF.out.versions)

    // Remove consensus variants for readability of minor vcf
    DEDUP_VCFS(
            ZIP_AND_INDEX_VCF.out.vcf
                .join(ch_con_vcf, by: [0])
    )
    ch_versions = ch_versions.mix(DEDUP_VCFS.out.versions)

    // Adjust filters and qual scores for viral minor variants
    FIX_VCF(
        DEDUP_VCFS.out.vcf
    )
    ch_complete_min_vcf = FIX_VCF.out.vcf
    ch_versions = ch_versions.mix(FIX_VCF.out.versions)

    // Publish a joined VCF with passing minor and major variants
    CAT_FULL_VCF(
        FIX_VCF.out.vcf
            .join(ch_con_vcf, by: [0])
    )

    PASS_VCF(
        CAT_FULL_VCF.out.vcf
    )
    ch_versions = ch_versions.mix(PASS_VCF.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    emit:
    vcf = ch_complete_min_vcf
    versions = ch_versions
}
