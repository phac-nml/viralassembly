/*
    Subworkflow to run pipeline steps for calling minor variants from nanopore data

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLAIRSTO_VARIANTS            } from '../../../modules/local/nanopore_minor_variants/main'
include { CAT_VCF                      } from '../../../modules/local/nanopore_minor_variants/main'
include { DEDUP_VCFS                   } from '../../../modules/local/nanopore_minor_variants/main'
include { FIX_VCF                      } from '../../../modules/local/nanopore_minor_variants/main'
include { CAT_PASS_VCF                      } from '../../../modules/local/nanopore_minor_variants/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow WF_NANOPORE_MINOR_VARIANTS {
    take:
    ch_bam       // channel: [ val(meta), file(bam) ]
    ch_reference    // channel: [ file(reference) ]
    ch_ref_fai      // channel: [ file(reference.fai) ]
    ch_con_vcf          // channel: [  val(meta), file(vcf) ]

    main:

    // Version tracking
    ch_versions = channel.empty()

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Variant Calling
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    CLAIRSTO_VARIANTS(
        ch_bam,
        ch_reference,
        ch_ref_fai
    )
    ch_versions = ch_versions.mix(CLAIRSTO_VARIANTS.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Vcf reformatting
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Combine indels and snvs (printed to different vcfs)
    CAT_VCF(
        CLAIRSTO_VARIANTS.out.vcf
    )

    // Remove consensus variants
    DEDUP_VCFS(
        CAT_VCF.out.vcf
           .join(ch_con_vcf, by: [0])
    )
    ch_primary_vcf = DEDUP_VCFS.out.vcf

    // Adjust filters and qual scores for viral minor variants
    FIX_VCF(
        DEDUP_VCFS.out.vcf
    )
    ch_complete_min_vcf = FIX_VCF.out.vcf

    // Publish a joined VCF with passing minor and major variants
    CAT_PASS_VCF(
        FIX_VCF.out.vcf
           .join(ch_con_vcf, by: [0])
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    emit:
    vcf = ch_complete_min_vcf
    versions = ch_versions
}
