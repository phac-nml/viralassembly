/*
    Subworkflow to run pipeline steps for subconsensus nanopore data

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLAIRSTO_VARIANTS            } from '../../../modules/local/nanopore_subconsensus/main'
include { CAT_VCF                      } from '../../../modules/local/nanopore_subconsensus/main'
include { DEDUP_VCFS                   } from '../../../modules/local/nanopore_subconsensus/main'
include { FIX_VCF                      } from '../../../modules/local/nanopore_subconsensus/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow WF_NANOPORE_SUBCONSENSUS {
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
    ch_complete_vcf = FIX_VCF.out.vcf

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    emit:
    vcf = ch_complete_vcf
    versions = ch_versions
}
