/*
    Subworkflow to run pipeline steps for subconsensus nanopore data

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Other tools
include { CLAIRSTO_VARIANTS            } from '../../../modules/local/nanopore_subconsensus/main'

// Artic subcommands steps
include { ZIP_AND_INDEX_VCF         } from '../../../modules/local/artic_subcommands/main'


// add the python scripts later
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
    ch_primary_vcf = CLAIRSTO_VARIANTS.out.vcf
    ch_versions = ch_versions.mix(CLAIRSTO_VARIANTS.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Vcf reformatting
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    //combine indels and snvs
    //CAT_VCF(
    //    CLAIRSTO_VARIANTS.out.primary_vcf
    //)
    //join to main vcf and remove duplicates
    //DEDUP_VCFS(
    //    CAT_VCF.out.vcf
    //       .join(ch_con_vcf, by: [0])
    //)
    //adjust filters and qual scores
    //FIX_VCF(
    //    DEDUP_VCFS.out.vcf
    //)
    //ch_complete_vcf = FIX_VCF.out.vcf

    // do i need this?

    //ZIP_AND_INDEX_VCF(
    //    ch_primary_vcf
    //)
    //ch_versions = ch_versions.mix(ZIP_AND_INDEX_VCF.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    emit:
    //vcf = ch_complete_vcf just produce the vcf for now
    vcf = ch_primary_vcf
    versions = ch_versions
}
