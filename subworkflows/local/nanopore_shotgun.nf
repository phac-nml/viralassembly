/*
    Subworkflow to run pipeline steps for non-amplicon nanopore data
        We will follow a similar process to the artic steps but not use any pools
        Those are their own spots

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Other tools
include { MINIMAP2_ALIGN            } from '../../modules/local/minimap2/main'
include { LONGSHOT                  } from '../../modules/local/longshot/main'
include { BCFTOOLS_NORM             } from '../../modules/local/bcftools/norm/main'
include { BCFTOOLS_CONSENSUS        } from '../../modules/local/bcftools/consensus/main'

// Variant calling tools
include { MEDAKA_CONSENSUS          } from '../../modules/local/nanopore_shotgun/main'
include { MEDAKA_VARIANT            } from '../../modules/local/nanopore_shotgun/main'
include { NANOPOLISH_VARIANTS       } from '../../modules/local/nanopore_shotgun/main'
include { CLAIR3_VARIANTS           } from '../../modules/local/nanopore_shotgun/main'

// Artic subcommands steps
include { ZIP_AND_INDEX_VCF         } from '../../modules/local/artic_subcommands/main'
include { CUSTOM_VCF_FILTER         } from '../../modules/local/artic_subcommands/main'
include { CUSTOM_MAKE_DEPTH_MASK    } from '../../modules/local/artic_subcommands/main'
include { ARTIC_MASK                } from '../../modules/local/artic_subcommands/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INITIALIZE CHANNELS FROM PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
ch_user_clair3_model = params.clair3_user_variant_model ? file(params.clair3_user_variant_model, checkIfExists: true) : []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow WF_NANOPORE_SHOTGUN {
    take:
    ch_fastqs       // channel: [ val(meta), file(fastq) ]
    ch_fast5s       // channel: [ file(fast5s) ]
    ch_seqSum       // channel: [ file(sequencing_summary) ]
    ch_reference    // channel: [ file(reference) ]
    ch_ref_fai      // channel: [ file(reference.fai) ]
    ch_refstats     // channel: [ file(refstats.txt) ]

    main:
    // Version tracking
    ch_versions = Channel.empty()

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Align
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    MINIMAP2_ALIGN(
        ch_fastqs,
        ch_reference
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    ch_bam = MINIMAP2_ALIGN.out.bam

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Variant Calling
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    if ( params.variant_caller == 'medaka' ) {
        MEDAKA_CONSENSUS(
            ch_bam
        )
        ch_versions = ch_versions.mix(MEDAKA_CONSENSUS.out.versions)
        MEDAKA_VARIANT(
            MEDAKA_CONSENSUS.out.hdf,
            ch_reference
        )
        ch_primary_vcf = MEDAKA_VARIANT.out.vcf
        ch_versions = ch_versions.mix(MEDAKA_VARIANT.out.versions)
    } else if ( params.variant_caller == 'nanopolish' ) {
        NANOPOLISH_VARIANTS(
            ch_fastqs
                .combine(ch_bam, by: [0]),
            ch_fast5s,
            ch_seqSum,
            ch_reference,
            ch_refstats
        )
        ch_primary_vcf = NANOPOLISH_VARIANTS.out.vcf
        ch_versions = ch_versions.mix(NANOPOLISH_VARIANTS.out.versions)
    } else {
        CLAIR3_VARIANTS(
            ch_bam,
            ch_reference,
            ch_ref_fai,
            ch_user_clair3_model
        )
        ch_primary_vcf = CLAIR3_VARIANTS.out.vcf
        ch_versions = ch_versions.mix(CLAIR3_VARIANTS.out.versions)
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Variant Handling
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Longshot for medaka only at the moment as that is how it is in artic minion!
    //  And I have yet to test it elsewhere
    if ( params.variant_caller == 'medaka' && ! params.skip_longshot ) {
        ZIP_AND_INDEX_VCF(
            ch_primary_vcf
        )
        ch_versions = ch_versions.mix(ZIP_AND_INDEX_VCF.out.versions)

        LONGSHOT(
            ZIP_AND_INDEX_VCF.out.vcf
                .join(ch_bam, by: [0]),
            ch_reference,
            ch_ref_fai
        )
        ch_intermediate_vcf = LONGSHOT.out.vcf
        ch_versions = ch_versions.mix(LONGSHOT.out.versions)
    } else {
        ch_intermediate_vcf = ch_primary_vcf
    }

    // Local variant filtering module based on artic_vcf_filter written to add clair3 in
    CUSTOM_VCF_FILTER(
        ch_intermediate_vcf
    )
    ch_versions = ch_versions.mix(CUSTOM_VCF_FILTER.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Consensus Generation
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Note that this step is slow and while I explored other methods
    //  Like bedtools or samtools depth I couldn't quite get the same formatting
    //  But will look more into it
    CUSTOM_MAKE_DEPTH_MASK(
        ch_bam,
        ch_reference
    )
    // ch_versions = ch_versions.mix(CUSTOM_MAKE_DEPTH_MASK.out.versions)

    ARTIC_MASK(
        CUSTOM_MAKE_DEPTH_MASK.out.coverage_mask
            .join(CUSTOM_VCF_FILTER.out.fail_vcf, by: [0]),
        ch_reference
    )
    ch_versions = ch_versions.mix(ARTIC_MASK.out.versions)

    // Norm to fix variants in both pass and fail VCF
    BCFTOOLS_NORM(
        ARTIC_MASK.out.preconsensus
            .join(CUSTOM_VCF_FILTER.out.pass_vcf, by: [0])
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    BCFTOOLS_CONSENSUS(
        ARTIC_MASK.out.preconsensus
            .join(CUSTOM_MAKE_DEPTH_MASK.out.coverage_mask, by: [0])
            .join(BCFTOOLS_NORM.out.vcf, by: [0])
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)

    // Remove tabix index from vcf as it is not needed and won't match the normal artic steps as output
    CUSTOM_VCF_FILTER.out.pass_vcf
        .map { it -> [ it[0], it[1] ] }
        .set { ch_pass_vcf }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    emit:
    consensus = BCFTOOLS_CONSENSUS.out.consensus
    bam = ch_bam
    vcf = ch_pass_vcf
    versions = ch_versions
}
