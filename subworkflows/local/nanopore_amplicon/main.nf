/*
    Subworkflow to run the artic pipeline steps as nextflow steps - see https://artic.readthedocs.io/en/latest/minion/#core-pipeline
        This is being tested to allow:
            - Better updating of dependencies
            - Clair3 or other variant callers to be used
            - Better flexibility in analysis

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Artic subcommands steps
include { ARTIC_ALIGN_TRIM as ARTIC_ALIGN_TRIM_START    } from '../../../modules/local/artic/align_trim/main'
include { ARTIC_ALIGN_TRIM as ARTIC_ALIGN_TRIM_PRIMERS  } from '../../../modules/local/artic/align_trim/main'
include { ARTIC_VCF_MERGE           } from '../../../modules/local/artic_subcommands/vcf_merge/main'
include { ZIP_AND_INDEX_VCF         } from '../../../modules/local/artic_subcommands/zip_and_index/main'
include { CUSTOM_VCF_FILTER         } from '../../../modules/local/artic_subcommands/vcf_filter/main'
include { ARTIC_MAKE_DEPTH_MASK     } from '../../../modules/local/artic_subcommands/make_depth_mask/main'
include { ARTIC_MASK                } from '../../../modules/local/artic_subcommands/mask/main'

// Variant calling tools
include { MEDAKA_CONSENSUS          } from '../../../modules/local/nanopore_amplicon/medaka/main'
include { MEDAKA_VARIANT            } from '../../../modules/local/nanopore_amplicon/medaka/main'
include { NANOPOLISH_VARIANTS       } from '../../../modules/local/nanopore_amplicon/nanopolish/main'
include { CLAIR3_VARIANTS           } from '../../../modules/local/nanopore_amplicon/clair3/main'

// Other tools
include { SPLIT_BED_BY_POOL         } from '../../../modules/local/custom/utils.nf'
include { CREATE_TILING_BED         } from '../../../modules/local/custom/utils.nf'
include { MINIMAP2_ALIGN            } from '../../../modules/local/minimap2/main'
include { LONGSHOT                  } from '../../../modules/local/longshot/main'
include { BCFTOOLS_NORM             } from '../../../modules/local/bcftools/norm/main'
include { BCFTOOLS_CONSENSUS as BCFTOOLS_CONSENSUS  } from '../../../modules/local/bcftools/consensus/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INITIALIZE CHANNELS FROM PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow WF_NANOPORE_AMPLICON {
    take:
    ch_fastqs       // channel: [ val(meta), file(fastq) ]
    ch_fast5s       // channel: [ file(fast5s) ]
    ch_seqSum       // channel: [ file(sequencing_summary) ]
    ch_reference    // channel: [ file(reference) ]
    ch_ref_fai      // channel: [ file(reference.fai) ]
    ch_refstats     // channel: [ file(refstats.txt) ]
    ch_primer_bed   // channel: [ file(primer_bed) ]
    ch_amplicon_bed // channel: [ file(amplicon_bed) ]
    ch_clair3_model // channel: [ file(clair3_model) ]

    main:
    // Version tracking
    ch_versions = channel.empty()

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Reference Align
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    MINIMAP2_ALIGN(
        ch_fastqs,
        ch_reference.collect{ _meta, ref -> ref }
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Alignment Processing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Start = softmask read alignments within their derived amplicon
    ARTIC_ALIGN_TRIM_START(
        MINIMAP2_ALIGN.out.bam,
        ch_primer_bed,
        'start',
        params.platform
    )
    ch_versions = ch_versions.mix(ARTIC_ALIGN_TRIM_START.out.versions)

    // Primers = softmask read alignments within their derived amplicon with additional softmasking to exclude primer sequences
    ARTIC_ALIGN_TRIM_PRIMERS (
        MINIMAP2_ALIGN.out.bam,
        ch_primer_bed,
        'primers',
        params.platform
    )
    // Setting the primertrimmed bams to a channel as its part of steps other than variant calling
    //  which needs to be combined with the primer-pool to work properly
    ch_primertrimmed_bams = ARTIC_ALIGN_TRIM_PRIMERS.out.bam
    ch_versions = ch_versions.mix(ARTIC_ALIGN_TRIM_PRIMERS.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Variant Calling
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    if ( params.variant_caller != 'clair3' ) {
        // Gives a channel with pools as values
        pools = ch_primer_bed
            .splitCsv(sep: '\t')
            .collect { row -> row[4] }
            .flatten()
            .unique()

        // Medaka or Nanopolish
        if ( params.variant_caller == 'medaka' ) {
            // This combines the pool name values with the bams for other steps
            //  Medaka uses the Primer trimmed bam according to docs
            ch_trimmed_bams_w_pool = ARTIC_ALIGN_TRIM_PRIMERS.out.bam.combine(pools) // channel: [ val(meta), file(bam), file(bai), val(pool) ]

            MEDAKA_CONSENSUS(
                ch_trimmed_bams_w_pool
            )
            ch_versions = ch_versions.mix(MEDAKA_CONSENSUS.out.versions)
            MEDAKA_VARIANT(
                MEDAKA_CONSENSUS.out.hdf,
                ch_reference.collect{ _meta, ref -> ref }
            )
            ch_tmp_vcfs = MEDAKA_VARIANT.out.vcf
            ch_versions = ch_versions.mix(MEDAKA_VARIANT.out.versions)
        } else {
            // This combines the pool name values with the bams for other steps
            //  Nanopolish uses the Start trimmed bam according to docs
            ch_trimmed_bams_w_pool = ARTIC_ALIGN_TRIM_START.out.bam.combine(pools) // channel: [ val(meta), file(bam), file(bai), val(pool) ]

            NANOPOLISH_VARIANTS(
                ch_fastqs
                    .combine(ch_trimmed_bams_w_pool, by: [0]),
                ch_fast5s,
                ch_seqSum,
                ch_reference.collect{ _meta, ref -> ref },
                ch_refstats
            )
            ch_tmp_vcfs = NANOPOLISH_VARIANTS.out.vcf
            ch_versions = ch_versions.mix(NANOPOLISH_VARIANTS.out.versions)
        }
    } else {
        // For clair3 need bed files for each amplicon pool named <POOL>.bed
        //  Clair3 doesn't seem to be dealing with the bed files as expected
        //  As such, add option to not split by pool and instead use the whole tiling region
        if ( ! params.clair3_no_pool_split ) {
            SPLIT_BED_BY_POOL(
                ch_amplicon_bed
            )
            SPLIT_BED_BY_POOL.out.bed
                .flatten()
                .map{ bed -> [ bed.baseName.replaceAll(~/\.bed$/, ''), file(bed) ] }
                .set { ch_bed_pools }
        } else {
            CREATE_TILING_BED(
                ch_amplicon_bed
            )
            CREATE_TILING_BED.out.bed
                .map { bed -> [ bed.baseName.replaceAll(~/\.bed$/, ''), file(bed) ] }
                .set { ch_bed_pools }
        }

        // Clair3 also uses the primer trimmed bams
        //  Based on testing, the way the pools and clair3 work, having the primer trimmed bams as input
        //  allows better and consistent calling in SNPs in primers
        ch_trimmed_bams_w_pool = ARTIC_ALIGN_TRIM_PRIMERS.out.bam.combine(ch_bed_pools) // channel: [ val(meta), path(bam), path(bai), val(pool), path(pool_bed) ]

        // Run clair3
        CLAIR3_VARIANTS(
            ch_trimmed_bams_w_pool,
            ch_reference.collect{ _meta, ref -> ref },
            ch_ref_fai,
            ch_clair3_model,
            params.clair3_no_pool_split
        )
        ch_tmp_vcfs = CLAIR3_VARIANTS.out.vcf
        ch_versions = ch_versions.mix(CLAIR3_VARIANTS.out.versions)
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Variant Handling
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Merge pools by merging the vcf files for each pool together
    ch_tmp_vcfs
        .map { meta, vcf, _pool -> [ meta, vcf ] }
        .groupTuple()
        .set { ch_pooled_vcfs } // channel: [ val(meta), [vcf1, vcf2, ... vcfx] ]
    // To merge vcfs, have to utilize the transformVCFList function based on how artic handles input data
    ARTIC_VCF_MERGE(
        ch_pooled_vcfs,
        ch_primer_bed
    )
    ch_versions = ch_versions.mix(ARTIC_VCF_MERGE.out.versions)

    //
    // Longshot for medaka only as that is how it is in artic minion!
    if ( params.variant_caller == 'medaka' && ! params.skip_longshot ) {
        ZIP_AND_INDEX_VCF(
            ARTIC_VCF_MERGE.out.vcf
        )
        ch_versions = ch_versions.mix(ZIP_AND_INDEX_VCF.out.versions)

        LONGSHOT(
            ZIP_AND_INDEX_VCF.out.vcf
                .join(ch_primertrimmed_bams, by: [0]),
            ch_reference.collect{ _meta, ref -> ref },
            ch_ref_fai
        )
        ch_merged_vcf = LONGSHOT.out.vcf
        ch_versions = ch_versions.mix(LONGSHOT.out.versions)
    } else {
        ch_merged_vcf = ARTIC_VCF_MERGE.out.vcf
    }

    // Local variant filtering module based on artic_vcf_filter written to add clair3 in
    CUSTOM_VCF_FILTER(
        ch_merged_vcf
    )
    ch_versions = ch_versions.mix(CUSTOM_VCF_FILTER.out.versions)

    ARTIC_MAKE_DEPTH_MASK(
        ch_primertrimmed_bams,
        ch_reference.collect{ _meta, ref -> ref }
    )
    ch_versions = ch_versions.mix(ARTIC_MAKE_DEPTH_MASK.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Consensus Generation
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    ARTIC_MASK(
        ARTIC_MAKE_DEPTH_MASK.out.coverage_mask
            .join(CUSTOM_VCF_FILTER.out.fail_vcf, by: [0]),
        ch_reference.collect{ _meta, ref -> ref }
    )
    ch_versions = ch_versions.mix(ARTIC_MASK.out.versions)

    // Norm to fix variants in both pass and fail VCF
    BCFTOOLS_NORM(
        ARTIC_MASK.out.preconsensus
            .join(CUSTOM_VCF_FILTER.out.pass_vcf, by: [0])
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    //
    // MODULE: Create final consensus sequence with all variants
    //
    BCFTOOLS_CONSENSUS(
        ARTIC_MASK.out.preconsensus
            .join(ARTIC_MAKE_DEPTH_MASK.out.coverage_mask, by: [0])
            .join(BCFTOOLS_NORM.out.vcf, by: [0])
            .map { meta, fasta, mask, vcf, tbi ->
                [ meta, vcf, tbi, fasta, mask ]
            }
    )
    ch_consensus = BCFTOOLS_CONSENSUS.out.consensus
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)

    // Remove tabix index from vcf as it is not needed and won't match the normal artic steps as output
    CUSTOM_VCF_FILTER.out.pass_vcf
        .map { meta, vcf, _tbi -> [ meta, vcf ] }
        .set { ch_pass_vcf }


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    emit:
    consensus   = ch_consensus          // channel: [ val(meta), file(consensus) ]
    bam         = ch_primertrimmed_bams // channel: [ val(meta), file(bam), file(bai) ]
    vcf         = ch_pass_vcf           // channel: [ val(meta), file(vcf) ]
    versions    = ch_versions           // channel: [ path(versions.yml) ]
}
