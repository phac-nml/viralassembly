/*
    Subworkflow for amplicon and non-amplicon consensus sequence generation for Nanopore data

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Clair3 model
include { ARTIC_GET_MODELS          } from '../../../modules/local/artic/get_model/main'

// Read QC
include { ARTIC_GUPPYPLEX           } from '../../../modules/local/artic/guppyplex/main'
include { CHOPPER                   } from '../../../modules/local/chopper/main'
include { NANOSTAT                  } from '../../../modules/local/nanostat/main'
include { RENAME_FASTQ              } from '../../../modules/local/custom/utils.nf'

// Alignment
include { MINIMAP2_ALIGN            } from '../../../modules/local/minimap2/main'

// Amplicon Specific Tools
include { ARTIC_ALIGN_TRIM          } from '../../../modules/local/artic/align_trim/main'
include { SPLIT_BED_BY_POOL         } from '../../../modules/local/custom/utils.nf'
include { CREATE_TILING_BED         } from '../../../modules/local/custom/utils.nf'
include { ARTIC_VCF_MERGE           } from '../../../modules/local/artic/vcf_merge/main'
include { ARTIC_MAKE_DEPTH_MASK     } from '../../../modules/local/artic/make_depth_mask/main'

// Variant Calling & Handling
include { CLAIR3_VARIANTS           } from '../../../modules/local/clair3/main'
include { ARTIC_VCF_FILTER         } from '../../../modules/local/artic/vcf_filter/main'
include { CUSTOM_MAKE_DEPTH_MASK    } from '../../../modules/local/artic/make_depth_mask/main'

// Consensus Generation
include { ARTIC_MASK                } from '../../../modules/local/artic/mask/main'
include { BCFTOOLS_NORM             } from '../../../modules/local/bcftools/norm/main'
include { BCFTOOLS_CONSENSUS        } from '../../../modules/local/bcftools/consensus/main'

// Artic Analysis Pipeline
include { ARTIC_MINION              } from '../../../modules/local/artic/minion/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow WF_NANOPORE_CONSENSUS {
    take:
    ch_fastqs       // channel: [ val(meta), file(fastqs) ]
    ch_reference    // channel: [ val(meta), file(reference) ]
    ch_fai          // channel: [ file(fai) ]
    ch_primer_bed   // channel: [ file(primer.bed) ]
    ch_amplicon_bed // channel: [ file(amplicon.bed) ]

    main:
    // Tool version tracking
    ch_versions = channel.empty()

    // Clair3 model
    ch_model = channel.empty()
    if ( params.local_model ) {
        ch_model = file(params.local_model, type: 'dir', checkIfExists: true)
    } else {
        ARTIC_GET_MODELS(params.model)
        ch_model = ARTIC_GET_MODELS.out.model
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Read QC and Statistics
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Filter reads by read lengths
    ARTIC_GUPPYPLEX(
        ch_fastqs
    )
    ch_fastqs = ARTIC_GUPPYPLEX.out.fastq
    ch_versions = ch_versions.mix(ARTIC_GUPPYPLEX.out.versions)


    // Filter reads by quality
    CHOPPER(
        ch_fastqs
    )
    ch_versions = ch_versions.mix(CHOPPER.out.versions)
    ch_fastqs = CHOPPER.out.fastq

    // Pass/fail reads based on count after length and quality filtering
    ch_fastqs
        .branch{ _meta, fastq ->
            pass: fastq.countFastq() >= params.min_reads
            empty: fastq.countFastq() < params.min_reads
        }.set{ ch_filtered_fastqs }

    // Stats on final reads
    NANOSTAT(
        ch_filtered_fastqs.pass
    )
    ch_versions = ch_versions.mix(NANOSTAT.out.versions)

    if ( !params.use_artic_tool ) {
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Alignment
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Align reads to the reference
        MINIMAP2_ALIGN(
            ch_filtered_fastqs.pass,
            ch_reference
        )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
        ch_bam = MINIMAP2_ALIGN.out.bam
        ch_clair3_input = ch_bam
            .map { meta, bam, bai ->
                tuple(meta, bam, bai, '', [])
            }

        if ( params.primer_bed ) {
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
            // Amplicon Specific Alignment Processing
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
            // Softmask read alignments within their derived amplicon with additional softmasking to exclude primer sequences
            ARTIC_ALIGN_TRIM (
                ch_bam,
                ch_primer_bed,
                'primers',
                params.platform
            )
            ch_bam = ARTIC_ALIGN_TRIM.out.bam
            ch_versions = ch_versions.mix(ARTIC_ALIGN_TRIM.out.versions)

            // For clair3 need bed files for each amplicon pool named <POOL>.bed
            //  Clair3 doesn't seem to be dealing with the bed files as expected
            //  As such, add option to not split by pool and instead use the whole tiling region
            if ( ! params.no_pool_split ) {
                SPLIT_BED_BY_POOL(
                    ch_amplicon_bed
                )
                ch_bed_pools = SPLIT_BED_BY_POOL.out.bed
                    .flatten()
                    .map{ bed -> [ bed.baseName.replaceAll(~/\.bed$/, ''), file(bed) ] }
            } else {
                CREATE_TILING_BED(
                    ch_amplicon_bed
                )
                ch_bed_pools = CREATE_TILING_BED.out.bed
                    .map { bed -> [ bed.baseName.replaceAll(~/\.bed$/, ''), file(bed) ] }
            }

            // Clair3 also uses the primer trimmed bams
            //  Based on testing, the way the pools and clair3 work, having the primer trimmed bams as input
            //  allows better and consistent calling in SNPs in primers
            ch_clair3_input = ch_bam
                .combine(ch_bed_pools) // channel: [ val(meta), path(bam), path(bai), val(pool), path(pool_bed) ]
        }
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Variant Calling
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        CLAIR3_VARIANTS(
            ch_clair3_input,
            ch_reference,
            ch_fai,
            ch_model,
            params.no_pool_split
        )
        ch_primary_vcf = CLAIR3_VARIANTS.out.vcf
        ch_versions = ch_versions.mix(CLAIR3_VARIANTS.out.versions)

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Variant Handling
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        if ( params.primer_bed ) {
            // Utilizing transformVCFList function to merge vcfs based on how artic handles input data
            ARTIC_VCF_MERGE(
                ch_primary_vcf
                    .groupTuple(),
                ch_primer_bed
            )
            ch_versions = ch_versions.mix(ARTIC_VCF_MERGE.out.versions)
            ch_primary_vcf = ARTIC_VCF_MERGE.out.vcf
        }

        // Filter VCF variants on depth, quality, and frequency
        ARTIC_VCF_FILTER(
            ch_primary_vcf
        )
        ch_versions = ch_versions.mix(ARTIC_VCF_FILTER.out.versions)

        // Make depth mask based on minimum depth to call position
        if ( params.primer_bed ) {
            ARTIC_MAKE_DEPTH_MASK(
                ch_bam,
                ch_reference
            )
            ch_depth_mask = ARTIC_MAKE_DEPTH_MASK.out.coverage_mask
            ch_versions = ch_versions.mix(ARTIC_MAKE_DEPTH_MASK.out.versions)
        } else {
            CUSTOM_MAKE_DEPTH_MASK(
                ch_bam,
                ch_reference
            )
            ch_depth_mask = CUSTOM_MAKE_DEPTH_MASK.out.coverage_mask
            ch_versions = ch_versions.mix(CUSTOM_MAKE_DEPTH_MASK.out.versions)
        }
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Consensus Generation
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Apply Masking
        ARTIC_MASK(
            ch_depth_mask
                .join(ARTIC_VCF_FILTER.out.fail_vcf, by: [0]),
            ch_reference
        )
        ch_versions = ch_versions.mix(ARTIC_MASK.out.versions)

        // Normalize variants for consensus generation
        BCFTOOLS_NORM(
            ARTIC_MASK.out.preconsensus
                .join(ARTIC_VCF_FILTER.out.pass_vcf, by: [0])
        )
        ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

        // Create final consensus sequence with all variants
        BCFTOOLS_CONSENSUS(
            ARTIC_MASK.out.preconsensus
                .join(ch_depth_mask, by: [0])
                .join(BCFTOOLS_NORM.out.vcf, by: [0])
                .map { meta, fasta, mask, vcf, tbi ->
                    [ meta, vcf, tbi, fasta, mask ]
                }
        )
        ch_consensus = BCFTOOLS_CONSENSUS.out.consensus
        ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)

        // Remove tabix index from vcf as it is not needed
        ch_vcf = ARTIC_VCF_FILTER.out.pass_vcf
            .map { meta, vcf, _tbi -> [ meta, vcf ] }
    } else {
        ARTIC_MINION(
            ch_filtered_fastqs.pass,
            ch_reference,
            ch_primer_bed,
            ch_model
        )
        ch_consensus = ARTIC_MINION.out.consensus
        ch_bam = ARTIC_MINION.out.bam
        ch_vcf = ARTIC_MINION.out.vcf
        ch_versions = ch_versions.mix(ARTIC_MINION.out.versions)
    }
    emit:
    consensus             = ch_consensus                // channel: [ val(meta), file(consensus) ]
    bam                   = ch_bam                      // channel: [ val(meta), file(bam), file(bai) ]
    vcf                   = ch_vcf                      // channel: [ val(meta), file(vcf) ]
    empty_filtered_fastqs = ch_filtered_fastqs.empty    // channel: [ val(meta), file(fastqs) ]
    stats                 = NANOSTAT.out.stats          // channel: [ val(meta), file(nanostat.txt) ]
    versions              = ch_versions                 // channel: [ path(versions.yml) ]
}
