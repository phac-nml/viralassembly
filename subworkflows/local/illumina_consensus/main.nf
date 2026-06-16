//
// Subworkflow for amplicon and non-amplicon consensus sequence generation for Illumina data
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Initial Steps
include { FASTP                         } from '../../../modules/nf-core/fastp/main'

// Amplicon Specific
include { ARTIC_ALIGN_TRIM              } from '../../../modules/local/artic/align_trim/main'

// Alignment
include { BOWTIE2_BUILD                 } from '../../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN                 } from '../../../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_SORT                 } from '../../../modules/nf-core/samtools/sort/main'

// Variant Calling (iVar) and Consensus Generation
include { IVAR_VARIANTS                 } from '../../../modules/local/ivar/variants/main'
include { IVAR_VARIANTS_TO_VCF          } from '../../../modules/local/ivar/variants_to_vcf/main'
include { BCFTOOLS_SORT                 } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_FILTER               } from '../../../modules/nf-core/bcftools/filter/main'
include { TABIX_TABIX                   } from '../../../modules/nf-core/tabix/tabix/main'
include { MAKE_BED_MASK                 } from '../../../modules/local/make_bed_mask/main'
include { BEDTOOLS_MERGE                } from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_MASKFASTA            } from '../../../modules/nf-core/bedtools/maskfasta/main'

// Variant Calling (Freebayes) and Consensus Generation
include { FREEBAYES                     } from '../../../modules/local/freebayes/main'
include { PROCESS_VCF                   } from '../../../modules/local/process_vcf/main'
include { CUSTOM_MAKE_DEPTH_MASK        } from '../../../modules/local/artic_subcommands/make_depth_mask/main'
include { BCFTOOLS_CONSENSUS as BCFTOOLS_CONSENSUS_AMBIGUOUS  } from '../../../modules/local/bcftools/consensus/main'

// Output Final Consensus and Adjust Sequence Header
include { BCFTOOLS_CONSENSUS as BCFTOOLS_CONSENSUS_FINAL      } from '../../../modules/local/bcftools/consensus/main'
include { ADJUST_FASTA_HEADER           } from '../../../modules/local/artic_subcommands/adjust_fasta_header/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO SETUP REFERENCE DATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow WF_ILLUMINA_CONSENSUS {

    take:
    ch_fastqs       // channel: [ val(meta), file(fastqs) ]
    ch_reference    // channel: [ val(meta), file(reference) ]
    ch_fai          // channel: [ file(fai) ]
    ch_primer_bed   // channel: [ file(primer.bed) ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Run fastp for read quality filtering
    //
    FASTP(
        ch_fastqs,
        params.fastp_adapter ? file(params.fastp_adapter, type: 'file', checkIfExists: true) : [],
        '',
        '',
        ''
    )
    ch_fastqs = FASTP.out.reads
    ch_versions = ch_versions.mix(FASTP.out.versions)

    // Pass/fail reads based on count after length and quality filtering
    ch_fastqs
        .branch{ _meta, fastq ->
            pass: fastq[0].countFastq() >= params.min_reads
            empty: fastq[0].countFastq() < params.min_reads
        }.set{ ch_filtered_fastqs }

    //
    // MODULE: Create index of reference
    //
    BOWTIE2_BUILD(
        ch_reference
    )
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)

    //
    // MODULE: Run BOWTIE2 to map to reference
    //
    BOWTIE2_ALIGN(
        ch_filtered_fastqs.pass,
        BOWTIE2_BUILD.out.index,
        ch_reference,
        '',
        ''
    )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

    //
    // MODULE: Sort and index output bam file from BWA/BOWTIE2
    //
    SAMTOOLS_SORT(
        BOWTIE2_ALIGN.out.bam,
        ch_reference,
        'bai'
    )
    ch_bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_SORT.out.bai, by: [0])
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    if ( params.primer_bed ) {
        ARTIC_ALIGN_TRIM(
            ch_bam_bai,
            ch_primer_bed,
            'primers'
        )
        ch_bam_bai = ARTIC_ALIGN_TRIM.out.bam
        ch_versions = ch_versions.mix(ARTIC_ALIGN_TRIM.out.versions)
    }

    if (params.use_ivar) {
        //
        // MODULE: Call variants with iVar
        //
        IVAR_VARIANTS(
            ch_bam_bai,
            ch_reference.collect{ _meta, ref -> ref },
            ch_fai
        )
        ch_versions = ch_versions.mix(IVAR_VARIANTS.out.versions)

        //
        // MODULE: Convert iVar output to vcf
        //
        IVAR_VARIANTS_TO_VCF(
            IVAR_VARIANTS.out.tsv,
            ch_reference.collect{ _meta, ref -> ref }
        )
        ch_versions = ch_versions.mix(IVAR_VARIANTS_TO_VCF.out.versions)

        //
        // MODULE: Sort variants with BCFTools
        //
        BCFTOOLS_SORT(
            IVAR_VARIANTS_TO_VCF.out.vcf
        )
        ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

        //
        // MODULE: Filter variants by alle frequency
        //
        BCFTOOLS_FILTER(
            BCFTOOLS_SORT.out.vcf.map{ meta, vcf -> [ meta, vcf, [] ] } // Just adding in blank tbi index
        )
        ch_vcf = BCFTOOLS_FILTER.out.vcf
        ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions)

        //
        // MODULE: Index vcf file with Tabix
        //
        TABIX_TABIX(
            BCFTOOLS_FILTER.out.vcf
        )
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

        ch_bam_vcf = ch_bam_bai
            .join(BCFTOOLS_FILTER.out.vcf, by: [0])

        //
        // MODULE: Create BED file to mask consensus regions
        //
        MAKE_BED_MASK(
            ch_bam_vcf,
            ch_reference.collect{ _meta, ref -> ref }
        )
        ch_versions = ch_versions.mix(MAKE_BED_MASK.out.versions)

        //
        // MODULE: Merge intervals with BEDTools
        //
        BEDTOOLS_MERGE(
            MAKE_BED_MASK.out.bed
        )
        ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions)

        //
        // MODULE: Mask consensus regions with BEDTools
        //
        BEDTOOLS_MASKFASTA(
            BEDTOOLS_MERGE.out.bed,
            ch_reference.collect{ _meta, ref -> ref }
        )
        ch_versions = ch_versions.mix(BEDTOOLS_MASKFASTA.out.versions)

        //
        // MODULE: Create final consensus sequence with all variants
        //
        ch_bcfcons_in = BEDTOOLS_MASKFASTA.out.fasta
            .join(BCFTOOLS_FILTER.out.vcf, by: [0])
            .join(TABIX_TABIX.out.tbi, by: [0])
            .map { meta, fasta, vcf, tbi ->
                [ meta, vcf, tbi, fasta, [] ]
            }
        BCFTOOLS_CONSENSUS_FINAL(
            ch_bcfcons_in
        )
        ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS_FINAL.out.versions)

        //
        // MODULE: Adjust final consensus sequence headers to make downstream processes easier
        //
        ADJUST_FASTA_HEADER(
            BCFTOOLS_CONSENSUS_FINAL.out.fasta,
            ch_reference,
            '.consensus',
            ''
        )
        ch_consensus = ADJUST_FASTA_HEADER.out.consensus
        ch_versions = ch_versions.mix(ADJUST_FASTA_HEADER.out.versions)
    } else {
        //
        // MODULE: Run Freebayes to call variants
        //
        FREEBAYES(
            ch_bam_bai,
            ch_reference.collect{ _meta, ref -> ref }
        )
        ch_versions = ch_versions.mix(FREEBAYES.out.versions)
        //
        // MODULE: Process freebayes variant calls with custom python script and bcftools norm
        //
        PROCESS_VCF(
            FREEBAYES.out.vcf,
            ch_reference.collect{ _meta, ref -> ref }
        )
        ch_versions = ch_versions.mix(PROCESS_VCF.out.versions)
        ch_vcf = PROCESS_VCF.out.consensus_vcf.map { meta, vcf, _tbi -> tuple(meta, vcf) }
        //
        // MODULE: Make a depth mask based on the minimum depth required to call a position
        //
        CUSTOM_MAKE_DEPTH_MASK(
            ch_bam_bai,
            ch_reference.collect{ _meta, ref -> ref }
        )
        ch_versions = ch_versions.mix(CUSTOM_MAKE_DEPTH_MASK.out.versions)
        //
        // MODULE: Create intermediate fasta file with IUPACs for ambiguous positions from freebayes
        //
        BCFTOOLS_CONSENSUS_AMBIGUOUS(
            PROCESS_VCF.out.ambiguous_vcf
                .combine(ch_reference)
                .map { meta, vcf, tbi, _meta_ref, ref ->
                    tuple(meta, vcf, tbi, ref, [])
                }
        )

        //
        // MODULE: Create final consensus sequence with all variants
        //
        BCFTOOLS_CONSENSUS_FINAL(
            BCFTOOLS_CONSENSUS_AMBIGUOUS.out.fasta
                .join(CUSTOM_MAKE_DEPTH_MASK.out.coverage_mask, by: [0])
                .join(PROCESS_VCF.out.consensus_vcf, by: [0])
                .map { meta, fasta, mask, vcf, tbi ->
                    [ meta, vcf, tbi, fasta, mask ]
                }
        )
        ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS_FINAL.out.versions)

        //
        // MODULE: Adjust final consensus sequence headers to contain sample id and reference info
        //
        ADJUST_FASTA_HEADER(
            BCFTOOLS_CONSENSUS_FINAL.out.fasta,
            ch_reference,
            '.consensus',
            ''
        )
        ch_consensus = ADJUST_FASTA_HEADER.out.consensus
        ch_versions = ch_versions.mix(ADJUST_FASTA_HEADER.out.versions)
    }

    emit:
    consensus               = ch_consensus              // channel: [ val(meta), file(consensus) ]
    bam                     = ch_bam_bai                // channel: [ val(meta), file(bam), file(bai) ]
    vcf                     = ch_vcf                    // channel: [ val(meta), file(vcf) ]
    empty_filtered_fastqs   = ch_filtered_fastqs.empty  // channel: [ val(meta), file(fastqs) ]
    versions                = ch_versions               // channel: [ path(versions.yml) ]
}
