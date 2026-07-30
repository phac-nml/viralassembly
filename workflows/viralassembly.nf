/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Functions
include { isSegmented               } from '../subworkflows/local/utils_nfcore_viralassembly_pipeline'
include { fastaHeaderId             } from '../subworkflows/local/utils_nfcore_viralassembly_pipeline'

// Utils / Custom checks / Primer Validate
include { GET_REF_STATS             } from '../modules/local/custom/utils.nf'
include { RENAME_FASTQ              } from '../modules/local/custom/utils.nf'
include { PRIMALBEDTOOLS_VALIDATE   } from '../modules/local/primalbedtools/validate/main.nf'
include { PRIMALBEDTOOLS_AMPLICON   } from '../modules/local/primalbedtools/amplicon/main.nf'
include { TRACK_FILTERED_SAMPLES as TRACK_INITIAL_FILTERED_SAMPLES } from '../modules/local/custom/filtering.nf'
include { TRACK_FILTERED_SAMPLES as TRACK_SIZE_FILTERED_SAMPLES    } from '../modules/local/custom/filtering.nf'

// QC related
include { SAMTOOLS_DEPTH            } from '../modules/local/samtools/depth/main'
include { MAKE_SAMPLE_QC_CSV        } from '../modules/local/qc/main'
include { FINAL_QC_CSV              } from '../modules/local/qc/main'

// SnpEff related
include { SNPEFF_DATABASE   } from '../modules/local/snpeff/database/main'

// Subworkflows
include { WF_NANOPORE_CONSENSUS     } from '../subworkflows/local/nanopore_consensus'
include { WF_ILLUMINA_CONSENSUS     } from '../subworkflows/local/illumina_consensus'
include { WF_NANOPORE_MINOR_VARIANTS       } from '../subworkflows/local/nanopore_minor_variants'
include { WF_SNPEFF_ANNOTATE        } from '../subworkflows/local/snpeff_annotate'
include { WF_SNPEFF_ANNOTATE as   WF_SNPEFF_ANNOTATE_MIN    } from '../subworkflows/local/snpeff_annotate'
include { WF_NEXTCLADE              } from '../subworkflows/local/nextclade'
include { WF_VIRUS_COVID            } from '../subworkflows/local/virus_specific/covid'
include { WF_CREATE_MULTIQC_REPORTS } from '../subworkflows/local/create_multiqc_reports'
include { WF_CREATE_CUSTOM_REPORT   } from '../subworkflows/local/create_custom_report'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow VIRALASSEMBLY {
    take:
    ch_fastqs           // channel: [ val(meta), file(fastq) ]
    ch_empty_fastqs     // channel: [ val(meta), file(fastq) ]

    main:
    // Optional value channel files from params
    ch_metadata = params.metadata ? file(params.metadata, type: 'file', checkIfExists: true) : []
    ch_pcr_primer_bed = params.pcr_primer_bed ? file(params.pcr_primer_bed, type: 'file', checkIfExists: true) : []

    // Tool version tracking
    ch_versions = channel.empty()

    // Rename if given metadata and not using input param, separate out fastqs with little data
    if ( params.platform == 'nanopore' && (ch_metadata && !params.input) ) {
        RENAME_FASTQ(
            ch_fastqs,
            ch_metadata
        )
        ch_versions = ch_versions.mix(RENAME_FASTQ.out.versions)
        // Remap the id based on the new name
        RENAME_FASTQ.out.fastq
            .map{ fastq -> [ [id: fastq.baseName.replaceAll(~/\.fastq.*$/, '')], file(fastq) ] }
            .set{ ch_fastqs }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Scheme and Reference
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

    // Check if the reference is a segmented virus
    def segmented = isSegmented(params.reference)

    // Create reference channel
    ch_reference = params.reference
        ? channel.value( file(params.reference, type: 'file', checkIfExists: true) )
        : []

    // Create amplicon channels
    ch_primer_bed = params.primer_bed   ? channel.value(file(params.primer_bed, type: 'file', checkIfExists: true)) : []
    ch_amplicon_bed = channel.empty()

    if ( params.primer_bed ) {
        // Amplicon information
        PRIMALBEDTOOLS_VALIDATE(
            ch_primer_bed,
            ch_reference
        )
        PRIMALBEDTOOLS_AMPLICON(
            ch_primer_bed
        )
        ch_amplicon_bed = PRIMALBEDTOOLS_AMPLICON.out.bed
        ch_versions = ch_versions.mix(PRIMALBEDTOOLS_AMPLICON.out.versions)
    }

    // Reference stats and files for various processes
    GET_REF_STATS(
        ch_reference
    )
    ch_fai = GET_REF_STATS.out.fai
    ch_versions = ch_versions.mix(GET_REF_STATS.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Consensus Generation
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // WORKFLOWS: Generate consensus and supporting files for either Nanopore or Illumina data
    //
    if( params.platform == 'nanopore' ) {
        //
        // WORKFLOW: Nanopore
        //
        WF_NANOPORE_CONSENSUS(
            ch_fastqs,
            ch_reference,
            ch_fai,
            ch_primer_bed,
            ch_amplicon_bed
        )
        ch_consensus = WF_NANOPORE_CONSENSUS.out.consensus
        ch_bam = WF_NANOPORE_CONSENSUS.out.bam
        ch_vcf = WF_NANOPORE_CONSENSUS.out.vcf
        ch_filtered_fastqs_empty = WF_NANOPORE_CONSENSUS.out.empty_filtered_fastqs
        ch_reads_stats = WF_NANOPORE_CONSENSUS.out.stats
        ch_versions = ch_versions.mix(WF_NANOPORE_CONSENSUS.out.versions)

    } else if( params.platform == 'illumina' ) {
        //
        // WORKFLOW: Illumina
        //
        WF_ILLUMINA_CONSENSUS(
            ch_fastqs,
            ch_reference,
            ch_fai,
            ch_primer_bed
        )
        ch_consensus = WF_ILLUMINA_CONSENSUS.out.consensus
        ch_bam = WF_ILLUMINA_CONSENSUS.out.bam
        ch_vcf = WF_ILLUMINA_CONSENSUS.out.vcf
        ch_filtered_fastqs_empty = WF_ILLUMINA_CONSENSUS.out.empty_filtered_fastqs
        ch_reads_stats = ch_consensus.map { meta, _fasta -> tuple(meta, []) }
        ch_versions = ch_versions.mix(WF_ILLUMINA_CONSENSUS.out.versions)

    } else {
        error "Please provide the --platform parameter with either 'nanopore' or 'illumina' to run"
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    //  Call Minor variants (i.e. AF below consensus level), optional
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    ch_min_vcf = channel.empty()
    // Only Nanopore for now
    if ( params.minor_variants && params.platform == 'nanopore') {
        ch_clairsto_model = Channel.value(params.clairsto_model)

        WF_NANOPORE_MINOR_VARIANTS(
            ch_bam,
            ch_reference,
            GET_REF_STATS.out.fai,
            ch_vcf, // major variants from the main pipeline required for deduplication of vcfs
            ch_clairsto_model
        )
        ch_min_vcf = WF_NANOPORE_MINOR_VARIANTS.out.vcf
        ch_versions = ch_versions.mix(WF_NANOPORE_MINOR_VARIANTS.out.versions)
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // SnpEff annotation
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

    ch_snpeff_csv = channel.empty()
    ch_snpeff_db = channel.empty()
    ch_snpeff_config = channel.empty()
    // Use gff if provided:
    ch_gff = params.gff ? file(params.gff, type: 'file', checkIfExists: true) : []

    if (! params.skip_snpeff) {

        // Store original VCF as fallback if SnpEff fails
        // Pipeline robust for typical snpeff issues but annotation files can be problematic
        ch_vcf_original = ch_vcf

        // Get reference ids
        ch_reference
            .map{ ref -> fastaHeaderId(ref) }
            .flatten()
            .collect() // To collect segmented and turn to a value channel
            .set{ ch_ref_ids }

        SNPEFF_DATABASE(
            ch_ref_ids,
            ch_reference,
            ch_gff,
            segmented
        )
        ch_versions = ch_versions.mix(SNPEFF_DATABASE.out.versions)
        ch_snpeff_db = SNPEFF_DATABASE.out.db
        ch_snpeff_config = SNPEFF_DATABASE.out.config

        WF_SNPEFF_ANNOTATE(
            ch_vcf,
            ch_snpeff_db,
            ch_snpeff_config,
            "Major"
        )
        // If SnpEff failed the original vcf will be used in reports, otherwise reports never run if SnpEff fails
        ch_vcf = WF_SNPEFF_ANNOTATE.out.vcf
            .mix(ch_vcf_original)
            .groupTuple()
            .map { meta, vcfs -> [meta, vcfs.find { it.name.endsWith('.ann.vcf.gz') } ?: vcfs[0]] }
        ch_snpeff_csv = WF_SNPEFF_ANNOTATE.out.csv
        ch_versions = ch_versions.mix(WF_SNPEFF_ANNOTATE.out.versions)

        if ( params.minor_variants ) {
            // Store original VCF as fallback if SnpEff fails, which is common
            ch_minvcf_original = ch_min_vcf
            WF_SNPEFF_ANNOTATE_MIN(
                ch_min_vcf,
                ch_snpeff_db,
                ch_snpeff_config,
                "Minor"
            )
            ch_min_vcf = WF_SNPEFF_ANNOTATE_MIN.out.vcf
                .mix(ch_minvcf_original)
                .groupTuple()
                .map { meta, vcfs -> [meta, vcfs.find { it.name.endsWith('.ann.vcf.gz') } ?: vcfs[0]] }
            ch_min_snpeff_csv = WF_SNPEFF_ANNOTATE_MIN.out.csv
        }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Nextclade
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    if ( ! params.skip_nextclade ) {
        WF_NEXTCLADE(
            ch_consensus,
            segmented
        )
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Virus specific tools
    //  More viruses to be added later
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    ch_pangolin_report = channel.empty()
    if ( params.virus_name == 'covid' ) {
        WF_VIRUS_COVID(
            ch_consensus
        )
        ch_pangolin_report = WF_VIRUS_COVID.out.pangolin_report
        ch_versions = ch_versions.mix(WF_VIRUS_COVID.out.versions)
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // QC and Tracking Workflow
    //  This is a stop for segmented viruses at the moment
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    if (! params.skip_qc) {
        //  Filtered out samples - might want to move this
        ch_filter_tracking = channel.empty()
        TRACK_INITIAL_FILTERED_SAMPLES(
            ch_empty_fastqs,
            ch_metadata,
            "TOO FEW INPUT READS"
        )
        ch_filter_tracking = ch_filter_tracking.mix(TRACK_INITIAL_FILTERED_SAMPLES.out.csv)
        TRACK_SIZE_FILTERED_SAMPLES(
            ch_filtered_fastqs_empty,
            ch_metadata,
            "TOO FEW SIZE SELECTED READS"
        )
        ch_filter_tracking = ch_filter_tracking.mix(TRACK_SIZE_FILTERED_SAMPLES.out.csv)

        //  All other samples sequence tracking
        SAMTOOLS_DEPTH(
            ch_bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions)

        // Pass minor vcf to qc or create dummy channel if not running minor_variants
        ch_min_vcf_for_qc = params.minor_variants ? ch_min_vcf
            : ch_consensus.map { meta, _consensus -> [meta, []] }

        MAKE_SAMPLE_QC_CSV(
            ch_consensus
                .join(ch_bam, by: [0])
                .join(SAMTOOLS_DEPTH.out.bed, by: [0])
                .join(ch_vcf, by: [0])
                .join(ch_min_vcf_for_qc, by: [0]),
            ch_primer_bed,
            ch_metadata,
            ch_pcr_primer_bed
        )
        ch_versions = ch_versions.mix(MAKE_SAMPLE_QC_CSV.out.versions)

        //  Combine QCs, check neg controls
        FINAL_QC_CSV(
            MAKE_SAMPLE_QC_CSV.out.csv
                .map{ _meta, csv -> csv }
                .collectFile(keepHeader: true, skip: 1, name: 'concat.qc.csv'),
            ch_filter_tracking
                .collectFile(keepHeader: true, skip: 1, name: 'filter_tracking.csv')
                .ifEmpty([]),
            ch_metadata,
            ch_reference,
            params.neg_control_threshold,
            params.neg_ctrl_substrings
        )
        ch_versions = ch_versions.mix(FINAL_QC_CSV.out.versions)

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        // Final reports workflow
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
        if ( params.multiqc_report ) {
            WF_CREATE_MULTIQC_REPORTS(
                ch_consensus,
                ch_bam,
                ch_vcf,
                MAKE_SAMPLE_QC_CSV.out.csv,
                ch_reads_stats,
                ch_snpeff_csv,
                ch_reference,
                ch_amplicon_bed,
                FINAL_QC_CSV.out.csv,
                ch_pangolin_report,
                ch_versions
            )
        } else {
            WF_CREATE_CUSTOM_REPORT(
                ch_consensus,
                ch_bam,
                ch_vcf,
                ch_reference,
                GET_REF_STATS.out.genome_bed,
                ch_amplicon_bed,
                FINAL_QC_CSV.out.csv,
                ch_versions
            )
        }
    }
}
