/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Utils / Custom checks
include { DOWNLOAD_SCHEME           } from '../modules/local/custom/utils.nf'
include { SIMPLE_SCHEME_VALIDATE    } from '../modules/local/custom/utils.nf'
include { GET_REF_STATS             } from '../modules/local/custom/utils.nf'
include { CREATE_AMPLICON_BED       } from '../modules/local/custom/utils.nf'
include { RENAME_FASTQ              } from '../modules/local/custom/utils.nf'
include { TRACK_FILTERED_SAMPLES as TRACK_INITIAL_FILTERED_SAMPLES } from '../modules/local/custom/filtering.nf'
include { TRACK_FILTERED_SAMPLES as TRACK_SIZE_FILTERED_SAMPLES    } from '../modules/local/custom/filtering.nf'

// Read QC
include { CHOPPER                   } from '../modules/local/chopper/main'
include { NANOSTAT                  } from '../modules/local/nanostat/main'

// Artic related
include { ARTIC_GUPPYPLEX           } from '../modules/local/artic/guppyplex/main'
include { ARTIC_MINION              } from '../modules/local/artic/minion/main'

// QC related
include { SAMTOOLS_DEPTH            } from '../modules/local/samtools/depth/main'
include { MAKE_SAMPLE_QC_CSV        } from '../modules/local/qc/main'
include { FINAL_QC_CSV              } from '../modules/local/qc/main'

// Subworkflows
include { WF_NANOPORE_AMPLICON      } from '../subworkflows/local/nanopore_amplicon.nf'
include { WF_NANOPORE_SHOTGUN       } from '../subworkflows/local/nanopore_shotgun.nf'
include { WF_SNPEFF_ANNOTATE        } from '../subworkflows/local/snpeff_annotate.nf'
include { WF_CREATE_MULTIQC_REPORTS } from '../subworkflows/local/create_multiqc_reports.nf'
include { WF_CREATE_CUSTOM_REPORT   } from '../subworkflows/local/create_custom_report.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INITIALIZE CHANNELS FROM PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Optional value channel files from params
ch_metadata = params.metadata ? file(params.metadata, type: 'file', checkIfExists: true) : []
ch_local_scheme = params.local_scheme ? file(params.local_scheme, type: 'dir', checkIfExists: true) : []
ch_pcr_primer_bed = params.pcr_primer_bed ? file(params.pcr_primer_bed, type: 'file', checkIfExists: true) : []

// Nanopolish required channels, will be ignored when running medaka but still passed to the process
ch_fast5s = params.fast5_pass ? file(params.fast5_pass, type: 'dir', checkIfExists: true) : []
ch_seqSum = params.sequencing_summary ? file(params.sequencing_summary, type: 'file', checkIfExists: true) : []

// Reference for if not using a scheme
ch_reference = params.reference ? Channel.value(file(params.reference, type: 'file', checkIfExists: true)) : []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow NANOPORE {
    take:
    ch_fastqs           // channel: [ val(meta), file(fastq) ]
    ch_empty_fastqs     // channel: [ val(meta), file(fastq) ]

    main:
    // Tool version tracking
    ch_versions = Channel.empty()

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Scheme and Reference
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    ch_amplicon_bed = Channel.empty()
    ch_primer_bed = Channel.value([]) // This has to be a value channel for qc creation to work
    if ( ! params.reference ) {
        if ( ! ch_local_scheme ) {
            DOWNLOAD_SCHEME()
            ch_local_scheme = DOWNLOAD_SCHEME.out.scheme
        }
        SIMPLE_SCHEME_VALIDATE(
            ch_local_scheme
        )
        ch_scheme = SIMPLE_SCHEME_VALIDATE.out.scheme
        ch_primer_bed = SIMPLE_SCHEME_VALIDATE.out.bed
        // Overwrite reference if one was given along with a scheme
        ch_reference = SIMPLE_SCHEME_VALIDATE.out.ref

        // Amplicon information
        CREATE_AMPLICON_BED(
            ch_primer_bed
        )
        ch_amplicon_bed = CREATE_AMPLICON_BED.out.amplicon_bed
        ch_versions = ch_versions.mix(CREATE_AMPLICON_BED.out.versions)
    }

    // Reference stats and files for various processes
    //  FAI, Ref-stats for nanopolish, genome.bed for bedtools
    GET_REF_STATS(
        ch_reference
    )
    ch_versions = ch_versions.mix(GET_REF_STATS.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Read QC and Statistics
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    ARTIC_GUPPYPLEX(
        ch_fastqs
    )
    ch_fastqs = ARTIC_GUPPYPLEX.out.fastq
    ch_versions = ch_versions.mix(ARTIC_GUPPYPLEX.out.versions)

    // Rename if given metadata and not using input param, separate out fastqs with little data
    if ( ch_metadata && ! params.input ) {
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

    // Chopper may be useless as we already filter based on length earlier
    //  But it also does add quality filtering
    CHOPPER(
        ch_fastqs
    )
    ch_versions = ch_versions.mix(CHOPPER.out.versions)
    ch_fastqs = CHOPPER.out.fastq

    // Pass/fail reads based on count after length and quality filtering
    ch_fastqs
        .branch{
            pass: it[1].countFastq() >= params.min_reads
            empty: it[1].countFastq() < params.min_reads
        }.set{ ch_filtered_fastqs }

    // Stats on final reads
    NANOSTAT(
        ch_filtered_fastqs.pass
    )
    ch_versions = ch_versions.mix(NANOSTAT.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Chose which pipeline to run based on input params
    //  The "proper" artic minion pipeline or re-implemented nextflow version
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    if ( params.reference ) {
        WF_NANOPORE_SHOTGUN(
            ch_filtered_fastqs.pass,
            ch_fast5s,
            ch_seqSum,
            ch_reference,
            GET_REF_STATS.out.fai,
            GET_REF_STATS.out.refstats
        )
        ch_consensus = WF_NANOPORE_SHOTGUN.out.consensus
        ch_bam = WF_NANOPORE_SHOTGUN.out.bam
        ch_vcf = WF_NANOPORE_SHOTGUN.out.vcf
        ch_versions = ch_versions.mix(WF_NANOPORE_SHOTGUN.out.versions)
    } else if ( ! params.use_artic_tool ) {
        WF_NANOPORE_AMPLICON(
            ch_filtered_fastqs.pass,
            ch_fast5s,
            ch_seqSum,
            ch_scheme,
            ch_reference,
            GET_REF_STATS.out.fai,
            GET_REF_STATS.out.refstats,
            ch_primer_bed,
            ch_amplicon_bed,
            CREATE_AMPLICON_BED.out.tiling_bed
        )
        ch_consensus = WF_NANOPORE_AMPLICON.out.consensus
        ch_bam = WF_NANOPORE_AMPLICON.out.bam
        ch_vcf = WF_NANOPORE_AMPLICON.out.vcf
        ch_versions = ch_versions.mix(WF_NANOPORE_AMPLICON.out.versions)
    } else {
        ARTIC_MINION(
            ch_filtered_fastqs.pass,
            ch_fast5s,
            ch_seqSum,
            ch_scheme
        )
        ch_consensus = ARTIC_MINION.out.consensus
        ch_bam = ARTIC_MINION.out.bam
        ch_vcf = ARTIC_MINION.out.vcf
        ch_versions = ch_versions.mix(ARTIC_MINION.out.versions)
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // SnpEff annotation
    //  Only run if we have one reference sequence for now
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    ch_snpeff_csv = Channel.empty()
    if ( (! params.skip_snpeff) || (ch_reference.countFasta() == 1) ) {
        WF_SNPEFF_ANNOTATE(
            ch_vcf,
            ch_reference
        )
        ch_vcf = WF_SNPEFF_ANNOTATE.out.vcf
        ch_snpeff_csv = WF_SNPEFF_ANNOTATE.out.csv
        ch_versions = ch_versions.mix(WF_SNPEFF_ANNOTATE.out.versions)
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // QC and Tracking Workflow
    //  This is a stop for segmented viruses at the moment
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    if ( (! params.skip_qc) || (ch_reference.countFasta() == 1) ) {
        //  Filtered out samples - might want to move this
        ch_filter_tracking = Channel.empty()
        TRACK_INITIAL_FILTERED_SAMPLES(
            ch_empty_fastqs,
            ch_metadata,
            "Too few found fastqs"
        )
        ch_filter_tracking = ch_filter_tracking.mix(TRACK_INITIAL_FILTERED_SAMPLES.out.csv)
        TRACK_SIZE_FILTERED_SAMPLES(
            ch_filtered_fastqs.empty,
            ch_metadata,
            "Too few size selected reads"
        )
        ch_filter_tracking = ch_filter_tracking.mix(TRACK_SIZE_FILTERED_SAMPLES.out.csv)

        //  All other samples sequence tracking
        SAMTOOLS_DEPTH(
            ch_bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions)

        MAKE_SAMPLE_QC_CSV(
            ch_consensus
                .join(ch_bam, by: [0])
                .join(SAMTOOLS_DEPTH.out.bed, by: [0])
                .join(ch_vcf, by: [0]),
            ch_primer_bed,
            ch_metadata,
            ch_pcr_primer_bed
        )
        ch_versions = ch_versions.mix(MAKE_SAMPLE_QC_CSV.out.versions)

        //  Combine QCs, check neg controls
        FINAL_QC_CSV(
            MAKE_SAMPLE_QC_CSV.out.csv
                .map{ it -> it[1] }
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
        if ( params.custom_report ) {
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
        } else {
            WF_CREATE_MULTIQC_REPORTS(
                ch_consensus,
                ch_bam,
                ch_vcf,
                MAKE_SAMPLE_QC_CSV.out.csv,
                NANOSTAT.out.stats,
                ch_snpeff_csv,
                ch_reference,
                ch_amplicon_bed,
                FINAL_QC_CSV.out.csv,
                ch_versions
            )
        }
    }
}
