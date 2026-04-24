/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Utils / Custom checks / Primer Validate
include { GET_REF_STATS             } from '../modules/local/custom/utils.nf'
include { RENAME_FASTQ              } from '../modules/local/custom/utils.nf'
include { PRIMALBEDTOOLS_VALIDATE   } from '../modules/local/primalbedtools/validate/main.nf'
include { PRIMALBEDTOOLS_AMPLICON   } from '../modules/local/primalbedtools/amplicon/main.nf'
include { TRACK_FILTERED_SAMPLES as TRACK_INITIAL_FILTERED_SAMPLES } from '../modules/local/custom/filtering.nf'
include { TRACK_FILTERED_SAMPLES as TRACK_SIZE_FILTERED_SAMPLES    } from '../modules/local/custom/filtering.nf'

// Read QC
include { CHOPPER                   } from '../modules/local/chopper/main'
include { NANOSTAT                  } from '../modules/local/nanostat/main'

// Clair3 model
include { GET_MODEL                 } from '../modules/local/get_model/main'

// Artic related
include { ARTIC_GUPPYPLEX           } from '../modules/local/artic/guppyplex/main'
include { ARTIC_MINION              } from '../modules/local/artic/minion/main'

// QC related
include { SAMTOOLS_DEPTH            } from '../modules/local/samtools/depth/main'
include { MAKE_SAMPLE_QC_CSV        } from '../modules/local/qc/main'
include { FINAL_QC_CSV              } from '../modules/local/qc/main'

// SnpEff related
include { SNPEFF_DATABASE   } from '../modules/local/snpeff/database/main'

// Subworkflows
include { WF_NANOPORE_AMPLICON      } from '../subworkflows/local/nanopore_amplicon'
include { WF_NANOPORE_SHOTGUN       } from '../subworkflows/local/nanopore_shotgun'
include { WF_NANOPORE_SUBCONSENSUS       } from '../subworkflows/local/nanopore_subconsensus'
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
workflow NANOPORE {
    take:
    ch_fastqs           // channel: [ val(meta), file(fastq) ]
    ch_empty_fastqs     // channel: [ val(meta), file(fastq) ]

    main:
    // Optional value channel files from params
    ch_metadata = params.metadata ? file(params.metadata, type: 'file', checkIfExists: true) : []
    ch_pcr_primer_bed = params.pcr_primer_bed ? file(params.pcr_primer_bed, type: 'file', checkIfExists: true) : []

    // Nanopolish required channels, will be ignored when running clair3 or medaka but still passed to the workflow
    ch_fast5s = params.fast5_pass ? file(params.fast5_pass, type: 'dir', checkIfExists: true) : []
    ch_seqsum = params.sequencing_summary ? file(params.sequencing_summary, type: 'file', checkIfExists: true) : []

    // Tool version tracking
    ch_versions = channel.empty()

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Scheme and Reference
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    ch_reference = params.reference     ? channel.value(file(params.reference, type: 'file', checkIfExists: true)) : []
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
    //  FAI, Ref-stats for nanopolish, genome.bed for bedtools
    GET_REF_STATS(
        ch_reference
    )
    ch_versions = ch_versions.mix(GET_REF_STATS.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Models (like me <3)
    //  Medaka just using params for now, should be
    //  in the container(?)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Clair3 model for when running with clair3
    ch_clair3_model = channel.empty()
    if ( params.clair3_local_model ) {
        ch_clair3_model = file(params.clair3_local_model, checkIfExists: true)
    } else {
        GET_MODEL(params.clair3_model)
        ch_clair3_model = GET_MODEL.out.model
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Read QC and Statistics
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    ARTIC_GUPPYPLEX(
        ch_fastqs
    )
    ch_fastqs = ARTIC_GUPPYPLEX.out.fastq
    ch_versions = ch_versions.mix(ARTIC_GUPPYPLEX.out.versions)

    // Rename if given metadata and not using input param, separate out fastqs with little data
    if ( ch_metadata && !params.input ) {
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
        .branch{ _meta, fastq ->
            pass: fastq.countFastq() >= params.min_reads
            empty: fastq.countFastq() < params.min_reads
        }.set{ ch_filtered_fastqs }

    // Stats on final reads
    NANOSTAT(
        ch_filtered_fastqs.pass
    )
    ch_versions = ch_versions.mix(NANOSTAT.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // Chose which pipeline to run based on input params
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    if ( !params.primer_bed ) {
        WF_NANOPORE_SHOTGUN(
            ch_filtered_fastqs.pass,
            ch_fast5s,
            ch_seqsum,
            ch_reference,
            GET_REF_STATS.out.fai,
            GET_REF_STATS.out.refstats,
            ch_clair3_model
        )
        ch_consensus = WF_NANOPORE_SHOTGUN.out.consensus
        ch_bam = WF_NANOPORE_SHOTGUN.out.bam
        ch_vcf = WF_NANOPORE_SHOTGUN.out.vcf
        ch_versions = ch_versions.mix(WF_NANOPORE_SHOTGUN.out.versions)
    } else if ( !params.use_artic_tool ) {
        WF_NANOPORE_AMPLICON(
            ch_filtered_fastqs.pass,
            ch_fast5s,
            ch_seqsum,
            ch_reference,
            GET_REF_STATS.out.fai,
            GET_REF_STATS.out.refstats,
            ch_primer_bed,
            ch_amplicon_bed,
            ch_clair3_model
        )
        ch_consensus = WF_NANOPORE_AMPLICON.out.consensus
        ch_bam = WF_NANOPORE_AMPLICON.out.bam
        ch_vcf = WF_NANOPORE_AMPLICON.out.vcf
        ch_versions = ch_versions.mix(WF_NANOPORE_AMPLICON.out.versions)
    } else {
        ARTIC_MINION(
            ch_filtered_fastqs.pass,
            ch_reference,
            ch_primer_bed,
            ch_clair3_model.ifEmpty([])
        )
        ch_consensus = ARTIC_MINION.out.consensus
        ch_bam = ARTIC_MINION.out.bam
        ch_vcf = ARTIC_MINION.out.vcf
        ch_versions = ch_versions.mix(ARTIC_MINION.out.versions)
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // subconsensus variants, optional
    //  Run before SnpEff to use unannotated VCF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    ch_min_vcf = channel.empty()
    if ( params.subconsensus ) {
        WF_NANOPORE_SUBCONSENSUS(
            ch_bam,
            ch_reference,
            GET_REF_STATS.out.fai,
            ch_vcf // major variants from the main pipeline
        )
        ch_min_vcf = WF_NANOPORE_SUBCONSENSUS.out.vcf
        ch_versions = ch_versions.mix(WF_NANOPORE_SUBCONSENSUS.out.versions)
    }        
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    // SnpEff annotation
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //    

    ch_snpeff_csv = channel.empty() 
    ch_snpeff_db = channel.empty()
    ch_snpeff_config = channel.empty()

    if (! params.skip_snpeff) {

        // Store original VCF as fallback if SnpEff fails, which is common
        ch_vcf_original = ch_vcf

        // Get reference id
        ch_reference.splitFasta( record: [ id: true ] )
        .map{ record -> record.id.toString() }
        .collect() // To collect segmented and turn to a value channel
        .set{ ch_ref_ids }

        // Use gff if provided:
        ch_gff = params.gff ? file(params.gff, type: 'file', checkIfExists: true) : []

        SNPEFF_DATABASE(
            ch_ref_ids,
            ch_reference,
            ch_gff
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
        // If SnpEff failed the original vcf will be used in reports, otherwise reports never run if SnpEff fails (common issue)
        ch_vcf = WF_SNPEFF_ANNOTATE.out.vcf
            .mix(ch_vcf_original)
            .groupTuple()
            .map { meta, vcfs -> [meta, vcfs.find { it.name.endsWith('.ann.vcf.gz') } ?: vcfs[0]] }
        ch_snpeff_csv = WF_SNPEFF_ANNOTATE.out.csv
        ch_versions = ch_versions.mix(WF_SNPEFF_ANNOTATE.out.versions)
    
        if ( params.subconsensus ) {
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
    if ( ! params.skip_nextclade && (params.nextclade_dataset_name || params.nextclade_dataset_dir) ) {
        WF_NEXTCLADE(
            ch_consensus
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
            ch_filtered_fastqs.empty,
            ch_metadata,
            "TOO FEW SIZE SELECTED READS"
        )
        ch_filter_tracking = ch_filter_tracking.mix(TRACK_SIZE_FILTERED_SAMPLES.out.csv)

        //  All other samples sequence tracking
        SAMTOOLS_DEPTH(
            ch_bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions)

        // Pass minor vcf to qc or create dummy channel if not running subconsensus
        ch_min_vcf_for_qc = params.subconsensus ? ch_min_vcf
            : ch_consensus.map { meta, consensus -> [meta, []] }

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
                NANOSTAT.out.stats,
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
