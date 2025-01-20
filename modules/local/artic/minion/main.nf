process ARTIC_MINION {
    label 'process_high'
    label 'error_retry'
    tag "$meta.id"
    publishDir "${params.outdir}/consensus", pattern: "${meta.id}.consensus.fasta", mode: "copy"
    publishDir "${params.outdir}/bam", pattern: "${meta.id}.*bam*", mode: "copy"
    publishDir "${params.outdir}/vcf", pattern: "${meta.id}.pass.vcf*", mode: "copy"
    publishDir "${params.outdir}/vcf", pattern: "${meta.id}.fail.vcf", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.2.4--pyh7cba7a3_1' :
        'biocontainers/artic:1.2.4--pyh7cba7a3_1' }"

    input:
    tuple val(meta), path(fastq)
    path fast5_dir
    path sequencing_summary
    path scheme

    output:
    tuple val(meta), path("${meta.id}.primertrimmed.rg.sorted.bam"), path("${meta.id}.primertrimmed.rg.sorted.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.pass.vcf.gz"), emit: vcf
    tuple val(meta), path("${meta.id}.consensus.fasta"), emit: consensus
    tuple val(meta), path("${meta.id}.fail.vcf"), emit: fail_vcf
    path "${meta.id}*", emit: all
    path "versions.yml", emit: versions

    script:
    // Setup args for medaka vs nanopolish
    def argsList = []
    if ( params.variant_caller == 'medaka' ) {
        argsList.add("--medaka")
        argsList.add("--medaka-model ${params.medaka_model}")
    } else {
        argsList.add("--fast5-directory $fast5_dir")
        argsList.add("--sequencing-summary $sequencing_summary")
    }
    if ( params.normalise ) {
        argsList.add("--normalise ${params.normalise}")
    }
    if ( params.no_frameshift ) {
        argsList.add("--no-frameshifts")
    }
    def argsConfig = argsList.join(" ")

    // Aligner
    def alignerArg = "--minimap2"
    if ( params.use_bwa ) {
        alignerArg = "--bwa"
    }

    // Cmd
    """
    artic minion \\
        ${argsConfig} \\
        ${alignerArg} \\
        --threads ${task.cpus} \\
        --read-file $fastq \\
        --scheme-version ${params.scheme_version} \\
        ${params.scheme} \\
        ${meta.id}

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.primertrimmed.rg.sorted.bam
    touch ${meta.id}.primertrimmed.rg.sorted.bam.bai
    touch ${meta.id}.pass.vcf.gz
    touch ${meta.id}.consensus.fasta

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
