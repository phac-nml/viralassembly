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
        'https://depot.galaxyproject.org/singularity/artic:1.6.2--pyhdfd78af_0' :
        'biocontainers/artic:1.6.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)
    path scheme

    output:
    tuple val(meta), path("${meta.id}.primertrimmed.rg.sorted.bam"), path("${meta.id}.primertrimmed.rg.sorted.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.pass.vcf.gz"), emit: vcf
    tuple val(meta), path("${meta.id}.consensus.fasta"), emit: consensus
    tuple val(meta), path("${meta.id}.fail.vcf"), emit: fail_vcf
    path "${meta.id}*", emit: all
    path "versions.yml", emit: versions

    script:
    // Clair3 model is added conditonally if it's been set
    // Setup args
    def argsList = []
    if ( params.normalise ) {
        argsList.add("--normalise ${params.normalise}")
    }
    if ( params.no_frameshift ) {
        argsList.add("--no-frameshifts")
    }
    if ( params.clair3_model && params.clair3_model != 'null') {
        argsList.add("--model ${params.clair3_model}")
    }
    def argsConfig = argsList.join(" ")

    // Cmd
    """
    artic minion \\
        ${argsConfig} \\
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
