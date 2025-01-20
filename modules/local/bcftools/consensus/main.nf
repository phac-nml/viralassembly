process BCFTOOLS_CONSENSUS {
    label 'process_single'
    tag "$meta.id"
    publishDir "${params.outdir}/consensus", pattern: "${meta.id}.consensus.fasta", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.19--h8b25389_0' :
        'biocontainers/bcftools:1.19--h8b25389_0' }"

    input:
    tuple val(meta), path(preconsensus), path(coverage_mask), path(pass_vcf), path(pass_vcf_tbi)

    output:
    tuple val(meta), path("${meta.id}.consensus.fasta"), emit: consensus
    path "versions.yml", emit: versions

    script:
    """
    # Command #
    bcftools consensus \\
        -f $preconsensus \\
        ${meta.id}.pass.norm.vcf.gz \\
        -m $coverage_mask \\
        -o ${meta.id}.consensus.fasta

    # Apply samplename as header but keep existing info #
    sed -i "s/>/>$meta.id /" ${meta.id}.consensus.fasta

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.consensus.fasta

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
