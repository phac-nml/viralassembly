process BCFTOOLS_CONSENSUS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : task.ext.override_configured_container_registry != false
            ? 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'
            : 'library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(fasta), path(mask)

    output:
    tuple val(meta), path('*.fasta'), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "${meta.id}-${task.ext.prefix}": "${meta.id}"
    def masking = mask ? "-m ${mask}" : ""
    """
    cat ${fasta} \\
        | bcftools \\
            consensus \\
            ${vcf} \\
            ${args} \\
            ${masking} \\
            > ${prefix}.fasta

    # Apply samplename as header but keep existing info #
    sed -i "s/>/>$meta.id /" ${prefix}.fasta

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ? "${meta.id}-${task.ext.prefix}": "${meta.id}"
    """
    touch ${prefix}.fasta

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
