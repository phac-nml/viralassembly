process BCFTOOLS_CONSENSUS_FINAL {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(preconsensus), path(coverage_mask), path(pass_vcf), path(pass_vcf_tbi)

    output:
    tuple val(meta), path("${meta.id}-consensus.fasta"), emit: consensus
    path "versions.yml", emit: versions

    script:
    def masking = coverage_mask ? "-m ${coverage_mask}" : ""
    """
    # Command #
    bcftools consensus \\
        -f $preconsensus \\
        ${meta.id}.consensus.norm.vcf.gz \\
        ${masking} \\
        -o ${meta.id}-consensus.fasta

    # Apply samplename as header but keep existing info #
    sed -i "s/>/>$meta.id /" ${meta.id}-consensus.fasta

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}-consensus.fasta

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}

process BCFTOOLS_CONSENSUS_AMBIGUOUS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(vcf), path(tbi), path(fasta), path(mask)

    output:
    tuple val(meta), path('*.fa'), emit: fasta
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def masking = mask ? "-m ${mask}" : ""
    """
    cat ${fasta} \\
        | bcftools \\
            consensus \\
            ${vcf} \\
            ${args} \\
            ${masking} \\
            > ${prefix}.fa
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa
    """
}
