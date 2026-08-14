process ARTIC_VCF_FILTER {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.11.1--pyhdfd78af_0' :
        'biocontainers/artic:1.11.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.pass.vcf.gz"), path("${meta.id}.pass.vcf.gz.tbi"), emit: pass_vcf
    tuple val(meta), path("${meta.id}.fail.vcf"), emit: fail_vcf
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    artic_vcf_filter \\
        $args \\
        $vcf \\
        ${meta.id}.pass.vcf \\
        ${meta.id}.fail.vcf \\
        ${meta.id}.ignore.vcf
    bgzip -f ${meta.id}.pass.vcf
    tabix -p vcf ${meta.id}.pass.vcf.gz

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.pass.vcf.gz
    touch ${meta.id}.pass.vcf.gz.tbi
    touch ${meta.id}.fail.vcf
    touch ${meta.id}.ignore.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
