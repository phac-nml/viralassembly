// Filter vcf to remove failed variants
process PASS_VCF {
    label 'process_single'
    label 'error_retry'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.19--h8b25389_0' :
        'biocontainers/bcftools:1.19--h8b25389_0' }"

    input:
    tuple val(meta), path(full_vcf)

    output:
    tuple val(meta), path("${meta.id}-full.vcf.gz"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    bcftools view --output-type z -f PASS $full_vcf > ${meta.id}-full.vcf.gz

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
