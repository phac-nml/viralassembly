// Filter vcf to remove failed variants
process PASS_VCF {
    label 'process_single'
    label 'error_retry'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

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
