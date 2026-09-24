// Combine two vcfs into one, especially useful for callers that output indels and snps separately
process CAT_VCF {
    label 'process_single'
    label 'error_retry'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.19--h8b25389_0' :
        'biocontainers/bcftools:1.19--h8b25389_0' }"

    input:
    tuple val(meta), path(snv_vcf), path(indel_vcf)

    output:
    tuple val(meta), path("${meta.id}-cat.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    real_snv=\$(readlink -f $snv_vcf)
    real_indel=\$(readlink -f $indel_vcf)

    # Need to use absolute paths in bcftools command because links were causing errors.
    bcftools concat \\
        -a -o ${meta.id}-cat.vcf \\
        \$real_snv \$real_indel

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
