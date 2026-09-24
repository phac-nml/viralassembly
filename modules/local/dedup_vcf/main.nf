// compares minor variants with consensus vcf and remove duplicates to output only minor variant vcf
process DEDUP_VCFS {
    label 'process_single'
    label 'error_retry'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.11.1--pyhdfd78af_0' :
        'biocontainers/artic:1.11.1--pyhdfd78af_0' }"

    input:
    tuple val(meta),
        path(cat_vcf),
        path(cat_tbi),
        path(pass_vcf)

    output:
    // dedup_vcfs/0000.vcf - records private to sample-cat.vcf.gz (unique to ClairS-TO)
    tuple val(meta), path("dedup_vcfs/0000.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    real_cat=\$(readlink -f $cat_vcf)
    real_pass=\$(readlink -f $pass_vcf)

    #Need to use absolute paths in bcftools command
    dedup_vcfs.py \\
        --consensus-vcf \$real_pass \\
        --clairSTO-vcf \$real_cat

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
