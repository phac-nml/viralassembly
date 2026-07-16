//Taken from https://github.com/nf-core/viralrecon/blob/master/modules/local/make_bed_mask.nf
process MAKE_BED_MASK {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-1a35167f7a491c7086c13835aaa74b39f1f43979:6b5cffa1187cfccf2dc983ed3b5359d49b999eb0-0'
        : 'biocontainers/mulled-v2-1a35167f7a491c7086c13835aaa74b39f1f43979:6b5cffa1187cfccf2dc983ed3b5359d49b999eb0-0'}"

    input:
    tuple val(meta), path(bam), path(bai), path(vcf)
    path reference

    output:
    tuple val(meta), path("*.bed")    , emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: 5
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        mpileup \\
        ${args} \\
        --reference ${reference} \\
        ${bam} \\
        | awk -v OFS='\\t' '{print \$1, \$2-1, \$2, \$4}' | awk '\$4 < ${args2}' > lowcov_positions.txt

    make_bed_mask.py \\
        ${vcf} \\
        lowcov_positions.txt \\
        ${prefix}.bed

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
