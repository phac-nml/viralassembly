process BEDTOOLS_COVERAGE_GENOME_BED {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path genome_bed

    output:
    tuple val(meta), path("${meta.id}.per_base_coverage.bed"), emit: cov_bed
    path "versions.yml", emit: versions

    script:
    """
    echo -e "reference_name	start	end	position	depth" \\
        > ${meta.id}.per_base_coverage.bed
    bedtools coverage \\
        -d \\
        -a $genome_bed \\
        -b $bam \\
        >> ${meta.id}.per_base_coverage.bed

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools v//' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.per_base_coverage.bed

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools v//' )
    END_VERSIONS
    """
}
process BEDTOOLS_COVERAGE_AMPLICON_BED {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"

    input:
    tuple val(meta), path(bam), path(bai)
    path amplicon_bed

    output:
    tuple val(meta), path("*.amplicon_coverage.bed"), emit: amplicon_coverage
    path "versions.yml", emit: versions

    script:
    """
    echo -e "reference_name	start	end	amplicon_id	pool	strand	read_count	covered_bases	amplicon_length	fraction_covered" \\
        > ${meta.id}.amplicon_coverage.bed
    bedtools coverage \\
        -a $amplicon_bed \\
        -b $bam \\
        -F 0.85 \\
        >> ${meta.id}.amplicon_coverage.bed

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools v//' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.amplicon_coverage.bed

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools v//' )
    END_VERSIONS
    """
}
