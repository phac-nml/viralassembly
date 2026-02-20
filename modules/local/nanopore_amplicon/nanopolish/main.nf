process NANOPOLISH_VARIANTS {
    label 'process_high'
    label 'error_retry'
    tag "${meta.id}-${pool}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanopolish:0.14.0--hd7c1219_0' :
        'biocontainers/nanopolish:0.14.0--hd7c1219_0' }"

    input:
    tuple val(meta), path(fastq), path(bam), path(bai), val(pool)
    path fast5s
    path seqsum
    path reference
    path reference_stats // Example: MN908947.3:1-29904

    output:
    tuple val(meta), path("${meta.id}.${pool}.vcf"), val(pool), emit: vcf
    path "versions.yml", emit: versions

    // Should look into if index can be a separate step per fastq file to speed up time
    //   As then instead of doing it 2x (for 2 pools) it'd only be 1x
    script:
    """
    refstats=\$(cat $reference_stats)
    nanopolish index \\
        -d $fast5s \\
        -s $seqsum \\
        $fastq
    nanopolish variants \\
        --min-flanking-sequence 10 \\
        -x 1000000 \\
        --progress \\
        -t ${task.cpus} \\
        --reads $fastq \\
        -b $bam \\
        -g $reference \\
        -w "\$refstats" \\
        --ploidy 1 \\
        -m 0.15 \\
        --read-group ${pool} \\
        -o ${meta.id}.${pool}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanopolish: \$(echo \$(nanopolish --version | grep nanopolish | sed 's/nanopolish version //'))
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.${pool}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanopolish: \$(echo \$(nanopolish --version | grep nanopolish | sed 's/nanopolish version //'))
    END_VERSIONS
    """
}
