process MINIMAP2_ALIGN {
    label 'process_medium'
    tag "$meta.id"
    publishDir "${params.outdir}/bam", pattern: "${sampleName}.sorted.bam*", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:365b17b986c1a60c1b82c6066a9345f38317b763-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:365b17b986c1a60c1b82c6066a9345f38317b763-0' }"

    input:
    tuple val(meta), path(fastq)
    path reference

    output:
    tuple val(meta), path("${sampleName}.sorted.bam"), path("${sampleName}.sorted.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    sampleName = "$meta.id"
    """
    minimap2 \\
        -a \\
        -x map-ont \\
        -t ${task.cpus} \\
        $reference \\
        $fastq \\
    | samtools view -bS -F 4 - \\
    | samtools sort -o ${sampleName}.sorted.bam

    samtools index ${sampleName}.sorted.bam

    # Versions from nf-core #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
