process COLLATE_CSVS {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/coreutils:8.31--h14c3975_0'
        : 'biocontainers/coreutils:8.31--h14c3975_0' }"

    input:
    tuple val(meta), path(csvs, stageAs: 'segment_*.csv')

    output:
    tuple val(meta), path("${meta.file_name.replace('.consensus.fasta','')}.csv"), emit: final_csv

    script:
    def filename = "${meta.file_name.replace('.consensus.fasta','')}.csv"
    """
    awk 'FNR==1 && NR!=1 { next } { print }' ${csvs.join(' ')} > ${filename}
    """

    stub:
    def filename = "${meta.file_name.replace('.consensus.fasta','')}.csv"
    """
    touch ${filename}
    """
}
