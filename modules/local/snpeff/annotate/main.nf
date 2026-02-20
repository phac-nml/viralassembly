process SNPEFF_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpeff:5.4.0a--hdfd78af_0' :
        'biocontainers/snpeff:5.4.0a--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)
    val genome
    path snpeff_db
    path config

    output:
    tuple val(meta), path("*.ann.vcf"), emit: vcf
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml", emit: versions

    script:
    // Memory
    def avail_mem = 6144
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    // Args for db and config
    def snpeff_db_command = snpeff_db ? "-dataDir ${snpeff_db}" : ""
    def config_command = config ? "-config ${config}" : ""
    """
    # Sporatic lock issue in tmp dir solution
    #  Partially from https://github.com/apache/arrow/pull/39115/files
    mkdir -p tmp
    export _JAVA_OPTIONS="-Djava.io.tmpdir=\$PWD/tmp -XX:-UsePerfData"

    # Run command
    snpEff \\
        -Xmx${avail_mem}M \\
        -csvStats ${meta.id}.csv \\
        $snpeff_db_command \\
        $config_command \\
        -no-intergenic \\
        -no-intron \\
        -hgvs1LetterAa \\
        ${genome} \\
        $vcf \\
        > ${meta.id}.ann.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.ann.vcf
    touch ${meta.id}.csv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
