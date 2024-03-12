process SNPEFF_DATABASE {
    label 'process_medium'
    label 'error_ignore' // If can't build we don't run snpeff
    publishDir "${params.outdir}/snpeff/database", pattern: "snpeff_db", mode: "copy"
    publishDir "${params.outdir}/snpeff/database", pattern: "snpEff.config", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpeff:5.2--hdfd78af_0' :
        'biocontainers/snpeff:5.2--hdfd78af_0' }"

    input:
    val genome
    path reference
    path gff

    output:
    path("snpeff_db"), emit: db
    path('snpEff.config'), optional: true, emit: config
    path "versions.yml", emit: versions

    script:
    def avail_mem = 6144
    // Memory
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    // Build with gff if that param is given
    if ( gff ) {
        """
        # Setup reference
        mkdir -p snpeff_db/genomes/
        cd snpeff_db/genomes/
        ln -s ../../$reference ${genome}.fa
        cd ../..

        # Setup gff
        mkdir -p snpeff_db/${genome}
        cd snpeff_db/${genome}
        ln -s ../../$gff genes.gff
        cd ../..

        # Create database
        echo "${genome}.genome : ${genome}" > snpEff.config
        snpEff \\
            -Xmx${avail_mem}M \\
            build \\
            -v \\
            -config snpeff.config \\
            -dataDir ./snpeff_db \\
            -gff3 \\
            ${genome}
        """
    } else {
        """
        # Check if we can find the reference name in the database
        if \$(snpEff databases | grep -q "$genome" ); then 
            echo "Found $genome in snpEff database"
            snpEff \\
                -Xmx${avail_mem}M \\
                download ${genome} \\
                -dataDir \${PWD}/snpeff_db

        # Otherwise try to make it from NCBI
        else
            # Set pathes
            echo "Attempting to make snpEff database for $genome from NCBI genbank file"
            DIR="snpeff_db/${genome}"
            GENE_FILE="\$DIR/genes.gbk"

            # Try to get gbk file
            mkdir -p "\$DIR"
            wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${genome}&rettype=gbwithparts&retmode=text" -O \$GENE_FILE

            # Create database
            echo "${genome}.genome : ${genome}" > snpEff.config
            snpEff \\
                -Xmx${avail_mem}M \\
                build \\
                -v \\
                -genbank \\
                -config snpEff.config \\
                -dataDir \$PWD/snpeff_db \\
                ${genome}
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
        END_VERSIONS
        """
    }
}
process SNPEFF_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/snpeff", pattern: "*.vcf", mode: "copy"
    publishDir "${params.outdir}/snpeff", pattern: "*.csv", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpeff:5.2--hdfd78af_0' :
        'biocontainers/snpeff:5.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)
    val genome
    path snpeff_db
    path config

    output:
    tuple val(meta), path("*.ann.vcf"), emit: vcf
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Memory
    def avail_mem = 6144
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    sampleName = "$meta.id"
    // Args for db and config
    def snpeff_db_command = snpeff_db ? "-dataDir \${PWD}/${snpeff_db}" : ""
    def config_command = config ? "-config \${PWD}/${config}" : ""
    """
    # Sporatic lock issue in tmp dir solution
    #  Partially from https://github.com/apache/arrow/pull/39115/files
    mkdir -p tmp
    export _JAVA_OPTIONS="-Djava.io.tmpdir=\$PWD/tmp -XX:-UsePerfData"

    # Run command
    snpEff \\
        -Xmx${avail_mem}M \\
        -csvStats ${sampleName}.csv \\
        $snpeff_db_command \\
        $config_command \\
        -no-intergenic \\
        -no-intron \\
        -hgvs1LetterAa \\
        ${genome} \\
        $vcf \\
        > ${sampleName}.ann.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
