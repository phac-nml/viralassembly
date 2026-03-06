process SNPEFF_DATABASE {
    label 'process_small'
    label 'error_ignore' // If can't build we don't run snpeff

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpeff:5.4.0a--hdfd78af_0' :
        'biocontainers/snpeff:5.4.0a--hdfd78af_0' }"

    input:
    val ref_ids
    path reference
    path gff

    output:
    tuple val(genome), path("snpeff_db"), emit: db
    path("snpeff.config"), optional: true, emit: config
    path "versions.yml", emit: versions

    script:
    def avail_mem = 6144
    // Memory
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    // Some setup based on segmented virus compared to non
    def segmented = ref_ids.size() > 1 ? true : false
    def str_ref_ids = ref_ids.join(' ')
    genome = str_ref_ids
    if (segmented) {
        genome = reference.name.split((/\./))[0]
    }

    // Build with gff if that param is given
    if ( gff ) {
        """
        # Setup reference
        mkdir -p snpeff_db/genomes/
        cd snpeff_db/genomes/
        ln -s ../../$reference ${genome}.fa
        cd ../../

        # Setup gff
        mkdir -p snpeff_db/${genome}/
        cd snpeff_db/${genome}/
        ln -s ../../$gff genes.gff
        cd ../../

        # Create config
        echo "${genome}.genome : ${genome}" > snpeff.config

        # Create database
        snpEff \\
            -Xmx${avail_mem}M \\
            build \\
            -config snpeff.config \\
            -dataDir ./snpeff_db \\
            -noCheckCds \\
            -noCheckProtein \\
            -gff3 \\
            -v \\
            ${genome}

        # Versions #
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
        END_VERSIONS
        """
    } else {
        """
        # Check if we can find the reference name in the database but only on non-segmented ones
        if \$(snpEff databases | grep -q "$genome" ) && [ "$segmented" = "false" ]; then
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
            for ref_id in $str_ref_ids ; do
                wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=\${ref_id}&rettype=gbwithparts&retmode=text" -O - >> \$GENE_FILE
            done

            # Create database
            echo "${genome}.genome : ${genome}" > snpeff.config
            snpEff \\
                -Xmx${avail_mem}M \\
                build \\
                -v \\
                -genbank \\
                -config snpeff.config \\
                -dataDir \$PWD/snpeff_db \\
                ${genome}
        fi

        # Versions #
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
        END_VERSIONS
        """
    }

    stub:
    """
    mkdir snpeff_db
    touch snpeff.config

    # Versions #
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
        END_VERSIONS
    """
}
