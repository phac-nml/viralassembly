/*
    Custom Utility Modules
        Modules focusing on quick single commands to make intermediate files including:
            * GET_REF_STATS     - Creates intermediate files from the reference
            * RENAME_FASTQ      - Renames barcodeXX fastqs to their sample name
            * SPLIT_BED_BY_POOL - Splits amplicon bed based on the primer pool
            * CREATE_TILING_BED - Creates bed file of the overall tiling region
*/
process GET_REF_STATS {
    label 'process_single'

    conda "bioconda::samtools=1.19.2 bioconda::htslib=1.19.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0':
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    path reference

    output:
    path "${reference}.fai", emit: fai
    path "genome.bed", emit: genome_bed
    path "versions.yml", emit: versions

    script:
    """
    samtools faidx $reference
    cat ${reference}.fai | awk '{ print \$1 "	0	" \$2 }' > genome.bed

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${reference}.fai
    touch genome.bed

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
process RENAME_FASTQ {
    label 'process_single'
    tag "$meta.id"

    conda "conda-forge::python=3.10.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.2' :
        'biocontainers/python:3.10.2' }"

    input:
    tuple val(meta), path(fastq)
    path metadata

    output:
    path "*.fastq", includeInputs: true, emit: fastq
    path "versions.yml", emit: versions

    script:
    """
    rename_fastq.py \\
        --fastq $fastq \\
        --metadata $metadata \\
        --barcode $meta.id

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.fastq

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
process SPLIT_BED_BY_POOL {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/coreutils:8.31--h14c3975_0'
        : 'biocontainers/coreutils:8.31--h14c3975_0' }"

    input:
    path bed

    output:
    path "*.bed", emit: bed

    script:
    """
    awk -F'\t' -v OFS='\t' 'NR>0{print \$1, \$2, \$3, \$4, \$5 > \$5".bed"}' $bed
    """

    stub:
    """
    touch 1.bed
    touch 2.bed
    """
}
process CREATE_TILING_BED {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/coreutils:8.31--h14c3975_0'
        : 'biocontainers/coreutils:8.31--h14c3975_0' }"

    input:
    path amplicon_bed

    output:
    path "tiling_region.bed", emit: bed

    script:
    """
    awk -F '\t' '
    {
        chr = \$1
        if (!(chr in min) || \$2 < min[chr]) min[chr] = \$2
        if (!(chr in max) || \$3 > max[chr]) max[chr] = \$3
    }
    END {
        for (chr in min) {
            print chr "\t" min[chr] "\t" max[chr]
        }
    }
    ' $amplicon_bed > tiling_region.bed
    """

    stub:
    """
    touch tiling_region.bed
    """
}
