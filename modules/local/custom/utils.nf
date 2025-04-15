// Custom Utility Modules
process CAT_FASTQ {
    label 'process_single'
    tag "$meta.id"

    conda "conda-forge::pigz=2.3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(gzipped_reads), path(reads)

    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: fastq

    script:
    // Check if input lists are empty or not
    def gzReadsExist = !gzipped_reads.empty
    def readsExist = !reads.empty

    def outName = "${meta.id}.merged.fastq.gz"
    """
    touch $outName
    if [ "$gzReadsExist" == "true" ]; then
        cat $gzipped_reads >> $outName
    fi
    if [ "$readsExist" == "true" ]; then
        cat $reads | pigz -ck >> $outName
    fi
    """

    stub:
    """
    touch ${meta.id}.merged.fastq.gz
    """
}
process DOWNLOAD_SCHEME {
    label 'process_single'
    tag { params.scheme_repo }
    publishDir "${params.outdir}/downloaded_scheme", pattern: "primer-schemes", mode: "copy"

    output:
    path "primer-schemes", emit: scheme

    script:
    """
    git clone ${params.scheme_repo} primer-schemes
    """

    stub:
    """
    mkdir -p primer-schemes/stub/V1
    touch primer-schemes/stub/V1/stub.reference.fasta
    touch primer-schemes/stub/V1/stub.scheme.bed
    """
}
process SIMPLE_SCHEME_VALIDATE {
    label 'process_single'

    input:
    path scheme

    output:
    path("primer-schemes/${params.scheme}/${params.scheme_version}/*reference.fasta"), emit: ref
    path("primer-schemes/${params.scheme}/${params.scheme_version}/*scheme.bed"), emit: bed
    path "primer-schemes", emit: scheme

    // No clue if this is the best way to validate but eh for now it works
    script:
    """
    if [[ "$scheme" != "primer-schemes" ]]; then
        mv $scheme primer-schemes
    fi
    if [ ! -f primer-schemes/${params.scheme}/${params.scheme_version}/*reference.fasta ]; then
        echo "ERROR: Reference Fasta not found in 'primer-schemes/${params.scheme}/${params.scheme_version}/*reference.fasta'"
        exit 1
    elif [ ! -f primer-schemes/${params.scheme}/${params.scheme_version}/*scheme.bed ]; then
        echo "ERROR: Scheme bed file not found in 'primer-schemes/${params.scheme}/${params.scheme_version}/*scheme.bed'"
        exit 1
    fi
    """
}
process GET_REF_STATS {
    label 'process_single'
    publishDir "${params.outdir}/reference", pattern: "${reference}*", mode: "copy"
    publishDir "${params.outdir}/reference", pattern: "refstats.txt", mode: "copy"
    publishDir "${params.outdir}/reference", pattern: "genome.bed", mode: "copy"

    conda "bioconda::samtools=1.19.2 bioconda::htslib=1.19.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0':
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    path reference

    output:
    path "${reference}.fai", emit: fai
    path "refstats.txt", emit: refstats
    path "genome.bed", emit: genome_bed
    path "versions.yml", emit: versions

    script:
    """
    samtools faidx $reference
    cat ${reference}.fai | awk '{print \$1 ":1-" \$2+1}' > refstats.txt
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
    touch refstats.txt
    touch genome.bed

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
process CREATE_AMPLICON_BED {
    label 'process_single'
    publishDir "${params.outdir}/bed", pattern: "amplicon.bed", mode: "copy"
    publishDir "${params.outdir}/bed", pattern: "tiling_region.bed", mode: "copy"

    conda "conda-forge::python=3.10.2"
    container "quay.io/biocontainers/python:3.10.2"

    input:
    path bed

    output:
    path "amplicon.bed", emit: amplicon_bed
    path "tiling_region.bed", emit: tiling_bed
    path "versions.yml", emit: versions

    script:
    """
    primers_to_amplicons.py \\
        --bed $bed

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch amplicon.bed
    touch tiling_region.bed

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
process RENAME_FASTQ {
    label 'process_single'
    tag "$meta.id"

    conda "conda-forge::python=3.10.2"
    container "quay.io/biocontainers/python:3.10.2"

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
    publishDir "${params.outdir}/bed", pattern: "*.split.bed", mode: "copy"

    input:
    path bed

    output:
    path "*.split.bed", emit: bed

    script:
    """
    awk -F'\t' -v OFS='\t' 'NR>0{print \$1, \$2, \$3, \$4, \$5, \$6 > \$5".split.bed"}' $bed
    """

    stub:
    """
    touch 1.split.bed
    touch 2.split.bed
    """
}
