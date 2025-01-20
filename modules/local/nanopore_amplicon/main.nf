/*
    Variant callers for nanopore amplicon data
        The parameters match those for shotgun data
        The main difference is that the amplicon processes
        are heavily coupled to use of the readgroup/bed file/regions
        while the shotgun processes are not
*/
process MEDAKA_CONSENSUS {
    label 'process_medium'
    label 'error_retry'
    tag "${meta.id}-${pool}"

    conda "${moduleDir}/env-medaka.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.11.3--py39h05d5c5e_0' :
        'biocontainers/medaka:1.11.3--py39h05d5c5e_0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(pool)

    output:
    tuple val(meta), path("${meta.id}.${pool}.hdf"), val(pool), emit: hdf
    path "versions.yml", emit: versions

    script:
    """
    medaka consensus \\
        --model ${params.medaka_model} \\
        --threads ${task.cpus} \\
        --chunk_len 800 \\
        --chunk_ovlp 400 \\
        --RG ${pool} \\
        $bam \\
        ${meta.id}.${pool}.hdf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.${pool}.hdf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
process MEDAKA_VARIANT {
    label 'process_medium'
    tag "${meta.id}-${pool}"
    // publishDir "${params.outdir}/articMinionNextflow", pattern: "${meta.id}.${pool}.vcf", mode: "copy"

    conda "${moduleDir}/env-medaka.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.11.3--py39h05d5c5e_0' :
        'biocontainers/medaka:1.11.3--py39h05d5c5e_0' }"

    input:
    tuple val(meta), path(hdf), val(pool)
    path reference

    output:
    tuple val(meta), path("${meta.id}.${pool}.vcf"), val(pool), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    medaka variant \\
        $reference \\
        $hdf \\
        ${meta.id}.${pool}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.${pool}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
process NANOPOLISH_VARIANTS {
    label 'process_high'
    label 'error_retry'
    tag "${meta.id}-${pool}"
    // publishDir "${params.outdir}/articMinionNextflow", pattern: "${meta.id}.${pool}.vcf", mode: "copy"

    conda "${moduleDir}/env-nanopolish.yml"
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
process CLAIR3_VARIANTS {
    label 'process_high'
    label 'error_retry'
    tag "${meta.id}-${pool}"
    // publishDir "${params.outdir}/articMinionNextflow", pattern: "${meta.id}.${pool}.vcf", mode: "copy"

    conda "${moduleDir}/env-clair3.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:1.0.6--py39hf5e1c6e_0' :
        'biocontainers/clair3:1.0.6--py39hf5e1c6e_0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(pool), path(pool_bed)
    path reference
    path fai
    // optional model_path
    path model_path

    output:
    tuple val(meta), path("${meta.id}.${pool}.vcf"), val(pool), emit: vcf
    path "versions.yml", emit: versions

    script:
    // Using some of the nf-flu work to get clair3 working
    model_suffix = "models/${params.clair3_model}"
    using_conda = (workflow.containerEngine == null || workflow.containerEngine == '')
    """
    CLAIR_BIN_DIR=\$(dirname \$(which run_clair3.sh))
    if [[ "${params.clair3_user_variant_model}" != "" ]] ; then
        MODEL_PATH=${model_path}
    else
        if [[ ${using_conda} = true ]] ; then
            MODEL_PATH="\$CLAIR_BIN_DIR/${model_suffix}"
        else
            MODEL_PATH="/usr/local/bin/models/${params.clair3_model}"
        fi
    fi

    run_clair3.sh \
        --bam_fn=$bam \\
        --bed_fn=$pool_bed \\
        --ref_fn=$reference \\
        --threads=${task.cpus} \\
        --platform='ont' \\
        --model_path="\$MODEL_PATH" \\
        --output="${meta.id}-out" \\
        --min_coverage=10 \\
        --haploid_precise \\
        --enable_long_indel \\
        --fast_mode \\
        --include_all_ctgs \\
        --no_phasing_for_fa

    gunzip ${meta.id}-out/merge_output.vcf.gz
    ln -s ${meta.id}-out/merge_output.vcf ${meta.id}.${pool}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(run_clair3.sh -v) | sed 's/Clair3 //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.${pool}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(run_clair3.sh -v) | sed 's/Clair3 //')
    END_VERSIONS
    """
}
