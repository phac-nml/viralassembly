// Visualization Modules
//  Custom scripts are versioned here
process CREATE_READ_VARIATION_CSV {
    label 'process_high_memory'
    tag "$meta.id"
    publishDir "${params.outdir}/variation_csvs", pattern: "*.csv", mode: "copy"

    conda "${moduleDir}/env-artic.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.6.2--pyhdfd78af_0' :
        'biocontainers/artic:1.6.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path reference

    output:
    tuple val(meta), path("${meta.id}_variation.csv"), optional: true, emit: csv
    path "versions.yml", emit: versions

    script:
    """
    calc_bam_variation.py \\
        --bam $bam \\
        --reference $reference \\
        --sample $meta.id

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        calc_bam_variation.py: 0.1.0
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_variation.csv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        calc_bam_variation.py: 0.1.0
    END_VERSIONS
    """
}
process CREATE_VARIANT_TSV {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/env-artic.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.6.2--pyhdfd78af_0' :
        'biocontainers/artic:1.6.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.vcf.tsv"), optional: true, emit: tsv
    path "versions.yml", emit: versions

    script:
    def annotated_arg = params.skip_snpeff ? "" : "--annotated"
    """
    vcf_to_tsv.py \\
        --sample $meta.id \\
        $annotated_arg \\
        --vcf $vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        vcf_to_tsv.py: 0.1.0
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.vcf.tsv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        vcf_to_tsv.py: 0.1.0
    END_VERSIONS
    """
}
process COMBINE_AMPLICON_COVERAGE {
    label 'process_single'

    conda "${moduleDir}/env-pandas.yml"
    container "quay.io/biocontainers/pandas:1.5.2"

    input:
    path(amplicon_beds)

    output:
    path "merged_amplicon_depth.csv", emit: amplicons

    script:
    """
    #!/usr/bin/env python3
    '''
    Simple script to merge together and transpose all of the amplicon coverage files for reporting later
    '''

    import glob
    import pandas as pd
    from pathlib import Path

    # Find files and take only the two needed columns
    tsv_files = glob.glob('*.amplicon_coverage.bed')
    df_list = []
    for f in tsv_files:
        f = Path(f)
        name = f.name.split('.amplicon_coverage.bed')[0]
        df = pd.read_csv(f, sep='\t')
        df = df[['amplicon_id', 'read_count']]
        df.rename(columns={'read_count': name}, inplace=True)
        df_list.append(df)

    # Transpose and Output
    df = pd.DataFrame().join([d.set_index('amplicon_id') for d in df_list], how='outer').dropna().reset_index()
    df.set_index('amplicon_id', inplace=True, drop=True)
    df = df.transpose()
    df.index.name = 'sample'
    df.sort_index(axis=0, inplace=True)
    df.to_csv('merged_amplicon_depth.csv', index=True)
    """

    stub:
    """
    touch merged_amplicon_depth.csv
    """
}
process CSVTK_SAMPLE_AMPLICON_DEPTH {
    // Just to get the two columns of the amplicon depth file with no header
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/env-csvtk.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.29.0--h9ee0642_0' :
        'biocontainers/csvtk:0.29.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("${meta.id}_ampdepth.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    """
    csvtk cut \\
        -tT \\
        -f amplicon_id,read_count \\
        $bed \\
        | csvtk replace \\
            -tTf read_count \\
            -p "^0\$" \\
            -r 0.1 \\
        | tail -n +2 \\
        > ${meta.id}_ampdepth.tsv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(csvtk version | sed 's/csvtk v//g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_ampdepth.tsv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(csvtk version | sed 's/csvtk v//g')
    END_VERSIONS
    """
}
process CREATE_AMPLICON_COMPLETENESS {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/env-artic.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.6.2--pyhdfd78af_0' :
        'biocontainers/artic:1.6.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(consensus)
    path amplicon_bed

    output:
    path "*_amplicon_completeness.csv", emit: amplicon_completeness
    path "versions.yml", emit: versions

    script:
    """
    calc_amplicon_completeness.py \\
        --sample $meta.id \\
        --amplicon_bed $amplicon_bed \\
        --consensus $consensus

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        calc_amplicon_completeness.py: 0.1.0
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_amplicon_completeness.csv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        calc_amplicon_completeness.py: 0.1.0
    END_VERSIONS
    """
}
process CREATE_ALL_SAMPLE_SUMMARY_REPORT {
    label 'process_medium'
    publishDir "${params.outdir}", pattern: "reportDashboard.html", mode: "copy"

    conda "${moduleDir}/env-custom-report.yml"
    //container "https://depot.galaxyproject.org/singularity/artic:1.6.2--pyhdfd78af_0" to-create

    input:
    path read_variation_tsvs
    path called_variant_tsvs
    path base_coverage_beds
    path merged_amplicon_depth_csv
    path merged_amplicon_comp_csv
    path qc_csv
    path versions_yml

    output:
    path "reportDashboard.html"

    script:
    def rmd_main = "$projectDir/assets/rmarkdown-reports/reportDashboard.Rmd"
    def rmd_sample = "$projectDir/assets/rmarkdown-reports/sampleSubpage.Rmd"
    def rmd_amp = "$projectDir/assets/rmarkdown-reports/sampleAmplicons.Rmd"
    def amp_arg = merged_amplicon_depth_csv ? "run_amplicons = TRUE" : "run_amplicons = FALSE"
    """
    # Setup sample files to be found by RMD#
    cp $rmd_main $rmd_sample $rmd_amp .
    mkdir -p all_variation_positions
    mkdir -p variant_tsvs
    mkdir -p base_coverages
    mv $read_variation_tsvs all_variation_positions
    mv $called_variant_tsvs variant_tsvs
    mv $base_coverage_beds base_coverages

    # Create RMD #
    Rscript \\
        -e "library(rmarkdown)" \\
        -e "library(flexdashboard)" \\
        -e "rmarkdown::render('reportDashboard.Rmd', params=list($amp_arg))"
    """

    stub:
    """
    touch reportDashboard.html
    """
}
