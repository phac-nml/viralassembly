process GET_MODEL {
    label 'process_single'
    tag "Download"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/requests:2.26.0' :
        'biocontainers/requests:2.26.0' }"

    input:
    val model

    output:
    path("clair3_models/${model}"), emit: model
    path("versions.yml"), emit: versions

    script:
    // Downloading specific clair3 model only if its not bundled with clair3
    //  If model exists in clair3, match the method and copy it over
    //  If not download
    //  Downloading each time is unideal, however, its only 1 model now
    """
    # Have to try this with conda to make sure it also works there
    CLAIR_BIN_DIR=\$(dirname "\$(which run_clair3.sh)")

    mkdir -p clair3_models
    if ls \$CLAIR_BIN_DIR/models/ | grep -q '$model' ; then
        cp -r \$CLAIR_BIN_DIR/models/$model ./clair3_models/
    else
        download_models.py --model-dir ./clair3_models --model $model
    fi

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        requests: 2.26.0
        get_model: 2026.01.09
    END_VERSIONS
    """

    stub:
    """
    mkdir -p clair3_models/$model

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        get_model: 2026.01.09
    END_VERSIONS
    """
}
