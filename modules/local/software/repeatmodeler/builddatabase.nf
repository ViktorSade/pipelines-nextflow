// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process REPEATMODELER_BUILDDATABASE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::repeatmodeler==2.0.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/repeatmodeler:2.0.1--pl526_0"
    } else {
        container "quay.io/biocontainers/repeatmodeler:2.0.1--pl526_0"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path "repeatmodeler_db"          , emit: db
    path '*.version.txt'                              , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    BuildDatabase -name ${meta.id} $options.args $fasta
    mkdir repeatmodeler_db
    mv ${meta.id}.* repeatmodeler_db
    RepeatModeler --version | sed -e 's/RepeatModeler version //' > ${software}.version.txt
    """
}
