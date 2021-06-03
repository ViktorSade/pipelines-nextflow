// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process REPEATMODELER_BUILDDATABASE {
    tag "$genome"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process),  meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::repeatmodeler==2.0.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/repeatmodeler:2.0.1--pl526_0"
    } else {
        container "quay.io/biocontainers/repeatmodeler:2.0.1--pl526_0"
    }

    input:
    // tuple val(meta), path(fasta)
    path genome
    val organism_name

    output:
    // tuple val(meta), path('*.blastn.txt'), emit: txt
    path "repeatmodeler_db"              , emit: db
    path '*.version.txt'                 , emit: version

    script:
    def software = getSoftwareName(task.process)
    // def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    # -engine ncbi into \$options.args
    BuildDatabase -name $organism_name $options.args $genome
    mkdir repeatmodeler_db
    mv ${organism_name}.* repeatmodeler_db
    # FIXME RepeatModeler version
    echo \$(RepeatModeler -version 2>&1) | sed 's/^.*blastn: //; s/ .*\$//' > ${software}.version.txt
    """
}
