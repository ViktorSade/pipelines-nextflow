// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = 1.2.0 // Help version text (`--help`) is not reliable for version

process GAAS_FILTERSEQ {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gaas=1.2.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gaas:1.2.0--pl526r35_0"
    } else {
        container "quay.io/biocontainers/gaas:1.2.0--pl526r35_0"
    }

    input:
    tuple val(meta), path(tophits)
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("proteins.filtered.fa"), emit: fasta
    path "*.version.txt"                         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    awk '{if(\$0 ~ /^[^\\/\\/.*]/) print \$5}' $tophits | sort -u > accessions.list
    gaas_fasta_removeSeqFromIDlist.pl -f $fasta -l accessions.list -o proteins.filtered.fa
    echo $VERSION > ${software}.version.txt
    """
}
