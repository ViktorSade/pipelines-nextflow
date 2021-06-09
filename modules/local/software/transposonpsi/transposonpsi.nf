// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = 1.0.0

process TRANSPOSONPSI {
    tag "$fasta"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::transposonpsi==1.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/transposonpsi:1.0.0--hdfd78af_2"
    } else {
        container "quay.io/biocontainers/transposonpsi:1.0.0--hdfd78af_2"
    }

    input:
    path(fasta)

    output:
    path "*.all.TPSI.allHits", emit: transposon_allhits
    path "*.all.TPSI.topHits", emit: transposon_tophits
    path "*.version.txt"     , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    transposonPSI.pl $fasta prot
    echo $VERSION > ${software}.version.txt
    """
}
