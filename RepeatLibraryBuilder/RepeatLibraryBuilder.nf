#! /usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.genome = ''                                // Genome to build Repeat library from
params.organism_name = ''                         // Organism scientific name
params.repeatmasker_library = ''                  // Path to repeat masker lib
params.repeatpeps_library = ''                    // Path to repeat peps lib
params.transposon_psi_lib = ''                    // Path to uniprot fasta
params.build_repeat_masker_lib = false            // True if build blast db on repeat masker lib
params.build_repeat_peps_lib = false              // True if build blast db on repeat peps lib
params.build_transposon_psi_lib = false           // True if build blast db on uniprot fasta

log.info """
NBIS
  _    _ ____ _____  _____
 | \\ | |  _ \\_  _|/ ____|
 |  \\| | |_) || | | (___
 | . `  |  _ < | |  \\___\\
 | |\\  | |_) || |_ ____) |
 |_| \\_|____/_____|_____/  Annotation Service

 Repeat Library Builder workflow
 ===================================

 General parameters

 """

workflow {

    main:
        REPEAT_LIBRARY_BUILDER()

}

workflow REPEAT_LIBRARY_BUILDER {

    take:
        genome

    main:
        // Analyses
        BLASTX_MAKEBLASTDB() // import as X, Y, Z // conditional execution / input.parameter?
        REPEATMODELER_BUILDDB()
        REPEATMODELER_RUN()
        TRANSPOSONPSI() // Check if Uniprot updated? Check input or parameter?
        GAAS_FILTERSEQ()
        BLAST()
        PROTEXCLUDER()

        // Report?

    emit:
        repeat_library = REPEATMODELER_RUN().out.repeat_library
        gene_filtered_repeat_library = PROTEXCLUDER().out.repeat_library
        filtered_proteins = GAAS_FILTERSEQ().out.filtered_proteins

}

process BLASTX_MAKEBLASTDB {  // Import as module, call 3x

    input:
    val type
    path library

    output:
    path "$library*", includeInputs: true, emit: blast_db

    script:
    """
    makeblastdb -dbtype $type -in $library
    """

}

process REPEATMODELER_BUILDDB {

    input:
    path genome
    val organism_name

    output:
    path "${organism_name}.*", emit: repeat_modeler_db

    script:
    """
    BuildDatabase -name $organism_name -engine ncbi $genome
    """

}

process REPEATMODELER_RUN {

    input:
    path dbfiles

    output:
    path "*/consensi.fa.classified", emit: repeat_library

    script:
    database_name = // generate prefix.
    """
    RepeatModeler –database $database_name -engine ncbi –pa ${task.cpus}
    """

}

process TRANSPOSONPSI {

    input:
    path fasta  // proteins.fasta

    output:
    path "${fasta}.all.TPSI.{allHits,topHits}", emit: transposon_hits

    script:
    """
    # TODO: check how this is distributed.
    gaas_transposonPSI2grid.pl -f $fasta -o transposonPSI
    """

}

process GAAS_FILTERSEQ {

    input:
    path fasta
    path tophits

    output:
    path "proteins.filtered.fa", emit: filtered_fasta

    script:
    """
    awk '{if(\$0 ~ /^[^\\/\\/.*]/) print \$5}' $tophits | sort -u > accessions.list
    gaas_fasta_removeSeqFromIDlist.pl -f $fasta -l accessions.list -o proteins.filtered.fa
    """

}

process BLAST {
    //versions ncbi-blast-2.2.28+ and ncbi-blast-2.4.0+

    input:
    path protein_db //  /projects/references/databases/uniprot/2020-12/transposon/proteins.filtered.fa
    path fasta // consensi.fa.classified

    output:
    path "blastx.out", emit: blast_hits

    script:
    """
    // makeblastdb –in proteins.filtered.fa –dbtype prot // do BLASTX_MAKEBLASTDB
    blastx -db $protein_db -query $fasta -num_threads ${task.cpus} -out blastx.out
    """

}

process PROTEXCLUDER {

    input:
    path blast_hits
    path fasta // consensi.fa.classified

    output:
    path "${fasta}noProtFinal", emit: repeat_library

    script:
    """
    ProtExcluder.pl $blast_hits $fasta
    """

}

workflow.onComplete {
    log.info ( workflow.success ? "\nRepeat Library Builder complete!\n" : "Oops .. something went wrong\n" )
}
