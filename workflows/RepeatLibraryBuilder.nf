#! /usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

/*
 Samplesheet column headers
 - name : The name to use for file naming. Often based on the scientific name.
 - fasta : The path to the fasta.
 */
params.input = ''                                 // Path to Samplesheet
params.repeatmasker_db = ''                       // Path to Repeat Masker db
params.transposible_element_db = ''               // Path to Transposible element db
params.uniprot_db = ''                            // Path to Uniprot fasta
params.uniprot_is_filtered = true                 // True if the Uniprot Fasta is filtered for transposible elements

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

include { get_sample_info } from './utilities.nf'

workflow {

    main:
        Channel.fromPath(params.input)
            .splitCsv(header:true, sep:',')
            .map { get_sample_info(it) }
            .set { input_ch }
        REPEAT_LIBRARY_BUILDER(input_ch,
            file(params.repeatmasker_db, checkIfExists:true),
            file(params.transposible_element_db, checkIfExists:true),
            file(params.uniprot_db, checkIfExists:true))

}

include { BLAST_MAKEBLASTDB as BUILD_REPEATMASKER_DB  } from '../modules/local/blast/makeblastdb'           addParams(options:modules['build_repeatmasker_db'])
include { BLAST_MAKEBLASTDB as BUILD_TRANSPOSIBLE_DB  } from '../modules/local/blast/makeblastdb'           addParams(options:modules['build_te_db'])
include { BLAST_MAKEBLASTDB as BUILD_UNIPROT_DB       } from '../modules/local/blast/makeblastdb'           addParams(options:modules['build_uniprot_db'])
include { REPEATMODELER_BUILDDATABASE                 } from '../modules/local/repeatmodeler/builddatabase' addParams(options:modules['repeatmodeler_builddb'])
include { REPEATMODELER_REPEATMODELER                 } from '../modules/local/repeatmodeler/repeatmodeler' addParams(options:modules['repeatmodeler'])
include { TRANSPOSONPSI                               } from '../modules/local/transposonpsi/transposonpsi' addParams(options:modules['transposonpsi'])
include { GAAS_FILTERSEQ                              } from '../modules/local/gaas/filterseq'              addParams(options:modules['gaas_filterseq'])
include { BLAST_BLASTX                                } from '../modules/local/blast/blastx'                addParams(options:modules['blastx'])
include { PROTEXCLUDER                                } from '../modules/local/protexcluder/protexcluder'   addParams(options:modules['protexcluder'])

workflow REPEAT_LIBRARY_BUILDER {

    take:
        genome                                  // [name, path to genome]
        rep_mask_db                             // path to repeat masker library
        rep_pep_db                              // path to transposible element library
        uniprot_db                              // path to Uniprot database

    main:
        // Analyses
        BUILD_REPEATMASKER_DB(rep_mask_db)       // uses `storeDir` to determine if build needed?
        BUILD_TRANSPOSIBLE_DB(rep_pep_db)        // uses `storeDir` to determine if build needed?
        REPEATMODELER_BUILDDATABASE(genome)      // Does this need the db's above?
        REPEATMODELER_REPEATMODELER(
            REPEATMODELER_BUILDDATABASE.out.db,
            BUILD_REPEATMASKER_DB.out.db,        // Which process needs the db's?
            BUILD_TRANSPOSIBLE_DB.out.db)
        if (params.uniprot_is_filtered){
            uniprot_db.set { filtered_uniprot_db }
        } else {
            TRANSPOSONPSI(uniprot_db)
            GAAS_FILTERSEQ(uniprot_db,TRANSPOSONPSI.out.tophits)
            BUILD_UNIPROT_DB(GAAS_FILTERSEQ.out.filtered_sequences)             // uses `storeDir` to determine if build needed
            BUILD_UNIPROT_DB.out.db.set { filtered_uniprot_db }
        }
        BLAST_BLASTX(REPEATMODELER_REPEATMODELER.out.repeat_sequences,
            filtered_uniprot_db)
        PROTEXCLUDER(REPEATMODELER_REPEATMODELER.out.repeat_sequences,
            BLAST_BLASTX.out.txt)

        // Report?

    emit:
        repeat_library = REPEATMODELER_REPEATMODELER.out.repeat_sequences
        gene_filtered_repeat_library = PROTEXCLUDER.out.repeat_sequences
        filtered_proteins = GAAS_FILTERSEQ.out.filtered_sequences.ifEmpty([])   // Or do you want the next stage?

}

/* process BLASTX_MAKEBLASTDB {  // Import as module, call 3x

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

process REPEATMODELER_REPEATMODELER {

    input:
    path dbfiles

    output:
    path "* /consensi.fa.classified", emit: repeat_library

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

process BLAST_BLASTX {
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
 */

workflow.onComplete {
    log.info ( workflow.success ? "\nRepeat Library Builder complete!\n" : "Oops .. something went wrong\n" )
}
