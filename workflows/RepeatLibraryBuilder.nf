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
params.fasta_chunk_size = 10                      // How many records to store per fasta, when distributing processing

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
            // Need to pass the whole library or move uniprot build db
            uniprot_db.set { filtered_uniprot }
        } else {
            TRANSPOSONPSI(uniprot_db.splitFasta( by: params.fasta_chunk_size))
            GAAS_FILTERSEQ(uniprot_db,TRANSPOSONPSI.out.tophits.collectFile())
            GAAS_FILTERSEQ.out.fasta.set { filtered_uniprot }
        }
        BUILD_UNIPROT_DB(filtered_uniprot)       // uses `storeDir` to determine if build needed?
        BLAST_BLASTX(REPEATMODELER_REPEATMODELER.out.repeat_sequences,
            BUILD_UNIPROT_DB.out.db)
        PROTEXCLUDER(REPEATMODELER_REPEATMODELER.out.repeat_sequences
            .join(BLAST_BLASTX.out.txt))

        // Report?
        // grep ">" consensi.fa.classifiednoProtFinal| sed 's/.*#//' | sed 's/ .*//' | sort -k1,1 | uniq -c  | awk 'BEGIN{OFS="\t";print "Repeat","Frequency"}{OFS="\t";print $2,$1}'
        /*
        Repeat	Frequency
        DNA/TcMar-Ant1	1
        DNA/TcMar-Fot1	3
        DNA/TcMar-Fot1	2
        LTR/Gypsy	3
        LTR/Gypsy	1
        Unknown	42
        Unknown	5
        */

    emit:
        unfiltered_repeat_library = REPEATMODELER_REPEATMODELER.out.repeat_sequences
        gene_filtered_repeat_library = PROTEXCLUDER.out.repeat_sequences
        filtered_proteins = GAAS_FILTERSEQ.out.filtered_sequences.ifEmpty([])   // Or do you want the next stage?

}

workflow.onComplete {
    log.info ( workflow.success ? "\nRepeat Library Builder complete!\n" : "Oops .. something went wrong\n" )
}
