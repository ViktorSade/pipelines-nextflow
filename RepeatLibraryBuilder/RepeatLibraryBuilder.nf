#! /usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

 params.genome = ''
 params.lib = ''

log.info """
NBIS
  _   _ ____ _____  _____
 | \\ | |  _ \\_   _|/ ____|
 |  \\| | |_) || | | (___
 | . ` |  _ < | |  \\___ \\
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

    main:
        // Analyses
        REPEAT_MODELER()
        REPEAT_FILTER()
        // Report?

    emit:

}

workflow REPEAT_MODELER {

    take:

    main:
        REPEATMODELER_PREPARELIBRARIES()
        REPEATMODELER_BUILDDB()
        REPEATMODELER_RUN()

    emit:
}

workflow REPEAT_FILTER {

    take:

    main:
        TRANSPOSONPSI()
        FILTERSEQ()
        BLAST()
        PROTEXCLUDER()

    emit:
}

process REPEATMODELER_PREPARELIBRARIES {

    input:

    output:

    script:
    """
    """

}

process REPEATMODELER_BUILDDB {

    input:

    output:

    script:
    """
    """

}

process REPEATMODELER_RUN{

    input:

    output:

    script:
    """
    """

}

process TRANSPOSONPSI{

    input:

    output:

    script:
    """
    """

}

process PROTEXCLUDER {

    input:

    output:

    script:
    """
    """

}

process BLAST {

    input:

    output:

    script:
    """
    """

}

process FILTERSEQ {

    input:

    output:

    script:
    """
    """

}

workflow.onComplete {
    log.info ( workflow.success ? "\nRepeat Library Builder complete!\n" : "Oops .. something went wrong\n" )
}
