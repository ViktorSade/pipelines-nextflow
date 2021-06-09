/*
 Formats the sample sheet row into the structure expected by
 module processes. Modified from the nf-core DSL2 workflow template

 Samplesheet expected headers
 name,fasta
 */
def get_sample_info(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.name
    // meta.single_end   = row.single_end.toBoolean()

    if (!file(row.fasta).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Fasta file does not exist!\n${row.fasta}"
    }
    def array = [ meta, [ file(row.fasta) ] ]
    /* if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    } */
    return array
}
