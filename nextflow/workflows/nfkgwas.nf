/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    k-mers-based GWAS workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
 * Pipeline input parameters
*/

params.samplesheet = "/home/u16/cedar/git/pollen_quantitative_genetics/nextflow/samplesheets/test_samplesheet.tsv"
params.outdir = "." // Defaults to where the script is run

log.info """\
    NF - K - M E R S - G W A S   P I P E L I N E
    ============================================
    samplesheet  : ${params.samplesheet}
    outdir       : ${params.outdir}
    """
    .stripIndent()


/*
 * Setting up the input channel
*/

Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true, sep: '\t', strip: true)
    .view { row -> "${row[0]}, ${row[1]}, ${row[2]}" }
    .set { samples_ch }


/*
 * Processes
*/

process TEST_PROCESS_ONE {
    tag "TEST_PROCESS_ONE on $accession_id"

    input:
    tuple val(accession_id), path(reads1), path(reads2)

    output:
    path "${sample_id}_test_process_one_logs"

    script:
    """
    mkdir ${sample_id}_test_process_one_logs
    echo ${sample_id} ${reads1} ${reads2} > ${sample_id}_output.txt
    """
}

process TEST_PROCESS_TWO {
    tag "TEST_PROCESS_TWO on $accession_id"

    input:
    path(reads1)

    output:
    path "test_process_two_logs"

    script:
    """
    echo ${reads1}
    """
}


/*
 * Workflow
*/

workflow {
    test_process_one_ch = TEST_PROCESS_ONE(samples_ch)
    TEST_PROCESS_TWO(test_process_one_ch)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nSuccess" : "Failure" )
}
