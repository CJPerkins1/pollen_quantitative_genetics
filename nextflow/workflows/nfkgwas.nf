/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    k-mers-based GWAS workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
 * Setting up the input channel with metadata
*/

Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true, sep: '\t', strip: true)
    .map { row -> 
        def meta = [ accession_id: row.accession_id ]
        tuple(meta, file(row.fq1), file(row.fq2))
    }
    .set { samples_meta_ch }

/*
 * Processes
*/

process TEST_PROCESS_ONE {
    conda "/home/u16/cedar/git/pollen_quantitative_genetics/nextflow/config/mamba/kmc.yaml"
    
    tag "TEST_PROCESS_ONE on ${meta.accession_id}"

    publishDir "${params.outdir}/results/process_one", mode: 'symlink'

    input:
    tuple val(meta), path(fq1), path(fq2)

    output:
    path "${meta.accession_id}_output.txt"

    script:
    """
    echo "${meta.accession_id} ${fq1} ${fq2}" > ${meta.accession_id}_output.txt
    kmc --version
    """
}

process TEST_PROCESS_TWO {
    tag "TEST_PROCESS_TWO on ${file}"

    publishDir "${params.outdir}/results/process_two", mode: 'symlink'

    input:
    path file

    output:
    path "concatenated_files.txt"

    script:
    """
    cat ${file} > concatenated_files.txt
    """
}

/*
 * Workflow
*/

workflow {
    test_process_one_ch = TEST_PROCESS_ONE(samples_meta_ch)
    TEST_PROCESS_TWO(test_process_one_ch.collect())
}

workflow.onComplete {
    log.info ( workflow.success ? "\nSuccess" : "Failure" )
}

