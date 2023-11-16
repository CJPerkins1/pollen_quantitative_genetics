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
        tuple(row.accession_id, [file: file(row.fq), paired_end: row.paired_end, srr_id: row.srr_id])
    }
    .groupTuple(by: [0])
//    .view { "Grouped data: $it" }  // Debugging statement
    .set { samples_meta_ch }

/*
 * Processes
*/

// First kmc run with canonization
process KMC_COUNT_ONE_CANONIZED {
    conda "/home/u16/cedar/git/pollen_quantitative_genetics/nextflow/config/mamba/kmc.yaml"

    tag "TEST_PROCESS_ONE on ${accession_id}"

    publishDir "${params.outdir}/results/kmc_count_one_canonized/${accession_id}", mode: 'symlink'

    input:
    tuple val(accession_id), val(file_maps)

    output:
    path "output_kmc_canon_${accession_id}*"

    script:
    // Concatenate all read paths into a single string
    def read_paths = file_maps.collect { it.file }.join(' ')

    println("Processing accession: ${accession_id}, read paths: ${read_paths}") // debugging


    // Run kmc command
    """
    kmc -t${task.cpus} -k31 -ci2 ${read_paths} output_kmc_canon_${accession_id} ./ 1> kmc_canon.1 2> kmc_canon.2
    """
}


/*
 * Workflow
*/

workflow {
    KMC_COUNT_ONE_CANONIZED(samples_meta_ch)
    // TEST_PROCESS_TWO(test_process_one_ch.collect())
}

workflow.onComplete {
    log.info ( workflow.success ? "\nSuccess" : "Failure" )
}

