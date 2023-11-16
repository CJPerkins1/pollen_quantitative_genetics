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

// Short format samplesheet version
// Channel
//     .fromPath(params.samplesheet)
//     .splitCsv(header: true, sep: '\t', strip: true)
//     .map { row -> 
//         def meta = [ accession_id: row.accession_id ]
//         tuple(meta, file(row.fq1), file(row.fq2))
//     }
//     .set { samples_meta_ch }

// Long format samplesheet version
Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true, sep: '\t', strip: true)
    .map { row ->
        def meta = [
            accession_id: row.accession,
            paired_end: row.paired_end,
            srr_id: row.SRR
        ]
        tuple(meta, file(row.fq)) // Create a tuple for each row
    }
    .groupTuple(by: [0, 2]) // Group by accession_id and srr_id
    .set { samples_meta_ch }

/*
 * Processes
*/

// First kmc run with canonization
process KMC_COUNT_ONE_CANONIZED {
    conda "/home/u16/cedar/git/pollen_quantitative_genetics/nextflow/config/mamba/kmc.yaml"

    tag "TEST_PROCESS_ONE on ${meta.accession_id}"

    publishDir "${params.outdir}/results/kmc_count_one_canonized", mode: 'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.accession_id}_kmers_one"

    script:
    // Concatenate all read paths into a single string
    def readPaths = reads.join(' ')

    // Run kmc command
    """
    kmc -t${task.cpus} -k31 -ci2 ${readPaths} output_kmc_canon ${meta.accession_id}_kmers_one 1> kmc_canon.1 2> kmc_canon.2
    """
}


// process TEST_PROCESS_ONE {
//     conda "/home/u16/cedar/git/pollen_quantitative_genetics/nextflow/config/mamba/kmc.yaml"
//     
//     tag "TEST_PROCESS_ONE on ${meta.accession_id}"
// 
//     publishDir "${params.outdir}/results/process_one", mode: 'symlink'
// 
//     input:
//     tuple val(meta), path(fq1), path(fq2)
// 
//     output:
//     path "${meta.accession_id}_output.txt"
// 
//     script:
//     """
//     echo "${meta.accession_id} ${fq1} ${fq2}" > ${meta.accession_id}_output.txt
//     kmc --version
//     """
// }
// 
// process TEST_PROCESS_TWO {
//     tag "TEST_PROCESS_TWO on ${file}"
// 
//     publishDir "${params.outdir}/results/process_two", mode: 'symlink'
// 
//     input:
//     path file
// 
//     output:
//     path "concatenated_files.txt"
// 
//     script:
//     """
//     cat ${file} > concatenated_files.txt
//     """
// }

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

