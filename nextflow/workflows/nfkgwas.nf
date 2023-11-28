/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    k-mers-based GWAS workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.samplesheet = "${projectDir}/../samplesheets/test_samplesheet.tsv"
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
 * Loading modules
*/

include { INSTALL_KMERS_GWAS       } from "../modules/install_kmers_gwas.nf"
include { MAKE_KMC_READ_PATHS_FILE } from "../modules/make_kmc_read_paths_file.nf"
include { KMC_COUNT_CANONIZED      } from "../modules/kmc_count_canonized.nf"
include { KMC_COUNT_ALL            } from "../modules/kmc_count_all.nf"
include { COMBINE_KMC_COUNT        } from "../modules/combine_kmc_count.nf"


/*
 * Workflow
*/

workflow {
    kmers_gwas_paths_ch = INSTALL_KMERS_GWAS()
    read_paths_ch = MAKE_KMC_READ_PATHS_FILE(samples_meta_ch)
    KMC_COUNT_CANONIZED(read_paths_ch)
    KMC_COUNT_ALL(read_paths_ch)
    kmc_count_ch = KMC_COUNT_CANONIZED.out.join(KMC_COUNT_ALL.out, by:[0])
    kmc_count_combined_ch = COMBINE_KMC_COUNT(
        kmers_gwas_paths_ch,
        kmc_count_ch
    )
}


workflow.onComplete {
    log.info ( workflow.success ? "\nSuccess" : "Failure" )
}

