/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    k-mers-based GWAS workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.samplesheet = "${projectDir}/../samplesheets/test_samplesheet.tsv"
params.phenotypes_samplesheet = "${projectDir}/../samplesheets/nf_test_pheno.tsv"
params.cleanup_fastq = false
params.outdir = "." // Defaults to where the script is run
params.only_unique_kmers = false
params.run_convert_kmers_table_to_plink = false

log.info """\
    NF - K - M E R S - G W A S   P I P E L I N E
    ============================================
    samplesheet  : ${params.samplesheet}
    phenotypes   : ${params.phenotypes_samplesheet}
    outdir       : ${params.outdir}
    cleanup      : ${params.cleanup_fastq}
    unique_kmers : ${params.only_unique_kmers}
    kmers_to_plink : ${params.run_convert_kmers_table_to_plink}
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
    .view()
    .set { samples_meta_ch }

Channel
    .fromPath(params.phenotypes_samplesheet)
    .splitCsv(sep: '\t', header: true)
    .map { row -> tuple(row.pheno_name, file(row.pheno_path)) }
    .view()
    .set { phenotype_channel }

/*
 * Loading modules
*/

include { INSTALL_KMERS_GWAS                   } from "../modules/install_kmers_gwas.nf"
include { MAKE_KMC_READ_PATHS_FILE             } from "../modules/make_kmc_read_paths_file.nf"
include { KMC_COUNT_CANONIZED                  } from "../modules/kmc_count_canonized.nf"
include { KMC_COUNT_ALL                        } from "../modules/kmc_count_all.nf"
include { CLEANUP_FASTQ                        } from "../modules/cleanup_fastq.nf"
include { DO_KMERS_STATS                       } from "../modules/do_kmers_stats.nf"
include { COMBINE_KMC_COUNT                    } from "../modules/combine_kmc_count.nf"
include { LIST_KMERS_FOUND_IN_MULTIPLE_SAMPLES } from "../modules/list_kmers_found_in_multiple_samples.nf"
include { BUILD_KMERS_TABLE                    } from "../modules/build_kmers_table.nf"
include { CONVERT_KMERS_TABLE_TO_PLINK         } from "../modules/convert_kmers_table_to_plink.nf"

/*
 * Workflow
*/

workflow {
    kmers_gwas_paths_ch = INSTALL_KMERS_GWAS()
    read_paths_ch = MAKE_KMC_READ_PATHS_FILE(samples_meta_ch)
    KMC_COUNT_CANONIZED(read_paths_ch)
    KMC_COUNT_ALL(read_paths_ch)

    kmc_count_ch = KMC_COUNT_CANONIZED.out.join(KMC_COUNT_ALL.out, by:[0])

    fastq_paths_for_cleanup_ch = samples_meta_ch
        .flatMap { accession_id, metas ->
            metas.collect { meta -> tuple(accession_id, meta.file) }
        }

    // Optional cleanup - now triggered by the completion of KMC processes
    log.info "params.cleanup_fastq is: ${params.cleanup_fastq}"
    if (params.cleanup_fastq) {

        samples_fastq_files_ch = samples_meta_ch
            .map { accession_id, metas ->
                def fastq_files = metas.collect { it.file }
                tuple(accession_id, fastq_files)
            }
        kmc_count_ch
            .join(samples_fastq_files_ch, by:[0])
            .map { accession_id, kmc_canon_pre, kmc_canon_suf, kmc_canon_out, kmc_all_pre, kmc_all_suf, kmc_all_out, fastq_files_list ->
                tuple(accession_id, fastq_files_list)
            }
            .set { fastq_cleanup_ch_triggered }

        CLEANUP_FASTQ(fastq_cleanup_ch_triggered)
    }

    DO_KMERS_STATS(
        kmc_count_ch.map{ tuple -> return [tuple[3], tuple[6]]}.collect()
    )
    kmc_count_combined_ch = COMBINE_KMC_COUNT(
        kmers_gwas_paths_ch,
        kmc_count_ch
    )
    kmc_count_combined_ch.view()
    kmers_list_ch = LIST_KMERS_FOUND_IN_MULTIPLE_SAMPLES(
        kmers_gwas_paths_ch,
        kmc_count_combined_ch.collect()
    )
    kmers_table_ch = BUILD_KMERS_TABLE(
        kmers_gwas_paths_ch,
        kmc_count_combined_ch.collect(),
        kmers_list_ch
    )
    kmers_table_ch.view()
    if (params.run_convert_kmers_table_to_plink) {
        (bed_ch, bim_ch, fam_ch, log_file_ch) = CONVERT_KMERS_TABLE_TO_PLINK(
            kmers_table_ch,
            phenotype_channel,
            kmers_gwas_paths_ch
        )
        bed_ch.view()
        bim_ch.view()
        fam_ch.view()
        log_file_ch.view()
    }
}


workflow.onComplete {
    log.info ( workflow.success ? "\nSuccess" : "Failure" )
}

