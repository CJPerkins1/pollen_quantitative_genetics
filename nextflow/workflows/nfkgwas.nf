/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    k-mers-based GWAS workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.samplesheet = "${projectDir}/../samplesheets/test_samplesheet.tsv"
params.cleanup_fastq = false
params.outdir = "." // Defaults to where the script is run

log.info """\
    NF - K - M E R S - G W A S   P I P E L I N E
    ============================================
    samplesheet  : ${params.samplesheet}
    outdir       : ${params.outdir}
    cleanup      : ${params.cleanup_fastq}
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


/*
 * Workflow
*/

workflow {
    kmers_gwas_paths_ch = INSTALL_KMERS_GWAS()
    read_paths_ch = MAKE_KMC_READ_PATHS_FILE(samples_meta_ch)
    KMC_COUNT_CANONIZED(read_paths_ch)
    KMC_COUNT_ALL(read_paths_ch)

    KMC_COUNT_CANONIZED.out.view() // Add view for Canonized output
    KMC_COUNT_ALL.out.view()       // Add view for All output

    // kmc_count_ch collects outputs after both KMC processes finish for a sample
    kmc_count_ch = KMC_COUNT_CANONIZED.out.join(KMC_COUNT_ALL.out, by:[0])

    //kmc_count_ch.view()

    // Create a channel with just accession_id and individual fastq file paths from samples_meta_ch
    // This channel will be joined with kmc_count_ch to trigger cleanup after KMC
    fastq_paths_for_cleanup_ch = samples_meta_ch
        .flatMap { accession_id, metas ->
            metas.collect { meta -> tuple(accession_id, meta.file) }
        }
        // .view() // Debug: see what this channel contains

    // Optional cleanup - now triggered by the completion of KMC processes
    log.info "params.cleanup_fastq is: ${params.cleanup_fastq}"
    if (params.cleanup_fastq) {

        // Create a channel emitting [accession_id, list_of_fastq_file_objects] per sample
        samples_fastq_files_ch = samples_meta_ch
            .map { accession_id, metas ->
                def fastq_files = metas.collect { it.file } // Extract file objects into a list
                tuple(accession_id, fastq_files)
            }
            // .view() // Debug: see what this channel contains

        // Join kmc_count_ch (emits when KMC for a sample is done)
        // with samples_fastq_files_ch (emits list of files for a sample) by accession_id
        // This ensures cleanup is triggered only after KMC processes complete for a sample
        kmc_count_ch
            .join(samples_fastq_files_ch, by:[0])
            .map { accession_id, kmc_canon_pre, kmc_canon_suf, kmc_canon_out, kmc_all_pre, kmc_all_suf, kmc_all_out, fastq_files_list ->
                // We only need accession_id and fastq_files_list for cleanup
                tuple(accession_id, fastq_files_list)
            }
            // .view() // Debug: see what is being passed to CLEANUP_FASTQ
            .set { fastq_cleanup_ch_triggered } // Channel name for the cleanup input

        // Feed the triggered channel to the CLEANUP_FASTQ process
        CLEANUP_FASTQ(fastq_cleanup_ch_triggered)
    }

    DO_KMERS_STATS(
        kmc_count_ch.map{ tuple -> return [tuple[3], tuple[6]]}.collect()
    )
    kmc_count_combined_ch = COMBINE_KMC_COUNT(
        kmers_gwas_paths_ch,
        kmc_count_ch
    )
    kmers_list_ch = LIST_KMERS_FOUND_IN_MULTIPLE_SAMPLES(
        kmers_gwas_paths_ch,
        kmc_count_combined_ch.collect()
    )
    kmers_table_ch = BUILD_KMERS_TABLE(
        kmers_gwas_paths_ch,
        kmc_count_combined_ch.collect(),
        kmers_list_ch
    )
}


workflow.onComplete {
    log.info ( workflow.success ? "\nSuccess" : "Failure" )
}

workflow.onComplete {
    log.info ( workflow.success ? "\nSuccess" : "Failure" )
}

