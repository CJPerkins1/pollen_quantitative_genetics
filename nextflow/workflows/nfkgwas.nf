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

// Install kmers-gwas
process INSTALL_KMERS_GWAS {
    conda "/home/u16/cedar/git/pollen_quantitative_genetics/nextflow/config/conda/kmers_gwas.yaml"

    tag "INSTALL_KMERS_GWAS"

    output:
    path "kmers_gwas_dir"

    script:
    """
    mkdir kmers_gwas_dir
    wget https://github.com/voichek/kmersGWAS/releases/download/v0.3-beta/v0_3_beta.zip
    unzip v0_3_beta.zip -d kmers_gwas_dir
    """
}

// Make files of file paths for kmc input
process MAKE_KMC_READ_PATHS_FILE {
    tag "MAKE_KMC_READ_PATHS_FILE on ${accession_id}"

    input:
    tuple val(accession_id), val(file_maps)

    output:
    tuple val(accession_id), val(file_maps), path("${task.workDir}/file_paths_${accession_id}.txt")

    exec:
    def file_list="${task.workDir}/file_paths_${accession_id}.txt"
    new File(file_list).withWriter { writer ->
        file_maps.each { writer.println(it.file) }
    }
}

// First kmc run with canonization
process KMC_COUNT_CANONIZED {
    conda "/home/u16/cedar/git/pollen_quantitative_genetics/nextflow/config/conda/kmc.yaml"

    tag "KMC_COUNT_CANONIZED on ${accession_id}"

    publishDir "${params.outdir}/results/kmc_count_one_canonized/${accession_id}", mode: 'symlink'

    input:
    tuple val(accession_id), val(file_maps), path(read_paths_file)

    output:
    path "${accession_id}_kmc_count_canonized"

    script:
    println("KMC count canonized on accession: ${accession_id}, file list: ${read_paths_file}") // debugging
    // Run kmc command
    """
    mkdir ${accession_id}_kmc_count_canonized && \
    kmc -t${task.cpus} -k31 -ci2 @${read_paths_file} output_kmc_canon_${accession_id} ./${accession_id}_kmc_count_canonized 1> ${accession_id}_kmc_canon.1 2> ${accession_id}_kmc_canon.2
    """
}

// Second kmc run without canonization
process KMC_COUNT_ALL {
    conda "/home/u16/cedar/git/pollen_quantitative_genetics/nextflow/config/conda/kmc.yaml"

    tag "KMC_COUNT_ALL on ${accession_id}"

    publishDir "${params.outdir}/results/kmc_count_all/${accession_id}", mode: 'symlink'

    input:
    tuple val(accession_id), val(file_maps), path(read_paths_file)

    output:
    path "${accession_id}_kmc_count_all"

    script:
    println("KMC count all on accession: ${accession_id}, file list: ${read_paths_file}") // debugging
    // Run kmc command
    """
    mkdir ${accession_id}_kmc_count_all && \
    kmc -t${task.cpus} -k31 -ci0 -b @${read_paths_file} output_kmc_all_${accession_id} ./${accession_id}_kmc_count_all 1> ${accession_id}_kmc_all.1 2> ${accession_id}_kmc_all.2
    """
}

// Combine the two kmc runs
process COMBINE_KMC_COUNT {
    conda "/home/u16/cedar/git/pollen_quantitative_genetics/nextflow/config/conda/kmers_gwas.yaml"

    tag "COMBINE_KMC_COUNT on ${accession_id}"

    publishDir "${params.outdir}/results/combine_kmc_count/${accession_id}", mode: 'symlink'

    input:
    tuple val(accession_id), val(file_maps), path(read_paths_file) // read_paths_ch
    path kmers_gwas_base_dir // kmers_gwas_paths_ch
    path kmc_count_canonized_dir // kmc_count_canonized_ch
    path kmc_count_all_dir // kmc_count_all_ch

    output:
    path "${accession_id}_kmc_count_combined"

    script:
    println("KMC count combine on accession: ${accession_id}, kmers_gwas_dir: ${kmers_gwas_base_dir}, kmc_count_canonized_dir: ${kmc_count_canonized_dir}, kmc_count_all_dir: ${kmc_count_all_dir}") // debugging
    """
    ${kmers_gwas_base_dir}/bin/kmers_add_strand_information -c ${kmc_count_canonized_dir} -n ${kmc_count_all_dir} -k 31 -o ${accession_id}_kmers_with_strand
    """
}

/*
 * Workflow
*/

workflow {
    kmers_gwas_paths_ch = INSTALL_KMERS_GWAS()
    read_paths_ch = MAKE_KMC_READ_PATHS_FILE(samples_meta_ch)
    kmc_count_canonized_ch = KMC_COUNT_CANONIZED(read_paths_ch)
    kmc_count_all_ch = KMC_COUNT_ALL(read_paths_ch)
    kmc_count_combined_ch = COMBINE_KMC_COUNT(
        read_paths_ch,
        kmers_gwas_paths_ch,
	kmc_count_canonized_ch,
	kmc_count_all_ch
    )
}

workflow.onComplete {
    log.info ( workflow.success ? "\nSuccess" : "Failure" )
}

