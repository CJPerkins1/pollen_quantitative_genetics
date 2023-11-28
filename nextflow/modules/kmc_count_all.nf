// Second kmc run without canonization

process KMC_COUNT_ALL {
    conda "${projectDir}/../config/conda/kmc.yaml"

    tag "KMC_COUNT_ALL on ${accession_id}"

    publishDir "${params.outdir}/results/kmc_count_all/${accession_id}", mode: 'symlink'

    input:
    tuple val(accession_id), val(file_maps), path(read_paths_file)

    output:
    tuple val(accession_id), path("${accession_id}_kmc_count_all.kmc_pre"), path("${accession_id}_kmc_count_all.kmc_suf")

    script:
    println("KMC count all on accession: ${accession_id}, file list: ${read_paths_file}") // debugging
    // Run kmc command
    """
    kmc -t${task.cpus} -k31 -ci0 -b @${read_paths_file} ${accession_id}_kmc_count_all ./ 1> ${accession_id}_kmc_all.1 2> ${accession_id}_kmc_all.2
    """
}
