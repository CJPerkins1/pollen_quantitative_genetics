// Combine the two kmc runs

process COMBINE_KMC_COUNT {
    conda "../config/conda/kmers_gwas_combine_kmc_count.yaml"

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
