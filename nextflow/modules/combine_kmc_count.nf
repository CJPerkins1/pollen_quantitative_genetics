// Combine the two kmc runs

process COMBINE_KMC_COUNT {
    conda "${projectDir}/../config/conda/kmers_gwas_combine_kmc_count.yaml"

    tag "COMBINE_KMC_COUNT on ${accession_id}"

    publishDir "${params.outdir}/results/combine_kmc_count/${accession_id}", mode: 'symlink'

    input:
    path kmers_gwas_base_dir // kmers_gwas_paths_ch
    tuple val(accession_id), path(canonized_pre), path(canonized_suf), path(all_pre), path(all_suf) // kmc_count_ch

    output:
    path "${accession_id}_kmers_with_strand"

    script:
    println("KMC count combine on accession: ${accession_id}, kmers_gwas_dir: ${kmers_gwas_base_dir}, kmc_count_canonized_pre: ${canonized_pre}, kmc_count_canonized_suf: ${canonized_suf}, kmc_count_all_pre: ${all_pre}, kmc_count_all_suf: ${all_suf}") // debugging
//    """
//    echo combine_count on kmc_count_canonized_pre: ${canonized_pre}, kmc_count_all_pre: ${all_pre}
//    """
    """
    export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
    ${kmers_gwas_base_dir}/bin/kmers_add_strand_information -c ${accession_id}_kmc_count_canonized -n ${accession_id}_kmc_count_all -k 31 -o ${accession_id}_kmers_with_strand
    """
}
