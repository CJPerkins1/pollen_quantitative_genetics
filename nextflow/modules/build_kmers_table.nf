// Make file of file paths for kmers lists

process BUILD_KMERS_TABLE{
    conda "${projectDir}/../config/conda/kmers_gwas_combine_kmc_count.yaml"

    tag "BUILD_KMERS_TABLE"

    publishDir "${params.outdir}/results/kmers_table/", mode: 'link'

    input:
    path kmers_gwas_base_dir // kmers_gwas_paths_ch
    path "*" // kmers_count_combined_ch.collect()
    tuple path(kmers_to_use), path(kmers_to_use_shareness), path(kmers_to_use_stats_both), path(kmers_to_use_stats_only_canonical), path(kmers_to_use_stats_only_non_canonical), path(kmers_list_paths) // kmers_list_ch

    output:
    tuple path("kmers_table.table"), path("kmers_table.names")

    script:
    """
    export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
    ${kmers_gwas_base_dir}/bin/build_kmers_table -l ${kmers_list_paths} -k 31 -a ${kmers_to_use} -o kmers_table
    """
}
