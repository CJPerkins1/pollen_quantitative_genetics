// Make file of file paths for kmers lists

process LIST_KMERS_FOUND_IN_MULTIPLE_SAMPLES {
    conda "${projectDir}/../config/conda/kmers_gwas_combine_kmc_count.yaml"

    tag "LIST_KMERS_FOUND_IN_MULTIPLE_SAMPLES "
    
    publishDir "${params.outdir}/results/kmers_found_in_multiple_samples", mode: 'symlink'

    input:
    path kmers_gwas_base_dir
    path "*"

    output:
    tuple path("kmers_to_use"), path("kmers_to_use.shareness"), path("kmers_to_use.stats.both"), path("kmers_to_use.stats.only_canonical"), path("kmers_to_use.stats.only_non_canonical"), path("kmers_list_paths.txt")

    script:
    """
    export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
    ls *_kmers_with_strand | awk -F'_' '{
        if (\$1 ~ /^[A-Z]/) {
            print \$0 "\\t" \$1
        } else if (\$1 ~ /^[a-z]/) {
            print \$0 "\\t" \$1"_"\$2
        }
    }' > kmers_list_paths.txt
    ${kmers_gwas_base_dir}/bin/list_kmers_found_in_multiple_samples -l kmers_list_paths.txt -k 31 --mac 2 -p 0.2 -o kmers_to_use
    """
}
