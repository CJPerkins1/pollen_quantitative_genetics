process KMERS_GWAS {
    conda "${projectDir}/../config/conda/kmers_gwas_py2.yaml"

    tag "KMERS_GWAS"

    publishDir "${params.outdir}/results/kmers_gwas/", mode: 'symlink'

    input:
    tuple path(kmers_table_file), path(kmers_names_file)
    path kmers_gwas_base_dir
    tuple val(pheno_id), path(phenotype_file)
    path kinship_matrix_file

    output:
    path "${pheno_id}"

    script:
    def kmers_tab_prefix = kmers_table_file.getSimpleName().replaceFirst(/\.table$/, "")
    def out_prefix = "${pheno_id}"
    def kmer_len = 31
    def min_data_points = 30
    def mac = 5
    def maf = 0.05
    def kmers_number = 10001
    def n_permutations = 100
    // def extra = can leave blank for now

    """
    echo "--- Environment Before Modification ---"
    echo "CONDA_PREFIX: \$CONDA_PREFIX"
    echo "Original LD_LIBRARY_PATH: \$LD_LIBRARY_PATH"
    ldd "${kmers_gwas_base_dir}/external_programs/gemma_0_96" || echo "ldd check failed, continuing..."

    export LD_LIBRARY_PATH="\$CONDA_PREFIX/lib:\$CONDA_PREFIX/lib64:\$LD_LIBRARY_PATH"

    echo "--- Environment After Modification ---"
    echo "New LD_LIBRARY_PATH: \$LD_LIBRARY_PATH"
    ldd "${kmers_gwas_base_dir}/external_programs/gemma_0_96" || echo "ldd check failed, continuing..."

    python2 ${kmers_gwas_base_dir}/kmers_gwas.py --min_data_points ${min_data_points} --pheno ${phenotype_file} --kmers_table ${kmers_tab_prefix} --kmers_number ${kmers_number} --permutations ${n_permutations} --maf ${maf} --mac ${mac} -l ${kmer_len} -p ${task.cpus} --outdir ${out_prefix} >${out_prefix}.log 2>&1
    """
}

