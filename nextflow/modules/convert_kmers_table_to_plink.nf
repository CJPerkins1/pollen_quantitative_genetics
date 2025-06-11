process CONVERT_KMERS_TABLE_TO_PLINK {
    conda "${projectDir}/../config/conda/kmers_gwas.yaml"

    tag "CONVERT_KMERS_TABLE_TO_PLINK_${pheno_id}"

    publishDir "${params.outdir}/results/kmers_table/plink/pheno_${pheno_id}/", mode: 'symlink'

    input:
    tuple path(kmers_table_file), path(kmers_names_file)  // <- receives both
    tuple val(pheno_id), path(phenotype_file)
    path kmers_gwas_base_dir

    output:
    path "kmers_table.presence_absence.0.0.bed", emit: bed
    path "kmers_table.presence_absence.0.0.bim", emit: bim
    path "kmers_table.presence_absence.0.0.fam", emit: fam
    path "kmers_table.presence_absence.0.log", emit: log_file

    script:
    def kmer_len = 31
    def max_num_var = 1000000
    def mac = 5
    def maf = 0.01
    def only_unique_flag = params.only_unique_kmers ? "-u" : ""

    // strip `.table` from the file name
    def input_table_prefix = kmers_table_file.getSimpleName().replaceFirst(/\.table$/, "")

    """
    echo "DEBUG CONVERT_KMERS_TABLE_TO_PLINK (Mirroring BUILD_KMERS_TABLE directly):"
    echo "  kmers_table_file: ${kmers_table_file}"
    echo "  input_table_prefix: ${input_table_prefix}"
    echo "  phenotype_file: ${phenotype_file}"
    echo "  kmers_gwas_base_dir: ${kmers_gwas_base_dir}"
    echo "------------------------------------"

    export LD_LIBRARY_PATH=\$CONDA_PREFIX/lib

    ${kmers_gwas_base_dir}/bin/kmers_table_to_bed \\
      -t ${input_table_prefix} \\
      -p ${phenotype_file} \\
      -o kmers_table.presence_absence.0 \\
      -k ${kmer_len} \\
      --maf ${maf} --mac ${mac} \\
      -b ${max_num_var} \\
      ${only_unique_flag} \\
      > kmers_table.presence_absence.0.log 2>&1
    """
}
