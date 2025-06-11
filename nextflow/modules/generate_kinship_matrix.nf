process GENERATE_KINSHIP_MATRIX {
    conda "${projectDir}/../config/conda/kmers_gwas.yaml"

    tag "GENERATE_KINSHIP_MATRIX"

    publishDir "${params.outdir}/results/kmers_table/", mode: 'symlink'

    input:
    tuple path(kmers_table_file), path(kmers_names_file)
    path kmers_gwas_base_dir

    output:
    path "kmers_table.kinship"

    script:
    def input_table_prefix = kmers_table_file.getSimpleName().replaceFirst(/\.table$/, "") // Need some form of a prefix: smk format (prefix = lambda wildcards, output: output[0][:-8])
    def kmer_len = 31
    def mac = 5
    def maf = 0.01

    """
    export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
    ${kmers_gwas_base_dir}/bin/emma_kinship_kmers -t ${input_table_prefix} -k ${kmer_len} --maf ${maf} > kmers_table.kinship 2> kmers_table_kinship.log
    """

}
