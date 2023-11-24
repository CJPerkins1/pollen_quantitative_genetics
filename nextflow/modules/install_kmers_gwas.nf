// Install kmers-gwas

process INSTALL_KMERS_GWAS {
    conda "${projectDir}/../config/conda/kmers_gwas.yaml"

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
