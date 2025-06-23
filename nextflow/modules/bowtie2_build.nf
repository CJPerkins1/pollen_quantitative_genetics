process BOWTIE2_BUILD {
    conda "${projectDir}/../config/conda/align_kmers.yaml"
    tag "bowtie2_build on ${reference_genome.baseName}"
    publishDir "${params.outdir}/results/bowtie2_index", mode: 'symlink'

    input:
    path reference_genome

    output:
    tuple val(reference_genome.baseName), path("index")

    script:
    def prefix = reference_genome.baseName
    """
    mkdir index
    bowtie2-build --threads ${task.cpus} ${reference_genome} index/${prefix}
    """
}
