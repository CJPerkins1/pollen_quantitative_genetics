params.bowtie2_extra = ""

process ALIGN_KMERS {
    conda "${projectDir}/../config/conda/align_kmers.yaml"
    tag "Aligning kmers for ${pheno}"
    publishDir "${params.outdir}/results/align_kmers/${pheno}", mode: 'symlink'

    input:
    tuple val(pheno), path(kmers_fa)
    tuple val(index_prefix), path(index_dir)

    output:
    tuple val(pheno), path("${pheno}_kmers_alignment.sam"), emit: sam

    script:
    """
    bowtie2 -p ${task.cpus} ${params.bowtie2_extra} \\
        -x ${index_dir}/${index_prefix} \\
        -f ${kmers_fa} \\
        -S ${pheno}_kmers_alignment.sam
    """

}

process SAM_TO_BAM {
    conda "${projectDir}/../config/conda/align_kmers.yaml"
    tag "Converting SAM to BAM for ${pheno}"
    publishDir "${params.outdir}/results/align_kmers/${pheno}", mode: 'symlink'

    input:
    tuple val(pheno), path(sam_file)

    output:
    tuple val(pheno), path("${pheno}_kmers_alignment.bam"), emit: bam

    script:
    """
    samtools view -@ ${task.cpus} -Sbh ${sam_file} > ${pheno}_kmers_alignment.bam
    """
}

process BAM_SORT {
    conda "${projectDir}/../config/conda/align_kmers.yaml"
    tag "Sorting BAM for ${pheno}"
    publishDir "${params.outdir}/results/align_kmers/${pheno}", mode: 'symlink'

    input:
    tuple val(pheno), path(bam_file)

    output:
    tuple val(pheno), path("${pheno}_kmers_alignment.sorted.bam"), emit: sorted_bam

    script:
    """
    samtools sort -@ ${task.cpus} ${bam_file} -o ${pheno}_kmers_alignment.sorted.bam
    """
}
