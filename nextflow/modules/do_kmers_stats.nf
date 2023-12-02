// Combine the two kmc runs

process DO_KMERS_STATS {
    //conda "${projectDir}/../config/conda/do_kmers_stats.yaml"

    tag "DO_KMERS_STATS"

    publishDir "${params.outdir}/results/do_kmers_stats/

    input:
    path "*" //

    output:
    path("kmers_stats_summary.tsv")

    script:
    """
    ${projectDir}/../scripts/make_kmers_stats_summary.sh
    """
}
