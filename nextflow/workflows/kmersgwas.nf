/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    k-mers-based GWAS workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.samplesheet = "$workflow.scriptDir/../samplesheets/test_samplesheet.tsv"

Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true, sep: '\t', strip: true)
    .map { row ->
        tuple(
            row.accession_id,
            path(row.fastq1, checkIfExists: true),
            path(row.fastq2, checkIfExists: true),
            path(row.fastqx, checkIfExists: true)
        )
    }
    .set { samples_ch }

process FilterSamples {
    input:
    tuple val(accession_id), path(fastq1), path(fastq2), path(fastqx) from samples_ch

    output:
    tuple val(accession_id), path(fastq1), path(fastq2), path(fastqx) into filtered_samples_ch

    script:
    """
    # Here you would have some condition to filter the samples
    # For the sake of this example, we'll just pass all samples through
    echo "$accession_id" > /dev/null
    """
}

process ProcessFilteredSamples {
    input:
    tuple val(accession_id), path(fastq1), path(fastq2), path(fastqx) from filtered_samples_ch

    script:
    """
    echo "Processing sample: $accession_id"
    echo "FASTQ1 file: $fastq1"
    echo "FASTQ2 file: $fastq2"
    echo "FASTQX file: $fastqx"
    """
}

workflow {
    take:
    ch_samples from samples_ch

    main:
    ch_samples | FilterSamples | ProcessFilteredSamples
}

