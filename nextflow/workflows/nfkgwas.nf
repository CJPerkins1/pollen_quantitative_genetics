/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    k-mers-based GWAS workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Inputs
// Check the input
println("samplesheet: ${params.samplesheet}")

params.samplesheet = "$workflow.scriptDir/../samplesheets/test_samplesheet.tsv"
println("samplesheet: ${params.samplesheet}")
params.outdir = "." // Defaults to where the script is run

Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true, sep: '\t', strip: true)
    .set { samples_ch }

process FilterSamples {
    input:
    tuple val(accession_id), path(fq1), path(fq2)

    output:
    tuple val(accession_id), path(fq1), path(fq2) into filtered_samples_ch

    script:
    """
    # Here you would have some condition to filter the samples
    # For the sake of this example, we'll just pass all samples through
    echo "$accession_id" > /dev/null
    """
}

process ProcessFilteredSamples {
    input:
    tuple val(accession_id), path(fq1), path(fq2) from filtered_samples_ch

    script:
    """
    echo "Processing sample: $accession_id"
    echo "FASTQ1 file: $fq1"
    echo "FASTQ2 file: $fq2"
    """
}

workflow {
    take:
    ch_samples from samples_ch

    main:
    ch_samples | FilterSamples | ProcessFilteredSamples
}

