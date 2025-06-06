//delete raw fastq files after they are used to conserve disk space
process CLEANUP_FASTQ {
    cache = false
    tag "${accession_id}"

    input:
    tuple val(accession_id), val(fq_files_list)

    script:
    """
    echo "Deleting FASTQ files for ${accession_id}"

    for fq_file in ${fq_files_list.join(' ')}; do
        echo "Deleting \$fq_file"
        rm -f "\$fq_file"
    done
    """
}
