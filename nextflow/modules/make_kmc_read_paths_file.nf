// Make files of file paths for kmc input

process MAKE_KMC_READ_PATHS_FILE {
    tag "MAKE_KMC_READ_PATHS_FILE on ${accession_id}"

    input:
    tuple val(accession_id), val(file_maps)

    output:
    tuple val(accession_id), val(file_maps), path("file_paths_${accession_id}.txt")

    exec:
    def file_list="${task.workDir}/file_paths_${accession_id}.txt"
    new File(file_list).withWriter { writer ->
        file_maps.each { writer.println(it.file) }
    }
}
