 // First kmc run with canonization

 process KMC_COUNT_CANONIZED {
     conda "${projectDir}/../config/conda/kmc.yaml"
 
     tag "KMC_COUNT_CANONIZED on ${accession_id}"
 
     publishDir "${params.outdir}/results/kmc_count_one_canonized/${accession_id}", mode: 'symlink'
 
     input:
     tuple val(accession_id), val(file_maps), path(read_paths_file)
 
     output:
     tuple val(accession_id), path("${accession_id}_kmc_count_canonized.kmc_pre"), path("${accession_id}_kmc_count_canonized.kmc_suf"), path("${accession_id}_kmc_canon.1")
 
     script:
     println("KMC count canonized on accession: ${accession_id}, file list: ${read_paths_file}") // debugging
     // Run kmc command
     """
     kmc -t${task.cpus} -k31 -ci2 @${read_paths_file} ${accession_id}_kmc_count_canonized ./ 1> ${accession_id}_kmc_canon.1 2> ${accession_id}_kmc_canon.2
     """
 }
