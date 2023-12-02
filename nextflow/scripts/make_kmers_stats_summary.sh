#!/bin/bash

output_file="kmers_stats_summary.tsv"

# Header 
echo -e "accession\ttotal_reads\tn_non_canon_kmers\tn_unique_canonized_kmers" > "$output_file"

# Associative array to hold data
declare -A data

# Get data from all kmers
for file in *_kmc_all.1; do
    accession=$(basename "$file" | cut -c 1-6)
    total_reads=$(grep "Total no. of reads" "$file" | awk -F': ' '{print $2}')
    n_non_canon_kmers=$(grep "No. of unique k-mers" "$file" | awk -F': ' '{print $2}')

    data["$accession"]="${total_reads}\t${n_non_canon_kmers}"
done

# Get data from canonized kmers
for file in *_kmc_canon.1; do
    accession=$(basename "$file" | cut -c 1-6)
    n_unique_canonized_kmers=$(grep "No. of unique counted k-mers" "$file" | awk -F': ' '{print $2}')

    # Append to the existing data in the associative array
    if [ -n "${data[$accession]}" ]; then
        data["$accession"]="${data[$accession]}\t${n_unique_canonized_kmers}"
    fi
done

# Output
for accession in "${!data[@]}"; do
    echo -e "${accession}\t${data[$accession]}" >> "$output_file"
done

echo "Extraction complete. Data saved in $output_file."

