#!/bin/bash

# find $PWD -name "*fastq.gz" | sort > ecoli_fastq_paths.txt
filepaths="/xdisk/rpalaniv/cedar/kmers-gwas/genomes/ecoli_test/fastq/ecoli_fastq_paths.txt"

samples="/home/u16/cedar/git/pollen_quantitative_genetics/kgwasflow/config/ecoli_test/samples.tsv"
temp_samples="/home/u16/cedar/git/pollen_quantitative_genetics/kgwasflow/config/ecoli_test/samples_temp.tsv"

# Check if filepaths file exists
if [ ! -f "$filepaths" ]; then
  echo "Error: Filepaths file $filepaths does not exist."
  exit 1
fi

# Check if samples TSV file exists
if [ ! -f "$samples" ]; then
  echo "Error: Samples TSV file $samples does not exist."
  exit 1
fi

# Read each file path from the filepaths file
while IFS= read -r filepath; do
  echo "Processing filepath: $filepath"
  
  # Extract the SRR identifier and the file type (either _1 or _2) from the file path
  srr=$(echo "$filepath" | grep -o "SRR[0-9]\+" | xargs)
  type=$(echo "$filepath" | grep -o "_[12]\." | tr -d '.' | xargs)
  
  # Check if the SRR identifier and file type were successfully extracted
  if [ -z "$srr" ] || [ -z "$type" ]; then
    echo "Could not extract SRR identifier and/or file type from $filepath"
    continue
  fi
  
  echo "Extracted SRR: $srr, Type: $type"
  
  # Update the samples TSV file
  awk -v srr="$srr" -v type="$type" -v filepath="$filepath" -F'\t' -v OFS='\t' '
    BEGIN { updated = 0 }
    { gsub(/^[ \t]+|[ \t]+$/, "", $5) }  # Trim whitespaces from SRR field
    $5 == srr {
      if (type == "_1") $3 = filepath;
      else if (type == "_2") $4 = filepath;
      updated = 1
    }
    { print }
    END { if (updated == 0) print "SRR not found in TSV: " srr > "/dev/stderr" }
  ' "$samples" > "$temp_samples" && mv "$temp_samples" "$samples"
  
done < "$filepaths"

echo "TSV file updated successfully."
