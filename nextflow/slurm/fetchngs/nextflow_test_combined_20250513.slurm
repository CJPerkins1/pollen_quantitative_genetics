#!/usr/bin/env bash
#SBATCH --job-name=nf_combined_read_test_20250512
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --time=10:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cjperkins1@arizona.edu

# Sourcing my .bashrc so that conda works and the path is how it should be
source /home/u28/cjperkins1/.bashrc

# Activating the conda environment
conda activate nf

cd '/xdisk/rpalaniv/cjperkins1/kmers-gwas/genomes/dp/test'

# Running Nextflow
nextflow run nf-core/fetchngs \
  -c '/xdisk/rpalaniv/cjperkins1/kmers-gwas/fetchngs/config/20250512_fetchngs_test.config' \
  -r dev \
  -profile singularity \
  -resume \
  --input '/xdisk/rpalaniv/cjperkins1/kmers-gwas/fetchngs/input/combined_reads.tsv' \
  --outdir '/xdisk/rpalaniv/cjperkins1/kmers-gwas/genomes/dp/test'
