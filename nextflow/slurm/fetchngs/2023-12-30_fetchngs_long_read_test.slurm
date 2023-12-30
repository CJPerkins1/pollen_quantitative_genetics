#!/usr/bin/env bash
#SBATCH --job-name=nf_fetchngs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --time=200:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# Sourcing my .bashrc so that conda works and the path is how it should be
source /home/u16/cedar/.bashrc

# Activating the conda environment
conda activate nf-core

cd '/xdisk/rpalaniv/cedar/kmers-gwas/genomes/alonge'

# Running Nextflow
nextflow run nf-core/fetchngs \
  -c '/home/u16/cedar/git/hpc/kmers-gwas/config/2023-10-09_fetchngs_large_files.config' \
  -r dev \
  -profile singularity \
  -resume \
  --input '/home/u16/cedar/git/pollen_quantitative_genetics/nextflow/samplesheets/fetchngs/2023-12-30_small_alonge_genomic.tsv' \
  --outdir '/xdisk/rpalaniv/cedar/kmers-gwas/genomes/alonge'
