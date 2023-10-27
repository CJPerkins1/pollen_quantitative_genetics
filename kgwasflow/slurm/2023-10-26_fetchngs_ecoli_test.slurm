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

cd '/xdisk/rpalaniv/cedar/kmers-gwas/genomes/varitome'

# Running Nextflow
nextflow run nf-core/fetchngs \
  -c '/home/u16/cedar/git/hpc/kmers-gwas/config/2023-05-08_fetchngs.config' \
  -r dev \
  -profile singularity \
  --input '/home/u16/cedar/git/hpc/kmers-gwas/sample_files/2023-10-23_ecoli_test.tsv' \
  --outdir '/xdisk/rpalaniv/cedar/kmers-gwas/genomes/ecoli_test'

