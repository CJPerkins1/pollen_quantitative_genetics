#!/usr/bin/env bash
#SBATCH --job-name=nfkgwas_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --time=2:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# Sourcing my .bashrc so that conda works and the path is how it should be
source /home/u16/cedar/.bashrc

# Activating the conda environment
conda activate nf-core

cd "/xdisk/rpalaniv/cedar/kmers-gwas/2023-11-13_nfkgwas_test"

# Running Nextflow
nextflow run /home/u16/cedar/git/pollen_quantitative_genetics/nextflow/workflows/nfkgwas.nf \
  -c "/home/u16/cedar/git/pollen_quantitative_genetics/nextflow/config/2023-11-13_nfkgwas.config" \
#  -r dev \
#  -profile singularity \
  --samplesheet "/home/u16/cedar/git/hpc/kmers-gwas/sample_files/2023-10-23_ecoli_test.tsv" \
  --outdir "/xdisk/rpalaniv/cedar/kmers-gwas/2023-11-13_nfkgwas_test"

