#!/usr/bin/env bash
#SBATCH --job-name=nfkgwas_kgwas_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10gb
#SBATCH --time=05:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cjperkins1@arizona.edu

# Sourcing my .bashrc so that conda works and the path is how it should be
source /home/u28/cjperkins1/.bashrc

# Activating the conda environment
conda activate nf-core

cd "/xdisk/rpalaniv/cjperkins1/kmers-gwas/nextflow/2025-06-16_nfkgwas_kgwas_test"

# Running Nextflow
nextflow run /xdisk/rpalaniv/cjperkins1/kmers-gwas/nextflow/workflows/nfkgwas.nf \
  -c "/home/u28/cjperkins1/git/pollen_quantitative_genetics/nextflow/config/2025-06-02_nfkgwas_combined_read.config" \
  --outdir "/xdisk/rpalaniv/cjperkins1/kmers-gwas/nextflow/2025-06-16_nfkgwas_kgwas_test" \
 -resume
