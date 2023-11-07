#!/usr/bin/env bash
#SBATCH --job-name=kgwasflow_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=300gb
#SBATCH --time=20:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# Sourcing my .bashrc so that conda works and the path is how it should be
source /home/u16/cedar/.bashrc

# Activating the conda environment
conda activate kgwasflow

cd "/xdisk/rpalaniv/cedar/kmers-gwas/2023-11-07_varitome_test_rerun"

# Running the test
kgwasflow run \
	-t 64 \
	-c /home/u16/cedar/git/pollen_quantitative_genetics/kgwasflow/config/flowers_and_varitome_test/config.yaml \
	--work-dir /xdisk/rpalaniv/cedar/kmers-gwas/2023-11-07_varitome_test_rerun

