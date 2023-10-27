#!/usr/bin/env bash
#SBATCH --job-name=kgwasflow_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40gb
#SBATCH --time=8:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# Sourcing my .bashrc so that conda works and the path is how it should be
source /home/u16/cedar/.bashrc

# Activating the conda environment
conda activate kgwasflow

cd '/xdisk/rpalaniv/cedar/kmers-gwas/ecoli_test_dir'

# Running the test
kgwasflow test -t 16 --work-dir /xdisk/rpalaniv/cedar/kmers-gwas/ecoli_test_dir
