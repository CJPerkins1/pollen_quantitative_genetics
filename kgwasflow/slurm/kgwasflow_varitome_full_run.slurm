#!/bin/bash
#SBATCH --job-name=kgwasflow_varitome_full_run
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400gb
#SBATCH --time=30:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cjperkins1@arizona.edu
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# Sourcing my .bashrc so that conda works and the path is how it should be
source /home/u28/cjperkins1/.bashrc

# Activating the conda environment
conda activate kgwasflow

# Getting to working directory
cd "/xdisk/rpalaniv/cjperkins1/kmers-gwas/kgwasflow/varitome"

# unlock the directory

# kgwasflow run -t 64 -c /xdisk/rpalaniv/cjperkins1/kmers-gwas/kgwasflow/varitome/config/config.yaml --unlock 

# Run kGWASflow with the default config file default arguments and 16 threads
kgwasflow run \
	-t 64 \
	-c /xdisk/rpalaniv/cjperkins1/kmers-gwas/kgwasflow/varitome/config/config.yaml \
	--work-dir /xdisk/rpalaniv/cjperkins1/kmers-gwas/kgwasflow/varitome

