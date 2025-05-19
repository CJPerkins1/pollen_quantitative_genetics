#!/bin/bash
#SBATCH --job-name=kgwasflow_resume
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400gb
#SBATCH --time=35:00:00
#SBATCH --account=rpalaniv
#SBATCH --partition=standard
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cjperkins1@arizona.edu
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# Source bashrc for conda
source /home/u28/cjperkins1/.bashrc

# Activate conda environment
conda activate kgwasflow

# Navigate to working directory
cd "/xdisk/rpalaniv/cjperkins1/kmers-gwas/kgwasflow/varitome"

# Unlock the directory
kgwasflow run -t 64 -c config/config.yaml --unlock

# Resume the workflow
kgwasflow run -t 64 -c config/config.yaml --work-dir /xdisk/rpalaniv/cjperkins1/kmers-gwas/kgwasflow/varitome --rerun-incomplete

