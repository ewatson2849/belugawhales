#!/bin/bash
#SBATCH --job-name=beluga_assembly
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --output=assembly.out
#SBATCH --error=assembly.err

module load miniconda
conda activate snakemake

snakemake -s assembly.smk -c4 --use-conda -conda-frontend mamba --rerun-incomplete