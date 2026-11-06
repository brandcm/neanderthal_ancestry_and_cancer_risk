#!/bin/bash
#$ -N introgressed_variants
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o logs/
#$ -e logs/
#$ -l h_rt=24:00:00
#$ -l mem_free=4G

# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

# load modules
module load CBI bedtools2/2.31.1

# per-job SGE logs
mkdir -p logs

# run
snakemake --profile ../profiles/sge "$@"