#!/bin/bash
#$ -N significant_GWAS_variants
#$ -cwd
#$ -o logs/
#$ -e logs/
#$ -l h_rt=24:00:00
#$ -l mem_free=4G

# load modules
module load CBI plink/1.90b6.26

# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

mkdir -p logs

snakemake --profile ../profiles/sge "$@"