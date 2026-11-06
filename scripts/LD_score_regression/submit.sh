#!/bin/bash
#$ -N LD_score_regression
#$ -cwd
#$ -o logs/
#$ -e logs/
#$ -l h_rt=24:00:00
#$ -l mem_free=4G

# load conda environment
# ldsc rather than snakemake: the LDSC scripts are python2-era code ported to
# python3 and pinned to this environment's pandas and scipy. -V in the profile's
# submit command exports it to the compute nodes, so the jobs inherit it.
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ldsc

mkdir -p logs

snakemake --profile ../profiles/sge "$@"
