#!/bin/bash
#$ -N GWAS_summary_stats_formatting_pipeline
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/scripts/GWAS_summary_stats/GWAS_summary_stats_formatting_pipeline.out
#$ -e /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/scripts/GWAS_summary_stats/GWAS_summary_stats_formatting_pipeline.err
#$ -l h_rt=8:00:00
#$ -l mem_free=20G

# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

# run
#snakemake --use-conda --cluster "qsub -l h_rt={resources.time} -l h_vmem={resources.mem}G -V" -j 21
snakemake --use-conda --executor cluster-generic --cluster-generic-submit-cmd "qsub -l h_rt={resources.time} -l h_vmem={resources.mem}G -V" -j 21