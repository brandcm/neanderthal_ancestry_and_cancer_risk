#!/bin/bash
#$ -N significant_GWAS_variants_pipeline
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/scripts/significant_GWAS_variants/significant_GWAS_variants_pipeline.out
#$ -e /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/scripts/significant_GWAS_variants/significant_GWAS_variants_pipeline.err
#$ -l h_rt=8:00:00
#$ -l mem_free=20G

# load modules
module load CBI
module load plink/1.90b6.26

# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

# run
snakemake --use-conda --executor cluster-generic --cluster-generic-submit-cmd "qsub -l h_rt={resources.time} -l h_vmem={resources.mem}G -V" -j 29

# remove LD clumping files
rm /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/data/significant_GWAS_variants/LD_clumped_significant_SNVs/*.log
rm /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/data/significant_GWAS_variants/LD_clumped_significant_SNVs/*.nosex