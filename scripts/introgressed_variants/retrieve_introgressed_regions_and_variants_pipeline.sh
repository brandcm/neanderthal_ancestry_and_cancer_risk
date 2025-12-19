#!/bin/bash
#$ -N retrieve_introgressed_regions_and_variants_pipeline
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/scripts/introgressed_variants/retrieve_introgressed_regions_and_variants_pipeline.out
#$ -e /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/scripts/introgressed_variants/retrieve_introgressed_regions_and_variants_pipeline.err
#$ -l h_rt=8:00:00
#$ -l mem_free=20G

# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

# load modules
module load CBI
module load bedtools2/2.31.1

# run
snakemake --use-conda --executor cluster-generic --cluster-generic-submit-cmd "qsub -l h_rt={resources.time} -l h_vmem={resources.mem}G -V" -j 22

# remove Skov et al. 2020 temp files
if [[ $? -eq 0 ]]; then
  echo "Snakemake finished successfully — cleaning external tmp files"

  rm -f \
    /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/data/introgressed_variants/Skov_et_al_2020_introgressed_variants_hg19.tmp \
    /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/data/introgressed_variants/EUR_variants/Skov_et_al_2020_EUR_introgressed_variants_hg19.tmp
else
  echo "Snakemake failed — not cleaning tmp files" >&2
  exit 1
fi