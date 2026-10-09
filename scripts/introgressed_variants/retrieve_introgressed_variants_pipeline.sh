#!/bin/bash
#$ -N retrieve_introgressed_variants_pipeline
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/scripts/introgressed_variants/retrieve_introgressed_variants_pipeline.out
#$ -e /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/scripts/introgressed_variants/retrieve_introgressed_variants_pipeline.err
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
snakemake \
	--executor cluster-generic \
	--cluster-generic-submit-cmd "qsub -terse -cwd -V -N {rule} -o logs/ -e logs/ -pe smp {threads} -l h_rt={resources.time} -l h_vmem={resources.mem}G -l mem_free={resources.mem}G" \
	--cluster-generic-status-cmd /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/scripts/utils/sge_status.sh \
	--cluster-generic-cancel-cmd qdel \
	--jobs 22 \
	--scheduler greedy \
	--latency-wait 60 \
	--rerun-incomplete