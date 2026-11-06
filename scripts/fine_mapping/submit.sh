#!/bin/bash
#$ -N fine_mapping
#$ -cwd
#$ -o logs/
#$ -e logs/
# the driver has to outlive the jobs it submits, and a refit starts at 24h and
# escalates to 48h then 72h on retry
#$ -l h_rt=72:00:00
#$ -l mem_free=4G

# load modules. The profile submits with -V, so the compute jobs inherit this PATH
module load CBI plink/1.90b6.26

# load conda environment. The R rules switch to the rgenomics environment
# themselves, through run_rscript.sh
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

# the output directories are created here rather than left to the compute jobs,
# which fail to create them on the group filesystem
fine_mapping_directory="/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/data/fine_mapping"
mkdir -p logs \
    "${fine_mapping_directory}/windows" \
    "${fine_mapping_directory}/results" \
    "${fine_mapping_directory}/refits"

# --keep-going so that one locus exhausting its memory or time does not stop the
# other 164, and --retries so that a locus that runs out of memory is retried with
# more, since fine_map_locus and refit_locus scale their requests with the attempt
# number. Snakemake resumes from what is already on disk, so if the driver's time
# runs out, resubmitting continues the run.
snakemake --profile ../profiles/sge --keep-going --retries 3 "$@"