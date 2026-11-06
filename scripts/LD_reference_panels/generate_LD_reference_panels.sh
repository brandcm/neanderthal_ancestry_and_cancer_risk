#!/bin/bash
#$ -N generate_LD_reference_panels
#$ -t 1-22
#$ -M colin.brand@ucsf.edu
#$ -m a
#$ -cwd
#$ -o logs/
#$ -e logs/
#$ -l h_rt=4:00:00
#$ -l mem_free=10G

set -euo pipefail

# load modules
module load CBI bcftools/1.21 plink2/2.0.0-a.6.9

# change into temporary directory
cd "$TMPDIR"

# assign variables
project_directory="/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/"
metadata_directory="${project_directory}data/metadata/"
output_directory="${project_directory}data/LD_reference_panels/"
chr="$SGE_TASK_ID"
vcf="/wynton/group/capra/data/1KG_phase_3/2025-11-05/1KG_phase_3_chr${chr}_genotypes_20130502.vcf.gz"
ancestry_groups=(EAS EUR)

# filter to biallelic SNVs
bcftools view -m2 -M2 -v snps "$vcf" -O z -o "chr${chr}_filtered.vcf.gz"

# generate one panel per ancestry group, keeping the first of any duplicated variants
for ancestry_group in "${ancestry_groups[@]}"; do
    mkdir -p "${output_directory}${ancestry_group}"

    # the variant ID template is single-quoted so that plink2, not the shell, expands $r and $a
    plink2 \
        --threads 1 \
        --vcf "chr${chr}_filtered.vcf.gz" \
        --keep "${metadata_directory}1000G_${ancestry_group}_sample_IDs.txt" \
        --set-missing-var-ids '@:#:$r:$a' \
        --rm-dup force-first \
        --make-bed \
        --out "${output_directory}${ancestry_group}/1000G_${ancestry_group}_chr${chr}"
done

