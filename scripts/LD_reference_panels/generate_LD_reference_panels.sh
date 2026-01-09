#!/bin/bash
#$ -N generate_LD_reference_panels
#$ -t 1-22
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/scripts/LD_reference_panels/generate_LD_reference_panels.out
#$ -e /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/scripts/LD_reference_panels/generate_LD_reference_panels.err
#$ -l h_rt=2:00:00
#$ -l mem_free=10G

# load modules
module load CBI
module load bcftools/1.21
#module load plink/1.90b6.26
module load plink2/2.0.0-a.6.9

# change into temporary directory
cd "$TMPDIR"

# assign variables
project_directory="/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/"
output_directory="${project_directory}data/LD_reference_panels/"
chr="$SGE_TASK_ID"
VCF="/wynton/group/capra/data/1KG_phase_3/2025-11-05/1KG_phase_3_chr${chr}_genotypes_20130502.vcf.gz"
EAS_samples="${output_directory}1000G_EAS_sample_IDs.txt"
EUR_samples="${output_directory}1000G_EUR_sample_IDs.txt"

# run
mkdir -p "${output_directory}EAS" "${output_directory}EUR"

# filter to biallelic SNPs
bcftools view -m2 -M2 -v snps "$VCF" -O z -o "chr${chr}_filtered.vcf.gz"
bcftools index -t "chr${chr}_filtered.vcf.gz"

# generate EAS panel and remove duplicates
plink2 \
  --threads 1 \
  --vcf "chr${chr}_filtered.vcf.gz" \
  --make-bed \
  --keep "$EAS_samples" \
  --set-missing-var-ids @:#:$r:$a \
  --keep-allele-order \
  --out "${output_directory}EAS/1000G_EAS_chr${chr}"

plink2 \
  --threads 1 \
  --bfile "${output_directory}EAS/1000G_EAS_chr${chr}" \
  --make-bed \
  --rm-dup force-first \
  --out "${output_directory}EAS/1000G_EAS_chr${chr}_nodups"

# generate EUR panel and remove duplicates
plink2 \
  --threads 1 \
  --vcf "chr${chr}_filtered.vcf.gz" \
  --make-bed \
  --keep "$EUR_samples" \
  --set-missing-var-ids @:#:$r:$a \
  --keep-allele-order \
  --out "${output_directory}EUR/1000G_EUR_chr${chr}"

plink2 \
  --threads 1 \
  --bfile "${output_directory}EUR/1000G_EUR_chr${chr}" \
  --make-bed \
  --rm-dup force-first \
  --out "${output_directory}EUR/1000G_EUR_chr${chr}_nodups"
