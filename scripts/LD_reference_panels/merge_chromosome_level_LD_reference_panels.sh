#!/bin/bash
#$ -N merge_chromosome_level_LD_reference_panels
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/scripts/LD_reference_panels/merge_chromosome_level_LD_reference_panels.out
#$ -e /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/scripts/LD_reference_panels/merge_chromosome_level_LD_reference_panels.err
#$ -l h_rt=2:00:00
#$ -l mem_free=10G

# load modules
module load CBI
module load plink/1.90b6.26

# assign variables
LD_reference_panels_directory="/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/data/LD_reference_panels/"
output_directory="/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/data/LD_reference_panels/merged/"

# run
mkdir -p "$output_directory"

for ancestry_group in EAS EUR; do
	anchor_prefix="${LD_reference_panels_directory}${ancestry_group}/1000G_${ancestry_group}_chr1_nodups"

	autosomes_merge_list="${output_directory}${ancestry_group}_autosomes_merge_list.txt" > "$autosomes_merge_list"

	for chr in {2..22}; do
		bed="${LD_reference_panels_directory}${ancestry_group}/1000G_${ancestry_group}_chr${chr}_nodups.bed"
		bim="${LD_reference_panels_directory}${ancestry_group}/1000G_${ancestry_group}_chr${chr}_nodups.bim"
		fam="${LD_reference_panels_directory}${ancestry_group}/1000G_${ancestry_group}_chr${chr}_nodups.fam"

		if [[ -f "$bed" && -f "$bim" && -f "$fam" ]]; then
			echo "${LD_reference_panels_directory}${ancestry_group}/1000G_${ancestry_group}_chr${chr}_nodups" >> "$autosomes_merge_list"
		else
			echo "Skipping $ancestry_group chr${chr}, missing files"
		fi
	done

	plink --bfile "$anchor_prefix" \
		  --merge-list "$autosomes_merge_list" \
		  --make-bed \
		  --out "${output_directory}1000G_${ancestry_group}_autosomes_merged"
done