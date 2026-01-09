#!/bin/bash

# assign variables
summary_stats_directory="/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/data/GWAS_summary_stats/"
pre_rsID_directory="${summary_stats_directory}/pre_rsIDs/"
urls="/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/data/metadata/cancer_GWAS_summary_stats_URLs.txt"

# run
mkdir -p "$summary_stats_directory"
mkdir -p "$pre_rsID_directory"

while read -r new_filename url; do
  target_file="${summary_stats_directory}${new_filename}"

  if [ ! -f "$target_file" ] && [ ! -f "${target_file}.gz" ]; then
    echo "Downloading $new_filename..."
    wget -O "$target_file" "$url"

  else
    echo "File ${new_filename} or ${new_filename}.gz already exists. Skipping download."
  fi
done < "$urls"

# move files that need rsID mapping to pre-rsID directory
mv "${summary_stats_directory}endometrial_cancer_EAS_Sakaue_GCST90018618.txt.gz" "$pre_rsID_directory"
mv "${summary_stats_directory}esophageal_cancer_EAS_Sakaue_GCST90018621.txt.gz" "$pre_rsID_directory"
mv "${summary_stats_directory}gastric_cancer_EAS_Sakaue_GCST90018629.txt.gz" "$pre_rsID_directory"
mv "${summary_stats_directory}liver_cancer_EAS_Sakaue_GCST90018583.txt.gz" "$pre_rsID_directory"
mv "${summary_stats_directory}lung_cancer_EAS_Sakaue_GCST90018655.txt.gz" "$pre_rsID_directory"
mv "${summary_stats_directory}ovarian_cancer_EAS_Sakaue_GCST90018668.txt.gz" "$pre_rsID_directory"
mv "${summary_stats_directory}pancreatic_cancer_EAS_Sakaue_GCST90018673.txt.gz" "$pre_rsID_directory"
mv "${summary_stats_directory}prostate_cancer_EAS_Wang_GCST90274716.txt" "$pre_rsID_directory"
mv "${summary_stats_directory}prostate_cancer_EUR_Wang_GCST90274714.txt" "$pre_rsID_directory"
mv "${summary_stats_directory}thyroid_cancer_EAS_Sakaue_GCST90018709.txt.gz" "$pre_rsID_directory"
