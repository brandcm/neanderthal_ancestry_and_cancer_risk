#!/bin/bash

# assign variables
downloads_directory="/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/data/GWAS_summary_stats/downloads/"
urls="/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/data/metadata/cancer_GWAS_summary_stats_URLs.txt"

# run
mkdir -p "$downloads_directory"

while read -r new_filename url; do
  target_file="${downloads_directory}${new_filename}"

  if [ ! -f "$target_file" ] && [ ! -f "${target_file}.gz" ]; then
    echo "Downloading $new_filename..."
    wget -O "$target_file" "$url"

    if [[ "$new_filename" != *.gz ]]; then
        echo "Compressing $new_filename..."
        gzip "$target_file"
    fi

  else
    echo "File ${new_filename} or ${new_filename}.gz already exists. Skipping download."
  fi
done < "$urls"
