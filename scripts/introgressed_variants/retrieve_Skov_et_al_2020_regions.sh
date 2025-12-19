#!/bin/bash

# load modules
module load CBI
module load bedtools2/2.31.1

# assign variables
liftOver="/wynton/group/capra/bin/liftOver/liftOver"
hg38_hg19_chain="/wynton/group/capra/bin/liftOver/hg38ToHg19.over.chain.gz"
out_directory="/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/data/introgressed_variants/Skov_et_al_2020"

# change directories
mkdir -p "$out_directory"
cd "$out_directory"

wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2225-9/MediaObjects/41586_2020_2225_MOESM3_ESM.txt
awk 'NR > 1 && ($4 == "DAV") && $11 ~ /^(Altai|ambiguous|Vindija)$/ { print "chr"$1,$2-1,$3}' OFS='\t' 41586_2020_2225_MOESM3_ESM.txt > Skov_et_al_2020_introgressed_regions_hg38.tmp
"$liftOver" Skov_et_al_2020_introgressed_regions_hg38.tmp "$hg38_hg19_chain" Skov_et_al_2020_introgressed_regions_hg19.tmp Skov_et_al_2020_introgressed_regions_hg38.unlifted
awk '{sub(/^chr/, "", $1); print $1,$2,$3}' OFS='\t' Skov_et_al_2020_introgressed_regions_hg19.tmp | sort -V -k1,1 -k2,2n | bedtools merge -i - > ../Skov_et_al_2020_introgressed_regions_hg19.bed
cp ../Skov_et_al_2020_introgressed_regions_hg19.bed ../EUR_regions/Skov_et_al_2020_EUR_introgressed_regions_hg19.bed
rm -R "$out_directory"