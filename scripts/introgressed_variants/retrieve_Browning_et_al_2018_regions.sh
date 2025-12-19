#!/bin/bash

# change directories
cd /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/data/introgressed_variants/

# load bedtools
module load CBI
module load bedtools2/2.31.1

# ensure EUR directory exists
mkdir -p EUR_regions/

# assign variables
pops=('BEB' 'CDX' 'CEU' 'CHB' 'CHS' 'CLM' 'FIN' 'GBR' 'GIH' 'IBS' 'ITU' 'JPT' 'KHV' 'MXL' 'PEL' 'PJL' 'PUR' 'Papuans' 'STU' 'TSI')
EUR_pops=('CEU' 'FIN' 'GBR' 'IBS' 'TSI')

# download and extract data
wget -O Sprime_results.zip https://prod-dcd-datasets-cache-zipfiles.s3.eu-west-1.amazonaws.com/y7hyt83vxr-1.zip
unzip Sprime_results.zip && rm Sprime_results.zip
mv "Sprime results for 1000 Genomes non-African populations and SGDP Papuans" Sprime_results
cd Sprime_results

for file in *.tar.gz; do
	gzip -d "$file"
	base="${file%.tar.gz}"
	tar -xf "${base}.tar"
done

# temporary files for segments
mkdir -p ../all_regions_tmp ../EUR_regions_tmp

for pop in "${pops[@]}"; do
	for f in mendeley_data/${pop}.chr*.ND_match; do
		[[ -e "$f" ]] || continue
		awk '
			NR>1 {
				key = $1 OFS $6
				if (!(key in min) || $2 < min[key]) min[key] = $2
				if (!(key in max) || $2 > max[key]) max[key] = $2
			}
			END {
				for (k in min) {
					split(k, a, OFS)
					print a[1], min[k]-1, max[k], a[2]
				}
			}' OFS='\t' "$f" >> ../all_regions_tmp/all_regions.tmp
	done
done

# collect per-population min/max segments for EUR populations
for pop in "${EUR_pops[@]}"; do
	for f in mendeley_data/${pop}.chr*.ND_match; do
		[[ -e "$f" ]] || continue
		awk '
			NR>1 {
				key = $1 OFS $6
				if (!(key in min) || $2 < min[key]) min[key] = $2
				if (!(key in max) || $2 > max[key]) max[key] = $2
			}
			END {
				for (k in min) {
					split(k, a, OFS)
					print a[1], min[k]-1, max[k], a[2]
				}
			}' OFS='\t' "$f" >> ../EUR_regions_tmp/EUR_regions.tmp
	done
done

cd ../

# merge into unified region sets
sort -k1,1 -k2,2n all_regions_tmp/all_regions.tmp | \
	bedtools merge -i - > Browning_et_al_2018_introgressed_regions_hg19.bed

sort -k1,1 -k2,2n EUR_regions_tmp/EUR_regions.tmp | \
	bedtools merge -i - > EUR_regions/Browning_et_al_2018_EUR_introgressed_regions_hg19.bed

# clean up
rm -r all_regions_tmp EUR_regions_tmp Sprime_results
