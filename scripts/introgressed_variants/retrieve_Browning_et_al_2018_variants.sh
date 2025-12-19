#!/bin/bash

# change directories
cd /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/data/introgressed_variants/

# ensure EUR directory exists
mkdir -p EUR_variants/

# assign variables
pops=('BEB' 'CDX' 'CEU' 'CHB' 'CHS' 'CLM' 'FIN' 'GBR' 'GIH' 'IBS' 'ITU' 'JPT' 'KHV' 'MXL' 'PEL' 'PJL' 'PUR' 'Papuans' 'STU' 'TSI' )
EUR_pops=('CEU' 'FIN' 'GBR' 'IBS' 'TSI' )

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

# concat variants across populations
for pop in "${pops[@]}"; do
	for f in mendeley_data/${pop}.chr*.ND_match; do
		[[ -e "$f" ]] || continue
		tail -n +2 "$f" | awk '{print $1,$2,$3,$4,$5,$7,$9 }' OFS='\t' | grep -v '^$' >> Sprime_introgressed_variants.tmp
	done
done

# concat variants across EUR populations
for pop in "${EUR_pops[@]}"; do
	for f in mendeley_data/${pop}.chr*.ND_match; do
		[[ -e "$f" ]] || continue
		tail -n +2 "$f" | awk '{print $1,$2,$3,$4,$5,$7,$9}' OFS='\t' | grep -v '^$' >> Sprime_EUR_introgressed_variants.tmp
	done
done

# generate files
awk 'NF {k=$1 OFS $2 OFS $4 OFS $5; if (!(k in seen) || ($3 ~ /^rs/ && seen[k] !~ /^rs/)) { seen[k]=$3; row[k]=$0 }} END {for (k in row) print row[k]}' OFS='\t' Sprime_introgressed_variants.tmp | sort -V -k1,1 -k2,2n > Browning_et_al_2018_introgressed_variants_hg19.txt
awk '{print $1,$2-1,$2,$3}' OFS='\t' Browning_et_al_2018_introgressed_variants_hg19.txt > Browning_et_al_2018_introgressed_variants_hg19.bed

awk 'NF {k=$1 OFS $2 OFS $4 OFS $5; if (!(k in seen) || ($3 ~ /^rs/ && seen[k] !~ /^rs/)) { seen[k]=$3; row[k]=$0 }} END {for (k in row) print row[k]}' OFS='\t' Sprime_EUR_introgressed_variants.tmp | sort -V -k1,1 -k2,2n > Browning_et_al_2018_EUR_introgressed_variants_hg19.txt
awk '{print $1,$2-1,$2,$3}' OFS='\t' Browning_et_al_2018_EUR_introgressed_variants_hg19.txt > Browning_et_al_2018_EUR_introgressed_variants_hg19.bed

mv Browning_et_al_2018_introgressed_variants_hg19.* ../
mv Browning_et_al_2018_EUR_introgressed_variants_hg19.* ../EUR_variants/

# clean up
cd ../ && rm -R Sprime_results