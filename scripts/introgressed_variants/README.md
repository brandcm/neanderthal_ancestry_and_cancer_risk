This directory contains scripts used to retrieve intrgrossed variant data. Please note that the Sankararaman et al. 2016 and Vernot et al. 2016 data used here are publicly available in Google Drive folders, which are difficult to curl/wget. Thus, these data were downloaded manually from https://drive.google.com/drive/folders/1DyhMw0E1mXQUDNeGQbQrcQbHe0uTyK2h and https://drive.google.com/drive/folders/0B9Pc7_zItMCVWUp6bWtXc2xJVkk?resourcekey=0-Cj8G4QYndXQLVIGPoWKUjQ, respectively, and then uploaded to the Capra Lab data directory: `/wynton/group/capra/data/`. Scripts are described below and listed in the order in which they should be run.

- `retrieve_Browning_et_al_2018_*.sh` retrieves the Browning et al. 2018 data.

- `retrieve_Skov_et_al_2020.sh` retrieves the Skov et al. 2020 data.

- `retrieve_introgressed_regions_and_variants_pipeline.sh` is a Snakemake pipeline, which implements `Snakefile` and `rsID_mapping.py` to retrieve the remaining data, perform rsID mapping, concatenate chromosome-level files, and generate an intersection and union of the four variant maps (Browning, Sankararaman, Skov, and Vernot).

These scripts should result in the following outputs in the `/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/data/introgressed_variants/` directory:
	Browning_et_al_2018_introgressed_regions_hg19.bed
	Browning_et_al_2018_introgressed_variants_hg19.bed
	Browning_et_al_2018_introgressed_variants_hg19.txt
	EUR_regions/
	EUR_variants/
	Skov_et_al_2020_introgressed_regions_hg19.bed
	Skov_et_al_2020_introgressed_variants_hg19.bed
	Skov_et_al_2020_introgressed_variants_hg19.txt
