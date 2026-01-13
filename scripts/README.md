This directory contains subdirectories that house the scripts used for this project. Subdirectories are described below in the order in which the scripts should be run.

- `introgressed_variants/` contains scripts used to retrieve and format the following introgressed variant sets: Browning et al. 2016, Sankararaman et al. 2016 (set 1), Skov et al. 2020, and Vernot et al. 2016. The scripts ultimately generate a .bed and .txt file per set and ensures that as many variants include rsIDs and have genomic positions using the GRCh37/hg19 reference.

- `GWAS_summary_stats/` contains scripts that retrieve and format GWAS summary stats for downstream analyses.

- `LD_reference_panels/` contains scripts that generate chromosome-level and merged LD reference panels from Thousand Genomes Phase 3 EAS and EUR samples.

- `significant_GWAS_variants/` contains scripts to extract and process significant GWAS variants and perform LD-clumping to identify independent loci.

- `significant_GWAS_loci_and_introgressed_variants_intersection/` contains a script to run an intersection between the identified significant GWAS loci and each of the four individual introgressed variants set as well as the intersection and union of those sets.