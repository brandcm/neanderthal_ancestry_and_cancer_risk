This directory contains scripts to retrieve and format cancer GWAS summary statistics. Scripts are described below and listed in the order in which they should be run.

- `retrieve_GWAS_summary_stats.sh` downloads summary statistics from the GWAS Catalog for cancer GWAS encompassing both East Asian and European ancestry individuals. Run this script on the command line rather than submitting to a compute node. Note that the Zhang et al. 2020 EUR breast cancer GWAS must be manually downloaded from [here](https://drive.google.com/file/d/1Co4OsQQlrNQuk7qf_K99Vc5HQuBF3Xph/view?usp=drive_link). It is not included in `cancer_GWAS_summary_stats_URLS.txt` because the data are housed in a Google Drive folder. Please note that the first 1017 data lines of this file are duplicated after the header. Run the following code block to fix:

```
head -n 1 icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt > header.txt
tail -n +1019 icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt > tail.txt
cat header.txt tail.txt > breast_cancer_EUR_Zhang_GCST001.txt
gzip breast_cancer_EUR_Zhang_GCST001.txt
rm header.txt icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt tail.txt
```

- `GWAS_summary_stats_formatting_pipeline.sh` is a Snakemake pipeline that runs the following scripts using `Snakefile`:

	- `format_EUR_breast_cancer_summary_stats.py` converts the Zhang et al. 2020 summary stats from a space-delimited to a tab-delimited file and maps rsIDs on to variant positions.

	- `prostate_cancer_summary_stats_rsID_mapping.py` reformats and maps rsIDs to genomic positions in both prostate cancer GWAS.

	- `Sakaue_et_al_2020_summary_stats_rsID_mapping.py` maps rsIDs to genomic positions for the the Sakaue et al. 2020 GWAS.