# GWAS summary statistics

Scripts to retrieve and format cancer GWAS summary statistics for East Asian (EAS) and European (EUR) ancestry studies. Run the steps below in order.

## 1. Download

`retrieve_GWAS_summary_stats.sh` downloads summary statistics from the GWAS Catalog (URLs in `cancer_GWAS_summary_stats_URLS.txt`) into `data/GWAS_summary_stats/downloads/`. Run it on the command line rather than submitting it to a compute node.

The Zhang et al. 2020 EUR breast cancer GWAS is not in the GWAS Catalog and must be downloaded manually from [Google Drive](https://drive.google.com/file/d/1Co4OsQQlrNQuk7qf_K99Vc5HQuBF3Xph/view?usp=drive_link). The first 1,017 data lines of this file are duplicated after the header; from `downloads/`, remove them with:

```bash
input=icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt
output=breast_cancer_EUR_Zhang_GCST001_initial.txt

{ head -n 1 "$input"; tail -n +1019 "$input"; } > "$output"
gzip "$output"
rm "$input"
```

## 2. Format, Harmonize and liftOver

`GWAS_summary_stats_formatting_pipeline.sh` runs the Snakemake workflow in `Snakefile` on the cluster. Submit it from this directory:

```bash
conda activate snakemake
snakemake -n                                   # dry run: lists the jobs without running them
qsub GWAS_summary_stats_formatting_pipeline.sh
```

The workflow reads each study's column names and genome build from `data/metadata/cancer_GWAS_metadata.yaml` and uses these scripts:

- `format_eur_breast_cancer_summary_stats.py` converts the Zhang et al. 2020 summary statistics from space- to tab-delimited and adds an rsID column taken from the variant ID.
- `summary_stats_rsid_mapping.py` adds rsIDs by chromosome and position (gnomAD v2.1.1, hg19) for the Wang et al. prostate cancer GWAS and the Sakaue et al. GWAS. The prostate cancer files are sorted by chromosome and position first.
- `harmonize_effect_sizes.py` standardizes all studies to a `beta` (log odds ratio, derived from odds ratios where needed) and `standard_error` (derived from beta and the p-value where missing), drops rows without a p-value, and removes columns that are empty in every row.

Studies on hg38 (`assembly: hg38` in the metadata) are then lifted to hg19 with UCSC liftOver, in the Snakefile itself, so that positions match the 1000 Genomes LD panel and the introgressed variant maps.

The job status script used by the launch script is shared across pipelines: `scripts/utils/sge_status.sh`.

## Outputs

In `data/GWAS_summary_stats/`:

| File | Contents |
|---|---|
| `{gwas}.txt.gz` | Harmonized summary statistics for every study |
| `{gwas}_hg19.txt.gz` | hg19 version of each hg38 study |
| `logs/{gwas}_hg38.unlifted` | liftOver's record of positions that could not be converted |
| `tmp/` | Intermediate files; emptied by Snakemake after a successful run |

Per-job cluster logs are written to `logs/` in this scripts directory.
