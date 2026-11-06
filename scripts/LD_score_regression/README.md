# LD score regression

Scripts to estimate total heritability and heritability enrichment in Neanderthal variants for each cancer GWAS with stratified LD score regression (LDSC).

## Inputs

- The formatted summary statistics from `GWAS_summary_stats/`.
- The four maps' introgressed variant BED files (hg19) from `introgressed_variants/`.
- `data/metadata/cancer_GWAS_ldsc_sumstats_formatting_metadata.yaml`: per-GWAS ancestry group, `munge_sumstats.py` flags and, where known, sample and population prevalence. `data/metadata/cancer_GWAS_metadata.yaml` supplies the case and control counts.
- LDSC and its reference data in the Capra Lab directory, `/wynton/group/capra/data/ldsc_new/`: the Python 3 port of LDSC (`munge_sumstats.py`, `make_annot.py`, `ldsc.py`), the ancestry-matched HapMap3 variant lists, 1000 Genomes PLINK filesets, baselineLD v2.2 LD scores and regression weights. These are not part of this repository.

## Pipeline

`submit.sh` runs the Snakemake workflow in `Snakefile` in the `ldsc` conda environment, which pins the pandas and scipy versions the LDSC port needs. Submit it from this directory:

```bash
mkdir -p logs
conda activate ldsc
snakemake -n                                   # dry run: lists the jobs without running them
qsub submit.sh
```

The workflow:

- formats each GWAS against its ancestry's HapMap3 variants with `munge_sumstats.py`;
- builds an annotation and LD scores for each map, ancestry group and autosome with `make_annot.py` and `ldsc.py --l2` (1 cM window);
- estimates each GWAS's total heritability, lambda GC, mean χ² and intercept under the baselineLD model alone (`ldsc.py --h2`), and collects them with `summarize_heritability_qc.py`;
- estimates the enrichment of each map, added as a single annotation to baselineLD v2.2, for each GWAS, and collects the results with `concat_heritability_results.py`.

Heritability is on the observed scale unless the metadata records both prevalences, in which case LDSC also reports the liability scale.

## Outputs

Written to `data/LD_score_regression/results/summaries/`:

- `GWAS_heritability_qc.txt`: one row per GWAS with the regression diagnostics and total heritability;
- `introgressed_variants_heritability_summary.txt`: one row per GWAS and map with the enrichment and its standard error.

The analysis notebook applies the QC thresholds (mean χ² ≥ 1.02, observed h² ≥ 0.05) and writes `heritability_qc.txt`, which the figures read.
