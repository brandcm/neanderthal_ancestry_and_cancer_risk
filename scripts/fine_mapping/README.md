# Fine-mapping

Scripts to fine-map cancer GWAS loci containing genome-wide significant Neanderthal variants with SuSiE-RSS, select one fit per locus, collapse duplicate credible sets, and annotate the results. All coordinates are hg19.

## Inputs

- `data/metadata/cancer_GWAS_metadata.yaml`: per-GWAS metadata (ancestry group, sample sizes, column names).
- The significant introgressed variants from `significant_GWAS_variants/`.
- The 1000 Genomes EAS and EUR LD reference panels from `LD_reference_panels/`.
- For the annotation steps: the GENCODE gene models and the Iasi et al. 2025 adaptively introgressed regions from `annotations/`, and the reintroduced ancestral allele calls in `data/reintroduced_ancestral_alleles/`.

## Pipeline

`submit.sh` runs the Snakemake workflow in `Snakefile` on the cluster. Submit it from this directory:

```bash
mkdir -p logs
conda activate snakemake
snakemake -n                                   # dry run: lists the jobs without running them
qsub submit.sh
```

R steps run through `run_rscript.sh`, which switches to the `rgenomics` conda environment so the pipeline and the notebooks use one R installation. The workflow:

1. identifies loci with `identify_fine_mapping_loci.py`: the most significant Neanderthal variant in each GWAS is the index of the first locus, every Neanderthal variant within 500 kb joins it, and the next most significant remaining variant starts the next (`fine_mapping_loci.txt`);
2. extracts each GWAS's locus windows once with `extract_locus_windows.R`, at twice the half-window so a locus can be rerun wider without re-reading the summary statistics;
3. fine-maps each locus with `fine_map_locus.R`, which calls `fine_mapping()` in `fine_mapping.R` (variant QC, allele alignment to the panel, PLINK LD, SuSiE-RSS and the kriging diagnostic);
4. collects the per-locus results with `summarize_fine_mapping_results.R`;
5. refits loci that did not converge with `refit_nonconverged.R` (refine = TRUE, L = 5 and L = 10), skipping loci listed in `refit_exclusions.txt`, and compares each refit with the original with `compare_refits.R`;
6. chooses one fit per locus with `select_fits.R` (`fit_status.txt` and the `_selected` tables);
7. merges credible sets found through overlapping index windows (Jaccard ≥ 0.5) with `deduplicate_credible_sets.R` (the `_deduplicated` tables);
8. names the genes overlapping each credible set and the gene nearest its top-PIP variant with `annotate_credible_sets.R` (`credible_set_genes.txt`);
9. intersects credible set variants with the adaptively introgressed regions with `intersect_adaptive_introgression.R`.

Per-locus settings that differ from the defaults (a wider window, or a different L) are read from `locus_overrides.txt`. `update_locus_overrides.R` rebuilds it after a round of fine-mapping from the loci whose credible sets sit near a window edge, holding loci in the MHC at 1 Mb.

The fine-mapping pipeline writes `<name>_current.<ext>` beside a table after a partial rerun; the notebooks prefer those files when they exist.

## Run separately

These are run by hand after the pipeline, through `run_rscript.sh`:

- `annotate_reintroduced_ancestral.R` flags credible set variants that are reintroduced ancestral alleles (`reintroduced_ancestral_variants.txt`). `intersect_reintroduced_ancestral.R` does the same for individual loci, with a check that the calls cover the window.
- `condition_locus.R` runs the approximate conditional analysis at one locus, conditioning every variant on one variant: z_cond = (z − r·z_c) / √(1 − r²). It writes `data/fine_mapping/conditional/<gwas>_<locus>_<ancestry>_conditioned_on_<rsID>.txt`, which the analysis notebook reads. `make_conditional_gwas.R` converts that output into the GWAS-shaped file the locuszoom plots read (`pthlh_breast_<ancestry>_conditioned.txt.gz`), so both come from the same run.
- `regenerate_plink_for_figures.R` rebuilds the PLINK LD files for the loci in the figures, since `fine_mapping()` deletes them after each fit. It needs `module load CBI plink/1.90b6.26`.
- `examine_locus.R` prints one locus's credible sets and the LD between the strongest Neanderthal variant and the strongest variant overall in each set.
- `simulate_fine_mapping.R` is a simulation test of `fine_mapping()`.

## Outputs

Written to `data/fine_mapping/`. The tables later steps read:

- `fine_mapping_loci.txt`, `fit_status.txt`, `locus_summary.txt`;
- `credible_set_summary` and `credible_set_variants`, each in raw, `_selected` and `_deduplicated` versions;
- `neanderthal_majority_signals.txt`, `kriging_flagged_variants.txt`;
- `credible_set_genes.txt`, `adaptive_introgression_overlap.txt`, `adaptive_introgression_variants.txt`, `reintroduced_ancestral_variants.txt`;
- `results/` (one `.rds` per locus), `refits/` and `conditional/`.
