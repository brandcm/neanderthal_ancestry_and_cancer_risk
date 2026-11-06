# Data

Small data files used in the analysis. Larger pipeline outputs (fine-mapping results, LD reference panels, introgressed variant sets, LD scores, GWAS summary statistics, and genomic annotations) are deposited on [Zenodo](https://zenodo.org/records/24015068).

## `metadata/`

- `cancer_GWAS_metadata.yaml` — Per-study metadata for all cancer GWAS: cancer type, ancestry, author, accession, sample sizes and effective sample size.
- `cancer_GWAS_ldsc_sumstats_formatting_metadata.yaml` — Column mappings and parameters for munging each GWAS with LDSC's `munge_sumstats.py`.
- `cancer_GWAS_summary_stats_URLs.txt` — Download URLs for all GWAS summary statistics from the GWAS Catalog and other repositories.
- `1000G_EAS_sample_IDs.txt` / `1000G_EUR_sample_IDs.txt` — 1000 Genomes Phase 3 sample IDs for the EAS and EUR superpopulations, used to build LD reference panels.

## `Open_Targets/`

Open Targets Genetics credible set data for the BNC2 and SLC24A4 loci, used for comparison with the SuSiE fine-mapping results.

## `ancient_genomes/`

- `rs2026804_AADR_genotypes.csv` — Ancient DNA genotypes at rs2026804 from the Allen Ancient DNA Resource (AADR), used for the derived allele frequency trajectory analysis.
- `rs2026804_binned_DAFs.txt` — Binned derived allele frequencies over time.

## `eQTL/gtex_v10/`

GTEx v10 eQTL results for the credible set variants at the PTHLH and SLC24A4 loci.

- `gtex_significant_eqtls.txt` — Significant eQTLs for credible set variants.
- `gtex_focus_gene_effects.txt` — Effect sizes for focus genes across tissues.
- `gtex_expression_by_genotype.txt.gz` — Expression by genotype data for the eQTL figures.

## `structural_variants/`

Structural variant linkage disequilibrium tables and LINE-1 insertion analysis results for the BNC2 and PTHLH loci.

- `Table_S11_BNC2_structural_variant_LD.tsv` — LD between structural variants and credible set SNPs at BNC2 (SVatalog; Supplementary Table S11).
- `Table_S14_PTHLH_structural_variant_LD.tsv` — LD between structural variants and credible set SNPs at PTHLH (SVatalog; Supplementary Table S14).
- `insertion_summary.tsv` — Summary statistics for the polymorphic LINE-1 insertion near PTHLH (allele frequency, size, orientation, TSD).
- `haplotype_classification.tsv` — Per-haplotype classification (578 haplotypes from long-read genome assemblies).
- `HG00512_h1_insertion.fa` — Extracted LINE-1 insertion sequence from exemplar haplotype HG00512_h1.

Input long-read assemblies (personal communication) are not redistributed; see `PTHLH_long_reads/` in `.gitignore`.
