# Scripts

Scripts for the Neanderthal Ancestry and Cancer Risk project, one subdirectory per stage. Each subdirectory has a README describing its inputs, how to run it and what it writes. All analyses use hg19 coordinates unless noted.

The scripts run on Wynton, UCSF's SGE cluster, with data in `/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/data/`, which is not part of this repository. Pipelines are Snakemake workflows submitted with each directory's `submit.sh`, using the cluster profile in `profiles/`. Steps that download data run on a dev node, since compute nodes have no internet access.

## Order

Run the stages in this order; each reads the outputs of the ones before it.

1. `introgressed_variants/` retrieves Neanderthal-introgressed variants and regions from four maps (Browning et al. 2018, Sankararaman et al. 2016, Skov 2024 and Vernot et al. 2016), merges them into one table and annotates them with gnomAD v4.1 minor allele frequencies.
2. `GWAS_summary_stats/` retrieves the cancer GWAS summary statistics, maps rsIDs, lifts coordinates to hg19 and harmonizes effect sizes.
3. `LD_reference_panels/` builds the 1000 Genomes Phase 3 EAS and EUR LD reference panels.
4. `significant_GWAS_variants/` extracts each GWAS's genome-wide significant variants, identifies the Neanderthal ones, and LD-clumps both sets into loci.
5. `annotations/` retrieves gene models, ENCODE cCREs and reintroduced ancestral allele calls, and lifts the cCREs to hg19.
6. `fine_mapping/` fine-maps the loci containing Neanderthal variants with SuSiE-RSS and annotates the credible sets.
7. `LD_score_regression/` estimates heritability and its enrichment in Neanderthal variants with stratified LDSC. It needs only stages 1 and 2.
8. `eQTL/` retrieves GTEx eQTL results for the *PTHLH* and *SLC24A4* credible set variants.
9. `notebooks/` builds the supplementary tables, runs the remaining analyses, checks every manuscript-quoted value, and draws the figures. Run `analysis.ipynb` before `figures.ipynb`.

## `structural_variants/`

- `analyze_PTHLH_insertion.py` — Characterizes a polymorphic LINE-1 insertion near *PTHLH* from long-read genome assemblies (578 haplotypes; personal communication). Classifies haplotypes by size, locates breakpoints with k-mer anchoring, extracts the inserted sequence, identifies target site duplication and orientation. Input data (personal communication) are not redistributed. Usage:

  ```
  python3 scripts/structural_variants/analyze_PTHLH_insertion.py \
      --fasta-dir data/structural_variants/PTHLH_long_reads/fasta \
      --reference data/structural_variants/PTHLH_long_reads/reference/PTHLH.t2t.fa \
      --outdir data/structural_variants \
      --exemplar HG00512_h1
  ```

`profiles/` holds the Snakemake profile for SGE that every pipeline uses.

## Environments

Pipelines run from conda environments activated in each `submit.sh`: `snakemake` for most workflows, `rgenomics` for R (entered through `fine_mapping/run_rscript.sh`) and `ldsc` for LD score regression. PLINK 1.9 is loaded with `module load CBI plink/1.90b6.26`.
