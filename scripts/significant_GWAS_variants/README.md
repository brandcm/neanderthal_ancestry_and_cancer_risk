# Significant GWAS variants

Scripts to retrieve the genome-wide significant SNVs from each cancer GWAS, identify which of them are Neanderthal introgressed, and LD-clump both sets into independent loci. All inputs and outputs are on hg19.

## Source data

| Input | Contents |
|---|---|
| `data/GWAS_summary_stats/{gwas}.txt.gz` and `{gwas}_hg19.txt.gz` | Harmonized summary statistics; the `_hg19` files are used for the studies the metadata records as hg38 |
| `data/introgressed_variants/merged_introgressed_variants_hg19.txt` | Merged introgressed variant table |
| `data/LD_reference_panels/{EAS,EUR}/1000G_{EAS,EUR}_chr{chr}` | PLINK filesets used for clumping |
| `data/metadata/cancer_GWAS_metadata.yaml` | Each study's cancer, ancestry group, genome assembly and column names |

## Pipeline

`submit.sh` runs the Snakemake workflow in `Snakefile` on the cluster. Submit it from this directory:

```bash
conda activate snakemake
mkdir -p logs
snakemake -n                                   # dry run: lists the jobs without running them
qsub submit.sh
```

Arguments are passed through to Snakemake, so `qsub submit.sh --forceall` reruns everything.

The workflow:

- retrieves each study's significant SNVs with `retrieve_significant_gwas_variants.py`, keeping variants on the autosomes whose effect and other alleles are both a single base and whose p-value is below 5e-8;
- counts each study's SNVs and significant SNVs into one summary table;
- matches the significant SNVs to the merged introgressed variant table with `retrieve_significant_introgressed_gwas_variants.py`, first by rsID and then, for the variants left unmatched, by hg19 position where the alleles agree with those the introgression maps report;
- LD-clumps the significant SNVs and the significant introgressed SNVs separately with PLINK, one chromosome at a time, at r² thresholds of 0.2, 0.5 and 0.8 within 250 kb, then concatenates the chromosomes for each study and threshold.

The pipeline is restricted to the autosomes, for consistency with LDSC, which does not model chromosome X. Variants on other chromosomes are dropped when significant variants are retrieved, so every count and output here describes the same autosomal variant set, and LD reference panels are built for chromosomes 1 to 22 only.

The cluster profile is shared by the pipelines in this repository: `scripts/profiles/sge/`.

## Outputs

In `data/significant_GWAS_variants/`:

| File | Contents |
|---|---|
| `summaries/cancer_GWAS_N_SNVs_and_N_significant_SNVs.txt` | Number of SNVs and significant SNVs per study |
| `all_variants/variants/{gwas}_hg19_significant_SNVs.txt.gz` | Significant SNVs for one study |
| `all_variants/LD_clumped/{gwas}_r2_{threshold}.clumped` | Clumped significant SNVs for one study and threshold |
| `introgressed_variants/variants/significant_introgressed_GWAS_variants.txt.gz` | Significant introgressed SNVs across all studies (see below) |
| `introgressed_variants/variants/significant_introgressed_GWAS_variants_summary.txt` | Per-study counts of significant SNVs called by each introgression map |
| `introgressed_variants/LD_clumped/{gwas}_r2_{threshold}.clumped` | Clumped significant introgressed SNVs for one study and threshold |
| `tmp/` | Intermediate files; emptied by Snakemake after a successful run |

Per-job cluster logs are written to `logs/` in this scripts directory.

### Significant introgressed variants

`significant_introgressed_GWAS_variants.txt.gz` has one row per significant introgressed SNV per study:

| Column | Description |
|---|---|
| `cancer`, `ancestry_group` | Study the variant is significant in |
| `rsID`, `chr`, `pos` | Variant (hg19) |
| `effect_allele`, `other_allele`, `beta`, `p_value` | Association in that study |
| `neanderthal_allele` | Neanderthal allele from the merged table; `NA` where no map reports one or the maps disagree |
| `introgression_maps`, `n_maps` | Maps that call the variant introgressed, and how many |
| `match` | Whether the variant matched the introgressed table by `rsID` or by `position` |
| `notes` | Conflicts between maps, carried from the merged table |
