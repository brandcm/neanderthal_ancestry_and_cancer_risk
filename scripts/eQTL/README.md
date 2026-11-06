# eQTL

Scripts to retrieve GTEx expression quantitative trait locus (eQTL) results for the *PTHLH* and *SLC24A4* credible set variants, for the eQTL figures and supplementary table.

## Configuration

`eqtl_config.yaml` lists, for each locus, the variants to query, the allele of interest at each (the Neanderthal allele, or the reintroduced ancestral allele for rs788464), and the focus genes whose effects are retrieved in every tissue. All effects in the outputs are oriented to the allele of interest: a positive normalized effect size (NES) means that allele increases expression. For *SLC24A4*, rs79342051 (Neanderthal allele T) is the variant plotted; rs12435365 is included only as a check, since its Neanderthal allele is inferred as the derived allele rather than called by the maps.

## Running

The GTEx API needs internet access, which Wynton compute nodes lack, so run the script on a dev node from this directory:

```bash
ssh dev1
cd /wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/scripts/eQTL
python3 retrieve_gtex_eqtls.py \
    --config eqtl_config.yaml \
    --dataset gtex_v10 \
    --output_directory ../../data/eQTL/gtex_v10
```

The figures use GTEx v10, whose sample sizes match the earlier rs788464 panel (n = 816 in skeletal muscle; v8 has 706). Write each release to its own folder.

It makes a few hundred API requests and takes a few minutes. Pairs that GTEx did not test (for example, a variant below its MAF filter in a tissue) are reported on stderr and left as NA.

## Outputs

Written to `data/eQTL/<release>/`:

- `gtex_significant_eqtls.txt`: every gene-tissue pair in which a queried variant is a significant single-tissue eQTL, with its NES, SE and p-value, the GTEx lead variant for that gene and tissue, and whether the queried variant is that lead. This is the basis of the supplementary table. For pairs where the queried variant is not the lead, check its LD with the lead (e.g. LDlink LDpair, EUR) before treating them as the same signal.
- `gtex_focus_gene_effects.txt`: the effect of each queried variant on each focus gene in every GTEx tissue with eQTL results, significant or not, for the tissue-specificity panel.
- `gtex_expression_by_genotype.txt.gz`: normalized expression per sample, by copies of the allele of interest, for every significant pair, for genotype-level plots.

GTEx does not report the standard error, so SE is computed as NES / t.
