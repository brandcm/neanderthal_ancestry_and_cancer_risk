# Introgressed variants

Scripts to retrieve Neanderthal-introgressed variants from four maps (Browning et al. 2018, Sankararaman et al. 2016, Skov 2024 and Vernot et al. 2016), combine them into one table, and annotate them with gnomAD minor allele frequencies. All outputs are on hg19.

## Source data

The [Sankararaman et al. 2016](https://drive.google.com/drive/folders/1DyhMw0E1mXQUDNeGQbQrcQbHe0uTyK2h) and [Vernot et al. 2016](https://drive.google.com/drive/folders/0B9Pc7_zItMCVWUp6bWtXc2xJVkk?resourcekey=0-Cj8G4QYndXQLVIGPoWKUjQ) data are in public Google Drive folders, which are difficult to download with curl or wget, so they were downloaded manually into the Capra Lab data directory, `/wynton/group/capra/data/`. The [Skov 2024](https://zenodo.org/records/14136628) data are also stored there.

## 1. Browning et al. 2018

`retrieve_Browning_et_al_2018_variants.sh` retrieves the Browning et al. 2018 variants (Neanderthal "match" calls only) into `data/introgressed_variants/`. Run it on a dev node before the pipeline; the pipeline uses its output but does not build it.

## 2. Pipeline

`retrieve_introgressed_regions_and_variants_pipeline.sh` runs the Snakemake workflow in `Snakefile` on the cluster. Submit it from this directory:

```bash
conda activate snakemake
snakemake -n                                   # dry run: lists the jobs without running them
qsub retrieve_introgressed_regions_and_variants_pipeline.sh
```

The workflow:

- extracts variants from each map: Sankararaman (posterior > 0.5), Skov (SNVs shared with a Neanderthal genome; lifted from hg38 with liftOver, carrying the ancestral and derived alleles) and Vernot (tag SNPs);
- builds Skov introgressed regions (Neanderthal or Both, posterior ≥ 0.8; lifted from hg38);
- adds rsIDs to the Sankararaman, Skov and Vernot variants with `rsid_mapping.py` (gnomAD v2.1.1, by hg19 position; when a position has several rsIDs the first is kept);
- merges the four maps into one table with `merge_introgressed_variants.py`;
- computes gnomAD v4.1.0 minor allele frequencies for all introgressed rsIDs with `quantify_introgressed_variant_mafs.py` (matched by rsID, since gnomAD v4.1 is GRCh38).

The job status script used by the launch script is shared across pipelines: `scripts/utils/sge_status.sh`.

## Outputs

In `data/introgressed_variants/`:

| File | Contents |
|---|---|
| `{map}_introgressed_variants_hg19.txt` / `.bed` | Variants from each map, with rsIDs |
| `intersection_introgressed_variants_hg19.bed` | Positions called by all four maps |
| `union_introgressed_variants_hg19.bed` | Positions called by any map (merged intervals) |
| `merged_introgressed_variants_hg19.txt` | One row per position across all maps (see below) |
| `Skov_2024_introgressed_regions_hg19.bed` | Skov 2024 introgressed regions |
| `gnomAD_v4.1.0_introgressed_variant_MAFs.txt.gz` | Minor allele frequencies overall and in AFR, EAS, FIN, NFE and SAS |
| `logs/*.unlifted` | liftOver's record of Skov positions that could not be converted |
| `tmp/` | Intermediate files; emptied by Snakemake after a successful run |

Per-job cluster logs are written to `logs/` in this scripts directory.

### Merged table

`merged_introgressed_variants_hg19.txt` has one row per hg19 position:

| Column | Description |
|---|---|
| `chr`, `pos_hg19` | Position (hg19) |
| `rsID` | rsID from any map that reports one; `NA` if none |
| `ref`, `alt` | hg19 reference and alternate alleles (from Browning or Vernot) |
| `alleles` | The variant's two alleles as an unordered pair, e.g. `A/G` (from any map that reports alleles) |
| `neanderthal_allele` | Neanderthal-derived allele, from Browning (archaic allele), Skov (derived allele) or Vernot (Neanderthal genotype; `NA` if heterozygous). `NA` if no map reports it or the maps disagree |
| `maps` | Maps that call the variant introgressed |
| `n_maps` | Number of maps |
| `notes` | Conflicts between maps (different rsIDs, alleles or Neanderthal alleles) |
