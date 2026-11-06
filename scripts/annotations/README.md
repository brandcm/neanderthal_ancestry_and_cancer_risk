# Annotations

Scripts to retrieve and prepare the annotation files used to interpret fine-mapped credible sets: gene models, ENCODE candidate cis-regulatory elements (cCREs), reintroduced ancestral alleles and adaptively introgressed regions. Run these before the annotation steps of `fine_mapping/` and before the notebooks. All outputs are on hg19.

## 1. Download

`retrieve_annotations.sh` downloads into `data/annotations/` and `data/reintroduced_ancestral_alleles/`. Run it on a dev node, since compute nodes have no internet access:

```bash
ssh dev1
./retrieve_annotations.sh
```

| File | Source | Used by |
|---|---|---|
| `gencode.v50lift37.basic.annotation.gff3.gz` | GENCODE release 50, basic set, GRCh37 mapping | gene tracks in the locuszoom figures; nearest and overlapping genes (`fine_mapping/annotate_credible_sets.R`) |
| `GRCh38-cCREs.bed` | ENCODE V3 cCREs, all cell types and tissues | cCRE overlap, after lifting |
| `Rinker_Suppl.File1.tar.gz` | Rinker et al. 2020, Supplementary File 1 | reintroduced ancestral alleles |

## 2. Lift the cCREs to hg19

`lift_regulatory_annotations.sh` lifts the cCREs to hg19 with liftOver for the chromosomes of the locus examples (chr1, chr9, chr12 and chr14), writing `GRCh37-cCREs_<chromosome>.bed` with every column of the original file. Extra columns are carried through liftOver in the name field and split back out afterwards. Intervals that fail to lift are written to `GRCh37-cCREs_<chromosome>.unlifted`.

```bash
./lift_regulatory_annotations.sh
```

## 3. Prepared by hand

- **Reintroduced ancestral alleles.** The fine-mapping scripts read one file per ancestry group, `data/reintroduced_ancestral_alleles/<EAS|EUR>_reintroduced_ancestral_alleles.txt.gz`, headerless, with two tab-separated columns (chromosome and hg19 position), built from the calls in the Rinker et al. 2020 supplement.
- **Adaptively introgressed regions.** `data/adaptive_introgression/Iasi_et_al_2025_high_frequency_ancient_modern_regions.bed` (chr, start, end; hg19) is built from Table S26 of Iasi et al. 2025, downloaded from the article's supplementary materials, which do not allow scripted download.
