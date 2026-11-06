# Neanderthal Ancestry and Cancer Risk

This repository houses data and scripts used in the Neanderthal Ancestry and Cancer Risk project.

## Data

- The `data/` directory contains small data files used in the analysis (metadata, Open Targets comparisons, ancient genomes, eQTL results, structural variant LD tables, and LINE-1 insertion analysis results). See `data/README.md` for details.

- Larger pipeline outputs — fine-mapping results, LD reference panels, introgressed variant sets, LD scores, GWAS summary statistics, and genomic annotations — are deposited on [Zenodo](https://zenodo.org/records/24015068).

- GWAS summary statistics are publicly available from the GWAS Catalog. Download URLs and accession numbers are in `data/metadata/cancer_GWAS_summary_stats_URLs.txt`.

- Analyses of UK Biobank individual-level data (GENIE heritability enrichment, application 226956) are not included, because UK Biobank's access terms do not permit sharing individual-level data. Summary estimates are in `scripts/notebooks/figures.ipynb`; code is available on request to researchers with approved UK Biobank access.

- Long-read assembly data for the PTHLH LINE-1 analysis were shared via personal communication and are not redistributed. The analysis script and results are included; see `scripts/structural_variants/analyze_PTHLH_insertion.py` and `data/structural_variants/`.

## Scripts

The `scripts/` directory contains all scripts used in the analysis for data generation, analysis, and visualization. Each subdirectory has a README describing its inputs, how to run it, and what it writes. See `scripts/README.md` for the order of operations.

## Contact

Colin M. Brand — colin.brand@ucsf.edu
John A. Capra — tony@capralab.org
