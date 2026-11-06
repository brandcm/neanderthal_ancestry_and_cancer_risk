# Notebooks

Two notebooks that read the pipeline outputs on Wynton, report every number quoted in the manuscript, and draw the figures. Their sections are numbered alike and follow the order of the manuscript.

- `analysis.ipynb` (Python) builds the supplementary tables from the pipeline outputs, runs the analyses that are not part of a pipeline (the map overlap, the *PTHLH* conditional analysis and reintroduced ancestral alleles, the ancient DNA allele frequencies and the heritability QC), and checks every manuscript-quoted value against the outputs in its final section. Tables it writes go to `data/supplementary_tables/` and are named by content rather than table number.
- `figures.ipynb` (R) draws the figures from the pipeline outputs and the tables `analysis.ipynb` writes, so run `analysis.ipynb` first. Main-figure components go to `figures/main_fig_components/` and supplementary figures to `figures/suppl_figs/`; panels marked as placeholders are completed in Inkscape.
- `utils.py` holds helpers shared by the cells of `analysis.ipynb`.

`figures.ipynb` also sources three helper scripts for the locuszoom and credible set PIP panels: `generate_locuszoom_plot.R`, `locuszoom_helpers.R` and `plot_credible_set_pips.R`. These are expected in this directory on Wynton.

## Environments

- `analysis.ipynb`: Python 3 with pandas, numpy, pybedtools (and bedtools on the PATH) and PyYAML.
- `figures.ipynb`: R with the tidyverse packages loaded at the top of the notebook, cowplot, patchwork, ComplexUpset and ggtext. Each notebook prints its package versions in its last section.

## Not included

The GENIE heritability enrichment estimates for breast cancer come from UK Biobank individual-level genotypes (application 226956). UK Biobank's access terms do not permit sharing individual-level data, so those analyses are not reproduced here; the summary estimates are entered directly in `figures.ipynb`.
