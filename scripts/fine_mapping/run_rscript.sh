#!/bin/bash
# Run an R script with the project's R.
#
# Snakemake runs from the snakemake environment, which has no R, and the notebook
# uses the rgenomics environment. Switching here keeps the pipeline and the notebook
# on one R installation and one package library.

set -euo pipefail

# conda's shell hook reads variables that may be unset, so relax -u around it
set +u
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rgenomics
set -u

exec Rscript "$@"
