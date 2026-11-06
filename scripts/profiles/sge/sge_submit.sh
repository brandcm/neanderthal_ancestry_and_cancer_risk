#!/bin/bash
# Snakemake cluster submit command for SGE.
# Usage (from the profile): sge_submit.sh <qsub options...> <jobscript>
# Runs the jobscript through sge_run.sh, which records the job's exit code so that
# sge_status.sh never needs qacct (slow on Wynton). Prints the SGE job ID.
profile_dir="$(cd "$(dirname "$0")" && pwd)"
jobscript="${@: -1}"
qsub_options=("${@:1:$#-1}")

qsub -terse "${qsub_options[@]}" "$profile_dir/sge_run.sh" "$jobscript"
