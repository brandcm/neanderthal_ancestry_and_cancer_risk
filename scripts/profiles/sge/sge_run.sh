#!/bin/bash
# Runs a Snakemake jobscript on the compute node and records its exit code in
# logs/job_status/<SGE job ID> (relative to the pipeline directory, via qsub -cwd).
bash "$1"
exit_code=$?

mkdir -p logs/job_status
echo "$exit_code" > "logs/job_status/$JOB_ID"
exit "$exit_code"
