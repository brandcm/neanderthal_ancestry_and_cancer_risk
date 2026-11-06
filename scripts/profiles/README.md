# Profiles

Snakemake configuration for running the pipelines on Wynton, which uses the SGE scheduler. Each pipeline's `submit.sh` passes `--profile ../profiles/sge`.

- `sge/config.yaml` submits each rule as an SGE job with the threads, memory (GB per thread) and run time the rule requests, runs up to 25 jobs at once, and reruns incomplete outputs.
- `sge/sge_submit.sh` submits a job through `sge_run.sh` and prints its job ID.
- `sge/sge_run.sh` runs the job and records its exit code.
- `sge/sge_status.sh` reports a job's state from that record, so Snakemake never calls `qacct`, which is slow on Wynton.

The profile's paths are absolute (`/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/scripts/profiles/sge/`), so a copy of the repository elsewhere needs them updated.
