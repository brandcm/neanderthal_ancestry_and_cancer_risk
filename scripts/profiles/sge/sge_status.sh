#!/bin/bash
# Snakemake cluster status command for SGE: prints running, success or failed for a job ID.
# Reads the exit code written by sge_run.sh; does not use qacct (slow on Wynton).
job_id="$1"
status_dir="logs/job_status"
status_file="$status_dir/$job_id"

# finished normally: exit code recorded
if [ -f "$status_file" ]; then
    if [ "$(cat "$status_file")" = "0" ]; then echo success; else echo failed; fi
    exit 0
fi

# still queued or running
if qstat -j "$job_id" &> /dev/null; then
    echo running
    exit 0
fi

# gone from the queue without an exit code: either the file is not visible yet
# (shared filesystem delay) or the job was killed (e.g. over h_rt or h_vmem).
# Wait 2 minutes from the first time it was seen gone, then call it failed.
mkdir -p "$status_dir"
gone_marker="$status_dir/$job_id.gone"
[ -f "$gone_marker" ] || date +%s > "$gone_marker"
if (( $(date +%s) - $(cat "$gone_marker") > 120 )); then
    echo failed
else
    echo running
fi
