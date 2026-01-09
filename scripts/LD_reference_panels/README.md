This directory contains scripts to generate LD reference panels used in downstream analyses. Scripts are described below and listed in the order in which they should be run.

- `generate_LD_reference_panels.sh` generates LD reference panels for both EAS and EUR populations from 1000 genomes Phase 3 VCFs. Missing rsIDs are replaced with "{chr}:{pos}:{ref}:{alt}" and the first entry for duplicate variants is retained.

- `merge_chromosome_level_LD_reference_panels.sh` merges the chromosome-level LD reference panels to use in downstream analyses.

