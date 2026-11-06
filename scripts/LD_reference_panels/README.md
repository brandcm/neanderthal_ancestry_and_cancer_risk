# LD reference panels

A script to generate linkage disequilibrium reference panels for East Asian (EAS) and European (EUR) ancestry groups from 1000 Genomes Phase 3 genotypes. The panels are on hg19 and are used for fine-mapping and other downstream analyses.

## Source data

1000 Genomes Phase 3 VCFs in the Capra Lab data directory, `/wynton/group/capra/data/1KG_phase_3/2025-11-05/`, and the sample IDs for each ancestry group in `data/metadata/1000G_EAS_sample_IDs.txt` (n = 504) and `data/metadata/1000G_EUR_sample_IDs.txt` (n = 502).

## Pipeline

`generate_LD_reference_panels.sh` runs on the cluster as an array job, one task per chromosome. Submit it from this directory:

```bash
mkdir -p logs
qsub generate_LD_reference_panels.sh
```

For each chromosome, the script:

- filters the VCF to biallelic SNVs with bcftools;
- writes one PLINK binary fileset per ancestry group, keeping only that group's samples, naming variants without an rsID `{chr}:{pos}:{ref}:{alt}` and retaining the first of any duplicated variants.

## Outputs

In `data/LD_reference_panels/`:

| File | Contents |
|---|---|
| `EAS/1000G_EAS_chr{chr}.bed` / `.bim` / `.fam` | EAS panel for one chromosome |
| `EUR/1000G_EUR_chr{chr}.bed` / `.bim` / `.fam` | EUR panel for one chromosome |
| `{EAS,EUR}/1000G_{EAS,EUR}_chr{chr}.log` | PLINK's record of the variants and samples retained |

Per-job cluster logs are written to `logs/` in this scripts directory.
