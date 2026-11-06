"""
Collect the per-GWAS regression diagnostics and heritability into one table.

Each ldsc.py --h2 run under the baselineLD model, without any introgression
annotation, writes a log containing the total heritability, lambda GC, the mean
chi-square and the regression intercept. Those four describe the GWAS rather than
any annotation, so they belong in one row per study rather than repeated across the
partitioned results, which re-estimate all four for every variant set.

The intercept is the reason to report the other three: lambda GC and the mean
chi-square say how inflated the statistics are, and the intercept says how much of
that inflation is confounding rather than polygenicity. An intercept near 1 with a
high mean chi-square is the signature of a well-behaved, highly polygenic trait; an
intercept well above 1 is the alternative explanation a reader will reach for.

Heritability is on the observed scale unless the GWAS metadata records both
prevalences, in which case LDSC also reports the liability scale and both are
carried through. Observed-scale estimates depend on each study's case fraction and
are not comparable between cancers.

Inputs:  the ldsc --h2 logs, the LDSC formatting metadata and the GWAS metadata.
Output:  tab-delimited, one row per GWAS.
"""

import argparse
import re

import pandas as pd
import yaml


COLUMNS = [
    "cancer", "ancestry_group", "study", "n_cases", "n_controls", "n_effective",
    "n_snps", "lambda_GC", "mean_chi2", "intercept", "intercept_SE",
    "h2_observed", "h2_observed_SE", "sample_prevalence", "population_prevalence",
    "h2_liability", "h2_liability_SE",
]

# "Total Observed scale h2: 0.0799 (0.0095)" and the liability-scale equivalent
ESTIMATE_AND_SE = re.compile(r":\s*(-?[\d.eE+-]+)\s*\(\s*(-?[\d.eE+-]+)\s*\)")
# "Lambda GC: 1.1428"
ESTIMATE_ONLY = re.compile(r":\s*(-?[\d.eE+-]+)")
# "After merging with regression SNP LD, 1022315 SNPs remain."
N_SNPS = re.compile(r"After merging with regression SNP LD,\s*([\d,]+)\s*SNPs remain")


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--logs", required=True, nargs="+", help="the ldsc --h2 logs, one per GWAS"
    )
    parser.add_argument(
        "--ldsc_metadata", required=True, help="path to the LDSC formatting metadata"
    )
    parser.add_argument(
        "--cancer_metadata", required=True, help="path to the GWAS metadata"
    )
    parser.add_argument("--output", required=True, help="path to output file")
    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.ldsc_metadata) as f:
        ldsc_metadata = yaml.safe_load(f)
    with open(args.cancer_metadata) as f:
        cancer_metadata = yaml.safe_load(f)

    rows = []

    for path in args.logs:
        gwas = path.split("/")[-1].removesuffix("_total_h2.log")
        with open(path) as f:
            log = f.read()

        values = parse_log(log)
        ldsc_entry = ldsc_metadata.get(gwas, {})
        # the two metadata files are keyed the same way, but a GWAS may be in one
        # and not the other, so neither is assumed to carry a field
        cancer_entry = cancer_metadata.get(gwas, {})

        rows.append(
            {
                "cancer": cancer_entry.get("cancer", ldsc_entry.get("cancer")),
                "ancestry_group": ldsc_entry.get(
                    "ancestry_group", cancer_entry.get("ancestry_group")
                ),
                "study": gwas,
                "n_cases": cancer_entry.get("n_cases"),
                "n_controls": cancer_entry.get("n_controls"),
                "n_effective": cancer_entry.get("n_effective"),
                "sample_prevalence": ldsc_entry.get("sample_prevalence"),
                "population_prevalence": ldsc_entry.get("population_prevalence"),
                **values,
            }
        )

    table = pd.DataFrame(rows)
    for column in COLUMNS:
        if column not in table.columns:
            table[column] = pd.NA
    table = table[COLUMNS].sort_values(["cancer", "ancestry_group"])
    table.to_csv(args.output, sep="\t", index=False)

    missing = table[table["h2_observed"].isna()]
    if not missing.empty:
        print(
            f"No heritability parsed for {len(missing)} GWAS: "
            f"{', '.join(missing['study'])}. Check those logs for an LDSC error.",
        )
    print(f"{len(table)} GWAS written to {args.output}")


def parse_log(log):
    """Pull the diagnostics out of one ldsc --h2 log.

    A log from a failed or truncated run is missing the lines rather than
    malformed, so every field is optional and absent ones stay missing instead of
    being silently reported as zero.
    """
    values = {}

    for line in log.splitlines():
        if line.startswith("Total Observed scale h2"):
            match = ESTIMATE_AND_SE.search(line)
            if match:
                values["h2_observed"] = float(match.group(1))
                values["h2_observed_SE"] = float(match.group(2))
        elif line.startswith("Total Liability scale h2"):
            match = ESTIMATE_AND_SE.search(line)
            if match:
                values["h2_liability"] = float(match.group(1))
                values["h2_liability_SE"] = float(match.group(2))
        elif line.startswith("Intercept"):
            match = ESTIMATE_AND_SE.search(line)
            if match:
                values["intercept"] = float(match.group(1))
                values["intercept_SE"] = float(match.group(2))
        elif line.startswith("Lambda GC"):
            match = ESTIMATE_ONLY.search(line)
            if match:
                values["lambda_GC"] = float(match.group(1))
        elif line.startswith("Mean Chi^2"):
            match = ESTIMATE_ONLY.search(line)
            if match:
                values["mean_chi2"] = float(match.group(1))

    match = N_SNPS.search(log)
    if match:
        values["n_snps"] = int(match.group(1).replace(",", ""))

    return values


if __name__ == "__main__":
    main()
