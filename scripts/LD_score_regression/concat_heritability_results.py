"""
Collect the introgressed-variant partition from every LDSC partitioned heritability run.

LDSC writes one .results file per GWAS and introgressed variant call set, holding a row
per annotation in the model: the ~97 baselineLD annotations followed by the introgressed
variant annotation. Only the last of these is of interest, and it is read together with
the model-level statistics LDSC reports in the .log rather than the .results file.

The h2, intercept and mean chi-square are needed to judge the partition: an enrichment
estimated on a trait with negligible or badly estimated heritability is not interpretable,
and that is visible only in the log.

Input:  .results files named <gwas>_<variant_set>.results, each beside its .log.
Output: tab-delimited, one row per input, with gwas and variant_set as separate columns.
"""

import argparse
import os
import re

import pandas as pd


# LDSC names the annotation added after the baseline model by its position in --ref-ld-chr
INTROGRESSED_CATEGORY = "L2_1"

LOG_STATISTICS = ["n_snps", "h2", "h2_se", "lambda_gc", "mean_chi2", "intercept", "intercept_se"]

# LDSC prints estimates in both decimal and scientific notation, and nan for a failed fit
NUMBER = r"[-+\d.eEnan]+"


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--results", required=True, nargs="+", help="paths to LDSC .results files"
    )
    parser.add_argument(
        "--variant_sets",
        required=True,
        nargs="+",
        help="introgressed variant call set names, used to split <gwas>_<variant_set>",
    )
    parser.add_argument(
        "--category",
        default=INTROGRESSED_CATEGORY,
        help=f"annotation to retain from each .results file (default: {INTROGRESSED_CATEGORY})",
    )
    parser.add_argument("--output", required=True, help="path to output file")
    return parser.parse_args()


def main():
    args = parse_args()

    rows = []
    for results_path in args.results:
        gwas, variant_set = split_file_stem(results_path, args.variant_sets)

        results = pd.read_csv(results_path, sep="\t")
        partition = results[results["Category"] == args.category]
        if partition.empty:
            raise ValueError(
                f"No category '{args.category}' in {results_path}. "
                f"Categories present: {', '.join(results['Category'].tail(5))}"
            )
        if len(partition) > 1:
            raise ValueError(f"Category '{args.category}' appears more than once in {results_path}")

        row = partition.drop(columns=["Category"]).iloc[0].to_dict()
        row.update(retrieve_ldsc_statistics(results_path.replace(".results", ".log")))
        rows.append({"gwas": gwas, "variant_set": variant_set, **row})

    summary = pd.DataFrame(rows).sort_values(["variant_set", "gwas"])
    summary.to_csv(args.output, sep="\t", index=False)

    print(f"Wrote {len(summary)} partitions to {args.output}")


def split_file_stem(results_path, variant_sets):
    """Split <gwas>_<variant_set>.results. Both names contain underscores, so the call set
    is matched as a suffix rather than by splitting on the separator."""
    stem = os.path.basename(results_path).replace(".results", "")
    for variant_set in sorted(variant_sets, key=len, reverse=True):
        if stem.endswith(f"_{variant_set}"):
            return stem[: -len(variant_set) - 1], variant_set
    raise ValueError(
        f"{stem} does not end in any of the given call sets: {', '.join(variant_sets)}"
    )


def retrieve_ldsc_statistics(log_file):
    """Model-level statistics from one LDSC log. Absent values are left empty rather than
    defaulted, so that a partition estimated on a failed fit is visible in the summary."""
    statistics = {name: None for name in LOG_STATISTICS}

    if not os.path.exists(log_file):
        print(f"Warning: no log beside {os.path.basename(log_file)}; its statistics are empty")
        return statistics

    with open(log_file) as f:
        content = f.read()

    # the number of variants is reported at several filtering steps; the last is the number used
    remaining = re.findall(r"(\d+) SNPs remain", content)
    if remaining:
        statistics["n_snps"] = int(remaining[-1])

    statistics["h2"], statistics["h2_se"] = match_estimate(
        rf"Total Observed scale h2: ({NUMBER}) \(({NUMBER})\)", content
    )
    statistics["intercept"], statistics["intercept_se"] = match_estimate(
        rf"Intercept: ({NUMBER}) \(({NUMBER})\)", content
    )
    statistics["lambda_gc"] = match_estimate(rf"Lambda GC: ({NUMBER})", content)[0]
    statistics["mean_chi2"] = match_estimate(rf"Mean Chi\^2: ({NUMBER})", content)[0]

    return statistics


def match_estimate(pattern, content):
    """Groups of the first match as numbers, or as many Nones as the pattern has groups."""
    match = re.search(pattern, content)
    if match is None:
        return (None,) * re.compile(pattern).groups
    return tuple(to_number(group) for group in match.groups())


def to_number(text):
    """LDSC writes nan for an estimate it could not produce, which float() accepts silently."""
    try:
        value = float(text)
    except ValueError:
        return None
    return None if pd.isna(value) else value


if __name__ == "__main__":
    main()
