"""
Harmonize effect sizes in GWAS summary statistics to a log-odds beta and its standard error.

- beta: copied from the beta column, or log(OR) when the beta column is an odds ratio
  ("OR" / "odds_ratio", any case). Non-positive odds ratios become missing.
- standard_error: copied from the standard-error column; where it is missing (or the
  column is absent), derived from beta and the two-sided p-value: SE = |beta| / z,
  z = Phi^-1(1 - p/2).
- Rows without a p-value are dropped, and columns that are empty in every row are removed.

Input:  tab-delimited summary statistics with a header (plain or gzipped).
Output: gzipped, tab-delimited; all original columns plus "beta" and "standard_error";
        missing values written as ".".
"""

import argparse
import gzip
import sys

import numpy as np
import pandas as pd
from scipy import stats

NA_VALUES = ["", " ", ".", "NA", "NaN"]
ODDS_RATIO_NAMES = {"or", "odds_ratio"}
MIN_P = 1e-300        # p-values below this (or parsed as 0) are clipped so z stays finite


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="path to GWAS summary statistics")
    parser.add_argument("--beta_column", required=True, help="beta (or odds ratio) column")
    parser.add_argument("--se_column", required=True, help="standard error column")
    parser.add_argument("--p_value_column", required=True, help="p-value column")
    parser.add_argument("--chunksize", type=int, default=100_000, help="rows processed at a time")
    parser.add_argument("--output", required=True, help="path to output file (.gz)")
    return parser.parse_args()


def main():
    args = parse_args()
    beta_column, se_column, p_column = (c.strip() for c in
                                        (args.beta_column, args.se_column, args.p_value_column))
    is_odds_ratio = beta_column.lower() in ODDS_RATIO_NAMES

    # pass 1: find columns that are empty in every row
    empty_columns = find_empty_columns(args.input, args.chunksize)
    drop_columns = [c for c in empty_columns if c not in {beta_column, se_column, p_column}]
    if drop_columns:
        print(f"dropping empty columns: {', '.join(drop_columns)}", file=sys.stderr)

    # pass 2: harmonize and write
    n_in = n_out = n_se_derived = 0
    with gzip.open(args.output, "wt") as out:
        for i, chunk in enumerate(read_chunks(args.input, args.chunksize)):
            n_in += len(chunk)
            chunk = chunk.drop(columns=drop_columns)

            # beta: odds ratios are kept and log(OR) is added as "beta"; otherwise renamed in place
            if is_odds_ratio:
                chunk["beta"] = harmonized_beta(chunk[beta_column], is_odds_ratio=True)
            else:
                chunk[beta_column] = harmonized_beta(chunk[beta_column], is_odds_ratio=False)
                chunk = chunk.rename(columns={beta_column: "beta"})

            # standard error: fill missing values (or a missing column) from beta and p
            p_value = pd.to_numeric(chunk[p_column], errors="coerce")
            if se_column not in chunk.columns:
                chunk[se_column] = np.nan
            se = pd.to_numeric(chunk[se_column], errors="coerce")
            missing_se = se.isna()
            se[missing_se] = se_from_p(chunk.loc[missing_se, "beta"], p_value[missing_se])
            n_se_derived += int((missing_se & se.notna()).sum())
            chunk[se_column] = se
            chunk = chunk.rename(columns={se_column: "standard_error"})

            chunk = chunk[p_value.notna()]
            n_out += len(chunk)
            chunk.to_csv(out, sep="\t", index=False, header=(i == 0), na_rep=".")

    print(f"{n_out} of {n_in} rows written ({n_in - n_out} without a p-value dropped; "
          f"SE derived from beta and p for {n_se_derived}); beta from "
          f"{'log(odds ratio)' if is_odds_ratio else 'beta column'}", file=sys.stderr)


def read_chunks(path, chunksize):
    """Read tab-delimited summary statistics in chunks, with column names stripped of whitespace."""
    reader = pd.read_csv(path, sep="\t", chunksize=chunksize, na_values=NA_VALUES,
                         keep_default_na=True, low_memory=False)
    for chunk in reader:
        chunk.columns = chunk.columns.str.strip()
        yield chunk


def find_empty_columns(path, chunksize):
    """Columns that are missing in every row of the file."""
    all_empty = None
    for chunk in read_chunks(path, chunksize):
        empty = chunk.isna().all()
        all_empty = empty if all_empty is None else (all_empty & empty)
    return [] if all_empty is None else all_empty[all_empty].index.tolist()


def harmonized_beta(values, is_odds_ratio):
    """Numeric beta; log of the odds ratio when is_odds_ratio (non-positive ORs -> NaN)."""
    values = pd.to_numeric(values, errors="coerce")
    if not is_odds_ratio:
        return values
    return np.log(values.where(values > 0))


def se_from_p(beta, p_value):
    """SE = |beta| / z with z = Phi^-1(1 - p/2); NaN where beta or p is missing or beta is 0."""
    p_value = p_value.clip(lower=MIN_P, upper=1.0)
    z = stats.norm.isf(p_value / 2)        # isf is accurate for very small p
    se = np.abs(beta) / z
    return se.where((beta != 0) & (z > 0))


if __name__ == "__main__":
    main()
