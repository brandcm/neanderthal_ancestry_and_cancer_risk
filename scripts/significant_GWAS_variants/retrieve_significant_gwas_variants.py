"""
Retrieve genome-wide significant SNVs from one GWAS's harmonized summary statistics.

A variant is retained when it lies on an autosome, its effect and other alleles are
both a single base, and its p-value is below the significance threshold. Analyses are
restricted to the autosomes for consistency with LDSC, which does not model
chromosome X.

Input:  harmonized, tab-delimited summary statistics with a header (plain or gzipped).
Output: gzipped, tab-delimited, with the columns chr, pos, rsID, effect_allele,
        other_allele, beta and p_value.
"""

import argparse
import sys

import pandas as pd


COLUMNS = ["chr", "pos", "rsID", "effect_allele", "other_allele", "beta", "p_value"]
AUTOSOMES = [str(chromosome) for chromosome in range(1, 23)]
CHUNK_SIZE = 2_000_000


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="path to summary statistics")
    parser.add_argument("--chr_column", required=True, help="name of the chromosome column")
    parser.add_argument("--pos_column", required=True, help="name of the position column")
    parser.add_argument("--rsid_column", required=True, help="name of the variant ID column")
    parser.add_argument(
        "--effect_allele_column", required=True, help="name of the effect allele column"
    )
    parser.add_argument(
        "--other_allele_column", required=True, help="name of the other allele column"
    )
    parser.add_argument("--p_value_column", required=True, help="name of the p-value column")
    parser.add_argument(
        "--p_value_threshold",
        type=float,
        default=5e-8,
        help="significance threshold (default: 5e-8)",
    )
    parser.add_argument("--output", required=True, help="path to output file (.gz)")
    return parser.parse_args()


def main():
    args = parse_args()

    column_names = {
        args.chr_column: "chr",
        args.pos_column: "pos",
        args.rsid_column: "rsID",
        args.effect_allele_column: "effect_allele",
        args.other_allele_column: "other_allele",
        args.p_value_column: "p_value"
    }

    n_snvs = n_autosomal = n_significant = 0
    header_written = False

    # the summary statistics are read in chunks because the largest studies have tens
    # of millions of rows, while the significant variants are at most a few thousand
    chunks = pd.read_csv(
        args.input,
        sep="\t",
        usecols=list(column_names) + ["beta"],
        dtype={args.chr_column: str, args.rsid_column: str,
               args.effect_allele_column: str, args.other_allele_column: str},
        chunksize=CHUNK_SIZE,
    )

    for chunk in chunks:
        chunk = chunk.rename(columns=column_names)

        # chromosomes are reported without a "chr" prefix, as in the introgressed
        # variant table and the LD reference panels
        chunk["chr"] = chunk["chr"].str.replace(r"^chr", "", regex=True)

        snvs = chunk[
            (chunk["effect_allele"].str.len() == 1) & (chunk["other_allele"].str.len() == 1)
        ]
        autosomal_snvs = snvs[snvs["chr"].isin(AUTOSOMES)]
        significant_snvs = autosomal_snvs[autosomal_snvs["p_value"] < args.p_value_threshold]

        n_snvs += len(snvs)
        n_autosomal += len(autosomal_snvs)
        n_significant += len(significant_snvs)

        significant_snvs[COLUMNS].to_csv(
            args.output,
            sep="\t",
            index=False,
            header=not header_written,
            mode="w" if not header_written else "a",
            compression="gzip",
        )
        header_written = True

    if not header_written:
        pd.DataFrame(columns=COLUMNS).to_csv(
            args.output, sep="\t", index=False, compression="gzip"
        )

    print(
        f"{n_autosomal} of {n_snvs} SNVs are autosomal; "
        f"{n_significant} are significant at p < {args.p_value_threshold}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
