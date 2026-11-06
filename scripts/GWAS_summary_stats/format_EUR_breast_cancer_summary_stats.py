"""
Format the EUR breast cancer summary statistics (Zhang, GCST001): convert the
space-delimited download to tab-delimited and add an rsID column taken from the
variant ID (e.g. "rs123:100:A:G" -> "rs123").

Input:  space-delimited summary statistics with a header (plain or gzipped).
Output: gzipped, tab-delimited; the original columns plus rsID ("." when the
        variant ID does not start with an rsID).
"""

import argparse
import gzip
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="path to summary statistics")
    parser.add_argument(
        "--variant_id_column_index",
        type=int,
        default=7,
        help="0-based index of the variant ID column (default: 7)",
    )
    parser.add_argument("--output", required=True, help="path to output file (.gz)")
    return parser.parse_args()


def main():
    args = parse_args()
    n_rows = n_with_rsid = 0

    with open_text(args.input) as infile, gzip.open(args.output, "wt") as outfile:
        header = to_fields(infile.readline())
        outfile.write("\t".join(header + ["rsID"]) + "\n")

        for line in infile:
            fields = to_fields(line)
            rsid = rsid_from_variant_id(fields[args.variant_id_column_index])
            n_rows += 1
            n_with_rsid += rsid != "."
            outfile.write("\t".join(fields + [rsid]) + "\n")

    print(f"{n_with_rsid} of {n_rows} rows have an rsID", file=sys.stderr)


def to_fields(line):
    """Split a space-delimited line into fields (each space is one separator, as with `tr ' ' '\\t'`)."""
    return line.rstrip("\n").replace(" ", "\t").split("\t")


def rsid_from_variant_id(variant_id):
    """'rs123:100:A:G' -> 'rs123'; '.' if the ID does not start with an rsID."""
    rsid = variant_id.split(":")[0]
    return rsid if rsid.startswith("rs") else "."


def open_text(path):
    """Open a plain or gzipped text file for reading."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


if __name__ == "__main__":
    main()
