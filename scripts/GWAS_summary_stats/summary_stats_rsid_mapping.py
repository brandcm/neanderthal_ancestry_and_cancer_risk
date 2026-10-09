"""
Add an rsID column to GWAS summary statistics by matching chromosome and position.

Input:  tab-delimited summary statistics with a header (plain or gzipped). Positions
        must be GRCh37/hg19, the build of the gnomAD v2.1.1 rsID dictionaries.
Output: the same rows with an rsID column appended (gzipped); "." where no rsID is found.

Dictionaries are loaded one chromosome at a time, so the input should be sorted (or at
least grouped) by chromosome. Ungrouped input still works but reloads dictionaries
repeatedly; the number of reloads is reported.
"""

import argparse
import gzip
import os
import pickle
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="path to GWAS summary statistics")
    parser.add_argument(
        "--rsid_mapping_directory",
        required=True,
        help="directory with gnomAD_v2.1.1_PASS_SNVs_chr{chr}_pos_rsIDs.pkl.gz files",
    )
    parser.add_argument("--chr_column_index", type=int, required=True, help="0-based chromosome column")
    parser.add_argument("--pos_column_index", type=int, required=True, help="0-based position column")
    parser.add_argument("--output", required=True, help="path to output file (.gz)")
    return parser.parse_args()


def main():
    args = parse_args()

    n_rows = n_mapped = n_multiple = n_unparsed = 0
    loaded_chromosomes = []          # in load order, to detect ungrouped input
    current_chr, rsid_dict = None, {}

    with open_text(args.input) as infile, gzip.open(args.output, "wt") as outfile:
        header = infile.readline().rstrip("\n")
        outfile.write(f"{header}\trsID\n")

        for line in infile:
            line = line.rstrip("\n")
            fields = line.split("\t")
            n_rows += 1

            try:
                chr_ = normalize_chromosome(fields[args.chr_column_index])
                pos = int(fields[args.pos_column_index])
            except (IndexError, ValueError):
                n_unparsed += 1
                outfile.write(f"{line}\t.\n")
                continue

            if chr_ != current_chr:
                rsid_dict = load_rsid_dict(args.rsid_mapping_directory, chr_)
                loaded_chromosomes.append(chr_)
                current_chr = chr_

            rsid = rsid_dict.get(pos, ".")
            if "," in rsid:
                n_multiple += 1
                rsid = rsid.split(",")[0]          # several rsIDs at one position: keep the first
            if rsid != ".":
                n_mapped += 1

            outfile.write(f"{line}\t{rsid}\n")

    n_reloads = len(loaded_chromosomes) - len(set(loaded_chromosomes))
    print(f"{n_mapped} of {n_rows} rows mapped to an rsID "
          f"({n_multiple} with multiple rsIDs, first kept; {n_unparsed} with unreadable chr/pos)",
          file=sys.stderr)
    if n_reloads:
        print(f"warning: input is not grouped by chromosome ({n_reloads} dictionary reloads); "
              "sorting it first is much faster", file=sys.stderr)


def normalize_chromosome(value):
    """'chr1' -> '1', '23' -> 'X'."""
    chr_ = str(value).removeprefix("chr")
    return "X" if chr_ == "23" else chr_


def load_rsid_dict(directory, chr_):
    """Position -> rsID dictionary for one chromosome; empty (with a warning) if the file is missing."""
    path = os.path.join(directory, f"gnomAD_v2.1.1_PASS_SNVs_chr{chr_}_pos_rsIDs.pkl.gz")
    if not os.path.exists(path):
        print(f"warning: no rsID dictionary for chr{chr_} ({path}); rows get '.'", file=sys.stderr)
        return {}

    print(f"loading rsID dictionary for chr{chr_}", file=sys.stderr)
    with gzip.open(path, "rb") as f:
        return pickle.load(f)


def open_text(path):
    """Open a plain or gzipped text file for reading."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


if __name__ == "__main__":
    main()
