"""
Add rsIDs to introgressed variants for one chromosome, matching on hg19 position.

Input (tab-delimited, no header):   chr  pos  [other columns...]
Output (tab-delimited, no header):  chr  pos  rsID  [other columns...]

rsIDs come from a gnomAD v2.1.1 position -> rsID dictionary. Variants with no rsID
get ".". When a position has several rsIDs, the first is kept and the rest are reported.
"""

import argparse
import gzip
import pickle
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="path to input variants")
    parser.add_argument("--dict", required=True, help="path to the position -> rsID dictionary (.pkl.gz)")
    parser.add_argument("--chr", required=True, help="chromosome to map (without 'chr')")
    parser.add_argument("--output", required=True, help="path to output file")
    return parser.parse_args()


def main():
    args = parse_args()
    chromosome = args.chr.removeprefix("chr")

    with gzip.open(args.dict, "rb") as f:
        rsid_dict = pickle.load(f)

    n_variants = n_mapped = n_multiple = 0
    with open(args.input) as infile, open(args.output, "w") as outfile:
        for line in infile:
            fields = line.rstrip("\n").split("\t")
            if not fields[0]:
                continue

            # strip "chr" before comparing so "chr1" and "1" both match
            chr_ = fields[0].removeprefix("chr")
            if chr_ != chromosome:
                continue

            pos = int(fields[1])
            rsid = rsid_dict.get(pos, ".")
            n_variants += 1

            if "," in rsid:
                n_multiple += 1
                print(f"multiple rsIDs at {chr_}:{pos} ({rsid}); keeping the first", file=sys.stderr)
                rsid = rsid.split(",")[0]
            if rsid != ".":
                n_mapped += 1

            outfile.write("\t".join([chr_, str(pos), rsid] + fields[2:]) + "\n")

    print(f"chr{chromosome}: {n_mapped} of {n_variants} variants mapped to an rsID "
          f"({n_multiple} positions with multiple rsIDs)", file=sys.stderr)


if __name__ == "__main__":
    main()