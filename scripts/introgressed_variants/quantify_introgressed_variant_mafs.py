"""
Minor allele frequencies (gnomAD v4.1.0) for introgressed variants on one chromosome.

Variants are the union of rsIDs across the introgression maps and are matched to
gnomAD by rsID (gnomAD v4.1 is GRCh38; the maps are hg19, so positions are not used).

Output (tab-delimited, with header): rsID  MAF  MAF_afr  MAF_eas  MAF_fin  MAF_nfe  MAF_sas
MAF is the frequency of the rarest allele at the site ("." if gnomAD reports no AF).
"""

import argparse
import gzip
import sys

AF_KEYS = ["AF", "AF_afr", "AF_eas", "AF_fin", "AF_nfe", "AF_sas"]
HEADER = ["rsID"] + [key.replace("AF", "MAF", 1) for key in AF_KEYS]


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--vcf", required=True, help="path to gnomAD v4.1.0 VCF")
    parser.add_argument("--chr", required=True, help="target chromosome")
    parser.add_argument(
        "--variants",
        required=True,
        nargs="+",
        help="introgressed variant files (chr pos rsID ...), one per map",
    )
    parser.add_argument("--output", required=True, help="path to output file")
    return parser.parse_args()


def main():
    args = parse_args()
    introgressed_rsids = load_introgressed_rsids(args.variants, args.chr)
    print(f"chr{args.chr}: {len(introgressed_rsids)} introgressed rsIDs", file=sys.stderr)
    n_written = write_mafs(args.vcf, introgressed_rsids, args.output)
    print(f"chr{args.chr}: wrote MAFs for {n_written} variants", file=sys.stderr)


def load_introgressed_rsids(paths, target_chr):
    """Union of introgressed rsIDs on one chromosome across all variant files."""
    target_chr = target_chr.removeprefix("chr")
    rsids = set()

    for path in paths:
        with open(path) as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue

                chr_, rsid = fields[0].removeprefix("chr"), fields[2]
                if chr_ == target_chr and rsid not in {".", ""}:
                    rsids.add(rsid)

    return rsids


def open_text(path):
    """Open a plain or gzipped text file for reading."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_info(info_field):
    """VCF INFO column as a dictionary (flags without "=" are skipped)."""
    info = {}
    for item in info_field.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
    return info


def minor_allele_frequency(af_value):
    """Frequency of the rarest allele given a comma-separated list of ALT AFs, or "."."""
    try:
        alt_afs = [float(x) for x in af_value.split(",") if x != "."]
    except ValueError:
        return "."
    if not alt_afs:
        return "."

    ref_af = 1.0 - sum(alt_afs)
    return str(round(min([ref_af] + alt_afs), 6))


def write_mafs(vcf_path, introgressed_rsids, output_path):
    """Write MAFs for VCF records whose ID matches an introgressed rsID; returns the count."""
    n_written = 0

    with open_text(vcf_path) as vcf, open(output_path, "w") as out:
        out.write("\t".join(HEADER) + "\n")

        for line in vcf:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if fields[2] == ".":
                continue

            # a record can list several rsIDs (rs1;rs2); report the first introgressed one
            matching = [rsid for rsid in fields[2].split(";") if rsid in introgressed_rsids]
            if not matching:
                continue

            info = parse_info(fields[7])
            mafs = [minor_allele_frequency(info.get(key, ".")) for key in AF_KEYS]
            out.write("\t".join([matching[0]] + mafs) + "\n")
            n_written += 1

    return n_written


if __name__ == "__main__":
    main()