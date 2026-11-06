"""
Identify the loci to fine-map from the significant introgressed GWAS variants.

Within each GWAS, the most strongly associated introgressed variant becomes the first
locus, every introgressed variant within the minimum distance of it is set aside, and
the next strongest remaining variant becomes the next locus. Each locus is therefore
centered on the strongest signal near it. Lead variants are at least the minimum
distance apart, so windows overlap when that distance is smaller than the window;
setting the minimum distance to the window size gives windows that never overlap.
Variants are not filtered on whether the introgression maps agree about the
Neanderthal allele.

Input:  the variant-level output of retrieve_significant_introgressed_gwas_variants.py
        and the GWAS metadata.
Output: tab-delimited, one row per locus. The columns config_key, ancestry_group,
        locus_id, n_samples and lead_rsid are named after the arguments of
        fine_mapping(), so the file can be read and passed straight to
        run_all_loci(); the remaining columns are context.
"""

import argparse
import sys

import pandas as pd
import yaml


COLUMNS = [
    "config_key", "cancer", "ancestry_group", "locus_id", "lead_rsid", "n_samples",
    "chr", "pos", "effect_allele", "other_allele", "beta", "p_value",
    "neanderthal_allele", "introgression_maps", "n_maps", "window_size",
    "window_start", "window_end", "n_introgressed_variants"
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", required=True, help="path to the significant introgressed variants"
    )
    parser.add_argument("--metadata", required=True, help="path to the GWAS metadata")
    parser.add_argument(
        "--window_size",
        type=int,
        default=1000,
        help="width of the window around each lead variant, in kb (default: 1000)",
    )
    parser.add_argument(
        "--minimum_distance",
        type=int,
        help="minimum distance between lead variants, in kb (default: half the window)",
    )
    parser.add_argument("--output", required=True, help="path to output file")
    return parser.parse_args()


def main():
    args = parse_args()
    half_window = args.window_size * 1000 // 2
    minimum_distance = (
        args.minimum_distance * 1000 if args.minimum_distance is not None else half_window
    )

    variants = pd.read_csv(
        args.input,
        sep="\t",
        dtype={"chr": str, "rsID": str, "effect_allele": str, "other_allele": str},
    )
    variants["pos"] = variants["pos"].astype(int)

    with open(args.metadata) as f:
        gwas_metadata = yaml.safe_load(f)

    # the variants are labeled by cancer and ancestry group, which identify one GWAS each
    studies = {
        (metadata["cancer"], metadata["ancestry_group"]): (gwas, metadata["n_effective"])
        for gwas, metadata in gwas_metadata.items()
    }
    # SuSiE needs effect sizes to turn the summary statistics into z scores, so a study
    # that reports none cannot be fine-mapped. Its variants are still in the input,
    # which counts significant introgressed associations; the filter belongs here
    # rather than upstream.
    no_effect_sizes = {
        (metadata["cancer"], metadata["ancestry_group"])
        for metadata in gwas_metadata.values()
        if not metadata["effect_sizes"]
    }

    loci = []

    for (cancer, ancestry_group), gwas_variants in variants.groupby(
        ["cancer", "ancestry_group"], sort=True
    ):
        if (cancer, ancestry_group) not in studies:
            raise ValueError(
                f"No GWAS in the metadata has cancer {cancer!r} and ancestry group "
                f"{ancestry_group!r}."
            )

        if (cancer, ancestry_group) in no_effect_sizes:
            print(
                f"Skipping {cancer} {ancestry_group}: the study reports no effect sizes, "
                f"so its {len(gwas_variants)} significant introgressed variants cannot "
                "be fine-mapped.",
                file=sys.stderr,
            )
            continue

        config_key, n_samples = studies[(cancer, ancestry_group)]

        for locus in identify_loci(gwas_variants, half_window, minimum_distance):
            locus["config_key"] = config_key
            locus["n_samples"] = n_samples
            locus["window_size"] = args.window_size
            loci.append(locus)

    loci = pd.DataFrame(loci)[COLUMNS]
    loci.to_csv(args.output, sep="\t", index=False)

    print(
        f"{len(loci)} loci across {loci['config_key'].nunique()} GWAS "
        f"from {len(variants)} significant introgressed variants",
        file=sys.stderr,
    )


def identify_loci(variants, half_window, minimum_distance):
    """Take the strongest signal, then the strongest beyond the minimum distance, and so on."""
    # ties are broken by position so that the same variants always give the same loci
    variants = variants.sort_values(["p_value", "chr", "pos"])
    loci = []

    for variant in variants.itertuples(index=False):
        if any(
            locus["chr"] == variant.chr
            and abs(locus["pos"] - variant.pos) <= minimum_distance
            for locus in loci
        ):
            continue

        window_start = max(variant.pos - half_window, 1)
        window_end = variant.pos + half_window
        in_window = variants[
            (variants["chr"] == variant.chr)
            & (variants["pos"] >= window_start)
            & (variants["pos"] <= window_end)
        ]

        loci.append(
            {
                "cancer": variant.cancer,
                "ancestry_group": variant.ancestry_group,
                "locus_id": f"{variant.chr}_{variant.pos}",
                "chr": variant.chr,
                "pos": variant.pos,
                "lead_rsid": variant.rsID,
                "effect_allele": variant.effect_allele,
                "other_allele": variant.other_allele,
                "beta": variant.beta,
                "p_value": variant.p_value,
                "neanderthal_allele": variant.neanderthal_allele,
                "introgression_maps": variant.introgression_maps,
                "n_maps": variant.n_maps,
                "window_start": window_start,
                "window_end": window_end,
                "n_introgressed_variants": len(in_window),
            }
        )

    return sorted(loci, key=lambda locus: (locus["chr"].zfill(2), locus["pos"]))


if __name__ == "__main__":
    main()
