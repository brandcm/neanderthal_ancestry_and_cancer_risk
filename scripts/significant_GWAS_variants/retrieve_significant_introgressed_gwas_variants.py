"""
Retrieve the genome-wide significant GWAS SNVs that are Neanderthal introgressed.

Each GWAS's significant SNVs are matched against the merged introgressed variant
table by rsID, and then, for the variants left unmatched, by hg19 position. A
position match is kept only when the two alleles agree with the alleles the
introgression maps report, or when no map reports alleles for that position.

Inputs:  the GWAS metadata, the merged introgressed variant table and the directory
         of per-GWAS significant SNVs.
Outputs: a variant-level table of significant introgressed SNVs and a per-GWAS
         summary of how many significant SNVs each introgression map calls.
"""

import argparse
import sys

import pandas as pd
import yaml


MAPS = ["Browning", "Sankararaman", "Skov", "Vernot"]
ANNOTATION_COLUMNS = ["alleles", "neanderthal_allele", "maps", "n_maps", "notes"]
VARIANT_COLUMNS = [
    "cancer", "ancestry_group", "rsID", "chr", "pos", "effect_allele", "other_allele",
    "beta", "p_value", "neanderthal_allele", "introgression_maps", "n_maps", "match", "notes"
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--metadata", required=True, help="path to the GWAS metadata")
    parser.add_argument(
        "--introgressed_variants",
        required=True,
        help="path to the merged introgressed variant table",
    )
    parser.add_argument(
        "--significant_variants_directory",
        required=True,
        help="directory of per-GWAS significant SNVs",
    )
    parser.add_argument("--summary_output", required=True, help="path to the summary output")
    parser.add_argument(
        "--variants_output", required=True, help="path to the variant-level output"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.metadata) as f:
        gwas_metadata = yaml.safe_load(f)

    introgressed = load_introgressed_variants(args.introgressed_variants)

    variants = []
    summary = []

    # every GWAS is matched, including those without effect sizes: whether a variant
    # is significant and introgressed is a question of position, rsID and p-value.
    # Effect sizes are a fine-mapping prerequisite, and identify_fine_mapping_loci.py
    # applies that filter where it belongs; excluding those studies here would drop
    # them from the association-level counts as well.
    for gwas, metadata in gwas_metadata.items():
        significant = load_significant_variants(args.significant_variants_directory, gwas)
        matched = match_introgressed(significant, introgressed)

        matched.insert(0, "ancestry_group", metadata["ancestry_group"])
        matched.insert(0, "cancer", metadata["cancer"])
        variants.append(matched[VARIANT_COLUMNS])

        summary.append(
            {
                "cancer": metadata["cancer"],
                "ancestry_group": metadata["ancestry_group"],
                "N_significant_GWAS_variants": len(significant),
                **{
                    map_name: int(matched["introgression_maps"].str.contains(map_name).sum())
                    for map_name in MAPS
                },
                "union": len(matched),
            }
        )

        print(
            f"{gwas}: {len(matched)} of {len(significant)} significant SNVs are introgressed",
            file=sys.stderr,
        )

    pd.concat(variants).to_csv(args.variants_output, sep="\t", index=False)
    pd.DataFrame(summary).to_csv(args.summary_output, sep="\t", index=False)


def load_introgressed_variants(path):
    """Load the merged introgressed variant table, keeping one row per position and rsID."""
    introgressed = pd.read_csv(path, sep="\t", dtype=str, na_values=["NA", ""])
    introgressed = introgressed.rename(columns={"pos_hg19": "pos"})

    return introgressed


def load_significant_variants(directory, gwas):
    """Load one GWAS's significant SNVs."""
    return pd.read_csv(
        f"{directory}{gwas}_hg19_significant_SNVs.txt.gz",
        sep="\t",
        dtype={"chr": str, "pos": str, "rsID": str,
               "effect_allele": str, "other_allele": str},
    )


def allele_pair(first_allele, second_allele):
    """Return two alleles as an unordered pair, e.g. 'A/G'."""
    return "/".join(sorted([str(first_allele).upper(), str(second_allele).upper()]))


def match_introgressed(significant, introgressed):
    """Match significant SNVs to introgressed variants by rsID, then by position."""
    significant = significant.reset_index(drop=True)
    significant["_row"] = significant.index

    by_rsid = introgressed.dropna(subset=["rsID"]).drop_duplicates(subset=["rsID"])
    matched_by_rsid = significant.merge(
        by_rsid[["rsID"] + ANNOTATION_COLUMNS], on="rsID", how="inner"
    )
    matched_by_rsid["match"] = "rsID"

    # variants whose rsID is absent from the table may still be called at their position
    unmatched = significant[~significant["_row"].isin(matched_by_rsid["_row"])]
    by_position = introgressed.drop_duplicates(subset=["chr", "pos"])
    matched_by_position = unmatched.merge(
        by_position[["chr", "pos"] + ANNOTATION_COLUMNS], on=["chr", "pos"], how="inner"
    )

    if not matched_by_position.empty:
        reported_pairs = matched_by_position["alleles"].map(
            lambda alleles: "/".join(sorted(alleles.split("/"))) if pd.notna(alleles) else None
        )
        gwas_pairs = matched_by_position.apply(
            lambda variant: allele_pair(variant["effect_allele"], variant["other_allele"]),
            axis=1,
        )
        matched_by_position = matched_by_position[
            reported_pairs.isna() | (reported_pairs == gwas_pairs)
        ]

    matched_by_position = matched_by_position.assign(match="position")

    matched = pd.concat([matched_by_rsid, matched_by_position]).sort_values("_row")
    matched = matched.rename(columns={"maps": "introgression_maps"})
    matched["introgression_maps"] = matched["introgression_maps"].fillna("")

    return matched.drop(columns=["_row", "alleles"])


if __name__ == "__main__":
    main()
