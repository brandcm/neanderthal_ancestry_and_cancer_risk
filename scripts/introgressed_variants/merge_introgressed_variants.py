#!/usr/bin/env python3
"""
Merge the per-map introgressed variant files (hg19) into one table with one row
per position: rsID, alleles, Neanderthal allele, and which maps call it.

Input formats (whitespace-delimited, no header):
  Browning:     chr pos rsID ref alt archaic_index(0=ref,1=alt) status   (already filtered to "match")
  Sankararaman: chr pos rsID                                            (no alleles)
  Skov:         chr pos rsID ancestral derived                          (derived = Neanderthal allele)
  Vernot:       chr pos rsID ref alt neanderthal_allele                 ("G/A" = Neanderthal heterozygous)

Output columns (tab-delimited):
  chr pos_hg19 rsID ref alt alleles neanderthal_allele maps n_maps notes
  - ref/alt: hg19 ref/alt from Browning or Vernot (NA otherwise; Skov's ancestral is not necessarily ref)
  - alleles: unordered allele pair from any map that reports two alleles, e.g. "A/G"
  - neanderthal_allele: NA if no map reports one or if maps disagree
"""

import argparse
import sys

import pandas as pd

MAP_NAMES = ["Browning", "Sankararaman", "Skov", "Vernot"]
STD_COLS = ["chr", "pos", "rsID", "ref", "alt", "allele_1", "allele_2", "neanderthal_allele", "map"]


def read_map(path, names, map_name):
    df = pd.read_csv(path, sep=r"\s+", header=None, names=names, dtype=str, engine="c")
    df["map"] = map_name
    print(f"{map_name}: {len(df):,} rows from {path}", file=sys.stderr)
    return df


def read_browning(path):
    df = read_map(path, ["chr", "pos", "rsID", "ref", "alt", "archaic_index", "status"], "Browning")
    print(df.groupby(["archaic_index", "status"]).size().to_string(), file=sys.stderr)
    df["neanderthal_allele"] = df["alt"].where(df["archaic_index"] == "1",
                                               df["ref"].where(df["archaic_index"] == "0"))
    df["allele_1"], df["allele_2"] = df["ref"], df["alt"]
    return df[STD_COLS]


def read_sankararaman(path):
    df = read_map(path, ["chr", "pos", "rsID"], "Sankararaman")
    for col in ["ref", "alt", "allele_1", "allele_2", "neanderthal_allele"]:
        df[col] = pd.NA
    return df[STD_COLS]


def read_skov(path):
    df = read_map(path, ["chr", "pos", "rsID", "ancestral", "derived"], "Skov")
    df["ref"] = df["alt"] = pd.NA                       # ancestral is not necessarily the hg19 reference
    df["allele_1"], df["allele_2"] = df["ancestral"], df["derived"]
    df["neanderthal_allele"] = df["derived"]           # Skov only calls Neanderthal-derived alleles
    return df[STD_COLS]


def read_vernot(path):
    df = read_map(path, ["chr", "pos", "rsID", "ref", "alt", "neanderthal_allele"], "Vernot")
    het = df["neanderthal_allele"].str.contains("/", na=False)
    print(f"Vernot: {het.sum():,} heterozygous Neanderthal genotypes set to NA", file=sys.stderr)
    df.loc[het, "neanderthal_allele"] = pd.NA
    df["allele_1"], df["allele_2"] = df["ref"], df["alt"]
    return df[STD_COLS]


def first_and_count(long, col, keys=("chr", "pos")):
    """First non-missing value of `col` per position and the number of distinct values."""
    sub = long.loc[long[col].notna(), [*keys, col]].drop_duplicates()
    out = sub.groupby(list(keys), sort=False)[col].agg(["first", "size"])
    return out.rename(columns={"first": col, "size": f"n_{col}"})


def merge_maps(browning, sankararaman, skov, vernot):
    long = pd.concat([read_browning(browning), read_sankararaman(sankararaman),
                      read_skov(skov), read_vernot(vernot)], ignore_index=True)

    # Clean
    long["chr"] = long["chr"].str.replace(r"^chr", "", case=False, regex=True)
    long["pos"] = long["pos"].astype("int64")
    long.loc[long["rsID"].isin([".", ""]), "rsID"] = pd.NA
    for col in ["ref", "alt", "allele_1", "allele_2", "neanderthal_allele"]:
        long[col] = long[col].str.upper()
    a1, a2 = long["allele_1"], long["allele_2"]
    both = a1.notna() & a2.notna()
    long["pair"] = pd.NA
    long.loc[both, "pair"] = (a1[both].where(a1[both] <= a2[both], a2[both]) + "/" +
                              a2[both].where(a1[both] <= a2[both], a1[both]))

    # Per-position summaries (drop_duplicates + groupby first/size: vectorized, fast)
    keys = ["chr", "pos"]
    rsids = first_and_count(long, "rsID")
    refalt = (long.loc[long["ref"].notna(), keys + ["ref", "alt"]].drop_duplicates()
              .groupby(keys, sort=False)[["ref", "alt"]].first())
    pairs = first_and_count(long, "pair").rename(columns={"pair": "alleles", "n_pair": "n_pairs"})
    nea = first_and_count(long, "neanderthal_allele")
    nea.loc[nea["n_neanderthal_allele"] > 1, "neanderthal_allele"] = pd.NA

    maps = long[keys + ["map"]].drop_duplicates()
    maps["bit"] = 2 ** maps["map"].map({m: i for i, m in enumerate(MAP_NAMES)})
    maps = maps.groupby(keys, sort=False)["bit"].agg(code="sum", n_maps="size")
    lookup = {c: ", ".join(m for i, m in enumerate(MAP_NAMES) if c & (1 << i)) for c in range(16)}
    maps["maps"] = maps["code"].map(lookup)

    out = (maps.drop(columns="code")
           .join([rsids, refalt, pairs, nea], how="left")
           .reset_index())
    for col in ["n_rsID", "n_pairs", "n_neanderthal_allele"]:
        out[col] = out[col].fillna(0).astype(int)

    # Notes (later notes overwrite earlier ones)
    out["notes"] = pd.NA
    out.loc[out["n_pairs"] > 1, "notes"] = "Maps report different alleles (multi-allelic site?)"
    out.loc[out["n_rsID"] > 1, "notes"] = "Maps report different rsIDs"
    out.loc[out["n_neanderthal_allele"] > 1, "notes"] = "Maps disagree on the Neanderthal allele"
    pair_a = out["alleles"].str.split("/").str[0]
    pair_b = out["alleles"].str.split("/").str[1]
    bad = (out["neanderthal_allele"].notna() & out["alleles"].notna() &
           (out["neanderthal_allele"] != pair_a) & (out["neanderthal_allele"] != pair_b))
    out.loc[bad, "notes"] = "Neanderthal allele is not one of the variant's alleles"

    out = out.rename(columns={"pos": "pos_hg19"})
    out["chr_order"] = pd.to_numeric(out["chr"], errors="coerce")
    out = out.sort_values(["chr_order", "chr", "pos_hg19"]).drop(columns="chr_order")
    out = out[["chr", "pos_hg19", "rsID", "ref", "alt", "alleles", "neanderthal_allele",
               "maps", "n_maps", "notes"]]

    # Report
    no_rs_any_map = long.loc[long["rsID"].isna(), keys].drop_duplicates().shape[0]
    print(f"\n{len(out):,} unique variants (chr + pos)", file=sys.stderr)
    print(f"  no rsID in any map:              {out['rsID'].isna().sum():,}", file=sys.stderr)
    print(f"  rsID filled in from another map: {no_rs_any_map - out['rsID'].isna().sum():,}", file=sys.stderr)
    print(f"  no alleles (Sankararaman-only):  {out['alleles'].isna().sum():,}", file=sys.stderr)
    print(f"  Neanderthal allele known:        {out['neanderthal_allele'].notna().sum():,}", file=sys.stderr)
    print("  variants per map:\n" + long[keys + ["map"]].drop_duplicates()["map"].value_counts().to_string(),
          file=sys.stderr)
    print("  notes:\n" + out["notes"].value_counts().to_string(), file=sys.stderr)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--browning", required=True)
    parser.add_argument("--sankararaman", required=True)
    parser.add_argument("--skov", required=True)
    parser.add_argument("--vernot", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    out = merge_maps(args.browning, args.sankararaman, args.skov, args.vernot)
    out.to_csv(args.output, sep="\t", index=False, na_rep="NA")
    print(f"Wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()