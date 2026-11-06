"""
Retrieve GTEx eQTL statistics and per-sample expression for credible set variants.

For the variants and genes in the config, three tables are written:

- gtex_significant_eqtls.txt: every gene-tissue pair in which a queried variant
  is a significant single-tissue eQTL, with the GTEx lead variant for that gene
  and tissue (the supplementary table and the all-eQTL panel);
- gtex_focus_gene_effects.txt: the effect of each queried variant on each focus
  gene in every GTEx tissue, significant or not (the tissue-specificity panel);
- gtex_expression_by_genotype.txt.gz: normalized expression per sample for every
  significant pair (the genotype-level panels).

GTEx reports effects for the alternate allele of its variant ID
(chr_pos_ref_alt_b38). Every effect here is oriented to the allele of interest
in the config instead: the sign is flipped when that allele is the reference
allele, and genotypes are recoded as copies of it. GTEx reports the t statistic
but not the standard error, so the SE is NES / t.

The GTEx API needs internet access, so run this on a Wynton dev node.

Input:  the YAML config (eqtl_config.yaml).
Output: the three tables above, in --output_directory.
"""

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd
import yaml


API = "https://gtexportal.org/api/v2"
PAGE_SIZE = 250
PAUSE_SECONDS = 0.1
RETRIES = 4
GENCODE_VERSIONS = {"gtex_v8": "v26", "gtex_v10": "v39"}


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", required=True, help="path to eqtl_config.yaml")
    parser.add_argument(
        "--dataset",
        default="gtex_v8",
        choices=sorted(GENCODE_VERSIONS),
        help="GTEx release (default: gtex_v8)",
    )
    parser.add_argument(
        "--output_directory", required=True, help="directory for the output tables"
    )
    return parser.parse_args()


def get(endpoint, **params):
    """One GTEx API request, retrying transient failures."""
    url = f"{API}/{endpoint}?{urllib.parse.urlencode(params)}"
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                result = json.load(response)
            time.sleep(PAUSE_SECONDS)
            return result
        except urllib.error.HTTPError as error:
            transient = error.code in (429, 500, 502, 503, 504)
            if transient and attempt < RETRIES - 1:
                time.sleep(2**attempt)
                continue
            body = error.read()[:300].decode(errors="replace")
            raise LookupError(f"HTTP {error.code} from {url}: {body}") from error
        except urllib.error.URLError:
            if attempt < RETRIES - 1:
                time.sleep(2**attempt)
                continue
            raise


def get_all(endpoint, **params):
    """Every row of a paginated GTEx endpoint."""
    rows, page = [], 0
    while True:
        result = get(endpoint, page=page, itemsPerPage=PAGE_SIZE, **params)
        rows.extend(result["data"])
        n_pages = result.get("paging_info", {}).get("numberOfPages", 1)
        page += 1
        if page >= n_pages:
            return rows


def resolve_variant(rsid, dataset):
    """GTEx variant ID (chr_pos_ref_alt_b38) for an rsID."""
    matches = get_all("dataset/variant", snpId=rsid, datasetId=dataset)
    if len(matches) != 1:
        raise LookupError(f"{rsid}: expected one GTEx variant, found {len(matches)}")
    return matches[0]["variantId"]


def resolve_gene(gene, dataset):
    """GENCODE ID and symbol for a gene symbol or (unversioned) Ensembl ID."""
    matches = get_all(
        "reference/gene",
        geneId=gene,
        gencodeVersion=GENCODE_VERSIONS[dataset],
        genomeBuild="GRCh38/hg38",
    )
    exact = [
        match for match in matches
        if match["geneSymbol"] == gene or match["gencodeId"].split(".")[0] == gene
    ]
    if len(exact) != 1:
        raise LookupError(f"{gene}: expected one GENCODE gene, found {len(exact)}")
    return exact[0]["gencodeId"], exact[0]["geneSymbol"]


def eqtl_tissues(dataset):
    """Tissues with eQTL results in this release."""
    tissues = get_all("dataset/tissueSiteDetail", datasetId=dataset)
    return sorted(
        tissue["tissueSiteDetailId"]
        for tissue in tissues
        if tissue.get("hasEGenes", True)
    )


def orientation(variant_id, allele):
    """+1 if the allele of interest is GTEx's alternate allele, -1 if reference."""
    _, _, ref, alt, _ = variant_id.split("_")
    if allele == alt:
        return 1
    if allele == ref:
        return -1
    raise ValueError(f"{allele} is neither allele of {variant_id}")


def dynamic_eqtl(gencode_id, variant_id, tissue, dataset):
    """Single-tissue statistics and per-sample expression for one pair, or None
    when GTEx did not test the pair (e.g. the variant fails its MAF filter)."""
    try:
        return get(
            "association/dyneqtl",
            gencodeId=gencode_id,
            variantId=variant_id,
            tissueSiteDetailId=tissue,
            datasetId=dataset,
        )
    except LookupError as error:
        print(f"  not tested: {error}", file=sys.stderr)
        return None


def oriented_statistics(result, sign):
    nes = result["nes"] * sign
    t_statistic = result.get("tStatistic")
    se = abs(result["nes"] / t_statistic) if t_statistic else None
    return {"nes": nes, "se": se, "p_value": result["pValue"]}


def gtex_lead(gencode_id, tissue, dataset, cache):
    """The most significant eQTL variant for a gene in a tissue."""
    key = (gencode_id, tissue)
    if key not in cache:
        hits = get_all(
            "association/singleTissueEqtl",
            gencodeId=gencode_id,
            tissueSiteDetailId=tissue,
            datasetId=dataset,
        )
        lead = min(hits, key=lambda hit: hit["pValue"])
        cache[key] = (lead["variantId"], lead.get("snpId"), lead["pValue"])
    return cache[key]


def main():
    args = parse_args()
    with open(args.config) as f:
        config = yaml.safe_load(f)
    os.makedirs(args.output_directory, exist_ok=True)

    tissues = eqtl_tissues(args.dataset)
    print(f"{len(tissues)} GTEx tissues with eQTL results in {args.dataset}")

    significant_rows, focus_rows, expression_frames = [], [], []
    lead_cache = {}

    for locus, entry in config.items():
        focus_genes = [resolve_gene(gene, args.dataset) for gene in entry["focus_genes"]]

        for variant in entry["variants"] or []:
            rsid, allele = variant["rsID"], variant["allele"]
            variant_id = resolve_variant(rsid, args.dataset)
            sign = orientation(variant_id, allele)
            labels = {
                "locus": locus,
                "signal": variant["signal"],
                "rsID": rsid,
                "variant_id": variant_id,
                "allele": allele,
            }
            print(f"{locus} {rsid} ({variant_id}), effects oriented to {allele}")

            hits = get_all(
                "association/singleTissueEqtl",
                variantId=variant_id,
                datasetId=args.dataset,
            )
            significant_pairs = set()
            for hit in hits:
                gencode_id, tissue = hit["gencodeId"], hit["tissueSiteDetailId"]
                significant_pairs.add((gencode_id, tissue))
                result = dynamic_eqtl(gencode_id, variant_id, tissue, args.dataset)
                if result is None:
                    continue
                lead_id, lead_rsid, lead_p = gtex_lead(
                    gencode_id, tissue, args.dataset, lead_cache
                )
                significant_rows.append({
                    **labels,
                    "gene_symbol": hit["geneSymbol"],
                    "gencode_id": gencode_id,
                    "tissue": tissue,
                    **oriented_statistics(result, sign),
                    "gtex_lead_variant": lead_id,
                    "gtex_lead_rsID": lead_rsid,
                    "gtex_lead_p_value": lead_p,
                    "is_gtex_lead": lead_id == variant_id,
                })

                genotypes = pd.Series(result["genotypes"], dtype=int)
                expression = pd.Series(result["data"], dtype=float)
                if len(genotypes) != len(expression):
                    raise ValueError(f"{rsid} {tissue}: genotypes and expression differ in length")
                expression_frames.append(pd.DataFrame({
                    "locus": locus,
                    "rsID": rsid,
                    "allele": allele,
                    "gene_symbol": hit["geneSymbol"],
                    "tissue": tissue,
                    # genotypes count alternate alleles
                    "copies_of_allele": genotypes if sign == 1 else 2 - genotypes,
                    "expression": expression,
                }))
            print(f"  {len(hits)} significant gene-tissue pairs")

            for gencode_id, gene_symbol in focus_genes:
                for tissue in tissues:
                    result = dynamic_eqtl(gencode_id, variant_id, tissue, args.dataset)
                    focus_rows.append({
                        **labels,
                        "gene_symbol": gene_symbol,
                        "gencode_id": gencode_id,
                        "tissue": tissue,
                        **(
                            oriented_statistics(result, sign) if result
                            else {"nes": None, "se": None, "p_value": None}
                        ),
                        "significant": (gencode_id, tissue) in significant_pairs,
                    })

    outputs = {
        "gtex_significant_eqtls.txt": pd.DataFrame(significant_rows),
        "gtex_focus_gene_effects.txt": pd.DataFrame(focus_rows),
        "gtex_expression_by_genotype.txt.gz": (
            pd.concat(expression_frames, ignore_index=True)
            if expression_frames else pd.DataFrame()
        ),
    }
    for file_name, table in outputs.items():
        if len(table):
            table.insert(0, "dataset", args.dataset)
        path = os.path.join(args.output_directory, file_name)
        table.to_csv(path, sep="\t", index=False, na_rep="NA")
        print(f"Wrote {len(table)} rows to {path}")


if __name__ == "__main__":
    main()
