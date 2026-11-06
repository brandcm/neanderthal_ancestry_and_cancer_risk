#!/usr/bin/env python3
"""
Analyze a polymorphic LINE-1 element near PTHLH from long-read genome assemblies.

Input data (not distributed; shared via personal communication):
    data/PTHLH_long_reads/fasta/       -- 578 FASTA files, one per phased haplotype
    data/PTHLH_long_reads/reference/   -- T2T reference FASTA for PTHLH region
    data/PTHLH_long_reads/summary/     -- per-haplotype summary TSVs
    data/PTHLH_long_reads/dotplots/    -- dot-plot images (optional)

Usage:
    python3 analyze_PTHLH_LINE_insertion.py \
        --fasta-dir data/PTHLH_long_reads/fasta \
        --reference data/PTHLH_long_reads/reference/PTHLH.t2t.fa \
        --outdir results/PTHLH_LINE_analysis
"""

import argparse
import glob
import os
import sys
from collections import defaultdict


def read_fasta(path):
    """Read a single-record FASTA file and return the uppercase sequence."""
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def classify_haplotypes(fasta_dir, ref_len):
    """Classify haplotypes by interval size relative to reference.

    Returns:
        dict mapping haplotype name to (sequence length, size difference)
        dict mapping category to list of haplotype names
    """
    sizes = {}
    categories = defaultdict(list)

    for ff in sorted(glob.glob(os.path.join(fasta_dir, "*.fa"))):
        name = os.path.basename(ff).replace(".PTHLH.fa", "")
        seq = read_fasta(ff)
        diff = len(seq) - ref_len
        sizes[name] = (len(seq), diff)

        if diff > 5000:
            categories["LINE_insertion"].append(name)
        elif diff > 2000:
            categories["medium_insertion"].append(name)
        elif diff > -2000:
            categories["reference_like"].append(name)
        elif diff > -5000:
            categories["deletion"].append(name)
        else:
            categories["other"].append(name)

    return sizes, categories


def find_breakpoint_kmer(ref, hap_seq, region_start=245000, region_end=245400,
                         k=20, coarse_step=20, fine_step=1):
    """Locate the LINE insertion breakpoint using k-mer anchoring.

    First pass: coarse scan with large step to find the approximate breakpoint.
    Second pass: fine scan with step=1 around the jump to get single-base resolution.

    Returns:
        (last_pre_ref, first_post_ref, pre_offset, post_offset)
    """
    def scan(start, end, step):
        pre_anchors = []
        post_anchors = []
        for rpos in range(start, end, step):
            kmer = ref[rpos:rpos + k]
            if len(kmer) < k:
                continue
            # Search pre-insertion range (small offset)
            idx = hap_seq.find(kmer, rpos, rpos + 1000)
            if idx >= 0 and (idx - rpos) < 1000:
                pre_anchors.append((rpos, idx, idx - rpos))
            # Search post-insertion range (large offset)
            idx2 = hap_seq.find(kmer, rpos + 5000, rpos + 8000)
            if idx2 >= 0:
                post_anchors.append((rpos, idx2, idx2 - rpos))
        return pre_anchors, post_anchors

    # Coarse pass
    pre_anchors, post_anchors = scan(region_start, region_end, coarse_step)

    if not pre_anchors or not post_anchors:
        return None

    # Narrow to the transition zone
    last_pre_coarse = max(pre_anchors, key=lambda x: x[0])
    first_post_coarse = min(post_anchors, key=lambda x: x[0])

    # Fine pass around the transition
    fine_start = max(region_start, last_pre_coarse[0] - coarse_step)
    fine_end = min(region_end, first_post_coarse[0] + coarse_step)
    pre_fine, post_fine = scan(fine_start, fine_end, fine_step)

    if pre_fine:
        last_pre = max(pre_fine, key=lambda x: x[0])
    else:
        last_pre = last_pre_coarse

    if post_fine:
        first_post = min(post_fine, key=lambda x: x[0])
    else:
        first_post = first_post_coarse

    return last_pre[0], first_post[0], last_pre[2], first_post[2]


def extract_insertion(hap_seq, last_pre_ref, first_post_ref, pre_offset, post_offset, k=20):
    """Extract the inserted sequence from the haplotype."""
    ins_start = last_pre_ref + pre_offset + k
    ins_end = first_post_ref + post_offset
    return hap_seq[ins_start:ins_end]


def find_tsd(hap_seq, ins_start, ins_end, max_tsd=20, search_margin=10):
    """Search for target site duplication flanking the insertion.

    Allows a small positional margin around the nominal breakpoints to
    account for k-mer step size in breakpoint localization.
    """
    for tsd_len in range(max_tsd, 4, -1):
        for offset_pre in range(-search_margin, search_margin + 1):
            for offset_post in range(-search_margin, search_margin + 1):
                pre_seq = hap_seq[ins_start + offset_pre - tsd_len:ins_start + offset_pre]
                post_seq = hap_seq[ins_end + offset_post:ins_end + offset_post + tsd_len]
                if pre_seq == post_seq and len(pre_seq) == tsd_len:
                    return tsd_len, pre_seq
    return 0, ""


def check_orientation(inserted_seq):
    """Determine LINE orientation from poly-A/T content at ends."""
    t_start = inserted_seq[:30].count("T")
    a_end = inserted_seq[-30:].count("A")

    if t_start > 20:
        return "antisense"
    elif a_end > 20:
        return "sense"
    else:
        return "unknown"


def genotype_summary(sizes, categories):
    """Compute diploid genotype counts from phased haplotype data.

    Only samples with both haplotypes assembled are genotyped.
    Samples with a single haplotype are counted separately.
    """
    line_haps = set(categories["LINE_insertion"])

    samples = defaultdict(list)
    for name, (seqlen, diff) in sizes.items():
        parts = name.rsplit("_", 1)
        if len(parts) == 2:
            sample, hap = parts
            samples[sample].append(name in line_haps)

    diploid = {s: haps for s, haps in samples.items() if len(haps) == 2}
    haploid = {s: haps for s, haps in samples.items() if len(haps) == 1}

    hom_ins = sum(1 for s, haps in diploid.items() if all(haps))
    het = sum(1 for s, haps in diploid.items() if any(haps) and not all(haps))
    hom_ref = sum(1 for s, haps in diploid.items() if not any(haps))

    return len(samples), len(diploid), hom_ins, het, hom_ref, len(haploid)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta-dir", required=True,
                        help="Directory containing per-haplotype FASTA files")
    parser.add_argument("--reference", required=True,
                        help="T2T reference FASTA for the PTHLH region")
    parser.add_argument("--outdir", default="results/PTHLH_LINE_analysis",
                        help="Output directory (default: results/PTHLH_LINE_analysis)")
    parser.add_argument("--exemplar", default=None,
                        help="Haplotype to use for breakpoint/insertion analysis "
                             "(default: first LINE-carrying haplotype)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Read reference
    print("Reading reference...", file=sys.stderr)
    ref = read_fasta(args.reference)
    ref_len = len(ref)
    print(f"  Reference length: {ref_len:,} bp", file=sys.stderr)

    # Classify haplotypes
    print("Classifying haplotypes...", file=sys.stderr)
    sizes, categories = classify_haplotypes(args.fasta_dir, ref_len)
    total = len(sizes)

    n_line = len(categories["LINE_insertion"])
    n_med = len(categories["medium_insertion"])
    n_ref = len(categories["reference_like"])
    n_del = len(categories["deletion"])
    n_other = len(categories["other"])

    print(f"  Total haplotypes: {total}", file=sys.stderr)
    print(f"  LINE insertion (>5 kb): {n_line} ({n_line/total:.1%})", file=sys.stderr)
    print(f"  Medium insertion (2-5 kb): {n_med}", file=sys.stderr)
    print(f"  Reference-like (±2 kb): {n_ref}", file=sys.stderr)
    print(f"  Deletion (2-5 kb shorter): {n_del}", file=sys.stderr)
    print(f"  Other: {n_other}", file=sys.stderr)

    # Insertion size stats
    line_sizes = [sizes[n][1] for n in categories["LINE_insertion"]]
    mean_size = sum(line_sizes) / len(line_sizes)
    min_size = min(line_sizes)
    max_size = max(line_sizes)

    # Diploid genotypes
    n_samples, n_diploid, hom_ins, het, hom_ref, n_haploid = genotype_summary(sizes, categories)
    af = n_line / total

    print(f"\n  Allele frequency (insertion): {af:.3f}", file=sys.stderr)
    print(f"  Total samples: {n_samples} ({n_diploid} diploid, {n_haploid} single-haplotype)", file=sys.stderr)
    print(f"  Diploid genotypes: {hom_ins} hom-ins, {het} het, {hom_ref} hom-ref", file=sys.stderr)

    # Breakpoint analysis on exemplar haplotype
    exemplar = args.exemplar or categories["LINE_insertion"][0]
    print(f"\nAnalyzing breakpoint in {exemplar}...", file=sys.stderr)
    exemplar_path = os.path.join(args.fasta_dir, f"{exemplar}.PTHLH.fa")
    hap_seq = read_fasta(exemplar_path)

    bp = find_breakpoint_kmer(ref, hap_seq)
    if bp is None:
        print("  WARNING: could not locate breakpoint", file=sys.stderr)
        tsd_len, tsd_seq, orientation, ins_len = 0, "", "unknown", 0
    else:
        last_pre_ref, first_post_ref, pre_offset, post_offset = bp
        net_insertion = post_offset - pre_offset - (first_post_ref - last_pre_ref)
        print(f"  Last pre-insertion anchor: ref {last_pre_ref}", file=sys.stderr)
        print(f"  First post-insertion anchor: ref {first_post_ref}", file=sys.stderr)
        print(f"  Net insertion size: ~{net_insertion:,} bp", file=sys.stderr)

        # Extract insertion
        k = 20
        ins_start = last_pre_ref + pre_offset + k
        ins_end = first_post_ref + post_offset
        inserted = hap_seq[ins_start:ins_end]
        ins_len = len(inserted)

        # TSD
        tsd_len, tsd_seq = find_tsd(hap_seq, ins_start, ins_end)
        if tsd_len:
            print(f"  Target site duplication: {tsd_len} bp ({tsd_seq})", file=sys.stderr)
        else:
            print("  No TSD found", file=sys.stderr)

        # Orientation
        orientation = check_orientation(inserted)
        print(f"  Orientation: {orientation}", file=sys.stderr)

        # Save insertion FASTA
        ins_fa = os.path.join(args.outdir, f"{exemplar}_LINE_insertion.fa")
        with open(ins_fa, "w") as f:
            f.write(f">{exemplar}_LINE_insertion_chr12_27945230_hg19\n")
            for i in range(0, len(inserted), 80):
                f.write(inserted[i:i + 80] + "\n")
        print(f"  Saved insertion sequence to {ins_fa}", file=sys.stderr)

    # Write summary TSV
    summary_path = os.path.join(args.outdir, "LINE_insertion_summary.tsv")
    with open(summary_path, "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"total_haplotypes\t{total}\n")
        f.write(f"LINE_insertion_count\t{n_line}\n")
        f.write(f"insertion_allele_frequency\t{af:.4f}\n")
        f.write(f"mean_insertion_size_bp\t{mean_size:.0f}\n")
        f.write(f"min_insertion_size_bp\t{min_size}\n")
        f.write(f"max_insertion_size_bp\t{max_size}\n")
        f.write(f"total_samples\t{n_samples}\n")
        f.write(f"diploid_samples\t{n_diploid}\n")
        f.write(f"single_haplotype_samples\t{n_haploid}\n")
        f.write(f"homozygous_insertion\t{hom_ins}\n")
        f.write(f"heterozygous\t{het}\n")
        f.write(f"homozygous_reference\t{hom_ref}\n")
        f.write(f"orientation\t{orientation}\n")
        f.write(f"TSD_length_bp\t{tsd_len}\n")
        f.write(f"TSD_sequence\t{tsd_seq}\n")
        f.write(f"extracted_insertion_length_bp\t{ins_len}\n")
        f.write(f"exemplar_haplotype\t{exemplar}\n")
    print(f"\nSummary written to {summary_path}", file=sys.stderr)

    # Write per-haplotype classification
    class_path = os.path.join(args.outdir, "haplotype_classification.tsv")
    with open(class_path, "w") as f:
        f.write("haplotype\tsequence_length\tsize_difference\tcategory\n")
        for name in sorted(sizes.keys()):
            seqlen, diff = sizes[name]
            cat = "LINE_insertion" if name in categories["LINE_insertion"] else \
                  "medium_insertion" if name in categories["medium_insertion"] else \
                  "reference_like" if name in categories["reference_like"] else \
                  "deletion" if name in categories["deletion"] else "other"
            f.write(f"{name}\t{seqlen}\t{diff}\t{cat}\n")
    print(f"Per-haplotype classification written to {class_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
