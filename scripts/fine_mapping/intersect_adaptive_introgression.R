#!/usr/bin/env Rscript
# Intersect credible set variants with high-frequency archaic ("adaptive
# introgression") regions, hg19.
#
#   ./run_rscript.sh intersect_adaptive_introgression.R
#
# Reports, for genome-wide significant, representative signals only:
#   - how many signals have at least one member inside an AI region
#   - how many loci that corresponds to
#   - the same restricted to Neanderthal-majority signals
#   - a first-pass comparison against all variants in the same windows
#
# The window comparison is NOT a null model: credible set variants are
# conditioned on being genome-wide significant and on LD structure, neither of
# which is matched here. It indicates whether a proper frequency- and
# LD-matched test is worth building, nothing more.
#
# Reads the deduplicated tables, so a signal recovered through two overlapping
# index windows is counted once.

suppressPackageStartupMessages({
    library(data.table)
    library(glue)
})

project <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/"
fm <- file.path(project, "data/fine_mapping")
ai_file <- file.path(project, "data/adaptive_introgression",
                     "Iasi_et_al_2025_high_frequency_ancient_modern_regions.bed")

regions <- fread(ai_file, header = FALSE, col.names = c("chr", "start", "end"))
regions[, chr := sub("^chr", "", chr)]
setkey(regions, chr, start, end)
message(glue("{nrow(regions)} adaptive introgression regions covering ",
             "{round(sum(regions$end - regions$start) / 1e6)} Mb."))

summary_file <- file.path(fm, "credible_set_summary_deduplicated.txt")
variants_file <- file.path(fm, "credible_set_variants_deduplicated.txt")
for (file in c(summary_file, variants_file)) {
    if (!file.exists(file)) {
        stop(glue("{basename(file)} does not exist. ",
                  "Run deduplicate_credible_sets.R first."))
    }
}

variants <- fread(variants_file, na.strings = c("NA", ""))
sets <- fread(summary_file, na.strings = c("NA", ""))

variants[, chr := as.character(chr)]
variants[, `:=`(start = pos, end = pos)]
setkey(variants, chr, start, end)

overlaps <- foverlaps(variants, regions, nomatch = NULL)
variants[, in_ai := FALSE]
if (nrow(overlaps) > 0) {
    key_overlap <- unique(overlaps[, .(cancer, ancestry_group, locus_id, CS, rsID)])
    variants[key_overlap, in_ai := TRUE,
             on = .(cancer, ancestry_group, locus_id, CS, rsID)]
}

# per-variant flags, for joining onto the variant-level supplementary table;
# per_set below aggregates these to one row per credible set
variant_file <- file.path(fm, "adaptive_introgression_variants.txt")
fwrite(
    variants[, .(cancer, ancestry_group, locus_id, CS, rsID, chr, pos,
                 is_neanderthal, in_ai)],
    variant_file, sep = "\t", na = "NA"
)
message(glue("Wrote {variant_file}"))

per_set <- variants[, .(
    n_variants = .N,
    n_in_ai = sum(in_ai),
    n_neanderthal_in_ai = sum(in_ai & is_neanderthal == TRUE, na.rm = TRUE)
), by = .(cancer, ancestry_group, locus_id, CS)]

per_set <- merge(
    per_set,
    sets[, .(cancer, ancestry_group, locus_id, CS, signal_id, is_representative,
             p_neanderthal, max_abs_z, purity, genome_wide_significant)],
    by = c("cancer", "ancestry_group", "locus_id", "CS"), all.x = TRUE
)

# one row per signal: the representative set, which is the strongest of its group
significant <- per_set[genome_wide_significant == TRUE & is_representative == TRUE]
archaic <- significant[p_neanderthal > 0.5]

count_loci <- function(dt) uniqueN(dt[, .(cancer, ancestry_group, locus_id)])

message("")
message("--- genome-wide significant signals ---")
message(glue("signals: {nrow(significant)} across {count_loci(significant)} loci"))
message(glue("signals with >=1 variant in an AI region: {nrow(significant[n_in_ai > 0])} ",
             "across {count_loci(significant[n_in_ai > 0])} loci"))

message("")
message("--- of those, Neanderthal-majority (p_neanderthal > 0.5) ---")
message(glue("signals: {nrow(archaic)} across {count_loci(archaic)} loci"))
message(glue("signals with >=1 variant in an AI region: {nrow(archaic[n_in_ai > 0])} ",
             "across {count_loci(archaic[n_in_ai > 0])} loci"))

message("")
message("--- first-pass comparison (not a null model) ---")
cs_rate <- variants[, mean(in_ai)]
message(glue("credible set variants inside an AI region: ",
             "{round(100 * cs_rate, 2)}% ({sum(variants$in_ai)} of {nrow(variants)})"))
genome_fraction <- sum(as.numeric(regions$end - regions$start)) / 3.1e9
message(glue("AI regions cover roughly {round(100 * genome_fraction, 2)}% of the genome"))

message("")
message("--- Neanderthal-majority signals overlapping an AI region ---")
hits <- archaic[n_in_ai > 0][order(-n_in_ai)]
if (nrow(hits) > 0) {
    print(hits[, .(cancer, anc = ancestry_group, locus_id, CS,
                   n = n_variants, in_ai = n_in_ai, nea_in_ai = n_neanderthal_in_ai,
                   p_nea = round(p_neanderthal, 2), maxZ = round(max_abs_z, 1))],
          nrows = 60)
} else {
    message("none")
}

out_file <- file.path(fm, "adaptive_introgression_overlap.txt")
fwrite(per_set, out_file, sep = "\t", na = "NA")
message("")
message(glue("Wrote {out_file}"))
