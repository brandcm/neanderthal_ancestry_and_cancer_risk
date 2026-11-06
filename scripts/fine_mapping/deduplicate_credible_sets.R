#!/usr/bin/env Rscript
# Collapse credible sets that are the same signal found through overlapping
# index windows, and write the Neanderthal subsets.
#
#   ./run_rscript.sh deduplicate_credible_sets.R [jaccard_threshold] [posterior_threshold]
#
# Index loci were defined around significant introgressed variants, so two
# variants a few hundred kb apart produce two windows that both contain the same
# signal. PTHLH is the worked example: 12_28153018 and 12_28658479 return the
# same L1/L2/L3 with the same lead variants, so every count that treats them as
# distinct is inflated.
#
# Sets are judged the same signal when they come from different index loci in
# the same GWAS and ancestry group and share most of their variants. Membership
# overlap is used rather than LD between leads because it needs no genotype
# panel and it is what downstream counts actually consume. The default threshold
# is deliberately lenient; a pair at 0.5 shares half its membership, which does
# not happen by chance between genuinely distinct signals at this scale.
#
# Writes, from credible_set_*_selected.txt:
#   credible_set_summary_deduplicated.txt   every set, with signal_id and
#                                           is_representative; nothing discarded
#   credible_set_variants_deduplicated.txt  the same columns carried onto the
#                                           variants, so nothing downstream has to
#                                           pair a summary with variants from a
#                                           different stage
#   neanderthal_credible_sets.txt           variants of the genome-wide
#                                           significant, representative signals
#                                           containing at least one Neanderthal
#                                           variant
#   neanderthal_majority_signals.txt        the same, restricted to signals whose
#                                           posterior is mostly on Neanderthal
#                                           variants

suppressPackageStartupMessages({
    library(data.table)
    library(glue)
})

project_directory <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/"
fine_mapping_directory <- file.path(project_directory, "data/fine_mapping")

arguments <- commandArgs(trailingOnly = TRUE)
threshold <- suppressWarnings(as.numeric(arguments[1]))
if (is.na(threshold)) threshold <- 0.5
posterior_threshold <- suppressWarnings(as.numeric(arguments[2]))
if (is.na(posterior_threshold)) posterior_threshold <- 0.5

# The summary and the variants must come from the same stage: pairing a refitted
# summary with the variants of the fits it replaced would compute overlaps between
# sets that were never fitted together, which is silently wrong rather than an
# error. select_fits.R writes both together.
summary_file <- file.path(fine_mapping_directory, "credible_set_summary_selected.txt")
variants_file <- file.path(fine_mapping_directory, "credible_set_variants_selected.txt")
for (file in c(summary_file, variants_file)) {
    if (!file.exists(file)) {
        stop(glue("{basename(file)} does not exist. Run select_fits.R first."))
    }
}
message(glue("Reading {basename(summary_file)} and {basename(variants_file)}"))

sets <- fread(summary_file, na.strings = c("NA", ""))
variants <- fread(variants_file, na.strings = c("NA", ""))

if (!"set_key" %in% names(sets)) {
    sets[, set_key := paste(cancer, ancestry_group, locus_id, CS, sep = "|")]
}
if (!"set_key" %in% names(variants)) {
    variants[, set_key := paste(cancer, ancestry_group, locus_id, CS, sep = "|")]
}

members <- split(variants$rsID, variants$set_key)

# --- union-find over sets that share most of their variants ---------------
keys <- sets$set_key
index_of <- setNames(seq_along(keys), keys)
parents <- seq_along(keys)

root_of <- function(i) {
    while (parents[i] != i) {
        parents[i] <<- parents[parents[i]]
        i <- parents[i]
    }
    i
}

join <- function(i, j) {
    ri <- root_of(i); rj <- root_of(j)
    if (ri != rj) parents[rj] <<- ri
}

merged_pairs <- list()

for (group_keys in split(keys, sub("\\|[^|]*\\|[^|]*$", "", keys))) {
    # group_keys: all sets within one cancer + ancestry group
    if (length(group_keys) < 2) next
    for (a in seq_len(length(group_keys) - 1)) {
        for (b in (a + 1):length(group_keys)) {
            key_a <- group_keys[a]; key_b <- group_keys[b]
            locus_a <- strsplit(key_a, "\\|")[[1]][3]
            locus_b <- strsplit(key_b, "\\|")[[1]][3]
            if (locus_a == locus_b) next          # same window, distinct effects

            variants_a <- members[[key_a]]; variants_b <- members[[key_b]]
            if (is.null(variants_a) || is.null(variants_b)) next
            shared <- length(intersect(variants_a, variants_b))
            if (shared == 0) next
            overlap <- shared / length(union(variants_a, variants_b))
            if (overlap >= threshold) {
                join(index_of[[key_a]], index_of[[key_b]])
                merged_pairs[[length(merged_pairs) + 1]] <- data.table(
                    set_a = key_a, set_b = key_b,
                    shared = shared, jaccard = round(overlap, 3)
                )
            }
        }
    }
}

sets[, signal_id := sapply(seq_along(keys), function(i) root_of(i))]
sets[, signal_id := as.integer(factor(signal_id))]

# keep the strongest set as the representative of each signal
setorder(sets, signal_id, -max_abs_z, -purity)
sets[, is_representative := seq_len(.N) == 1, by = signal_id]

# carry both onto the variants, so a variant-level consumer never has to join
# back to the summary to know which signal a row belongs to
variants <- merge(
    variants,
    sets[, .(set_key, signal_id, is_representative)],
    by = "set_key", all.x = TRUE
)

# --- report ----------------------------------------------------------------
message(glue("Jaccard threshold: {threshold}"))
message(glue("credible sets: {nrow(sets)} -> {uniqueN(sets$signal_id)} distinct signals"))

significant <- sets[genome_wide_significant == TRUE]
message(glue("genome-wide significant: {nrow(significant)} sets -> ",
             "{uniqueN(significant$signal_id)} signals"))

representative <- significant[is_representative == TRUE]
introgressed <- representative[n_neanderthal > 0]
majority <- representative[p_neanderthal > posterior_threshold]

message(glue("containing a Neanderthal variant: {nrow(introgressed)} signals across ",
             "{uniqueN(introgressed[, .(cancer, ancestry_group)])} GWAS"))
message(glue("Neanderthal-majority (p_neanderthal > {posterior_threshold}): ",
             "{nrow(majority)} signals across ",
             "{uniqueN(majority[, .(cancer, ancestry_group)])} GWAS"))

message("")
message("=== Neanderthal-majority signals by cancer ===")
print(majority[, .(signals = .N, loci = uniqueN(locus_id)),
               by = .(cancer, ancestry_group)][order(-signals)])

if (length(merged_pairs) > 0) {
    merged <- rbindlist(merged_pairs)
    message("")
    message(glue("{nrow(merged)} set pairs merged across index windows:"))
    print(merged[order(-jaccard)], nrows = 40)
} else {
    message("")
    message("No sets merged at this threshold.")
}

# --- write -----------------------------------------------------------------
write_table <- function(table, name) {
    file <- file.path(fine_mapping_directory, name)
    fwrite(table, file, sep = "\t", na = "NA")
    message(glue("Wrote {name} ({nrow(table)} rows)"))
}

message("")
write_table(sets, "credible_set_summary_deduplicated.txt")
write_table(variants, "credible_set_variants_deduplicated.txt")
write_table(variants[set_key %in% introgressed$set_key],
            "neanderthal_credible_sets.txt")
write_table(variants[set_key %in% majority$set_key],
            "neanderthal_majority_signals.txt")
