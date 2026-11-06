#!/usr/bin/env Rscript
# Intersect a locus's credible set variants with the reintroduced ancestral
# allele calls, which are ancestry-specific and in hg19.
#
#   ./run_rscript.sh intersect_reintroduced_ancestral.R <result.rds> [<result.rds> ...]
#
# The ancestry group is read from each result, so EAS and EUR results can be
# passed together and each is matched against its own call set.
#
# Coordinates are matched on chromosome and position, so this is only valid if
# the fine-mapping results are also hg19. The script prints the window
# coordinates and the overall match rate across the whole window so that a
# build mismatch is obvious: a window with hundreds of introgressed variants
# and zero reintroduced-ancestral matches anywhere in it is the signature.

suppressPackageStartupMessages({
    library(data.table)
    library(glue)
})

project_directory <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/"
reintroduced_directory <- file.path(project_directory, "data/reintroduced_ancestral_alleles")

result_files <- commandArgs(trailingOnly = TRUE)
if (length(result_files) == 0) stop("Pass one or more result .rds files.")

load_reintroduced <- function(ancestry_group) {
    path <- file.path(reintroduced_directory,
                      glue("{ancestry_group}_reintroduced_ancestral_alleles.txt.gz"))
    if (!file.exists(path)) stop(glue("No reintroduced ancestral calls for {ancestry_group}: {path}"))
    calls <- fread(path, header = FALSE, col.names = c("chr", "pos"))
    calls[, chr := sub("^chr", "", chr)]
    calls[, pos := as.integer(pos)]
    unique(calls)
}

reintroduced_cache <- list()

for (result_file in result_files) {
    result <- readRDS(result_file)
    summary_row <- result$locus_summary
    ancestry_group <- summary_row$ancestry_group
    chromosome <- as.character(summary_row$chr)

    if (is.null(reintroduced_cache[[ancestry_group]])) {
        reintroduced_cache[[ancestry_group]] <- load_reintroduced(ancestry_group)
    }
    calls <- reintroduced_cache[[ancestry_group]][chr == sub("^chr", "", chromosome)]

    message("")
    message(strrep("-", 72))
    message(glue("{summary_row$cancer} | {ancestry_group} | {summary_row$locus_id}"))
    message(glue("Window: chr{chromosome}:{summary_row$window_start}-{summary_row$window_end} ",
                 "({summary_row$window_size} kb); {nrow(calls)} reintroduced ancestral calls on this chromosome"))

    # Build sanity check across the entire window, not just credible sets.
    variants <- copy(result$variants)
    variants[, pos := as.integer(pos)]
    in_window <- calls[pos >= summary_row$window_start & pos <= summary_row$window_end]
    window_hits <- variants[pos %in% in_window$pos]
    message(glue("Window check: {nrow(in_window)} reintroduced calls fall inside these coordinates; ",
                 "{nrow(window_hits)} of {nrow(variants)} window variants match."))
    if (nrow(in_window) == 0) {
        message("  No calls anywhere in this window. Either there are genuinely none here, ",
                "or the results are not hg19 - check the window coordinates against the ",
                "expected hg19 position of the gene before reading anything into this.")
    }

    cs_table <- result$cs_table
    if (is.null(cs_table) || nrow(cs_table) == 0) {
        message("No credible sets at this locus.")
        next
    }
    cs_table <- copy(cs_table)
    cs_table[, pos := as.integer(pos)]
    cs_table[, is_reintroduced_ancestral := pos %in% calls$pos]

    per_set <- cs_table[, .(
        n_variants = .N,
        n_neanderthal = sum(is_neanderthal, na.rm = TRUE),
        n_reintroduced = sum(is_reintroduced_ancestral),
        max_abs_z = round(max(abs(Z)), 2),
        purity = round(unique(purity), 3)
    ), by = CS]
    per_set <- merge(per_set,
                     result$cs_summary[, .(CS, p_neanderthal = round(p_neanderthal, 3))],
                     by = "CS", all.x = TRUE)
    print(per_set[order(CS)])

    hits <- cs_table[is_reintroduced_ancestral == TRUE,
                     .(CS, rsID, pos, effect_allele, other_allele,
                       PIP = round(PIP, 4), Z = round(Z, 2),
                       is_neanderthal, neanderthal_allele,
                       introgression_maps)]
    if (nrow(hits) > 0) {
        message("")
        message("Credible set members that are reintroduced ancestral alleles:")
        print(hits[order(CS, -PIP)])
    } else {
        message("No credible set member is a reintroduced ancestral allele.")
    }
}

message("")
