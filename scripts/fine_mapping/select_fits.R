#!/usr/bin/env Rscript
# Choose one fit per locus and write the tables every later step reads.
#
#   ./run_rscript.sh select_fits.R
#
# One fit is chosen per locus:
#   - the original fit, if it converged
#   - otherwise a refine = TRUE refit that converged (L = 5 preferred over L = 10)
#   - otherwise the locus is unresolved and is excluded from every later table
#
# A locus whose refits disagree with the original about which variants carry the
# signal is also treated as unresolved, even if a refit converged: convergence to
# an answer that depends on the starting point is not a resolved locus.
#
# Credible sets are flagged genome-wide significant here, since SuSiE will build
# sets out of correlated noise and every later count is restricted to sets
# containing a significant variant.
#
# Signals found through two overlapping index windows are still counted twice in
# these tables; deduplicate_credible_sets.R collapses them.
#
# Writes fit_status.txt, credible_set_summary_selected.txt and
# credible_set_variants_selected.txt.

suppressPackageStartupMessages({
    library(data.table)
    library(glue)
})

project_directory <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/"
fine_mapping_directory <- file.path(project_directory, "data/fine_mapping")
refit_directory <- file.path(fine_mapping_directory, "refits")

GENOME_WIDE_Z <- 5.45
UNSTABLE_AGREEMENT <- 0.3
L_PREFERENCE <- c(5, 10)


# --- gather one fit per locus --------------------------------------------
loci <- fread(file.path(fine_mapping_directory, "locus_summary.txt"),
              na.strings = c("NA", ""))
original_sets <- fread(file.path(fine_mapping_directory, "credible_set_summary.txt"),
                       na.strings = c("NA", ""))
original_variants <- fread(file.path(fine_mapping_directory, "credible_set_variants.txt"),
                           na.strings = c("NA", ""))

jaccard <- function(a, b) {
    if (length(a) == 0 || length(b) == 0) return(0)
    length(intersect(a, b)) / length(union(a, b))
}

significant_variants <- function(sets, variants) {
    if (is.null(sets) || nrow(sets) == 0) return(character(0))
    keep <- sets[max_abs_z >= GENOME_WIDE_Z, CS]
    variants[CS %in% keep, unique(rsID)]
}

chosen_sets <- list()
chosen_variants <- list()
status <- list()

for (index in seq_len(nrow(loci))) {
    locus <- loci[index]
    key <- list(cancer = locus$cancer, ancestry_group = locus$ancestry_group,
                locus_id = locus$locus_id)

    these_sets <- original_sets[cancer == key$cancer &
                                ancestry_group == key$ancestry_group &
                                locus_id == key$locus_id]
    these_variants <- original_variants[cancer == key$cancer &
                                        ancestry_group == key$ancestry_group &
                                        locus_id == key$locus_id]

    if (isTRUE(locus$converged)) {
        chosen_sets[[length(chosen_sets) + 1]] <- these_sets
        chosen_variants[[length(chosen_variants) + 1]] <- these_variants
        status[[length(status) + 1]] <- data.table(
            cancer = key$cancer, ancestry_group = key$ancestry_group,
            locus_id = key$locus_id, fit = "original", state = "converged")
        next
    }

    # non-converged: look for a refit
    refit <- NULL
    refit_label <- NA_character_
    refits_present <- 0L
    for (L in L_PREFERENCE) {
        file <- file.path(refit_directory, key$cancer,
                          glue("{key$locus_id}_{key$ancestry_group}_L{L}.rds"))
        if (!file.exists(file)) next
        refits_present <- refits_present + 1L
        candidate <- readRDS(file)
        if (isTRUE(candidate$converged)) {
            refit <- candidate
            refit_label <- as.character(glue("refit L={L}"))
            break
        }
    }

    if (is.null(refit)) {
        # a refit that has not run yet is not the same as one that failed, and
        # reporting them together understates how much is still in flight
        state <- if (refits_present == 0) {
            "pending: no refit yet"
        } else if (refits_present < length(L_PREFERENCE)) {
            "pending: refits incomplete, none converged so far"
        } else {
            "unresolved: no refit converged"
        }
        status[[length(status) + 1]] <- data.table(
            cancer = key$cancer, ancestry_group = key$ancestry_group,
            locus_id = key$locus_id, fit = NA_character_, state = state)
        next
    }

    agreement <- jaccard(
        significant_variants(these_sets, these_variants),
        significant_variants(refit$cs_summary, refit$cs_table)
    )
    if (agreement < UNSTABLE_AGREEMENT) {
        status[[length(status) + 1]] <- data.table(
            cancer = key$cancer, ancestry_group = key$ancestry_group,
            locus_id = key$locus_id, fit = refit_label,
            state = as.character(glue("unresolved: fits disagree (Jaccard {round(agreement, 2)})")))
        next
    }

    chosen_sets[[length(chosen_sets) + 1]] <- refit$cs_summary
    chosen_variants[[length(chosen_variants) + 1]] <- refit$cs_table
    status[[length(status) + 1]] <- data.table(
        cancer = key$cancer, ancestry_group = key$ancestry_group,
        locus_id = key$locus_id, fit = refit_label,
        state = as.character(glue("rescued (Jaccard {round(agreement, 2)})")))
}

status <- rbindlist(status)
sets <- rbindlist(chosen_sets, fill = TRUE)
variants <- rbindlist(chosen_variants, fill = TRUE)

# --- significance ----------------------------------------------------------
sets[, genome_wide_significant := max_abs_z >= GENOME_WIDE_Z]
sets[, set_key := paste(cancer, ancestry_group, locus_id, CS, sep = "|")]
variants[, set_key := paste(cancer, ancestry_group, locus_id, CS, sep = "|")]

# Carry the significance flag onto the variants, so that anything variant-level
# reads the same fits these counts came from rather than the original tables,
# which still contain the unresolved loci.
variants <- merge(
    variants,
    sets[, .(set_key, genome_wide_significant)],
    by = "set_key", all.x = TRUE
)

# --- report ----------------------------------------------------------------
usable <- sets[genome_wide_significant == TRUE]

message("")
message("=== fits used ===")
print(status[, .N, by = state][order(-N)])

message("")
message("=== credible sets ===")
message(glue("loci with a usable fit: {uniqueN(sets[, .(cancer, ancestry_group, locus_id)])} ",
             "of {nrow(loci)}"))
message(glue("credible sets: {nrow(sets)}; genome-wide significant: {nrow(usable)}"))
message("Signal-level counts are reported by deduplicate_credible_sets.R, which ",
        "collapses sets found through overlapping index windows.")

pending <- status[state %like% "pending"]
if (nrow(pending) > 0) {
    message("")
    message(glue("=== {nrow(pending)} loci still awaiting refits; counts above are provisional ==="))
    print(pending[, .(cancer = substr(cancer, 1, 30), ancestry_group, locus_id, state)])
}

unresolved <- status[state %like% "unresolved"]
if (nrow(unresolved) > 0) {
    message("")
    message(glue("=== {nrow(unresolved)} unresolved loci, excluded from all counts ==="))
    print(unresolved[, .(cancer = substr(cancer, 1, 30), ancestry_group, locus_id, state)])
}

# --- write -----------------------------------------------------------------
fwrite(status, file.path(fine_mapping_directory, "fit_status.txt"),
       sep = "\t", na = "NA")
fwrite(sets, file.path(fine_mapping_directory, "credible_set_summary_selected.txt"),
       sep = "\t", na = "NA")
fwrite(variants, file.path(fine_mapping_directory, "credible_set_variants_selected.txt"),
       sep = "\t", na = "NA")
message("")
message(glue("Wrote fit_status.txt, credible_set_summary_selected.txt and ",
             "credible_set_variants_selected.txt to {fine_mapping_directory}"))
