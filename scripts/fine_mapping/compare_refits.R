#!/usr/bin/env Rscript
# Compare the original fit of each non-converged locus against its refine = TRUE
# refits at L = 5 and L = 10.
#
#   ./run_rscript.sh compare_refits.R
#
# The question is not "did it converge" but "do independent fits agree". Where
# the two refits select the same variants, the locus has settled on a stable
# answer and can be used. Where they disagree, the locus is genuinely ambiguous
# and should be reported as such rather than refit until it behaves.
#
# Agreement is measured as the Jaccard index between the variants in
# genome-wide significant credible sets, which is what downstream analysis
# actually consumes. Missing refits are reported as NA, so this can be run while
# the array job is still going.

suppressPackageStartupMessages({
    library(data.table)
    library(glue)
})

project_directory <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/"
fine_mapping_directory <- file.path(project_directory, "data/fine_mapping")
refit_directory <- file.path(fine_mapping_directory, "refits")

GENOME_WIDE_Z <- 5.45
L_VALUES <- c(5, 10)


describe_fit <- function(result) {
    if (is.null(result)) {
        # variants = NULL, not character(0): a missing fit must be
        # distinguishable from a fit that found nothing, or agreement with it
        # reads as 0 rather than unknown.
        return(list(converged = NA, sets = NA_integer_, significant = NA_integer_,
                    archaic = NA_integer_, single_sig = NA_integer_,
                    variants = NULL))
    }
    sets <- result$cs_summary
    if (is.null(sets) || nrow(sets) == 0) {
        return(list(converged = isTRUE(result$converged), sets = 0L,
                    significant = 0L, archaic = 0L, single_sig = 0L,
                    variants = character(0)))
    }
    sets <- copy(sets)
    sets[, significant := max_abs_z >= GENOME_WIDE_Z]
    keep <- sets[significant == TRUE, CS]
    variants <- result$cs_table[CS %in% keep, unique(rsID)]

    list(
        converged   = isTRUE(result$converged),
        sets        = nrow(sets),
        significant = sum(sets$significant),
        archaic     = nrow(sets[significant == TRUE & p_neanderthal > 0.5]),
        single_sig  = nrow(sets[significant == TRUE & n_variants == 1]),
        variants    = variants
    )
}


jaccard <- function(a, b) {
    if (is.null(a) || is.null(b)) return(NA_real_)       # a fit that has not run
    if (length(a) == 0 && length(b) == 0) return(NA_real_)
    union_size <- length(union(a, b))
    if (union_size == 0) return(NA_real_)
    round(length(intersect(a, b)) / union_size, 3)
}


summary_table <- fread(file.path(fine_mapping_directory, "locus_summary.txt"),
                       na.strings = c("NA", ""))
unconverged <- summary_table[converged == FALSE,
                             .(config_key = cancer, ancestry_group, locus_id, n_variants)]
setorder(unconverged, n_variants)

rows <- vector("list", nrow(unconverged))

for (index in seq_len(nrow(unconverged))) {
    locus <- unconverged[index]

    original_file <- file.path(fine_mapping_directory, "results", locus$config_key,
                               glue("{locus$locus_id}_{locus$ancestry_group}.rds"))
    original <- if (file.exists(original_file)) readRDS(original_file) else NULL
    fits <- list(original = describe_fit(original))

    for (L in L_VALUES) {
        refit_file <- file.path(refit_directory, locus$config_key,
                                glue("{locus$locus_id}_{locus$ancestry_group}_L{L}.rds"))
        refit <- if (file.exists(refit_file)) readRDS(refit_file) else NULL
        fits[[glue("L{L}")]] <- describe_fit(refit)
    }

    rows[[index]] <- data.table(
        cancer = locus$config_key,
        anc = locus$ancestry_group,
        locus_id = locus$locus_id,
        nvar = locus$n_variants,

        orig_sets = fits$original$sets,
        orig_sig  = fits$original$significant,
        orig_arch = fits$original$archaic,

        L5_conv = fits$L5$converged,
        L5_sets = fits$L5$sets,
        L5_sig  = fits$L5$significant,
        L5_arch = fits$L5$archaic,

        L10_conv = fits$L10$converged,
        L10_sets = fits$L10$sets,
        L10_sig  = fits$L10$significant,
        L10_arch = fits$L10$archaic,

        agree_L5_L10 = jaccard(fits$L5$variants, fits$L10$variants),
        agree_orig_L5 = jaccard(fits$original$variants, fits$L5$variants)
    )
}

comparison <- rbindlist(rows)
print(comparison, nrows = 40)

done <- comparison[!is.na(L5_sets) & !is.na(L10_sets)]
message("")
message(glue("{nrow(done)} of {nrow(comparison)} loci have both refits."))
if (nrow(done) > 0) {
    message(glue("both refits converged: {sum(done$L5_conv & done$L10_conv, na.rm = TRUE)}"))
    message(glue("refits agree closely (Jaccard >= 0.8): ",
                 "{sum(done$agree_L5_L10 >= 0.8, na.rm = TRUE)}"))
    message(glue("credible sets at L = 5 vs original: ",
                 "{sum(done$L5_sets)} vs {sum(done$orig_sets)}"))
    message(glue("introgressed-majority sets at L = 5 vs original: ",
                 "{sum(done$L5_arch)} vs {sum(done$orig_arch)}"))
}

out_file <- file.path(fine_mapping_directory, "refit_comparison.txt")
fwrite(comparison, out_file, sep = "\t", na = "NA")
message(glue("\nWrote {out_file}"))
