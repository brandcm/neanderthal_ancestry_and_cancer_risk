# Collect the per-locus fine-mapping results into flat tables.
#
# Each locus's result object carries its SuSiE fit and its kriging diagnostic, so
# reading all of them at once is heavy enough to lose a session. This reads one at
# a time and keeps only the tables.
#
# Outputs:
#   locus_summary.txt            one row per locus, including loci with no credible sets
#   credible_set_summary.txt     one row per credible set, with genome_wide_significant
#   credible_set_variants.txt    one row per variant in a credible set
#   kriging_flagged_variants.txt   one row per variant the kriging diagnostic flagged
#   neanderthal_credible_sets.txt  the credible set variants at loci where any
#                                  credible set contains an introgressed variant

suppressPackageStartupMessages({
    library(data.table)
    library(glue)
})

# |z| at p = 5e-8. A credible set with nothing above this contains no
# genome-wide significant variant: SuSiE fits L single-effect components whether
# or not L signals exist, and susie_get_cs filters on purity -- a property of the
# LD between set members, not of their evidence -- so sets form readily around
# correlated noise. Such sets inflate every count downstream and, through
# near_window_edge, can trigger expensive window widening for no reason.
GENOME_WIDE_Z <- 5.45

NEANDERTHAL_COLUMNS <- c(
    "cancer", "ancestry_group", "locus", "credible_set", "rsid", "chr", "pos",
    "effect_allele", "other_allele", "PIP", "beta", "z", "purity", "p_neanderthal",
    "genome_wide_significant",
    "is_neanderthal", "introgression_maps", "neanderthal_allele", "neanderthal_beta", "notes"
)


parse_args <- function(arguments = commandArgs(trailingOnly = TRUE)) {
    if (length(arguments) %% 2 != 0) stop("Arguments must be --name value pairs.")
    names_and_values <- split(arguments, rep(c("name", "value"), length(arguments) / 2))
    values <- as.list(names_and_values$value)
    names(values) <- sub("^--", "", names_and_values$name)
    values
}


read_results <- function(files) {
    locus_summaries <- vector("list", length(files))
    credible_set_summaries <- vector("list", length(files))
    credible_set_variants <- vector("list", length(files))
    kriging_flagged <- vector("list", length(files))

    for (index in seq_along(files)) {
        result <- readRDS(files[index])
        locus_summaries[[index]] <- result$locus_summary
        credible_set_summaries[[index]] <- result$cs_summary
        credible_set_variants[[index]] <- result$cs_table

        # the flag table carries no locus columns of its own
        if (!is.null(result$kriging_flagged) && nrow(result$kriging_flagged) > 0) {
            ids <- result$locus_summary[, .(cancer, ancestry_group, locus_id)]
            kriging_flagged[[index]] <- cbind(ids, result$kriging_flagged)
        }

        # the fit and the kriging plot are the bulk of each object
        rm(result)
        if (index %% 25 == 0) {
            gc()
            message(glue("  {index} of {length(files)}"))
        }
    }

    list(
        locus_summary = rbindlist(locus_summaries, fill = TRUE),
        credible_set_summary = rbindlist(credible_set_summaries, fill = TRUE),
        credible_set_variants = rbindlist(credible_set_variants, fill = TRUE),
        kriging_flagged_variants = rbindlist(kriging_flagged, fill = TRUE)
    )
}


mark_genome_wide_significance <- function(tables) {
    if (nrow(tables$credible_set_summary) == 0) return(tables)

    tables$credible_set_summary[, genome_wide_significant := max_abs_z >= GENOME_WIDE_Z]

    if (nrow(tables$credible_set_variants) > 0) {
        flags <- tables$credible_set_summary[
            , .(cancer, ancestry_group, locus_id, CS, genome_wide_significant)
        ]
        tables$credible_set_variants <- flags[
            tables$credible_set_variants, on = .(cancer, ancestry_group, locus_id, CS)
        ]
    }
    tables
}


# a flagged variant inside a credible set is the one that matters: it means the
# summary statistics and the LD panel disagree about a variant the model is relying on
mark_credible_set_membership <- function(flagged, credible_set_variants) {
    if (nrow(flagged) == 0) return(flagged)

    if (nrow(credible_set_variants) == 0) {
        flagged[, `:=`(in_credible_set = FALSE, CS = NA_character_)]
        return(flagged)
    }

    members <- unique(credible_set_variants[, .(cancer, ancestry_group, locus_id, rsID, CS)])
    flagged <- members[flagged, on = .(cancer, ancestry_group, locus_id, rsID)]
    flagged[, in_credible_set := !is.na(CS)]
    setorder(flagged, -in_credible_set, -log_lr)
    flagged
}


# every variant in every credible set at a locus where any credible set contains an
# introgressed variant, so that a Neanderthal signal is always reported alongside
# whatever else that locus found
neanderthal_credible_sets <- function(credible_set_variants) {
    if (nrow(credible_set_variants) == 0) return(data.table())

    neanderthal_loci <- unique(
        credible_set_variants[is_neanderthal == TRUE, .(cancer, ancestry_group, locus_id)]
    )
    variants <- credible_set_variants[
        neanderthal_loci, on = .(cancer, ancestry_group, locus_id), nomatch = NULL
    ]

    variants <- variants[, .(
        cancer, ancestry_group, locus = locus_id, credible_set = CS, rsid = rsID, chr, pos,
        effect_allele, other_allele, PIP, beta, z = Z, purity, p_neanderthal,
        genome_wide_significant,
        is_neanderthal, introgression_maps, neanderthal_allele, neanderthal_beta, notes
    )]
    setcolorder(variants, intersect(NEANDERTHAL_COLUMNS, names(variants)))
    setorder(variants, cancer, ancestry_group, locus, credible_set, -PIP)
    variants
}


main <- function() {
    args <- parse_args()

    files <- list.files(args$results_directory, pattern = "\\.rds$",
                        recursive = TRUE, full.names = TRUE)
    if (length(files) == 0) stop(glue("No .rds results found under {args$results_directory}."))
    message(glue("Reading {length(files)} loci..."))

    tables <- read_results(files)
    tables <- mark_genome_wide_significance(tables)
    tables$kriging_flagged_variants <- mark_credible_set_membership(
        tables$kriging_flagged_variants, tables$credible_set_variants
    )
    tables$neanderthal_credible_sets <- neanderthal_credible_sets(tables$credible_set_variants)

    dir.create(args$output_directory, recursive = TRUE, showWarnings = FALSE)
    for (name in names(tables)) {
        out_file <- file.path(args$output_directory, glue("{name}.txt"))
        fwrite(tables[[name]], out_file, sep = "\t", na = "NA")
        message(glue("Wrote {nrow(tables[[name]])} rows to {out_file}"))
    }

    summary <- tables$locus_summary
    message(glue(
        "\n{nrow(summary)} loci; {sum(!summary$converged)} did not converge; ",
        "{sum(summary$n_kriging_severe > 0)} with severe kriging flags; ",
        "{sum(summary$n_credible_sets == 0)} with no credible sets; ",
        "{sum(summary$n_neanderthal_cs > 0)} with an introgressed credible set"
    ))

    credible_sets <- tables$credible_set_summary
    if (nrow(credible_sets) > 0) {
        significant <- sum(credible_sets$genome_wide_significant)
        message(glue(
            "{nrow(credible_sets)} credible sets; {significant} contain a ",
            "genome-wide significant variant (|z| >= {GENOME_WIDE_Z}); ",
            "{nrow(credible_sets) - significant} do not and should be excluded ",
            "before interpreting p_neanderthal"
        ))
        neanderthal_majority <- credible_sets[p_neanderthal > 0.5]
        message(glue(
            "{nrow(neanderthal_majority)} sets with p_neanderthal > 0.5, of which ",
            "{sum(neanderthal_majority$genome_wide_significant)} are genome-wide significant"
        ))
        edge_on_noise <- credible_sets[near_window_edge == TRUE & !genome_wide_significant]
        if (nrow(edge_on_noise) > 0) {
            message(glue(
                "{nrow(edge_on_noise)} near-edge flags come from sets with no ",
                "genome-wide significant variant, across ",
                "{uniqueN(edge_on_noise[, .(cancer, ancestry_group, locus_id)])} loci; ",
                "those loci did not need widening"
            ))
        }
    }

    flagged <- tables$kriging_flagged_variants
    if (nrow(flagged) > 0) {
        message(glue(
            "{nrow(flagged)} kriging-flagged variants, ",
            "{sum(flagged$in_credible_set)} of them in a credible set"
        ))
    }
}

main()
