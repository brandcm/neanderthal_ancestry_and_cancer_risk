#!/usr/bin/env Rscript
# Flag credible set variants that are reintroduced ancestral alleles (RAAs).
#
#   ./run_rscript.sh annotate_reintroduced_ancestral.R
#
# intersect_reintroduced_ancestral.R answers the same question for one locus at
# a time from its .rds result, for interactive inspection. This runs across the
# deduplicated credible set variants instead and writes two tables, so the calls
# can be joined onto the supplementary tables.
#
# The calls are ancestry-specific and in hg19, matched on chromosome and
# position, so this is only valid if the credible set coordinates are also hg19.
# The per-ancestry match rate is printed as a build check: a set of windows with
# thousands of variants and almost no matches anywhere is the signature of a
# build mismatch rather than a biological result.

suppressPackageStartupMessages({
    library(data.table)
    library(glue)
})

project <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/"
fm <- file.path(project, "data/fine_mapping")
reintroduced_directory <- file.path(project, "data/reintroduced_ancestral_alleles")

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

load_reintroduced <- function(ancestry_group) {
    path <- file.path(
        reintroduced_directory,
        glue("{ancestry_group}_reintroduced_ancestral_alleles.txt.gz")
    )
    if (!file.exists(path)) {
        stop(glue("No reintroduced ancestral calls for {ancestry_group}: {path}"))
    }
    calls <- fread(path, header = FALSE, col.names = c("chr", "pos"))
    calls[, chr := sub("^chr", "", chr)]
    calls[, pos := as.integer(pos)]
    unique(calls)
}

# the calls differ by ancestry group, so each group is matched against its own
variants[, is_reintroduced_ancestral := FALSE]
for (group in sort(unique(variants$ancestry_group))) {
    calls <- load_reintroduced(group)
    matched <- variants[ancestry_group == group][
        calls, on = .(chr, pos), nomatch = NULL
    ]
    if (nrow(matched) > 0) {
        variants[unique(matched[, .(cancer, ancestry_group, locus_id, CS, rsID)]),
                 is_reintroduced_ancestral := TRUE,
                 on = .(cancer, ancestry_group, locus_id, CS, rsID)]
    }
    rate <- variants[ancestry_group == group, mean(is_reintroduced_ancestral)]
    message(glue("{group}: {round(100 * rate, 2)}% of credible set variants are ",
                 "reintroduced ancestral alleles ",
                 "({variants[ancestry_group == group, sum(is_reintroduced_ancestral)]} ",
                 "of {variants[ancestry_group == group, .N]})"))
}

variant_file <- file.path(fm, "reintroduced_ancestral_variants.txt")
fwrite(
    variants[, .(cancer, ancestry_group, locus_id, CS, rsID, chr, pos,
                 is_neanderthal, is_reintroduced_ancestral)],
    variant_file, sep = "\t", na = "NA"
)
message(glue("Wrote {variant_file}"))

per_set <- variants[, .(
    n_variants = .N,
    n_reintroduced_ancestral = sum(is_reintroduced_ancestral),
    n_neanderthal_reintroduced = sum(
        is_reintroduced_ancestral & is_neanderthal == TRUE, na.rm = TRUE
    )
), by = .(cancer, ancestry_group, locus_id, CS)]

per_set <- merge(
    per_set,
    sets[, .(cancer, ancestry_group, locus_id, CS, signal_id, is_representative,
             p_neanderthal, max_abs_z, purity, genome_wide_significant)],
    by = c("cancer", "ancestry_group", "locus_id", "CS"), all.x = TRUE
)

significant <- per_set[genome_wide_significant == TRUE & is_representative == TRUE]
archaic <- significant[p_neanderthal > 0.5]
count_loci <- function(dt) uniqueN(dt[, .(cancer, ancestry_group, locus_id)])

message("")
message("--- genome-wide significant signals ---")
message(glue("signals: {nrow(significant)} across {count_loci(significant)} loci"))
message(glue("signals with >=1 reintroduced ancestral allele: ",
             "{nrow(significant[n_reintroduced_ancestral > 0])} across ",
             "{count_loci(significant[n_reintroduced_ancestral > 0])} loci"))

message("")
message("--- of those, Neanderthal-majority (p_neanderthal > 0.5) ---")
message(glue("signals: {nrow(archaic)} across {count_loci(archaic)} loci"))
message(glue("signals with >=1 reintroduced ancestral allele: ",
             "{nrow(archaic[n_reintroduced_ancestral > 0])} across ",
             "{count_loci(archaic[n_reintroduced_ancestral > 0])} loci"))

out_file <- file.path(fm, "reintroduced_ancestral_overlap.txt")
fwrite(per_set, out_file, sep = "\t", na = "NA")
message("")
message(glue("Wrote {out_file}"))
