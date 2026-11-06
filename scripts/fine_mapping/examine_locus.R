# Look at one fine-mapped locus in detail.
#
# Prints the credible sets, the variants in them, and the LD between the strongest
# introgressed variant and the strongest variant overall in each set. That last part
# is the one that answers "is the introgressed signal real, or are those variants
# just in LD with a modern causal variant?" -- a credible set is a statement that the
# causal variant is one of its members, and the LD says how distinguishable they are.
#
# Usage: run_rscript.sh examine_locus.R <config_key> <locus_id> <ancestry_group>
 
suppressPackageStartupMessages({
    library(data.table)
    library(glue)
})
 
project_directory <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/"
source(file.path(project_directory, "scripts/fine_mapping/fine_mapping.R"))
 
fine_mapping_directory <- file.path(project_directory, "data/fine_mapping")
 
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3) stop("Usage: examine_locus.R <config_key> <locus_id> <ancestry_group>")
config_key <- arguments[1]
locus_id <- arguments[2]
ancestry_group <- arguments[3]
 
 
# r between the credible set's variants, from the same panel the fit used
credible_set_ld <- function(rsids, chromosome, ancestry_group) {
    prefix <- tempfile(pattern = "credible_set_ld_")
    dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(paste0(prefix, "*")), add = TRUE)
 
    panel <- credible_set_ld_panel(chromosome, ancestry_group)
    writeLines(rsids, glue("{prefix}_variants.txt"))
    run_plink(glue("--bfile {panel} --extract {prefix}_variants.txt --keep-allele-order ",
                   "--make-bed --out {prefix}"))
    run_plink(glue("--bfile {prefix} --keep-allele-order --r square --out {prefix}_ld"))
 
    bim <- fread(glue("{prefix}.bim"), header = FALSE,
                 col.names = c("chr", "rsID", "cm", "pos", "A1", "A2"),
                 colClasses = list(character = c(1, 2, 5, 6)))
    ld <- as.matrix(fread(glue("{prefix}_ld.ld"), header = FALSE))
    dimnames(ld) <- list(bim$rsID, bim$rsID)
    ld
}
 
 
credible_set_ld_panel <- function(chromosome, ancestry_group) {
    file.path(LD_PANEL_DIR, ancestry_group,
              glue("1000G_{ancestry_group}_chr{chromosome}"))
}
 
 
main <- function() {
    result_file <- file.path(fine_mapping_directory, "results", config_key,
                             glue("{locus_id}_{ancestry_group}.rds"))
    if (!file.exists(result_file)) stop(glue("No result at {result_file}"))
    result <- readRDS(result_file)
 
    message(glue("\n=== {config_key} | {locus_id} | {ancestry_group} ===\n"))
    print(result$locus_summary)
 
    if (is.null(result$cs_summary)) {
        message("No credible sets.")
        return(invisible(NULL))
    }
 
    message("\n--- credible sets ---")
    print(result$cs_summary[, .(CS, n_variants, n_neanderthal, p_neanderthal,
                                max_pip, max_abs_z, purity, near_window_edge)])
 
    variants <- as.data.table(result$cs_table)
 
    for (set in unique(variants$CS)) {
        members <- variants[CS == set][order(-PIP)]
        message(glue("\n--- {set}: {nrow(members)} variants, ",
                     "p_neanderthal {signif(members$p_neanderthal[1], 3)} ---"))
        print(members[, .(rsID, pos, effect_allele, other_allele, PIP, beta, Z,
                          is_neanderthal, introgression_maps, neanderthal_allele,
                          neanderthal_beta)])
 
        if (!any(members$is_neanderthal)) next
 
        # how distinguishable is the best introgressed variant from the best overall?
        ld <- tryCatch(
            credible_set_ld(members$rsID, result$locus_summary$chr, ancestry_group),
            error = function(e) {
                message(glue("LD unavailable: {conditionMessage(e)}"))
                NULL
            }
        )
        if (is.null(ld)) next
 
        top_overall <- members$rsID[1]
        top_neanderthal <- members[is_neanderthal == TRUE][order(-PIP)]$rsID[1]
        shared <- intersect(c(top_overall, top_neanderthal), rownames(ld))
 
        if (length(shared) == 2) {
            correlation <- ld[top_overall, top_neanderthal]
            message(glue(
                "\ntop variant {top_overall} (PIP {signif(members$PIP[1], 3)}) vs ",
                "top introgressed {top_neanderthal} ",
                "(PIP {signif(members[rsID == top_neanderthal]$PIP[1], 3)}): ",
                "r = {signif(correlation, 3)}, r2 = {signif(correlation^2, 3)}"
            ))
        }
 
        message("\nr between the set's variants:")
        print(round(ld[members$rsID[members$rsID %in% rownames(ld)],
                       members$rsID[members$rsID %in% rownames(ld)]], 2))
    }
 
    flagged <- as.data.table(result$kriging_flagged)
    message(glue("\n--- kriging: {nrow(flagged)} flagged, ",
                 "{result$locus_summary$n_kriging_severe} severe ---"))
    if (nrow(flagged) > 0) print(flagged)
}
 
main()
