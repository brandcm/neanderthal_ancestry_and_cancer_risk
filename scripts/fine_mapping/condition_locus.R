# Approximate conditional analysis at one locus.
#
# Recomputes each variant's association after conditioning on one variant, using
# only summary statistics and the LD panel:
#
#     z_conditional = (z_j - r_jc * z_c) / sqrt(1 - r_jc^2)
#
# This is the standard COJO approximation for a single conditioning variant. It
# answers the question a credible set cannot answer on its own: does a second
# signal survive once the first is accounted for? Because it uses nothing but
# marginal z scores and one correlation, it does not depend on susie_rss.
#
# Two assumptions worth stating when reporting it: every variant is taken to have
# the same sample size, and the panel's LD is taken to match the GWAS cohort's.
#
# Usage: run_rscript.sh condition_locus.R <config_key> <locus_id> <ancestry_group> <rsID>

suppressPackageStartupMessages({
    library(data.table)
    library(glue)
})

project_directory <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/"
source(file.path(project_directory, "scripts/fine_mapping/fine_mapping.R"))

fine_mapping_directory <- file.path(project_directory, "data/fine_mapping")

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 4) {
    stop("Usage: condition_locus.R <config_key> <locus_id> <ancestry_group> <rsID>")
}
config_key <- arguments[1]
locus_id <- arguments[2]
ancestry_group <- arguments[3]
conditioning_rsid <- arguments[4]


# signed r between the conditioning variant and every other variant in the window
correlations_with <- function(rsids, rsid, chromosome, ancestry_group) {
    prefix <- tempfile(pattern = "conditional_")
    dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(paste0(prefix, "*")), add = TRUE)

    panel <- file.path(LD_PANEL_DIR, ancestry_group,
                       glue("1000G_{ancestry_group}_chr{chromosome}"))
    writeLines(union(rsids, rsid), glue("{prefix}_variants.txt"))

    run_plink(glue("--bfile {panel} --extract {prefix}_variants.txt --keep-allele-order ",
                   "--make-bed --out {prefix}"))
    # --r gives signed r; --keep-allele-order keeps those signs tied to the .bim A1,
    # which is the allele the betas were aligned to
    run_plink(glue("--bfile {prefix} --keep-allele-order --r --ld-snp {rsid} ",
                   "--ld-window-kb 10000 --ld-window 999999 --ld-window-r2 0 ",
                   "--out {prefix}_r"))

    correlations <- fread(glue("{prefix}_r.ld"))
    correlations[, .(rsID = SNP_B, r = R)]
}


main <- function() {
    result_file <- file.path(fine_mapping_directory, "results", config_key,
                             glue("{locus_id}_{ancestry_group}.rds"))
    if (!file.exists(result_file)) stop(glue("No result at {result_file}"))
    result <- readRDS(result_file)

    variants <- as.data.table(result$variants)
    chromosome <- result$locus_summary$chr

    if (!(conditioning_rsid %in% variants$rsID)) {
        stop(glue("{conditioning_rsid} is not among the locus's fitted variants."))
    }

    conditioning_z <- variants[rsID == conditioning_rsid, Z]
    message(glue("Conditioning on {conditioning_rsid} (z = {round(conditioning_z, 2)}) ",
                 "at {config_key} | {locus_id} | {ancestry_group}"))

    correlations <- correlations_with(variants$rsID, conditioning_rsid,
                                      chromosome, ancestry_group)
    variants <- correlations[variants, on = "rsID"]

    # a variant in perfect LD with the conditioning variant has no information left
    variants[, r := fifelse(rsID == conditioning_rsid, 1, r)]
    variants[, z_conditional := fifelse(
        is.na(r) | abs(r) >= 0.999, NA_real_,
        (Z - r * conditioning_z) / sqrt(1 - r^2)
    )]

    variants <- variants[order(-abs(z_conditional), na.last = TRUE)]

    message("\n--- strongest remaining signals after conditioning ---")
    print(head(variants[, .(rsID, pos, Z, r = round(r, 3),
                            z_conditional = round(z_conditional, 2),
                            is_neanderthal)], 15))

    introgressed <- variants[is_neanderthal == TRUE & !is.na(z_conditional)]
    if (nrow(introgressed) > 0) {
        message("\n--- strongest introgressed variants after conditioning ---")
        print(head(introgressed[, .(rsID, pos, Z, r = round(r, 3),
                                    z_conditional = round(z_conditional, 2),
                                    introgression_maps)], 10))
        message(glue(
            "\nmax |z| among introgressed variants: {round(max(abs(introgressed$Z)), 2)} ",
            "before conditioning, {round(max(abs(introgressed$z_conditional)), 2)} after"
        ))
    }

    out_directory <- file.path(fine_mapping_directory, "conditional")
    dir.create(out_directory, recursive = TRUE, showWarnings = FALSE)
    out_file <- file.path(out_directory,
                          glue("{config_key}_{locus_id}_{ancestry_group}_",
                               "conditioned_on_{conditioning_rsid}.txt"))
    fwrite(variants[, .(rsID, chr, pos, beta, se, Z, r, z_conditional,
                        is_neanderthal, introgression_maps, neanderthal_allele)],
           out_file, sep = "\t", na = "NA")
    message(glue("\nWrote {out_file}"))
}

main()
