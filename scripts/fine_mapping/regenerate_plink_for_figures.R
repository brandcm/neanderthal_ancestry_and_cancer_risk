#!/usr/bin/env Rscript
# Regenerate the PLINK LD panel and .bim for the loci used in figures.
#
#   module load CBI plink/1.90b6.26
#   ./run_rscript.sh regenerate_plink_for_figures.R
#
# fine_mapping() writes its PLINK intermediates to a temp directory and deletes
# them on exit, so the LD matrix behind each fit no longer exists. Regenerating
# through fine_mapping() rather than calling PLINK directly keeps the variant set
# and allele alignment identical to the fit the credible sets came from.
#
# Files land in data/fine_mapping/plink/<config_key>/:
#   <config_key>_<locus_id>_<ancestry>_1KG_filtered.bim
#   <config_key>_<locus_id>_<ancestry>_1KG_filtered_ld.ld
#
# The fit itself is recomputed and discarded; nothing in results/ is touched.

suppressPackageStartupMessages({
    library(data.table)
    library(glue)
    library(susieR)
    library(yaml)
})

project_directory <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/"
source(file.path(project_directory, "scripts/fine_mapping/fine_mapping.R"))

loci <- fread(file.path(project_directory, "data/fine_mapping/fine_mapping_loci.txt"),
              na.strings = c("NA", ""))
overrides_file <- file.path(project_directory, "scripts/fine_mapping/locus_overrides.txt")
overrides <- if (file.exists(overrides_file)) {
    fread(overrides_file, na.strings = c("NA", ""))
} else NULL

introgressed <- load_introgressed_variants()


regenerate <- function(gwas, ancestry, locus) {
    row <- loci[config_key == gwas & ancestry_group == ancestry & locus_id == locus]
    if (nrow(row) != 1) stop(glue("Expected one row for {gwas} {ancestry} {locus}, got {nrow(row)}."))

    # use the window the fit actually used, which an override may have changed
    window <- row$window_size
    if (!is.null(overrides)) {
        override <- overrides[config_key == gwas & locus_id == locus]
        if (nrow(override) == 1 && !is.na(override$window_size)) window <- override$window_size
    }

    message(glue("\n=== {gwas} | {ancestry} | {locus} | window {window} kb ==="))
    fine_mapping(
        introgressed,
        config_key       = gwas,
        ancestry_group   = ancestry,
        locus_id         = locus,
        n_samples        = row$n_samples,
        lead_rsid        = row$lead_rsid,
        window_size      = as.integer(window),
        keep_plink_files = TRUE,
        save_results     = FALSE,
        main_file        = NULL
    )
    invisible(NULL)
}


# BNC2
regenerate("basal_cell_carcinoma_EUR_Seviiri_GCST90137411", "EUR", "9_16800341")
regenerate("squamous_cell_carcinoma_EUR_Seviiri_GCST90137412", "EUR", "9_16800341")

# KRT5
regenerate("basal_cell_carcinoma_EUR_Seviiri_GCST90137411", "EUR", "12_52908254")

# PABPC4
regenerate("squamous_cell_carcinoma_EUR_Seviiri_GCST90137412", "EUR", "1_40035928")

# SLC24A4
regenerate("squamous_cell_carcinoma_EUR_Seviiri_GCST90137412", "EUR", "14_92773295")

# PTHLH
regenerate("breast_cancer_EAS_Jia_GCST90429846", "EAS", "12_28143787")
regenerate("breast_cancer_EUR_Zhang_GCST001", "EUR", "12_28153018")

message(glue("\nPLINK files are under {file.path(project_directory, 'data/fine_mapping/plink')}"))
