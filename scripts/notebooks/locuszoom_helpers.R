# Convenience layer over generate_locuszoom_plot_with_annotations().
#
#   source("scripts/fine_mapping/locuszoom_helpers.R")
#
#   p <- locuszoom_for_credible_set(
#       cancer         = "squamous_cell_carcinoma_EUR_Seviiri_GCST90137412",
#       ancestry_group = "EUR",
#       locus_id       = "9_16800341",
#       credible_set   = "L2"
#   )
#
# Everything else is derived: the summary statistics path and its column names
# from the GWAS metadata, the LD matrix and .bim from the PLINK files
# fine_mapping() wrote, the credible set members and index variant from the
# credible set tables, and the cCRE file from the chromosome.
#
# Credible set members are NOT circled by default. Fill already carries r², so a
# second mark on the same points competes with it, and the PIP panel drawn below
# a locuszoom shares its x axis and already distinguishes the sets by colour: a
# reader reads membership off that panel by position. Pass
# mark_credible_sets = TRUE for a locuszoom shown on its own, where there is no
# PIP panel to carry it.
#
# Every resolved value is printed, and anything that cannot be resolved raises an
# error naming the missing file rather than failing later inside awk or a join.
# Any argument can still be passed explicitly to override what is resolved.

suppressPackageStartupMessages({
    library(data.table)
    library(yaml)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

PROJECT_DIRECTORY <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk"


# --- find the GWAS column names ---------------------------------------------
# The metadata keys have varied, so each field is looked up under several
# plausible names and what was found is reported.
gwas_columns <- function(config_key,
                         metadata_file = file.path(PROJECT_DIRECTORY,
                                                   "data/metadata/cancer_GWAS_metadata.yaml")) {
    metadata <- yaml::read_yaml(metadata_file)
    config <- metadata[[config_key]]
    if (is.null(config)) {
        stop(sprintf("No entry for '%s' in %s. Available: %s",
                     config_key, basename(metadata_file),
                     paste(head(names(metadata), 5), collapse = ", ")))
    }

    first_present <- function(candidates, label) {
        for (candidate in candidates) {
            if (!is.null(config[[candidate]])) return(config[[candidate]])
        }
        stop(sprintf("Could not find the %s column for %s. Keys present: %s",
                     label, config_key, paste(names(config), collapse = ", ")))
    }

    list(
        input    = first_present(c("input", "path", "sumstats"), "summary statistics"),
        chr_col  = first_present(c("chr_column", "chromosome_column", "chr"), "chromosome"),
        pos_col  = first_present(c("pos_column", "position_column", "bp_column"), "position"),
        pval_col = first_present(c("p_value_column", "pval_column", "p_column"), "p-value"),
        rsid_col = first_present(c("rsID_column", "rsid_column", "snp_column",
                                   "variant_column"), "rsID"),
        build    = config$assembly %||% config$positions_build
    )
}


# --- gather every argument, without plotting ---------------------------------
resolve_locus_inputs <- function(cancer, ancestry_group, locus_id, credible_set,
                                 index_snv = NULL,
                                 mark_credible_sets = FALSE,
                                 # override what the metadata would supply. All five
                                 # together describe a file that is not a registered
                                 # GWAS, such as the output of make_conditional_gwas.R,
                                 # and skip the metadata lookup entirely.
                                 gwas_path = NULL,
                                 chr_col = NULL,
                                 pos_col = NULL,
                                 pval_col = NULL,
                                 rsid_col = NULL,
                                 variants_file = NULL,
                                 ccres_path = NULL,
                                 adaptive_introgression_file = file.path(
                                     PROJECT_DIRECTORY,
                                     "data/adaptive_introgression",
                                     "Iasi_et_al_2025_high_frequency_ancient_modern_regions.bed")) {

    fine_mapping_directory <- file.path(PROJECT_DIRECTORY, "data/fine_mapping")

    if (is.null(variants_file)) {
        current <- file.path(fine_mapping_directory, "credible_set_variants_current.txt")
        variants_file <- if (file.exists(current)) current else
            file.path(fine_mapping_directory, "credible_set_variants.txt")
    }
    variants <- fread(variants_file, na.strings = c("NA", ""))

    # Copied to differently named locals because the arguments share names with
    # columns of `variants`. The `..` prefix only works in j for column
    # selection, not in an i expression.
    target_cancer   <- cancer
    target_ancestry <- ancestry_group
    target_locus    <- locus_id
    target_set      <- credible_set

    # the table holds only credible set members, so every row at this locus
    # belongs to one of its sets
    at_locus <- variants[cancer == target_cancer & ancestry_group == target_ancestry &
                         locus_id == target_locus]
    members <- at_locus[CS == target_set]
    if (nrow(members) == 0) {
        stop(sprintf("No credible set '%s' at %s %s %s. Available: %s",
                     credible_set, cancer, ancestry_group, locus_id,
                     if (nrow(at_locus) == 0) "none"
                     else paste(unique(at_locus$CS), collapse = ", ")))
    }

    others <- at_locus[CS != target_set]

    if (is.null(index_snv)) index_snv <- members[which.max(PIP), rsID]
    chromosome <- as.character(members$chr[1])

    plink_stem <- file.path(fine_mapping_directory, "plink", cancer,
                            sprintf("%s_%s_%s_1KG_filtered", cancer, locus_id, ancestry_group))
    ld_path  <- paste0(plink_stem, "_ld.ld")
    bim_path <- paste0(plink_stem, ".bim")
    for (required in c(ld_path, bim_path)) {
        if (!file.exists(required)) {
            stop(sprintf(paste("Missing %s.\nRegenerate it by running fine_mapping() for this",
                               "locus with keep_plink_files = TRUE (see",
                               "regenerate_plink_for_figures.R)."), required))
        }
    }

    if (is.null(ccres_path)) {
        ccres_path <- file.path(PROJECT_DIRECTORY, "data/annotations",
                                sprintf("GRCh37-cCREs_chr%s.bed", chromosome))
        if (!file.exists(ccres_path)) {
            warning(sprintf("No cCRE file at %s; the track will be omitted.", ccres_path))
            ccres_path <- NULL
        }
    }

    adaptive <- NULL
    if (!is.null(adaptive_introgression_file) && file.exists(adaptive_introgression_file)) {
        adaptive <- fread(adaptive_introgression_file, header = FALSE,
                          col.names = c("chr", "start", "end"))
    }

    # the metadata is only consulted for what was not supplied, so a file with no
    # metadata entry can be plotted by giving all five
    supplied <- c(gwas_path, chr_col, pos_col, pval_col, rsid_col)
    columns <- if (length(supplied) < 5) gwas_columns(cancer) else list()

    gwas_path <- gwas_path %||% columns$input
    chr_col   <- chr_col   %||% columns$chr_col
    pos_col   <- pos_col   %||% columns$pos_col
    pval_col  <- pval_col  %||% columns$pval_col
    rsid_col  <- rsid_col  %||% columns$rsid_col

    # NULL suppresses the circle layers, which are guarded on it. The index
    # variant carries its own symbol and is never circled as well.
    inputs <- list(
        gwas_path = gwas_path,
        ld_path   = ld_path,
        bim_path  = bim_path,
        index_snv = index_snv,
        cs_snvs   = if (mark_credible_sets) setdiff(members$rsID, index_snv) else NULL,
        other_cs_snvs = if (mark_credible_sets && nrow(others) > 0) others$rsID else NULL,
        chr_col   = chr_col,
        pos_col   = pos_col,
        pval_col  = pval_col,
        rsid_col  = rsid_col,
        ccres_path = ccres_path,
        adaptive_introgressed_regions = adaptive
    )

    message(sprintf("%s | %s | %s | %s", cancer, ancestry_group, locus_id, credible_set))
    message(sprintf("  index variant: %s (chr%s:%s, PIP %.3f)",
                    index_snv, chromosome,
                    format(members[rsID == index_snv, pos], big.mark = ","),
                    members[rsID == index_snv, PIP]))
    message(sprintf("  credible set:  %d variants, %d introgressed, spanning %.0f kb",
                    nrow(members), sum(members$is_neanderthal, na.rm = TRUE),
                    (max(members$pos) - min(members$pos)) / 1000))

    # the other sets are reported because they decide the plotting window: one of
    # them may carry a stronger peak than the set being shown
    if (nrow(others) > 0) {
        for (other_set in unique(others$CS)) {
            other_members <- others[CS == other_set]
            lead <- other_members[which.max(PIP)]
            message(sprintf("  also at locus: %s, %d variants, lead %s (chr%s:%s, PIP %.3f)",
                            other_set, nrow(other_members), lead$rsID, chromosome,
                            format(lead$pos, big.mark = ","), lead$PIP))
        }
    }

    message(sprintf("  sumstats:      %s", basename(gwas_path)))
    message(sprintf("  columns:       %s / %s / %s / %s",
                    chr_col, pos_col, pval_col, rsid_col))
    if (!is.null(columns$build) && !columns$build %in% c("hg19", "GRCh37", "grch37")) {
        warning(sprintf(paste("%s positions are %s, not hg19; the gene track, cCREs and",
                              "introgression regions are hg19."), cancer, columns$build))
    }
    if (!file.exists(gwas_path)) {
        stop(sprintf("Summary statistics not found: %s", gwas_path))
    }

    invisible(inputs)
}


# --- resolve, then plot ------------------------------------------------------
# window_kb is the total width, centred on the index variant. upstream_kb and
# downstream_kb are distances from the index variant and replace it; they are
# required together, because supplying one alone would silently take the other
# from window_kb / 2 and produce a window of neither the requested width nor the
# requested centre.
locuszoom_for_credible_set <- function(cancer, ancestry_group, locus_id, credible_set,
                                       index_snv = NULL,
                                       window_kb = 200,
                                       upstream_kb = NULL,
                                       downstream_kb = NULL,
                                       mark_credible_sets = FALSE,
                                       gwas_path = NULL,
                                       chr_col = NULL,
                                       pos_col = NULL,
                                       pval_col = NULL,
                                       rsid_col = NULL,
                                       ccres_path = NULL, ...) {

    if (xor(is.null(upstream_kb), is.null(downstream_kb))) {
        stop("Pass both upstream_kb and downstream_kb, or neither. ",
             "window_kb is the total width of a window centred on the index variant.")
    }

    inputs <- resolve_locus_inputs(cancer, ancestry_group, locus_id, credible_set,
                                   index_snv = index_snv,
                                   mark_credible_sets = mark_credible_sets,
                                   gwas_path = gwas_path,
                                   chr_col = chr_col, pos_col = pos_col,
                                   pval_col = pval_col, rsid_col = rsid_col,
                                   ccres_path = ccres_path)

    window <- if (is.null(upstream_kb)) list(window_kb = window_kb)
              else list(upstream_kb = upstream_kb, downstream_kb = downstream_kb)

    do.call(generate_locuszoom_plot_with_annotations, c(inputs, window, list(...)))
}
