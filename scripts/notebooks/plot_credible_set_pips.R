# Plot posterior inclusion probabilities across a locus, highlighting one
# credible set.
#
#   p <- plot_credible_set_pips(
#       result_file  = "data/fine_mapping/results/<gwas>/<locus>_<anc>.rds",
#       highlight_cs = "L2"
#   )
#
# Reads the fine-mapping result rather than the GWAS summary statistics. The
# credible set tables hold only set members, so a plot built from them has no
# PIPs for the background variants; the result object's `variants` table has a
# PIP for every variant in the fit. It also avoids re-deriving the window from
# the GWAS file, which has to agree with the fit's variant set to be meaningful.
#
# Only the named credible set is highlighted. Variants in other credible sets
# keep their true PIPs and are drawn in a third style: zeroing them would
# misrepresent the fine-mapping, which assigns posterior across all effects at
# once, and would contradict a locuszoom drawn over the same window, where those
# variants carry the association they were given posterior for.

plot_credible_set_pips <- function(
    result_file,
    highlight_cs,
    index_snv     = NULL,     # default: the highlighted set's top-PIP variant
    window_kb     = 200,
    upstream_kb   = NULL,
    downstream_kb = NULL,
    highlight_colour = "#007c83",
    other_set_colour = "grey55",
    show_other_sets  = TRUE,
    label_index      = TRUE,
    label_other_sets = TRUE
) {

    result <- readRDS(result_file)

    variants <- as.data.table(result$variants)[, .(SNP = rsID, BP = pos, PIP,
                                                   is_neanderthal)]
    cs_table <- as.data.table(result$cs_table)
    members <- cs_table[CS == highlight_cs, unique(rsID)]
    if (length(members) == 0) {
        stop(sprintf("No credible set named '%s' in %s. Available: %s",
                     highlight_cs, basename(result_file),
                     paste(unique(cs_table$CS), collapse = ", ")))
    }

    variants[, in_set := SNP %in% members]
    variants[, in_other_set := SNP %in% setdiff(cs_table$rsID, members)]

    if (is.null(index_snv)) {
        index_snv <- variants[in_set == TRUE][which.max(PIP), SNP]
    }
    if (!index_snv %in% variants$SNP) {
        stop(sprintf("Index variant %s is not in this locus.", index_snv))
    }
    center_bp <- variants[SNP == index_snv, BP][1]
    chrom <- as.character(result$locus_summary$chr)

    if (is.null(upstream_kb) && is.null(downstream_kb)) {
        start_bp <- center_bp - window_kb * 1000 / 2
        end_bp   <- center_bp + window_kb * 1000 / 2
    } else {
        if (is.null(upstream_kb))   upstream_kb   <- window_kb / 2
        if (is.null(downstream_kb)) downstream_kb <- window_kb / 2
        start_bp <- center_bp - upstream_kb * 1000
        end_bp   <- center_bp + downstream_kb * 1000
    }

    window <- variants[BP >= start_bp & BP <= end_bp]
    if (nrow(window) == 0) stop("No variants in the plotting window.")

    outside <- variants[in_set == TRUE & (BP < start_bp | BP > end_bp)]
    if (nrow(outside) > 0) {
        warning(sprintf("%d members of %s fall outside the window and are not shown.",
                        nrow(outside), highlight_cs))
    }

    ceiling_pip <- max(window$PIP, na.rm = TRUE)
    upper <- min(1, ceiling(ceiling_pip * 10) / 10)

    # drawn back to front: background, then the other credible sets, then the
    # highlighted set, so no highlighted member is hidden by a neighbour
    p <- ggplot(window, aes(x = BP / 1e6, y = PIP)) +
        geom_point(
            data = window[in_set == FALSE & (in_other_set == FALSE | !show_other_sets)],
            size = 1.5, colour = "grey78"
        ) +
        geom_point(
            data = if (show_other_sets) window[in_other_set == TRUE & in_set == FALSE]
                   else window[0],
            size = 2.8, shape = 21, fill = other_set_colour, colour = "black",
            stroke = 0.5
        ) +
        geom_point(
            data = window[in_set == TRUE & is_neanderthal == FALSE],
            size = 2.8, shape = 21, fill = "black", colour = "black", stroke = 0.5
        ) +
        geom_point(
            data = window[in_set == TRUE & is_neanderthal == TRUE],
            size = 2.8, shape = 21, fill = highlight_colour, colour = "black",
            stroke = 0.5
        ) +
        scale_y_continuous(
            limits = c(0, upper),
            breaks = scales::pretty_breaks(3),
            # headroom so a label on the top-PIP variant is not clipped
            expand = expansion(mult = c(0.02, 0.15))
        ) +
        labs(
            x = paste("Chromosome", chrom, "Position (Mb)"),
            y = "PIP"
        ) +
        expand_limits(x = c(start_bp / 1e6, end_bp / 1e6)) +
        fig_theme

    if (label_index) {
        p <- p + geom_text(
            data = window[SNP == index_snv],
            aes(label = SNP),
            vjust = -1, size = 9 / .pt, family = "Arial"
        )
    }

    # the top-PIP variant of each other credible set, so that a peak carried by a
    # different signal is named rather than left for the reader to wonder about
    if (show_other_sets && label_other_sets) {
        # PIPs are taken from the variants table rather than cs_table, which is a
        # membership listing and need not carry one
        other_sets <- unique(cs_table[!rsID %in% members, CS])
        other_leads <- vapply(other_sets, function(set) {
            in_set_variants <- variants[SNP %in% cs_table[CS == set, rsID]]
            if (nrow(in_set_variants) == 0) NA_character_
            else in_set_variants[which.max(PIP), SNP]
        }, character(1))
        other_leads <- other_leads[!is.na(other_leads)]
        p <- p + geom_text(
            data = window[SNP %in% other_leads],
            aes(label = SNP),
            vjust = -1, size = 9 / .pt, family = "Arial", colour = "grey35"
        )
    }

    p
}
