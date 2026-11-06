# LocusZoom-style regional plot with gene, cCRE and adaptive-introgression tracks.
#
# All layout lives in this file. The figures are not to be re-spaced in
# Inkscape: if a gap is wrong, change `track_heights`, `track_bar_fraction` or
# `track_label_width` here and every figure moves together.
#
# Layout model
# ------------
#   * `track_heights` fixes the four row heights. Boundaries sit midway between
#     adjacent tracks, so the visible gaps are a consequence of the row heights
#     and of how much of its row each bar fills -- not of plot margins. Every
#     margin is zero.
#   * `track_bar_fraction` is the share of its row that the cCRE and adaptive
#     introgression bars occupy; the rest is the gap above and below them.
#   * Each track carries a horizontal y-axis title ("Genes", "cCREs",
#     "Adaptive\nIntrogression"). Giving all four plots a y-axis title is also
#     what makes patchwork align their panels with the scatter.
#   * `title` draws the panel title (14 pt, centred on the full plot width,
#     `plot.title.position = "plot"`).
#
# Defaults are measured off figures/main_figs/fig_4.svg at 1152 px:
#   title band 0-17; genes 50-122; cCREs 140-151; adaptive introgression
#   164-188; scatter 211-497; x-axis labels to 511. Row allocations
#   106 / 31 / 38 / 311 px give c(1, 0.29, 0.36, 2.93); the cCRE bar is 8 px of
#   its 31 px row and the introgression bar 9 px of its 38 px row, hence 0.25.
#
# Fixes retained from the previous version:
#   1. ccres_path = NULL no longer errors (ccres was referenced outside the if)
#   2. an empty GWAS window fails with a readable message instead of a type
#      error forty frames deep inside left_join()
#   3. gff_path is an argument with an absolute default, so the function does
#      not depend on the working directory
#   4. the index variant must be present in the GWAS region, not just the .bim

generate_locuszoom_plot_with_annotations <- function(gwas_path, ld_path, bim_path, index_snv,
                                                     cs_snvs = NULL,
                                                     other_cs_snvs = NULL,
                                                     window_kb = 200,
                                                     upstream_kb = NULL,
                                                     downstream_kb = NULL,
                                                     start_bp_override = NULL,
                                                     end_bp_override = NULL,
                                                     chr_col = "CHR",
                                                     pos_col = "BP",
                                                     pval_col = "P",
                                                     rsid_col = "SNP",
                                                     ccres_path = NULL,
                                                     adaptive_introgressed_regions = NULL,
                                                     custom_features = NULL,
                                                     title = NULL,
                                                     max_gene_tracks = 2,
                                                     track_heights = c(
                                                       genes = 1,
                                                       ccres = 0.4,
                                                       introgression = 0.4,
                                                       main = 2.93
                                                     ),
                                                     track_bar_fraction = 0.25,
                                                     track_labels = c(
                                                       genes = "Genes",
                                                       ccres = "cCREs",
                                                       introgression = "Adaptive\nIntrogression"
                                                     ),
                                                     show_track_labels = TRUE,
                                                     arrows_per_window = 65,
                                                     exon_height = 0.2,
                                                     chevron_height = 0.07,
                                                     chevron_width = 0.15,
                                                     gff_path = file.path(
                                                       "/wynton/group/capra/projects",
                                                       "neanderthal_ancestry_and_cancer_risk",
                                                       "data/annotations/gencode.v50lift37.basic.annotation.gff3.gz"
                                                     )
                                                    ) {

  hairline <- 0.125 / 0.377

  bim <- fread(bim_path, col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  index_snv_info <- bim[SNP == index_snv]
  if (nrow(index_snv_info) == 0) stop(paste("Target SNP", index_snv, "not found in .bim."))

  chrom <- as.character(index_snv_info$CHR)
  center_bp <- index_snv_info$BP

  if (is.null(upstream_kb) && is.null(downstream_kb)) {

    window_bp <- window_kb * 1000

    start_bp <- center_bp - window_bp / 2
    end_bp   <- center_bp + window_bp / 2

  } else {

    if (is.null(upstream_kb))
      upstream_kb <- window_kb / 2

    if (is.null(downstream_kb))
      downstream_kb <- window_kb / 2

    start_bp <- center_bp - upstream_kb * 1000
    end_bp   <- center_bp + downstream_kb * 1000
  }

  # Two panels of the same locus in different ancestry groups have different
  # index variants, so a window derived from the index cannot be shared. Setting
  # the bounds directly is what lets both panels span identical coordinates.
  if (!is.null(start_bp_override)) start_bp <- start_bp_override
  if (!is.null(end_bp_override))   end_bp   <- end_bp_override

  if (end_bp <= start_bp) {
    stop(sprintf("Window end (%s) is not greater than window start (%s).",
                 format(end_bp, scientific = FALSE),
                 format(start_bp, scientific = FALSE)))
  }

  ld_mat <- fread(ld_path, header = FALSE)
  index_snv_idx <- which(bim$SNP == index_snv)
  ld_lookup <- tibble(SNP = bim$SNP, R2 = as.numeric(ld_mat[[index_snv_idx]])^2)

  header <- names(fread(gwas_path, nrows = 0))
  c_idx <- which(header == chr_col)
  p_idx <- which(header == pos_col)

  # FIX 2a: check the column names before building the awk command. A missing
  # name produces an empty index, sprintf writes "$ ==" and awk exits 1 with a
  # message about the temporary directory being full.
  for (pair in list(c(chr_col, length(c_idx)), c(pos_col, length(p_idx)),
                    c(pval_col, sum(header == pval_col)),
                    c(rsid_col, sum(header == rsid_col)))) {
    if (as.integer(pair[2]) != 1) {
      stop(sprintf("Column '%s' is not present exactly once in %s. Columns are: %s",
                   pair[1], basename(gwas_path), paste(header, collapse = ", ")))
    }
  }

  gwas_cmd <- sprintf("zcat %s | awk -F'\\t' 'NR==1 || ($%d == \"%s\" && $%d >= %d && $%d <= %d)'",
                      gwas_path, c_idx, chrom, p_idx, start_bp, p_idx, end_bp)

  gwas_region <- fread(cmd = gwas_cmd)

  # FIX 2b: an empty region types every column logical, and the failure then
  # surfaces as an incompatible-type error inside left_join(). Say what is
  # actually wrong.
  if (nrow(gwas_region) == 0) {
    stop(sprintf(paste("No variants in chr%s:%s-%s of %s.",
                       "Check that the file is the same genome build as the .bim",
                       "and that the chromosome is written the same way in both."),
                 chrom, format(start_bp, scientific = FALSE),
                 format(end_bp, scientific = FALSE), basename(gwas_path)))
  }

  gwas_region <- gwas_region %>%
    rename(CHR = !!chr_col, BP = !!pos_col, P = !!pval_col, SNP = !!rsid_col) %>%
    left_join(ld_lookup, by = "SNP") %>%
    replace_na(list(R2 = 0)) %>%
    mutate(R2_bin = cut(R2,
                        breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
                        labels = c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")))

  # FIX 4: the index variant is drawn from gwas_region, so being in the .bim is
  # not enough; without this it silently disappears from the plot.
  if (!index_snv %in% gwas_region$SNP) {
    stop(sprintf("%s is in the .bim but not in the GWAS window; it cannot be drawn as the index variant.",
                 index_snv))
  }

  gff_cmd <- sprintf("zcat %s | awk -F'\\t' '!/^#/ && ($1 == \"%s\" || $1 == \"chr%s\") && $5 >= %d && $4 <= %d'",
                     gff_path, chrom, chrom, start_bp, end_bp)

  gff <- fread(cmd = gff_cmd, header = FALSE)
  if (nrow(gff) == 0) stop("No gene features found in this window.")

  gff <- gff %>%
    rename(
      seqid = V1,
      type = V3,
      start = V4,
      end = V5,
      strand = V7,
      attr = V9
    ) %>%
    filter(type %in% c("transcript", "exon"))

  gff[, gene_name := str_extract(attr, "(?<=gene_name=)[^;]+")]
  gff[, trans_id := str_extract(attr, "(?<=transcript_id=)[^;]+")]
  gff[, gene_type := str_extract(attr, "(?<=gene_type=)[^;]+")]

  gene_list <- gff[type == "transcript"] %>%
    mutate(len = end - start) %>%
    group_by(gene_name) %>%
    slice_max(len, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      midpoint = (start + end) / 2
    )

  # Track packing. Genes are considered protein-coding first and longest first,
  # so that when the window holds more genes than `max_gene_tracks` lanes can
  # take, what is dropped is the lncRNA / pseudogene / snRNA clutter rather than
  # the gene the panel is about. This is what keeps RP11-472N19.3 and RNU6-366P
  # out of the SLC24A4 window.
  gene_list <- gene_list %>%
    arrange(
      desc(!is.na(gene_type) & gene_type == "protein_coding"),
      desc(len)
    )

  gene_list$track <- NA_integer_

  for (i in seq_len(nrow(gene_list))) {
    for (t in seq_len(max_gene_tracks)) {
      existing <- which(gene_list$track == t)

      if (
        length(existing) == 0 ||
        all(
          gene_list$start[i] > gene_list$end[existing] |
          gene_list$end[i] < gene_list$start[existing]
        )
      ) {
        gene_list$track[i] <- t
        break
      }
    }
  }

  dropped_genes <- gene_list$gene_name[is.na(gene_list$track)]

  if (length(dropped_genes) > 0) {
    message(sprintf(
      "%d gene(s) did not fit %d track(s) and were dropped: %s",
      length(dropped_genes), max_gene_tracks,
      paste(dropped_genes, collapse = ", ")
    ))
    gene_list <- gene_list %>% filter(!is.na(track))
  }

  gene_list <- gene_list %>% arrange(start, end)

  gene_features <- gff %>%
    filter(trans_id %in% gene_list$trans_id) %>%
    select(-gene_name) %>%
    left_join(
      gene_list %>%
        select(trans_id, gene_name, track, midpoint),
      by = "trans_id"
    )

  if (!is.null(custom_features)) {

    custom_gene_list <- custom_features %>%
      filter(type == "gene") %>%
      mutate(
        midpoint = (start + end) / 2
      )

    custom_gene_list$track <- NA_integer_

    for (i in seq_len(nrow(custom_gene_list))) {
      for (t in seq_len(max_gene_tracks)) {
        existing <- which(gene_list$track == t)

        if (
          length(existing) == 0 ||
          all(
            custom_gene_list$start[i] > gene_list$end[existing] |
            custom_gene_list$end[i] < gene_list$start[existing]
          )
        ) {
          custom_gene_list$track[i] <- t
          break
        }
      }
    }

    if (any(is.na(custom_gene_list$track))) {
      message(sprintf(
        "%d custom feature(s) did not fit %d track(s) and were placed in track %d.",
        sum(is.na(custom_gene_list$track)), max_gene_tracks, max_gene_tracks
      ))
      custom_gene_list$track[is.na(custom_gene_list$track)] <- max_gene_tracks
    }

    gene_list <- bind_rows(
      gene_list,
      custom_gene_list
    )

    custom_features <- custom_features %>%
      left_join(
        custom_gene_list %>%
          select(
            gene_name,
            track,
            midpoint
          ),
        by = "gene_name"
      )

    gene_features <- bind_rows(
      gene_features,
      custom_features
    )
  }

  # Strand arrows. Open chevrons at a fixed spacing in window units, so their
  # density does not change with the window width. fig_4 carries roughly 65
  # across the panel.
  arrow_step <- (end_bp - start_bp) / arrows_per_window

  # coord_cartesian keeps ggplot's default 5% expansion, so the panel extends
  # past [start_bp, end_bp] and a gene body running off the edge is drawn into
  # that margin. Clamp the chevrons to the panel, not the window, or they stop
  # 5% short of where the line does.
  x_pad <- (end_bp - start_bp) * 0.05
  panel_start <- start_bp - x_pad
  panel_end   <- end_bp + x_pad

  gene_strand_marks <- NULL

  if (nrow(gene_list) > 0 && "strand" %in% names(gene_list)) {

    marks <- lapply(seq_len(nrow(gene_list)), function(i) {

      s <- gene_list$strand[i]

      if (is.na(s) || !s %in% c("+", "-")) return(NULL)

      visible_start <- max(gene_list$start[i], panel_start)
      visible_end   <- min(gene_list$end[i], panel_end)

      half_w <- arrow_step * chevron_width

      span <- visible_end - visible_start

      if (span <= 2 * half_w) return(NULL)

      # length.out rather than by, so a chevron lands at both ends of the gene
      # instead of leaving up to a full step of bare line at the far end
      n <- max(2, round(span / arrow_step))

      centres <- seq(visible_start + half_w, visible_end - half_w, length.out = n)

      direction <- if (s == "+") 1 else -1
      tr <- gene_list$track[i]

      data.frame(
        id = paste(i, rep(seq_len(n), each = 3), sep = "_"),
        x = as.vector(rbind(
          centres - direction * half_w,
          centres + direction * half_w,
          centres - direction * half_w
        )),
        y = rep(c(tr + chevron_height, tr, tr - chevron_height), times = n)
      )
    })

    gene_strand_marks <- do.call(rbind, marks)
  }

  # cCREs Data
  # FIX 1: ccres is now always defined, so a NULL ccres_path draws an empty
  # track instead of erroring on a missing object.
  ccres <- NULL

  if (!is.null(ccres_path)) {
    ccres <- fread(ccres_path, header = FALSE)

    if (ncol(ccres) == 6) {
      setnames(
        ccres,
        c("chr", "start", "end", "EH38D", "EH38E", "class")
      )
    } else if (ncol(ccres) == 11) {
      setnames(
        ccres,
        c(
          "chr", "start", "end", "EH38E", "score", "strand",
          "thickStart", "thickEnd", "RGB", "class", "description"
        )
      )
    } else {
      stop(
        paste0(
          "Unexpected number of columns in cCRE file: ",
          ncol(ccres),
          ". Expected 6 or 11."
        )
      )
    }

    ccres <- ccres %>%
      filter(
        chr %in% c(chrom, paste0("chr", chrom)),
        end >= start_bp,
        start <= end_bp,
        class != "Unclassified"
      )
  }

  # Adaptive Introgressed Region Data
  if (!is.null(adaptive_introgressed_regions)) {

    adaptive_introgressed_regions <- adaptive_introgressed_regions %>%
      filter(
        chr %in% c(chrom, paste0("chr", chrom)),
        end >= start_bp,
        start <= end_bp
      )

  } else {
    adaptive_introgressed_regions <- NULL
  }

  # ---------------------------------------------------------------------------
  # Shared track geometry
  # ---------------------------------------------------------------------------

  # A bar drawn between y = 0 and y = 1 should fill `track_bar_fraction` of its
  # row, centred, with the remainder as the gap above and below it.
  bar_half_span <- 0.5 / track_bar_fraction
  bar_ylim <- c(0.5 - bar_half_span, 0.5 + bar_half_span)

  track_label <- function(label) {
    if (!show_track_labels || is.null(label) || is.na(label)) {
      return(theme(axis.title.y = element_blank()))
    }

    theme(
      axis.title.y = element_text(
        family = "Arial",
        size = 11,
        angle = 0,
        hjust = 1,
        vjust = 0.5,
        margin = margin(r = 4)
      )
    )
  }

  # Every track carries the same zero margins; all spacing comes from
  # `track_heights` and `bar_ylim`.
  track_theme <- function(label) {
    list(
      theme_void(),
      theme(
        plot.margin = margin(t = 0, r = 5, b = 0, l = 5),
        axis.title.x = element_blank()
      ),
      track_label(label),
      ylab(label)
    )
  }

  # An empty track, used wherever an annotation is absent, so that the row
  # heights stay the same across figures.
  empty_track <- function(label) {
    ggplot() +
      track_theme(label) +
      coord_cartesian(
        xlim = c(start_bp / 1e6, end_bp / 1e6),
        ylim = bar_ylim,
        clip = "on"
      )
  }

  # ---------------------------------------------------------------------------
  # Gene track
  # ---------------------------------------------------------------------------

  p_genes <- ggplot() +
    geom_segment(
      data = gene_list,
      aes(
        x = start / 1e6,
        xend = end / 1e6,
        y = track,
        yend = track
      ),
      colour = "grey40",
      linewidth = 0.5
    ) +
    { if (!is.null(gene_strand_marks) && nrow(gene_strand_marks) > 0)
      geom_path(
        data = gene_strand_marks,
        aes(
          x = x / 1e6,
          y = y,
          group = id
        ),
        colour = "grey40",
        linewidth = hairline,
        lineend = "butt",
        linejoin = "mitre"
      )
    } +
    geom_rect(
      data = filter(gene_features, type == "exon"),
      aes(
        xmin = start / 1e6,
        xmax = end / 1e6,
        ymin = track - exon_height,
        ymax = track + exon_height
      ),
      fill = "black",
      colour = "black",
      linewidth = 0.1
    ) +
    geom_text(
      data = gene_list,
      aes(
        x = (
          pmax(start, start_bp) +
          pmin(end, end_bp)
        ) / 2e6,
        y = track + 0.25,
        label = gene_name
      ),
      fontface = "italic",
      family = "Arial",
      size = 11 / .pt,
      hjust = 0.5,
      vjust = 0
    ) +
    track_theme(track_labels[["genes"]]) +
    coord_cartesian(
      xlim = c(start_bp / 1e6, end_bp / 1e6),
      ylim = c(0.5, max_gene_tracks + 0.8),
      clip = "on"
    )

  locuszoom_colors <- c(
    "0.8-1.0" = "red",
    "0.6-0.8" = "orange",
    "0.4-0.6" = "#1ec41e",
    "0.2-0.4" = "#57daed",
    "0.0-0.2" = "navy"
  )

  # ---------------------------------------------------------------------------
  # cCRE track
  # ---------------------------------------------------------------------------

  if (!is.null(ccres) && nrow(ccres) > 0) {
    p_ccres <- ggplot() +
      geom_rect(
        data = ccres,
        aes(
          xmin = start / 1e6,
          xmax = end / 1e6,
          ymin = 0,
          ymax = 1
        ),
        colour = "black",
        fill = "black",
        linewidth = 0.1
      ) +
      track_theme(track_labels[["ccres"]]) +
      coord_cartesian(
        xlim = c(start_bp / 1e6, end_bp / 1e6),
        ylim = bar_ylim,
        clip = "on"
      )
  } else {
    p_ccres <- empty_track(track_labels[["ccres"]])
  }

  # ---------------------------------------------------------------------------
  # Adaptive introgression track
  # ---------------------------------------------------------------------------

  if (!is.null(adaptive_introgressed_regions) && nrow(adaptive_introgressed_regions) > 0) {
    p_adaptive_introgressed_region <- ggplot() +
      geom_rect(
        data = adaptive_introgressed_regions,
        aes(
          xmin = start / 1e6,
          xmax = end / 1e6,
          ymin = 0,
          ymax = 1
        ),
        colour = "black",
        fill = "black",
        linewidth = 0.1
      ) +
      track_theme(track_labels[["introgression"]]) +
      coord_cartesian(
        xlim = c(start_bp / 1e6, end_bp / 1e6),
        ylim = bar_ylim,
        clip = "on"
      )
  } else {
    p_adaptive_introgressed_region <- empty_track(track_labels[["introgression"]])
  }

  # ---------------------------------------------------------------------------
  # LocusZoom scatter
  # ---------------------------------------------------------------------------

  p_main <- ggplot(gwas_region, aes(x = BP / 1e6, y = -log10(P))) +
    geom_point(
      aes(fill = R2_bin),
      shape = 21,
      colour = "black",
      stroke = 0.2,
      size = 1.8
    ) +
    # members of the locus's other credible sets, drawn first and more faintly, so
    # that a second peak in the window is marked as a signal in its own right
    # without competing with the credible set the panel is about
    { if (!is.null(other_cs_snvs))
      geom_point(
        data = filter(gwas_region, SNP %in% setdiff(other_cs_snvs, cs_snvs)),
        shape = 21,
        size = 3.5,
        stroke = 0.5,
        fill = NA,
        colour = "grey55"
      )
    } +
    { if (!is.null(cs_snvs))
      geom_point(
        data = filter(gwas_region, SNP %in% cs_snvs),
        shape = 21,
        size = 3.5,
        stroke = 0.5,
        fill = NA,
        colour = "black"
      )
    } +
    geom_point(
      data = filter(gwas_region, SNP == index_snv),
      colour = "purple",
      shape = 18,
      size = 5.5
    ) +
    geom_text(
      data = filter(gwas_region, SNP == index_snv),
      aes(label = SNP),
      vjust = -0.8,
      colour = "black",
      family = "Arial",
      size = 9 / .pt
    ) +
    scale_fill_manual(
      values = locuszoom_colors,
      # plotmath sets a bare r upright; italic(r) is what matches the figures
      name = expression(italic(r)^2),
      guide = guide_legend(reverse = TRUE, title.hjust = 0.5)
    ) +
    theme_bw() +
    xlab(paste("Chromosome", chrom, "Position (Mb)")) +
    ylab(expression(-log[10](italic(P)))) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    expand_limits(x = c(start_bp / 1e6, end_bp / 1e6)) +
    fig_theme +
    theme(
      plot.margin = margin(t = 0, r = 5, b = 0, l = 5),
      axis.line.x.bottom = element_line(linewidth = hairline),
      axis.line.y.left = element_line(linewidth = hairline),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.ticks = element_line(linewidth = hairline),
      axis.title = element_text(size = 12),
      axis.title.x = element_text(margin = margin(b = 0, l = 0, r = 0, t = 15)),
      axis.title.y = element_text(margin = margin(b = 0, l = 0, r = 15, t = 0)),
      legend.position = c(0.97, 0.6),
      legend.justification = c(1, 0.5),
      legend.background = element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.3
      ),
      legend.key = element_blank(),
      legend.spacing.y = unit(0.05, "cm"),
      legend.text = element_text(
        family = "Arial",
        size = 11
      ),
      legend.title = element_text(
        family = "Arial",
        size = 11,
        hjust = 0.5
      )
    )

  # ---------------------------------------------------------------------------
  # Combine
  # ---------------------------------------------------------------------------

  combined <- p_genes /
    p_ccres /
    p_adaptive_introgressed_region /
    p_main +
    plot_layout(
      heights = track_heights
    )

  if (!is.null(title)) {
    combined <- combined +
      plot_annotation(
        title = title,
        theme = theme(
          plot.title = element_text(
            family = "Arial",
            size = 14,
            hjust = 0.5,
            margin = margin(b = 6)
          ),
          # "plot" centres on the full width rather than on the panel; in
          # fig_4.svg the title's centre is 575.5 of 1152, i.e. the canvas
          # centre, not the panel centre at 612
          plot.title.position = "plot"
        )
      )
  }

  return(combined)
}