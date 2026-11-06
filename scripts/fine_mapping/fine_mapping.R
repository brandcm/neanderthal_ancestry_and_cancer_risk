# =============================================================================
# Fine-mapping of cancer GWAS loci with SuSiE-RSS
#
# For one lead variant: pull a window of summary statistics, build an
# ancestry-matched LD matrix from 1000 Genomes, align alleles to the LD panel,
# check z-score/LD consistency (kriging), run susie_rss, and summarize credible
# sets with respect to Neanderthal-introgressed variants. Credible-set variants are
# annotated with the introgression maps, Neanderthal allele and Neanderthal-allele
# beta, and loci with any Neanderthal credible-set variant are added to a main TSV.
#
# Requires: data.table, glue, yaml, susieR; PLINK 1.9 available via `module load`.
# =============================================================================

library(data.table)
library(glue)
library(susieR)

# --- PROJECT PATHS ---
PROJECT_DIR   <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk"
METADATA_YAML <- file.path(PROJECT_DIR, "data/metadata/cancer_GWAS_metadata.yaml")
LD_PANEL_DIR  <- file.path(PROJECT_DIR, "data/LD_reference_panels")
FM_DIR        <- file.path(PROJECT_DIR, "data/fine_mapping")
PLINK_MODULE  <- "module load CBI plink/1.90b6.26"
INTROGRESSED_FILE <- file.path(PROJECT_DIR, "data/introgressed_variants/merged_introgressed_variants_hg19.txt")

PALINDROMIC_PAIRS <- c("AT", "TA", "CG", "GC")


# =============================================================================
# Helpers
# =============================================================================

# --- CONFIG ---
load_gwas_config <- function(config_key, metadata_yaml = METADATA_YAML, quiet = FALSE) {
  metadata <- yaml::read_yaml(metadata_yaml)
  config <- metadata[[config_key]]
  if (is.null(config)) stop(glue("Key '{config_key}' not found in {metadata_yaml}."))

  # LD panel and introgression maps are hg19, so use the lifted-over sumstats when the GWAS is hg38.
  # positions_build records which build the positions in config$input are on.
  config$positions_build <- if (is.null(config$assembly)) NA_character_ else config$assembly
  if (identical(config$assembly, "hg38")) {
    hg19_path <- file.path(dirname(config$input), paste0(config_key, "_hg19.txt.gz"))
    if (file.exists(hg19_path)) {
      if (!quiet) message(glue("Assembly is hg38. Switching input to hg19 file: {hg19_path}"))
      config$input <- hg19_path
      config$positions_build <- "hg19"
    } else {
      warning(glue("hg38 detected but {hg19_path} not found. Proceeding with original (hg38) input."))
    }
  }
  config
}

# --- SUMMARY STATISTICS ---
load_sumstats <- function(config) {
  source_cols <- c(config$rsID_column, config$chr_column, config$pos_column,
                   config$effect_allele_column, config$other_allele_column,
                   "beta", "standard_error")
  std_cols <- c("rsID", "chr", "pos", "gwas_ea", "gwas_oa", "beta", "se")

  # a sentinel such as "." in a numeric column makes fread type the whole column as
  # character, which fails later at beta / se rather than here
  sumstats <- fread(config$input, select = source_cols,
                    na.strings = c("NA", "na", ".", ""),
                    colClasses = list(character = c(config$effect_allele_column, config$other_allele_column)))
  setnames(sumstats, old = source_cols, new = std_cols)

  for (column in c("pos", "beta", "se")) {
    if (!is.numeric(sumstats[[column]])) {
      coerced <- suppressWarnings(as.numeric(sumstats[[column]]))
      n_unreadable <- sum(is.na(coerced) & !is.na(sumstats[[column]]))
      message(glue("{column}: {n_unreadable} values could not be read as a number; set to NA."))
      sumstats[, (column) := coerced]
    }
  }

  # Same cleaning as the Neanderthal list so exact matching is fair
  sumstats[, rsID := normalize_rsids(rsID)]
  sumstats
}

extract_window <- function(sumstats, lead_rsid, window_size, drop_variants = NULL) {
  if (!(lead_rsid %in% sumstats$rsID)) {
    stop(glue("Lead rsID {lead_rsid} not found in summary statistics."))
  }

  lead <- sumstats[rsID == lead_rsid][1]
  half_window_bp <- (window_size * 1000) / 2
  window_start <- lead$pos - half_window_bp
  window_end   <- lead$pos + half_window_bp

  in_window <- sumstats[chr == lead$chr & pos >= window_start & pos <= window_end]
  in_window <- in_window[!is.na(rsID) & rsID != "."]
  in_window <- in_window[!duplicated(rsID)]
  in_window <- in_window[!is.na(se) & se > 0]

  if (!is.null(drop_variants)) {
    message(glue("Manually removing {length(drop_variants)} requested variants from the window..."))
    in_window <- in_window[!(rsID %in% drop_variants)]
  }

  list(variants = in_window,
       chr = gsub("chr", "", as.character(lead$chr)),
       start = window_start,
       end = window_end)
}

# --- RSID CLEANING ---
# Exact matching only (match / %in%), never grep: grepl("rs123", ...) would also match rs1234.
normalize_rsids <- function(ids) {
  ids <- trimws(as.character(ids))
  ids <- sub("^RS", "rs", ids)
  sub("^(rs[0-9]+):.*$", "\\1", ids)                # rs123:A:G -> rs123
}

# --- INTROGRESSED VARIANTS ---
# Table built by build_introgressed_variants.R: one row per hg19 position with
# rsID, alleles (unordered pair, e.g. "A/G"), neanderthal_allele, maps, notes.
load_introgressed_variants <- function(path = INTROGRESSED_FILE) {
  intro <- fread(path, sep = "\t", colClasses = list(character = c(
    "chr", "rsID", "ref", "alt", "alleles", "neanderthal_allele", "maps", "notes")))
  intro[, rsID := normalize_rsids(rsID)]
  intro[, pos_key := paste0(chr, "_", pos_hg19)]
  message(glue("Loaded {nrow(intro)} introgressed variants from {path} ",
               "({intro[is.na(rsID), .N]} without rsID; {intro[!is.na(neanderthal_allele), .N]} with a Neanderthal allele)"))
  intro
}

allele_pair <- function(a, b) {
  a <- toupper(a); b <- toupper(b)
  fifelse(is.na(a) | is.na(b), NA_character_, paste(pmin(a, b), pmax(a, b), sep = "/"))
}

# Row of `intro` for each variant (NA if not introgressed).
# 1. rsID match.  2. Otherwise hg19 chr + pos, but only if the GWAS positions are hg19,
#    and only if the alleles agree (either strand) when the map reports alleles.
#    This catches variants without an rsID, or with a different rsID across dbSNP builds,
#    without matching a different variant at the same position.
match_introgressed <- function(rsID, chr, pos, a1, a2, intro, use_position = TRUE) {
  a1 <- toupper(a1); a2 <- toupper(a2)
  idx <- match(rsID, intro$rsID, incomparables = NA)
  if (use_position) {
    need <- which(is.na(idx))
    pos_key <- paste0(sub("^chr", "", as.character(chr[need]), ignore.case = TRUE), "_", pos[need])
    pidx <- match(pos_key, intro$pos_key)
    map_pair <- intro$alleles[pidx]
    ok <- !is.na(pidx) & (is.na(map_pair) |
                          map_pair == allele_pair(a1[need], a2[need]) |
                          map_pair == allele_pair(complement(a1[need]), complement(a2[need])))
    ok <- ok %in% TRUE                                  # NA (e.g. missing GWAS allele) -> no match
    idx[need[ok]] <- pidx[ok]
  }
  idx
}

# Adds is_neanderthal, introgression_maps, neanderthal_allele, neanderthal_beta and
# notes to aligned variants (beta oriented to ref_A1).
annotate_introgressed <- function(variants, intro, use_position) {
  idx <- variants[, match_introgressed(rsID, chr, pos, gwas_ea, gwas_oa, intro, use_position)]
  variants[, `:=`(
    is_neanderthal     = !is.na(idx),
    introgression_maps = intro$maps[idx],
    neanderthal_allele = intro$neanderthal_allele[idx],
    notes              = intro$notes[idx]
  )]
  variants[, neanderthal_beta := NA_real_]
  known <- variants[, !is.na(neanderthal_allele)]
  variants[known, neanderthal_beta := fcase(
    neanderthal_allele == ref_A1, beta,
    neanderthal_allele == ref_A2, -beta,
    !palindromic & complement(neanderthal_allele) == ref_A1, beta,
    !palindromic & complement(neanderthal_allele) == ref_A2, -beta,
    default = NA_real_
  )]
  variants[!is.na(neanderthal_allele) & is.na(neanderthal_beta),
           notes := "Neanderthal allele matches neither GWAS allele"]
  invisible(variants)
}

# --- PLINK / LD ---
run_plink <- function(args) {
  cmd <- glue("{PLINK_MODULE}; plink {args}")
  status <- system(cmd)
  if (status != 0) stop(glue("PLINK failed (exit status {status}):\n{cmd}"))
  invisible(status)
}

compute_ld <- function(rsids, chr, ancestry_group, out_prefix) {
  ld_ref <- file.path(LD_PANEL_DIR, ancestry_group, glue("1000G_{ancestry_group}_chr{chr}"))
  rsid_file <- glue("{out_prefix}_variants.txt")
  writeLines(rsids, rsid_file)

  message(glue("Running PLINK for chr{chr}..."))
  # --mac 1 drops variants monomorphic in this ancestry group and keeps singletons.
  # Equivalent to MAF >= 1/(2N) for any panel size N, without floating-point edge cases.
  n_ref <- nrow(fread(glue("{ld_ref}.fam"), header = FALSE))
  message(glue("LD panel {ancestry_group}: N = {n_ref}; min MAF = 1/(2N) = {signif(1 / (2 * n_ref), 3)}"))
  run_plink(glue("--bfile {ld_ref} --extract {rsid_file} --mac 1 ",
                 "--keep-allele-order --make-bed --out {out_prefix}"))

  # --keep-allele-order is required here too: without it PLINK 1.9 swaps A2 to
  # the major allele before --r, so LD signs would not match the .bim A1 used
  # to align betas
  run_plink(glue("--bfile {out_prefix} --keep-allele-order --r square --out {out_prefix}_ld"))

  bim <- fread(glue("{out_prefix}.bim"), header = FALSE,
               col.names = c("chr", "rsID", "cm", "pos", "ref_A1", "ref_A2"),
               colClasses = list(character = c(1, 2, 5, 6)))           # alleles as character ("T" != TRUE)
  ld <- as.matrix(fread(glue("{out_prefix}_ld.ld"), header = FALSE))
  dimnames(ld) <- list(bim$rsID, bim$rsID)

  list(bim = bim, ld = ld, n_ref = n_ref)
}

# --- ALLELE ALIGNMENT ---
complement <- function(alleles) chartr("ATCG", "TAGC", alleles)

align_alleles <- function(window_variants, bim, drop_palindromic = TRUE) {
  aligned <- merge(bim[, .(rsID, ref_A1, ref_A2)], window_variants, by = "rsID", sort = FALSE)
  allele_cols <- c("gwas_ea", "gwas_oa", "ref_A1", "ref_A2")
  aligned[, (allele_cols) := lapply(.SD, toupper), .SDcols = allele_cols]
  aligned[, palindromic := paste0(ref_A1, ref_A2) %in% PALINDROMIC_PAIRS]

  n_palindromic <- 0L
  if (drop_palindromic) {
    n_palindromic <- sum(aligned$palindromic)
    if (n_palindromic > 0) message(glue("Dropping {n_palindromic} strand-ambiguous (A/T, C/G) variants."))
    aligned <- aligned[palindromic == FALSE]
  }

  # Direct matches take priority; complement (strand) matching only for
  # non-palindromic variants, where it is unambiguous
  aligned[, `:=`(
    match_direct = gwas_ea == ref_A1 & gwas_oa == ref_A2,
    flip_direct  = gwas_ea == ref_A2 & gwas_oa == ref_A1,
    match_comp   = !palindromic & complement(gwas_ea) == ref_A1 & complement(gwas_oa) == ref_A2,
    flip_comp    = !palindromic & complement(gwas_ea) == ref_A2 & complement(gwas_oa) == ref_A1
  )]
  aligned[, keep_orient := match_direct | (!flip_direct & match_comp)]
  aligned[, do_flip := !keep_orient & (flip_direct | flip_comp)]

  n_matched   <- sum(aligned$keep_orient)
  n_flipped   <- sum(aligned$do_flip)
  n_unmatched <- sum(!(aligned$keep_orient | aligned$do_flip))
  message(glue("Allele alignment: {n_matched} matched, {n_flipped} flipped, {n_unmatched} dropped."))

  aligned[do_flip == TRUE, beta := -beta]
  aligned <- aligned[keep_orient | do_flip]
  aligned[, c("match_direct", "flip_direct", "match_comp", "flip_comp", "keep_orient", "do_flip") := NULL]

  list(variants = aligned, n_palindromic = n_palindromic, n_unmatched = n_unmatched)
}

prepare_ld_matrix <- function(ld_full, rsids, lambda_shrink = 0) {
  ld <- ld_full[rsids, rsids, drop = FALSE]
  ld[is.na(ld)] <- 0
  ld <- (ld + t(ld)) / 2
  diag(ld) <- diag(ld) + 1e-4
  ld <- stats::cov2cor(ld)

  if (lambda_shrink > 0) {
    message(glue("Applying {lambda_shrink * 100}% LD shrinkage toward the identity."))
    ld <- (1 - lambda_shrink) * ld + lambda_shrink * diag(nrow(ld))
  }
  ld
}

# --- DIAGNOSTICS ---
run_kriging <- function(z_scores, ld, n_samples, variants) {
  message("Running consistency diagnostic (Kriging)...")
  kriging <- susieR::kriging_rss(z = z_scores, R = ld, n = n_samples)
  cd <- kriging$conditional_dist

  # susieR vignette: flag logLR > 2 & |z| > 2. logLR > 10 reported as severe.
  flagged <- which(cd$logLR > 2 & abs(cd$z) > 2)
  n_severe <- sum(cd$logLR > 10)

  flag_table <- data.table(
    rsID           = variants$rsID[flagged],
    pos            = variants$pos[flagged],
    z_score        = cd$z[flagged],
    log_lr         = cd$logLR[flagged],
    is_neanderthal = variants$is_neanderthal[flagged]
  )[order(-log_lr)]

  if (nrow(flag_table) > 0) {
    message(glue("{nrow(flag_table)} SNPs flagged (logLR > 2 & |z| > 2); {n_severe} with logLR > 10:"))
    print(flag_table)
  } else {
    message("Kriging: 0 SNPs flagged (logLR > 2 & |z| > 2).")
  }

  list(result = kriging, flagged = flag_table, n_severe = n_severe)
}

# --- CREDIBLE SET SUMMARIES ---
summarize_credible_sets <- function(fitted, cs_info, variants, z_scores,
                                    window_start, window_end, edge_buffer_bp = 1e5) {
  if (is.null(cs_info$cs)) return(list(cs_table = NULL, cs_summary = NULL))

  is_nea <- variants$is_neanderthal
  cs_names <- names(cs_info$cs)
  effect_idx <- as.integer(sub("L", "", cs_names))

  cs_table <- rbindlist(lapply(seq_along(cs_info$cs), function(k) {
    idx <- cs_info$cs[[k]]
    data.table(
      CS                 = cs_names[k],
      rsID               = variants$rsID[idx],
      pos                = variants$pos[idx],
      effect_allele      = variants$ref_A1[idx],   # beta is oriented to this allele after alignment
      other_allele       = variants$ref_A2[idx],
      PIP                = fitted$pip[idx],
      alpha              = fitted$alpha[effect_idx[k], idx],
      beta               = variants$beta[idx],
      Z                  = z_scores[idx],
      purity             = cs_info$purity[k, "min.abs.corr"],
      is_neanderthal     = is_nea[idx],
      introgression_maps = variants$introgression_maps[idx],
      neanderthal_allele = variants$neanderthal_allele[idx],
      neanderthal_beta   = variants$neanderthal_beta[idx],
      notes              = variants$notes[idx]
    )
  }))

  # p_neanderthal: posterior probability that the causal variant for this effect
  # is a Neanderthal allele (sum of this effect's alpha over Neanderthal variants)
  cs_summary <- data.table(
    CS             = cs_names,
    n_variants     = lengths(cs_info$cs),
    n_neanderthal  = sapply(cs_info$cs, function(idx) sum(is_nea[idx])),
    p_neanderthal  = rowSums(fitted$alpha[effect_idx, is_nea, drop = FALSE]),
    max_pip        = sapply(cs_info$cs, function(idx) max(fitted$pip[idx])),
    max_abs_z      = sapply(cs_info$cs, function(idx) max(abs(z_scores[idx]))),
    purity         = cs_info$purity[, "min.abs.corr"],
    prior_variance = fitted$V[effect_idx],
    near_window_edge = sapply(cs_info$cs, function(idx) {
      any(variants$pos[idx] < window_start + edge_buffer_bp |
          variants$pos[idx] > window_end - edge_buffer_bp)
    })
  )

  cs_table <- cs_summary[, .(CS, p_neanderthal)][cs_table, on = "CS"]
  list(cs_table = cs_table[order(CS, -PIP)], cs_summary = cs_summary)
}


# =============================================================================
# Main function
# =============================================================================

fine_mapping <- function(introgressed,                    # from load_introgressed_variants()
                         config_key,
                         ancestry_group,
                         locus_id,
                         n_samples,
                         lead_rsid,
                         window_size = 1000,               # total window in kb (1000 = +/-500 kb)
                         L = 10,
                         max_iter = 1000,
                         min_lambda = 0.01,
                         estimate_prior_variance = TRUE,
                         refine = FALSE,
                         lambda_shrink = 0,
                         drop_variants = NULL,
                         drop_palindromic = TRUE,
                         keep_plink_files = FALSE,
                         save_results = TRUE,
                         main_file = MAIN_FILE,           # NULL to skip appending
                         sumstats = NULL) {               # optional pre-loaded sumstats (faster for many loci)

  # PLINK intermediates (.bed/.bim/.fam/.ld/.log/variant list) go to a per-run
  # temp directory and are deleted on exit. Set keep_plink_files = TRUE to keep
  # them under FM_DIR/plink/<config_key>/ for debugging.
  run_name <- glue("{config_key}_{locus_id}_{ancestry_group}_1KG_filtered")
  if (keep_plink_files) {
    plink_dir <- file.path(FM_DIR, "plink", config_key)
  } else {
    plink_dir <- tempfile(pattern = paste0(run_name, "_"))
    on.exit(unlink(plink_dir, recursive = TRUE), add = TRUE)
  }
  dir.create(plink_dir, recursive = TRUE, showWarnings = FALSE)
  out_prefix <- file.path(plink_dir, run_name)

  # --- LOAD AND WINDOW ---
  sumstats_provided <- !is.null(sumstats)
  config <- load_gwas_config(config_key, quiet = sumstats_provided)
  if (!sumstats_provided) sumstats <- load_sumstats(config)

  # Introgression maps are hg19; match by position only when GWAS positions are hg19 too
  use_position <- identical(config$positions_build, "hg19")
  if (!use_position) {
    warning(glue("{config_key}: GWAS positions are on '{config$positions_build}', not hg19. ",
                 "Matching introgressed variants by rsID only."))
  }
  count_nea <- function(dt) {
    sum(!is.na(match_introgressed(dt$rsID, dt$chr, dt$pos, dt$gwas_ea, dt$gwas_oa, introgressed, use_position)))
  }
  # Genome-wide count is slow on a full GWAS; run_all_loci reports it once per GWAS instead
  n_nea_in_gwas <- NA_integer_
  if (!sumstats_provided) {
    n_nea_in_gwas <- count_nea(sumstats)
    message(glue("Introgressed variants found in GWAS: {n_nea_in_gwas} of {nrow(introgressed)} ",
                 "({round(100 * n_nea_in_gwas / nrow(introgressed), 1)}%)."))
  }
  locus    <- extract_window(sumstats, lead_rsid, window_size, drop_variants)
  rm(sumstats)

  # --- LD REFERENCE ---
  ld_panel <- compute_ld(locus$variants$rsID, locus$chr, ancestry_group, out_prefix)

  # --- ALLELE ALIGNMENT ---
  alignment <- align_alleles(locus$variants, ld_panel$bim, drop_palindromic)
  variants  <- alignment$variants
  ld_mat    <- prepare_ld_matrix(ld_panel$ld, variants$rsID, lambda_shrink)
  z_scores  <- variants$beta / variants$se
  annotate_introgressed(variants, introgressed, use_position)

  attrition <- list(in_gwas     = n_nea_in_gwas,
                    window      = count_nea(locus$variants),
                    in_ld_panel = count_nea(locus$variants[rsID %in% ld_panel$bim$rsID]),
                    aligned     = sum(variants$is_neanderthal))

  # --- DIAGNOSTICS ---
  kriging <- run_kriging(z_scores, ld_mat, n_samples, variants)

  # --- SUSIE ---
  message("Executing susie_rss...")
  fitted <- susie_rss(
    bhat = variants$beta,
    shat = variants$se,
    R    = ld_mat,
    n    = n_samples,
    L    = L,
    max_iter = max_iter,
    min_lambda = min_lambda,
    estimate_prior_variance = estimate_prior_variance,
    refine = refine,
    check_ridge = TRUE,
    verbose = FALSE,
    tol = 1e-3
  )
  if (!isTRUE(fitted$converged)) warning(glue("susie_rss did not converge for {locus_id} ({config_key})."))
  message(glue("Model ELBO for L = {L}: {round(tail(fitted$elbo, 1), 2)}"))

  # --- CREDIBLE SETS ---
  cs_info <- susie_get_cs(fitted, Xcorr = ld_mat)
  cs <- summarize_credible_sets(fitted, cs_info, variants, z_scores,
                                locus$start, locus$end)

  if (is.null(cs$cs_summary)) {
    message("No credible sets detected.")
  } else {
    message(glue("Found {nrow(cs$cs_summary)} credible sets."))
    print(cs$cs_summary)
    print(cs$cs_table)
    sign_check <- cs$cs_table[!is.na(neanderthal_beta), .(n_signs = uniqueN(sign(neanderthal_beta))), by = CS]
    if (any(sign_check$n_signs > 1)) {
      warning(glue("Neanderthal-allele betas change sign within a credible set at {locus_id}; check allele coding."))
    }
    if (any(cs$cs_summary$near_window_edge)) {
      warning(glue("Credible set within 100 kb of the window edge for {locus_id}; rerun with a larger window_size."))
    }
  }

  message(glue("Neanderthal variants retained: window {attrition$window} -> ",
               "LD panel {attrition$in_ld_panel} -> aligned {attrition$aligned}"))

  # --- OUTPUT ---
  ids <- data.table(cancer = config_key, ancestry_group = ancestry_group,
                    locus_id = locus_id, lead_rsid = lead_rsid, chr = locus$chr)

  locus_summary <- cbind(ids, data.table(
    window_start      = locus$start,
    window_end        = locus$end,
    window_size       = window_size,
    L                 = L,
    n_variants        = nrow(variants),
    n_neanderthal     = attrition$aligned,
    n_credible_sets   = if (is.null(cs$cs_summary)) 0L else nrow(cs$cs_summary),
    n_neanderthal_cs  = if (is.null(cs$cs_summary)) 0L else sum(cs$cs_summary$n_neanderthal > 0),
    converged         = isTRUE(fitted$converged),
    elbo              = tail(fitted$elbo, 1),
    n_kriging_flagged = nrow(kriging$flagged),
    n_kriging_severe  = kriging$n_severe
  ))

  if (!is.null(cs$cs_summary)) {
    cs$cs_summary <- cbind(ids, cs$cs_summary)
    cs$cs_table   <- cbind(ids, cs$cs_table)
  }

  variant_table <- variants[, .(rsID, chr, pos, effect_allele = ref_A1, other_allele = ref_A2, beta, se)]
  variant_table[, `:=`(Z = z_scores, PIP = fitted$pip)]
  variant_table <- cbind(variant_table, variants[, .(is_neanderthal, introgression_maps, neanderthal_allele,
                                                     neanderthal_beta, notes)])

  result <- list(
    locus_summary         = locus_summary,
    model                 = fitted,
    variants              = variant_table,
    cs                    = cs_info$cs,
    cs_table              = cs$cs_table,
    cs_summary            = cs$cs_summary,
    pip                   = fitted$pip,
    ld                    = ld_mat,
    kriging               = kriging$result,
    kriging_flagged       = kriging$flagged,
    converged             = isTRUE(fitted$converged),
    neanderthal_attrition = attrition,
    settings              = list(window_size = window_size, L = L,
                                 estimate_prior_variance = estimate_prior_variance,
                                 refine = refine,
                                 lambda_shrink = lambda_shrink,
                                 drop_palindromic = drop_palindromic,
                                 drop_variants = drop_variants,
                                 n_ld_panel = ld_panel$n_ref)
  )

  if (save_results) save_locus_result(result)
  if (!is.null(main_file)) append_to_main_file(result, main_file)
  invisible(result)
}


# =============================================================================
# Saving and collecting results
# =============================================================================

RESULTS_DIR <- file.path(FM_DIR, "results")
MAIN_FILE   <- file.path(FM_DIR, "neanderthal_credible_sets.tsv")

MAIN_FILE_COLUMNS <- c("cancer", "ancestry_group", "locus", "credible_set", "rsid", "chr", "pos",
                       "effect_allele", "other_allele", "PIP", "beta", "z", "purity", "p_neanderthal",
                       "is_neanderthal", "introgression_maps", "neanderthal_allele", "neanderthal_beta", "notes")

# Adds this locus's rows to the main tab-delimited file. If ANY credible set at the
# locus contains a Neanderthal variant, every variant in every credible set at that
# locus is written (Neanderthal or not). Loci with no Neanderthal variant in any
# credible set are not written.
# Rerunning a locus replaces its old rows (no duplicates), and a rerun that no
# longer finds a Neanderthal credible set removes them. A lock file keeps
# parallel jobs from overwriting each other.
append_to_main_file <- function(result, main_file = MAIN_FILE) {
  ids <- result$locus_summary
  this_locus <- function(dt) dt$cancer == ids$cancer & dt$ancestry_group == ids$ancestry_group &
                             dt$locus == ids$locus_id

  new_rows <- NULL
  if (!is.null(result$cs_table) && any(result$cs_table$is_neanderthal)) {
    new_rows <- result$cs_table[, .(
      cancer, ancestry_group, locus = locus_id, credible_set = CS, rsid = rsID, chr, pos,
      effect_allele, other_allele, PIP, beta, z = Z, purity, p_neanderthal,
      is_neanderthal, introgression_maps, neanderthal_allele, neanderthal_beta, notes
    )]
    setorder(new_rows, credible_set, -PIP)
  }
  if (is.null(new_rows) && !file.exists(main_file)) {
    message(glue("Main file: no Neanderthal credible sets at {ids$locus_id}; nothing to write."))
    return(invisible(main_file))
  }

  dir.create(dirname(main_file), recursive = TRUE, showWarnings = FALSE)
  lock <- filelock::lock(paste0(main_file, ".lock"), timeout = 60000)
  if (is.null(lock)) stop(glue("Could not lock {main_file} within 60 s."))
  on.exit(filelock::unlock(lock), add = TRUE)

  if (file.exists(main_file)) {
    # Alleles must be read as character: fread would otherwise parse "T" as TRUE
    main <- fread(main_file, sep = "\t", colClasses = list(character = c(
      "chr", "credible_set", "effect_allele", "other_allele",
      "introgression_maps", "neanderthal_allele", "notes")))
    n_old <- sum(this_locus(main))
    main <- main[!this_locus(main)]
  } else {
    main <- NULL
    n_old <- 0L
  }

  main <- rbindlist(list(main, new_rows), use.names = TRUE, fill = TRUE)
  if (nrow(main) > 0) {
    setcolorder(main, intersect(MAIN_FILE_COLUMNS, names(main)))
    setorder(main, cancer, ancestry_group, locus, credible_set, -PIP)
  }
  fwrite(main, main_file, sep = "\t", na = "NA")

  n_new <- if (is.null(new_rows)) 0L else nrow(new_rows)
  message(glue("Main file: {ids$cancer} | {ids$locus_id} | {ids$ancestry_group}: ",
               "replaced {n_old} rows with {n_new} ({main_file})"))
  invisible(main_file)
}

# One .rds per locus: results/<cancer>/<locus_id>_<ancestry>.rds
# The LD matrix is dropped to keep files small (it can be rebuilt from the panel).
save_locus_result <- function(result, results_dir = RESULTS_DIR) {
  ids <- result$locus_summary
  out_dir <- file.path(results_dir, ids$cancer)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, glue("{ids$locus_id}_{ids$ancestry_group}.rds"))
  result$ld <- NULL
  saveRDS(result, out_file)
  message(glue("Saved results to {out_file}"))
  invisible(out_file)
}

# Combine every saved locus into three tables:
#   locus_summary: one row per locus (including loci with no credible sets)
#   cs_summary:    one row per credible set
#   cs_variants:   one row per variant in a credible set
collect_fine_mapping_results <- function(results_dir = RESULTS_DIR) {
  files <- list.files(results_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) stop(glue("No .rds results found under {results_dir}."))
  results <- lapply(files, readRDS)

  list(
    locus_summary = rbindlist(lapply(results, `[[`, "locus_summary"), fill = TRUE),
    cs_summary    = rbindlist(lapply(results, `[[`, "cs_summary"), fill = TRUE),
    cs_variants   = rbindlist(lapply(results, `[[`, "cs_table"), fill = TRUE)
  )
}

# Run fine_mapping() over a table of loci. Needs columns:
#   config_key, ancestry_group, locus_id, n_samples, lead_rsid
# Each GWAS is loaded once; a failure at one locus is logged and the rest continue.
run_all_loci <- function(loci, introgressed = load_introgressed_variants(), ...) {
  loci <- as.data.table(loci)
  failures <- list()

  for (key in unique(loci$config_key)) {
    config <- load_gwas_config(key)
    sumstats <- load_sumstats(config)
    n_nea <- sum(!is.na(match_introgressed(sumstats$rsID, sumstats$chr, sumstats$pos, sumstats$gwas_ea,
                                           sumstats$gwas_oa, introgressed,
                                           identical(config$positions_build, "hg19"))))
    message(glue("{key}: {n_nea} of {nrow(introgressed)} introgressed variants found in GWAS ",
                 "({round(100 * n_nea / nrow(introgressed), 1)}%)."))
    for (i in which(loci$config_key == key)) {
      row <- loci[i]
      message(glue("\n=== {row$config_key} | {row$locus_id} | {row$ancestry_group} ==="))
      tryCatch(
        fine_mapping(introgressed,
                     config_key     = row$config_key,
                     ancestry_group = row$ancestry_group,
                     locus_id       = row$locus_id,
                     n_samples      = row$n_samples,
                     lead_rsid      = row$lead_rsid,
                     sumstats       = sumstats,
                     ...),
        error = function(e) {
          message(glue("FAILED: {conditionMessage(e)}"))
          failures[[length(failures) + 1]] <<- cbind(row, error = conditionMessage(e))
        }
      )
    }
    rm(sumstats); gc()
  }

  failures <- rbindlist(failures, fill = TRUE)
  if (nrow(failures) > 0) warning(glue("{nrow(failures)} loci failed; see returned table."))
  invisible(failures)
}
