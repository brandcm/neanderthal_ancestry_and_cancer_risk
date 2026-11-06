# =============================================================================
# Simulation test for fine_mapping()
#
# Simulates loci with known causal variants, writes the files fine_mapping() reads
# (a GWAS metadata entry, summary statistics, a PLINK LD panel and an introgressed
# variant table), fine-maps them, and checks the results against what was simulated.
# Covers the cases that cannot be verified on real data: that a causal variant is
# recovered, that a variant coded with the opposite effect allele to the LD panel has
# its beta flipped rather than inverted, that palindromic variants are dropped, and
# that the introgression annotation and p_neanderthal are right.
#
# The last section reruns one simulated locus the way the earlier code treated it, to
# show what the two allele-handling bugs did to the results, and injects a known
# error to show that the kriging diagnostic now detects one.
#
# Usage: Rscript simulate_fine_mapping.R [output directory] [replicates]
# Requires: data.table, glue, yaml, susieR; PLINK 1.9 via `module load`.
# =============================================================================

library(data.table)
library(glue)

FINE_MAPPING_R <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/scripts/fine_mapping/fine_mapping.R"

source(FINE_MAPPING_R)

args <- commandArgs(trailingOnly = TRUE)
OUTPUT_DIR <- if (length(args) > 0) args[1] else file.path(tempdir(), "fine_mapping_simulation")
N_REPLICATES <- if (length(args) > 1) as.integer(args[2]) else 5L

# --- SIMULATION SETTINGS ---
N_GWAS        <- 8000      # individuals in the simulated GWAS
N_PANEL       <- 500       # individuals in the simulated LD panel
N_VARIANTS    <- 120       # variants in the locus
SPACING       <- 5000      # bp between variants
FIRST_POS     <- 1000000
CHROMOSOME    <- "22"
ANCESTRY      <- "SIM"
CONFIG_KEY    <- "simulated_cancer_SIM_Sim_GCST000"
LD_RHO        <- 0.97      # correlation between neighbouring variants' latent values
CAUSAL        <- c(40, 80) # indices of the causal variants
CAUSAL_BETA   <- c(0.12, -0.10)
CAUSAL_FREQ   <- 0.3       # the causal variants' allele frequency, fixed so that the
                           # checks test the code rather than the luck of the draw
# measured in the 1000 Genomes EUR panel, among variants with a minor allele
# frequency of at least 1%: the share whose alternate allele is the major allele, and
# the share that are palindromic. These set how often each bug could bite.
ALT_MAJOR_FRACTION   <- 0.182
PALINDROMIC_FRACTION <- 0.145
ALWAYS_PALINDROMIC   <- c(41, 81) # beside each causal variant, so the flip can matter
SWAPPED       <- c(30, 70) # variants whose summary statistics use the panel's A2 as effect allele


# =============================================================================
# Simulating a locus
# =============================================================================

# a haplotype is drawn as a correlated series of latent values, each thresholded at
# its variant's own allele frequency. Correlation between neighbouring variants gives
# the LD, and the threshold gives the frequency, so one does not constrain the other.
simulate_genotypes <- function(n_individuals, frequencies) {
  n_haplotypes <- 2 * n_individuals
  latent <- matrix(0, nrow = n_haplotypes, ncol = length(frequencies))
  latent[, 1] <- rnorm(n_haplotypes)

  for (variant in 2:length(frequencies)) {
    latent[, variant] <- LD_RHO * latent[, variant - 1] +
      sqrt(1 - LD_RHO^2) * rnorm(n_haplotypes)
  }

  thresholds <- qnorm(frequencies)
  haplotypes <- matrix(0L, nrow = n_haplotypes, ncol = length(frequencies))

  for (variant in seq_along(frequencies)) {
    haplotypes[, variant] <- as.integer(latent[, variant] < thresholds[variant])
  }

  haplotypes[seq(1, n_haplotypes, by = 2), ] + haplotypes[seq(2, n_haplotypes, by = 2), ]
}


# marginal association of each variant with the phenotype, by simple linear regression
simulate_summary_statistics <- function(genotypes, effects) {
  phenotype <- as.vector(genotypes %*% effects) + rnorm(nrow(genotypes))
  centered_genotypes <- scale(genotypes, center = TRUE, scale = FALSE)
  centered_phenotype <- phenotype - mean(phenotype)

  sum_squares <- colSums(centered_genotypes^2)
  beta <- as.vector(crossprod(centered_genotypes, centered_phenotype)) / sum_squares
  residual_sum_squares <- sum(centered_phenotype^2) - beta^2 * sum_squares
  se <- sqrt(residual_sum_squares / ((nrow(genotypes) - 2) * sum_squares))

  data.table(true_beta = beta, true_se = se, true_z = beta / se)
}


simulate_locus <- function(seed, directory) {
  set.seed(seed)
  dir.create(file.path(directory, "panel", ANCESTRY), recursive = TRUE, showWarnings = FALSE)

  # the alternate allele is the major allele at ALT_MAJOR_FRACTION of variants, which
  # is where PLINK assigns A1 differently without --keep-allele-order
  alt_is_major <- runif(N_VARIANTS) < ALT_MAJOR_FRACTION
  frequencies <- ifelse(alt_is_major, runif(N_VARIANTS, 0.5, 0.95),
                        runif(N_VARIANTS, 0.05, 0.5))
  # a causal variant's power depends on its frequency, so it is held at one common
  # value; a rare draw would leave its signal undetectable and fail the checks for
  # reasons that have nothing to do with the code under test
  frequencies[CAUSAL] <- CAUSAL_FREQ

  gwas_genotypes <- simulate_genotypes(N_GWAS, frequencies)
  panel_genotypes <- simulate_genotypes(N_PANEL, frequencies)

  variants <- data.table(
    chr = CHROMOSOME,
    pos = FIRST_POS + (seq_len(N_VARIANTS) - 1) * SPACING,
    rsID = sprintf("rs%06d", seq_len(N_VARIANTS)),
    ref = "A",
    alt = "G"
  )
  variants[seq_len(N_VARIANTS) %% 2 == 0, `:=`(ref = "C", alt = "T")]

  # palindromic variants at the measured rate, always including one beside each causal
  # variant; a causal variant is never palindromic, since the fixed code would drop it
  n_palindromic <- round(PALINDROMIC_FRACTION * N_VARIANTS) - length(ALWAYS_PALINDROMIC)
  palindromic <- c(ALWAYS_PALINDROMIC,
                   sample(setdiff(seq_len(N_VARIANTS),
                                  c(CAUSAL, SWAPPED, ALWAYS_PALINDROMIC)),
                          max(n_palindromic, 0)))
  variants[palindromic, `:=`(ref = "A", alt = "T")]
  variants[, is_palindromic := seq_len(.N) %in% palindromic]
  variants[, is_swapped := seq_len(.N) %in% SWAPPED]

  effects <- rep(0, N_VARIANTS)
  effects[CAUSAL] <- CAUSAL_BETA
  variants <- cbind(variants, simulate_summary_statistics(gwas_genotypes, effects))

  # the summary statistics report the alternate allele as the effect allele, except
  # for the swapped variants, which report the reference allele and the opposite sign
  sumstats <- variants[, .(
    chromosome = chr,
    base_pair_location = pos,
    rsid = rsID,
    effect_allele = ifelse(is_swapped, ref, alt),
    other_allele = ifelse(is_swapped, alt, ref),
    beta = ifelse(is_swapped, -true_beta, true_beta),
    standard_error = true_se,
    p_value = 2 * pnorm(-abs(true_z))
  )]

  sumstats_file <- file.path(directory, glue("{CONFIG_KEY}.txt"))
  fwrite(sumstats, sumstats_file, sep = "\t")

  metadata_file <- file.path(directory, "metadata.yaml")
  writeLines(c(
    glue("{CONFIG_KEY}:"),
    "  cancer: simulated",
    glue("  ancestry_group: {ANCESTRY}"),
    "  assembly: hg19",
    "  chr_column: chromosome",
    "  pos_column: base_pair_location",
    "  rsID_column: rsid",
    "  effect_allele_column: effect_allele",
    "  other_allele_column: other_allele",
    "  beta_column: beta",
    "  standard_error_column: standard_error",
    "  p_value_column: p_value",
    glue("  input: {sumstats_file}")
  ), metadata_file)

  # the panel as a VCF, so that PLINK takes REF and ALT from the file rather than
  # from the allele frequencies, then as a PLINK fileset
  vcf_file <- file.path(directory, "panel.vcf")
  genotype_strings <- apply(panel_genotypes, 2, function(genotypes) {
    c("0/0", "0/1", "1/1")[genotypes + 1]
  })

  vcf_body <- data.table(
    CHROM = CHROMOSOME, POS = variants$pos, ID = variants$rsID,
    REF = variants$ref, ALT = variants$alt, QUAL = ".", FILTER = "PASS",
    INFO = ".", FORMAT = "GT"
  )
  vcf_body <- cbind(vcf_body, as.data.table(t(genotype_strings)))

  writeLines(c(
    "##fileformat=VCFv4.2",
    glue("##contig=<ID={CHROMOSOME},length=60000000>"),
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
            sprintf("sample%03d", seq_len(N_PANEL))), collapse = "\t")
  ), vcf_file)
  fwrite(vcf_body, vcf_file, sep = "\t", col.names = FALSE, append = TRUE)

  panel_prefix <- file.path(directory, "panel", ANCESTRY,
                            glue("1000G_{ANCESTRY}_chr{CHROMOSOME}"))
  run_plink(glue("--vcf {vcf_file} --keep-allele-order --double-id --make-bed ",
                 "--out {panel_prefix}"))

  introgressed_indices <- c(CAUSAL[1], 20, 45, 100, ALWAYS_PALINDROMIC[1])
  introgressed_file <- file.path(directory, "introgressed_variants.txt")
  fwrite(variants[introgressed_indices, .(
    chr, pos_hg19 = pos, rsID, ref, alt,
    alleles = paste(pmin(ref, alt), pmax(ref, alt), sep = "/"),
    neanderthal_allele = alt, maps = "Browning, Skov", n_maps = 2L, notes = NA_character_
  )], introgressed_file, sep = "\t", na = "NA")

  list(
    directory = directory,
    variants = variants,
    metadata_file = metadata_file,
    panel_directory = file.path(directory, "panel"),
    panel_prefix = panel_prefix,
    introgressed_file = introgressed_file,
    introgressed_indices = introgressed_indices,
    causal_rsids = variants$rsID[CAUSAL],
    palindromic_rsids = variants$rsID[palindromic],
    swapped_rsids = variants$rsID[SWAPPED]
  )
}


# =============================================================================
# Checking a result against the simulation
# =============================================================================

check_result <- function(result, simulation) {
  checks <- list()
  fitted_variants <- as.data.table(result$variants)
  cs_table <- if (is.null(result$cs_table)) NULL else as.data.table(result$cs_table)
  causal_rsids <- simulation$causal_rsids

  checks[["causal variants in credible sets"]] <-
    !is.null(cs_table) && all(causal_rsids %in% cs_table$rsID)

  # a causal variant need not top its own set, since a variant in tight LD with it can
  # carry more of the posterior; what matters is that no set is spurious
  checks[["every credible set contains a causal variant"]] <-
    !is.null(cs_table) && all(vapply(split(cs_table, by = "CS"), function(credible_set) {
      any(credible_set$rsID %in% causal_rsids)
    }, logical(1)))

  checks[["each causal variant's set holds most of the posterior"]] <-
    !is.null(cs_table) && all(vapply(causal_rsids, function(rsid) {
      sets <- cs_table[rsID == rsid, unique(CS)]
      length(sets) == 1 && sum(cs_table[CS %in% sets, PIP]) > 0.5
    }, logical(1)))

  expected <- merge(fitted_variants[, .(rsID, beta, effect_allele)],
                    simulation$variants[, .(rsID, alt, true_beta, true_z)], by = "rsID")
  strong <- expected[abs(true_z) > 4]

  checks[["beta signs match the simulation"]] <-
    nrow(strong) > 0 && all(sign(strong$beta) == sign(strong$true_beta))
  checks[["effect allele is the panel's A1"]] <- all(expected$effect_allele == expected$alt)

  swapped <- expected[rsID %in% simulation$swapped_rsids]
  checks[["swapped variants were flipped"]] <-
    nrow(swapped) == length(simulation$swapped_rsids) &&
    all(sign(swapped$beta) == sign(swapped$true_beta))

  checks[["palindromic variants dropped"]] <-
    !any(simulation$palindromic_rsids %in% fitted_variants$rsID)

  expected_introgressed <- setdiff(
    simulation$variants$rsID[simulation$introgressed_indices], simulation$palindromic_rsids
  )
  checks[["introgressed variants annotated"]] <-
    setequal(fitted_variants[is_neanderthal == TRUE, rsID], expected_introgressed)

  if (!is.null(result$cs_summary)) {
    causal_set <- cs_table[rsID == causal_rsids[1], unique(CS)]
    p_neanderthal <- as.data.table(result$cs_summary)[CS %in% causal_set, p_neanderthal]
    checks[["p_neanderthal high for the introgressed causal set"]] <-
      length(p_neanderthal) == 1 && p_neanderthal > 0.5
  } else {
    checks[["p_neanderthal high for the introgressed causal set"]] <- FALSE
  }

  checks[["model converged"]] <- isTRUE(result$converged)
  checks[["no severe kriging flags"]] <- result$locus_summary$n_kriging_severe == 0

  vapply(checks, isTRUE, logical(1))
}


fine_map_simulation <- function(simulation) {
  METADATA_YAML     <<- simulation$metadata_file
  LD_PANEL_DIR      <<- simulation$panel_directory
  FM_DIR            <<- simulation$directory
  INTROGRESSED_FILE <<- simulation$introgressed_file

  fine_mapping(
    load_introgressed_variants(),
    config_key = CONFIG_KEY,
    ancestry_group = ANCESTRY,
    locus_id = "simulated_locus",
    n_samples = N_GWAS,
    lead_rsid = simulation$causal_rsids[1],
    save_results = FALSE,
    main_file = NULL
  )
}


# =============================================================================
# Demonstrating the effect of the bugs that were fixed
# =============================================================================

# label each variant with the credible set it belongs to, or NA
label_credible_sets <- function(n_variants, credible_sets) {
  labels <- rep(NA_character_, n_variants)

  if (length(credible_sets) > 0) {
    set_names <- names(credible_sets)
    if (is.null(set_names)) set_names <- paste0("L", seq_along(credible_sets))

    for (set in seq_along(credible_sets)) labels[credible_sets[[set]]] <- set_names[set]
  }

  labels
}

# Before the fix, compute_ld ran `--r square` without --keep-allele-order, so PLINK
# set A1 to the major allele for the LD matrix while the betas had been aligned to the
# A1 of the .bim. This recomputes the LD that way and refits the same summary
# statistics, so the difference is the bug's effect and nothing else.
demonstrate_ld_allele_order <- function(result, simulation) {
  fitted_variants <- as.data.table(result$variants)
  prefix <- file.path(simulation$directory, "unordered_ld")
  writeLines(fitted_variants$rsID, glue("{prefix}_variants.txt"))

  run_plink(glue("--bfile {simulation$panel_prefix} --extract {prefix}_variants.txt ",
                 "--r square --out {prefix}"))
  run_plink(glue("--bfile {simulation$panel_prefix} --extract {prefix}_variants.txt ",
                 "--make-bed --out {prefix}_bed"))

  bim <- fread(glue("{prefix}_bed.bim"), header = FALSE,
               col.names = c("chr", "rsID", "cm", "pos", "A1", "A2"),
               colClasses = list(character = c(1, 2, 5, 6)))
  ld <- as.matrix(fread(glue("{prefix}.ld"), header = FALSE))
  dimnames(ld) <- list(bim$rsID, bim$rsID)

  variants <- fitted_variants[match(bim$rsID, rsID)]
  fitted <- susie_rss(bhat = variants$beta, shat = variants$se, R = ld, n = N_GWAS,
                      L = 10, estimate_prior_variance = TRUE, refine = FALSE,
                      verbose = FALSE)
  credible_sets <- susie_get_cs(fitted, Xcorr = ld)

  n_disagreeing <- sum(bim$A1 != fitted_variants$effect_allele[match(bim$rsID, fitted_variants$rsID)])

  list(
    n_disagreeing_alleles = n_disagreeing,
    n_credible_sets = length(credible_sets$cs),
    causal_recovered = sum(vapply(credible_sets$cs, function(indices) {
      any(variants$rsID[indices] %in% simulation$causal_rsids)
    }, logical(1))),
    top_pip_variant = variants$rsID[which.max(fitted$pip)],
    max_pip = max(fitted$pip),
    pips = data.table(rsID = variants$rsID, pos = variants$pos, pip = fitted$pip,
                      credible_set = label_credible_sets(nrow(variants), credible_sets$cs))
  )
}


# Before the fix, align_alleles flipped every palindromic variant rather than dropping
# it, which inverts the beta whenever the strands already agreed. This refits with
# those variants kept and their betas inverted, as the earlier code left them.
demonstrate_palindromic_flip <- function(simulation) {
  METADATA_YAML     <<- simulation$metadata_file
  LD_PANEL_DIR      <<- simulation$panel_directory
  FM_DIR            <<- simulation$directory
  INTROGRESSED_FILE <<- simulation$introgressed_file

  # the fixed code drops palindromic variants, so they are first kept, and then their
  # betas are inverted the way the earlier align_alleles inverted every one of them
  kept <- suppressMessages(fine_mapping(
    load_introgressed_variants(),
    config_key = CONFIG_KEY,
    ancestry_group = ANCESTRY,
    locus_id = "simulated_locus",
    n_samples = N_GWAS,
    lead_rsid = simulation$causal_rsids[1],
    drop_palindromic = FALSE,
    save_results = FALSE,
    main_file = NULL
  ))

  inverted <- as.data.table(kept$variants)
  ld <- kept$ld

  palindromic_pairs <- c("AT", "TA", "CG", "GC")
  is_palindromic <- paste0(inverted$effect_allele, inverted$other_allele) %in% palindromic_pairs
  inverted[is_palindromic, beta := -beta]

  fitted <- susie_rss(bhat = inverted$beta, shat = inverted$se, R = ld, n = N_GWAS,
                      L = 10, estimate_prior_variance = TRUE, refine = FALSE,
                      verbose = FALSE)
  credible_sets <- susie_get_cs(fitted, Xcorr = ld)

  list(
    n_inverted = sum(is_palindromic),
    n_credible_sets = length(credible_sets$cs),
    causal_recovered = sum(vapply(credible_sets$cs, function(indices) {
      any(inverted$rsID[indices] %in% simulation$causal_rsids)
    }, logical(1))),
    pips = data.table(rsID = inverted$rsID, pos = inverted$pos, pip = fitted$pip,
                      credible_set = label_credible_sets(nrow(inverted), credible_sets$cs))
  )
}


# Before the fix, run_kriging read a field that kriging_rss does not return, so it
# never flagged anything. This injects a sign error into the strongest variant and
# checks that the diagnostic now reports it.
demonstrate_kriging <- function(result, simulation) {
  fitted_variants <- as.data.table(result$variants)
  ld <- result$ld

  z_scores <- fitted_variants$Z
  worst <- which.max(abs(z_scores))
  corrupted <- z_scores
  corrupted[worst] <- -corrupted[worst]

  clean <- run_kriging(z_scores, ld, N_GWAS, fitted_variants)
  flagged <- run_kriging(corrupted, ld, N_GWAS, fitted_variants)

  list(
    corrupted_variant = fitted_variants$rsID[worst],
    n_flagged_clean = nrow(clean$flagged),
    n_flagged_corrupted = nrow(flagged$flagged),
    detected = fitted_variants$rsID[worst] %in% flagged$flagged$rsID
  )
}


# =============================================================================
# Figure
# =============================================================================

# one panel per version of the code, each showing the posterior inclusion probability
# of every variant in the locus against its position, with the simulated causal
# variants marked. Under the fixed code the causal variants carry the posterior;
# under either bug the posterior moves elsewhere.
plot_comparison <- function(panels, causal_positions, file) {
  set_colours <- c("#B2182B", "#2166AC", "#1B7837", "#E08214",
                   "#762A83", "#4393C3", "#5AAE61", "#D6604D")

  pdf(file, width = 6.5, height = 7.2)
  par(mfrow = c(length(panels), 1), mar = c(2.4, 4.2, 3.2, 1.2), oma = c(5.2, 0, 0, 0),
      mgp = c(2.6, 0.6, 0), las = 1, bty = "n", family = "sans", cex = 0.95)

  for (index in seq_along(panels)) {
    panel <- panels[[index]]
    pips <- panel$pips
    sets <- sort(unique(stats::na.omit(pips$credible_set)))
    colours <- setNames(rep(set_colours, length.out = max(length(sets), 1))[seq_along(sets)], sets)

    plot(pips$pos / 1e6, pips$pip, type = "n", ylim = c(0, 1.08), axes = FALSE,
         xlab = "", ylab = "")

    # The simulated causal variants. The dashed line alone is hidden wherever a
    # credible set sits on top of it, which is exactly where the reader needs to
    # see it, so each panel also carries a tick above the plot: a coloured stem
    # with no tick over it is a set containing no causal variant.
    abline(v = causal_positions / 1e6, col = "grey60", lty = 2, lwd = 0.9)
    segments(causal_positions / 1e6, 1.01, causal_positions / 1e6, 1.08,
             col = "grey25", lwd = 1.4, xpd = NA)
    if (index == 1) {
      text(causal_positions / 1e6, 1.09, "causal", col = "grey25", cex = 0.72, pos = 3,
           offset = 0.05, xpd = NA)
    }

    axis(1, lwd = 0, lwd.ticks = 0.8, tck = -0.025, col.axis = "grey25")
    axis(2, at = seq(0, 1, 0.5), lwd = 0, lwd.ticks = 0.8, tck = -0.025, col.axis = "grey25")
    mtext("PIP", side = 2, line = 2.6, las = 0, cex = 0.8, col = "grey25")

    # every variant as a stem, so the credible sets stand out against a flat locus
    in_set <- !is.na(pips$credible_set)
    segments(pips$pos[!in_set] / 1e6, 0, pips$pos[!in_set] / 1e6, pips$pip[!in_set],
             col = "grey88", lwd = 0.8)
    points(pips$pos[!in_set] / 1e6, pips$pip[!in_set], pch = 16, cex = 0.45, col = "grey72")

    if (any(in_set)) {
      members <- pips[in_set]
      segments(members$pos / 1e6, 0, members$pos / 1e6, members$pip,
               col = colours[members$credible_set], lwd = 1.6)
      points(members$pos / 1e6, members$pip, pch = 16, cex = 1.1,
             col = colours[members$credible_set])
    }

    mtext(panel$title, side = 3, line = 1.7, adj = 0, font = 2, cex = 0.92, col = "grey15")
    mtext(panel$subtitle, side = 3, line = 0.5, adj = 0, cex = 0.78, col = "grey40")
  }

  mtext("Position (Mb)", side = 1, line = 1.4, cex = 0.8, col = "grey25")

  # the caption is wrapped by hand, since mtext draws a single line and clips the rest
  caption <- strwrap(
    glue("Coloured stems are credible sets; dashed lines are the simulated causal ",
         "variants. {round(100 * ALT_MAJOR_FRACTION)}% of variants have a major ",
         "alternate allele and {round(100 * PALINDROMIC_FRACTION)}% are palindromic, ",
         "as measured in the 1000 Genomes EUR panel."),
    width = 92
  )

  # line increases downward in the outer margin, so the first line of the caption
  # takes the smallest value; starting below 0 would overlap the axis label
  for (line_number in seq_along(caption)) {
    mtext(caption[line_number], side = 1, line = 0.9 + 0.85 * (line_number - 1),
          outer = TRUE, adj = 0.5, cex = 0.66, col = "grey45")
  }

  dev.off()
}


# a one-line description of what a panel shows, built from the result rather than written
describe_panel <- function(pips, causal_rsids) {
  sets <- sort(unique(stats::na.omit(pips$credible_set)))
  recovered <- sum(vapply(sets, function(set) {
    any(pips$rsID[which(pips$credible_set == set)] %in% causal_rsids)
  }, logical(1)))

  spurious <- length(sets) - recovered
  glue("{length(sets)} credible {ifelse(length(sets) == 1, 'set', 'sets')}",
       ifelse(spurious > 0, glue(", {spurious} containing no causal variant"), ""),
       "; {recovered} of {length(causal_rsids)} causal variants recovered")
}


# =============================================================================
# Run
# =============================================================================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

results <- list()
check_matrix <- NULL

for (seed in seq_len(N_REPLICATES)) {
  message(glue("\n=== Replicate {seed} of {N_REPLICATES} ==="))
  directory <- file.path(OUTPUT_DIR, glue("seed_{seed}"))
  simulation <- simulate_locus(seed, directory)

  result <- tryCatch(
    suppressMessages(fine_map_simulation(simulation)),
    error = function(e) {
      message(glue("FAILED: {conditionMessage(e)}"))
      NULL
    }
  )

  if (is.null(result)) next

  checks <- check_result(result, simulation)
  check_matrix <- rbind(check_matrix, checks)
  results[[seed]] <- list(result = result, simulation = simulation)

  message(glue("{sum(checks)} of {length(checks)} checks passed"))
  if (any(!checks)) {
    message("failed: ", paste(names(checks)[!checks], collapse = "; "))
  }
  message(glue("causal z scores: ",
               "{paste(sprintf('%.2f', simulation$variants$true_z[CAUSAL]), collapse = ', ')}"))
}

message("\n=== Checks across replicates ===")
for (check in colnames(check_matrix)) {
  passed <- sum(check_matrix[, check])
  message(sprintf("%-55s %d/%d", check, passed, nrow(check_matrix)))
}
message(glue("\n{sum(check_matrix)} of {length(check_matrix)} checks passed across ",
             "{nrow(check_matrix)} replicates"))

# --- the bugs, on the first replicate ---
first <- results[[1]]

message("\n=== What the fixed code produces ===")
message(glue("Credible sets: {first$result$locus_summary$n_credible_sets}; ",
             "variants: {first$result$locus_summary$n_variants}; ",
             "introgressed: {first$result$locus_summary$n_neanderthal}"))
print(first$result$cs_summary[, .(CS, n_variants, p_neanderthal, max_pip, purity)])

message("\n=== LD computed without --keep-allele-order (the earlier behavior) ===")
unordered_ld <- demonstrate_ld_allele_order(first$result, first$simulation)
print(unordered_ld[setdiff(names(unordered_ld), "pips")])

message("\n=== Palindromic variants flipped rather than dropped (the earlier behavior) ===")
palindromic <- demonstrate_palindromic_flip(first$simulation)
print(palindromic[setdiff(names(palindromic), "pips")])

message("\n=== Kriging diagnostic, with an injected sign error ===")
print(demonstrate_kriging(first$result, first$simulation))

cs_table <- as.data.table(first$result$cs_table)
fixed_pips <- as.data.table(first$result$variants)[, .(rsID, pos, pip = PIP)]
fixed_pips[, credible_set := cs_table$CS[match(rsID, cs_table$rsID)]]

causal_rsids <- first$simulation$causal_rsids

panels <- list(
  list(pips = fixed_pips,
       title = "Current code",
       subtitle = describe_panel(fixed_pips, causal_rsids)),
  list(pips = unordered_ld$pips,
       title = "Earlier code: LD computed without --keep-allele-order",
       subtitle = describe_panel(unordered_ld$pips, causal_rsids)),
  list(pips = palindromic$pips,
       title = "Earlier code: palindromic variants flipped rather than dropped",
       subtitle = describe_panel(palindromic$pips, causal_rsids))
)

figure_file <- file.path(OUTPUT_DIR, "fine_mapping_simulation.pdf")
plot_comparison(panels, first$simulation$variants$pos[CAUSAL], figure_file)

# the plotted values, for redrawing the figure elsewhere
plotted <- rbindlist(lapply(panels, function(panel) {
  cbind(version = panel$title, panel$pips)
}))
plotted[, is_causal := rsID %in% causal_rsids]
pips_file <- file.path(OUTPUT_DIR, "fine_mapping_simulation_pips.txt")
fwrite(plotted, pips_file, sep = "\t", na = "NA")

message(glue("\nFigure written to {figure_file}"))
message(glue("Plotted values written to {pips_file}"))

message(glue("\nSimulation files are in {OUTPUT_DIR}"))
