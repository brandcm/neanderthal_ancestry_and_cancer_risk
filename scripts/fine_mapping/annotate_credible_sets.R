#!/usr/bin/env Rscript
# Attach gene names to credible sets.
#
#   ./run_rscript.sh annotate_credible_sets.R [--all]
#
# By default annotates the genome-wide significant, Neanderthal-majority signals
# -- the table that goes to co-authors. With --all, every genome-wide significant
# signal.
#
# Each set is described by the span of its variants and by its top-PIP variant.
# Genes are reported two ways, because they answer different questions:
#   overlapping - genes whose body intersects the credible set span
#   nearest     - the closest gene to the top-PIP variant, with the distance
# A set inside a gene gets both; a set in a gene desert gets only the second,
# and the distance is what tells you whether to take the name seriously.
#
# Only protein-coding and lncRNA genes are considered, so that the name is not
# some overlapping pseudogene or small RNA. Gene models are GENCODE v50 (basic
# set, GRCh37 mapping), the release used for the gene tracks in the figures.
# GENCODE merged lincRNA into the broader lncRNA type in release 31, so lncRNA
# here also includes antisense and other non-coding transcripts.
#
# Reads the deduplicated tables, so a signal recovered through two overlapping
# index windows is annotated once.

suppressPackageStartupMessages({
    library(data.table)
    library(glue)
})

project_directory <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk"
fine_mapping_directory <- file.path(project_directory, "data/fine_mapping")
gff_file <- file.path(project_directory, "data/annotations/gencode.v50lift37.basic.annotation.gff3.gz")

annotate_all <- "--all" %in% commandArgs(trailingOnly = TRUE)

summary_file <- file.path(fine_mapping_directory,
                          "credible_set_summary_deduplicated.txt")
variants_file <- file.path(fine_mapping_directory,
                           "credible_set_variants_deduplicated.txt")
for (file in c(summary_file, variants_file)) {
    if (!file.exists(file)) {
        stop(glue("{basename(file)} does not exist. ",
                  "Run deduplicate_credible_sets.R first."))
    }
}

sets <- fread(summary_file, na.strings = c("NA", ""))
variants <- fread(variants_file, na.strings = c("NA", ""))

# is_representative already keeps one set per signal, so there is nothing further
# to collapse here
selected <- sets[genome_wide_significant == TRUE & is_representative == TRUE]
if (!annotate_all) selected <- selected[p_neanderthal > 0.5]
message(glue("Annotating {nrow(selected)} signals."))

# --- genes -----------------------------------------------------------------
message("Reading gene models...")
genes <- fread(
    cmd = glue("zcat {gff_file} | awk -F'\\t' '!/^#/ && $3 == \"gene\"'"),
    header = FALSE, select = c(1, 4, 5, 7, 9),
    col.names = c("chr", "start", "end", "strand", "attributes")
)
# GENCODE ships both GFF3 (key=value;) and GTF (key "value";) attribute syntax
# and the file extension does not always say which. sub() returns its input
# unchanged when the pattern misses, so a wrong guess yields a gene_type of the
# entire attribute string and an empty table rather than an error -- hence the
# explicit check below.
extract_attribute <- function(attributes, key) {
    gff3 <- sprintf('.*%s=([^;]+).*', key)
    gtf  <- sprintf('.*%s "([^"]+)".*', key)
    value <- ifelse(grepl(gff3, attributes), sub(gff3, "\\1", attributes), NA_character_)
    unmatched <- is.na(value)
    value[unmatched] <- ifelse(
        grepl(gtf, attributes[unmatched]),
        sub(gtf, "\\1", attributes[unmatched]),
        NA_character_
    )
    value
}

genes[, gene_name := extract_attribute(attributes, "gene_name")]
genes[, gene_type := extract_attribute(attributes, "gene_type")]

if (all(is.na(genes$gene_name))) {
    stop(glue("Could not parse gene names from {basename(gff_file)}. ",
              "First attribute string:\n{genes$attributes[1]}"))
}
genes <- genes[gene_type %in% c("protein_coding", "lncRNA")]
genes[, chr := sub("^chr", "", chr)]
genes[, attributes := NULL]
message(glue("{nrow(genes)} protein-coding and lncRNA genes."))

# --- annotate ---------------------------------------------------------------
rows <- vector("list", nrow(selected))

for (index in seq_len(nrow(selected))) {
    row <- selected[index]
    members <- variants[cancer == row$cancer & ancestry_group == row$ancestry_group &
                        locus_id == row$locus_id & CS == row$CS]
    if (nrow(members) == 0) next

    chromosome <- as.character(members$chr[1])
    span_start <- min(members$pos)
    span_end   <- max(members$pos)
    lead <- members[which.max(PIP)]

    on_chromosome <- genes[chr == chromosome]

    overlapping <- on_chromosome[start <= span_end & end >= span_start, gene_name]
    overlapping <- unique(overlapping)

    distances <- pmax(0, pmax(on_chromosome$start - lead$pos, lead$pos - on_chromosome$end))
    nearest_index <- which.min(distances)
    nearest_gene <- on_chromosome$gene_name[nearest_index]
    nearest_distance <- distances[nearest_index]

    rows[[index]] <- data.table(
        cancer = row$cancer,
        ancestry_group = row$ancestry_group,
        locus_id = row$locus_id,
        CS = row$CS,
        signal_id = row$signal_id,
        lead_rsid = lead$rsID,
        lead_pos = lead$pos,
        lead_is_neanderthal = lead$is_neanderthal,
        span_kb = round((span_end - span_start) / 1000, 1),
        n_variants = row$n_variants,
        n_neanderthal = row$n_neanderthal,
        p_neanderthal = round(row$p_neanderthal, 2),
        max_abs_z = round(row$max_abs_z, 1),
        purity = round(row$purity, 2),
        overlapping_genes = if (length(overlapping) == 0) NA_character_
                            else paste(overlapping, collapse = ","),
        nearest_gene = nearest_gene,
        nearest_kb = round(nearest_distance / 1000, 1)
    )
}

annotated <- rbindlist(rows)
setorder(annotated, -p_neanderthal)

print(annotated[, .(cancer = substr(cancer, 1, 20), anc = ancestry_group,
                    locus_id, lead = lead_rsid,
                    genes = fifelse(is.na(overlapping_genes),
                                    paste0("[", nearest_gene, " +", nearest_kb, "kb]"),
                                    overlapping_genes),
                    p_nea = p_neanderthal, maxZ = max_abs_z, pur = purity)],
      nrows = 60)

out_file <- file.path(fine_mapping_directory, "credible_set_genes.txt")
fwrite(annotated, out_file, sep = "\t", na = "NA")
message(glue("\nWrote {out_file}"))
message("Genes in [brackets] do not overlap the credible set; the distance is to the lead variant.")
