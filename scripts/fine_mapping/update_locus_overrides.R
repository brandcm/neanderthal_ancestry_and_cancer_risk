# Record which loci need settings other than the defaults.
#
# Run after a round of fine mapping. Reads the credible set summaries, finds the
# loci whose credible sets sit near a window edge, and writes them into
# locus_overrides.txt with a wider window -- except in the MHC, which is held at
# 1 Mb deliberately.
#
# The MHC exception is the reason this is a script rather than a one-off command.
# A 2 Mb window there overlaps the neighbouring MHC loci, so the same signal would
# be reported under several locus IDs, and the variant count makes the LD matrix too
# large to decompose. Those loci are reported with the window-edge caveat instead,
# and this script re-applies that decision every time the overrides are rebuilt.
#
# Existing rows are kept: a locus already in the file is left alone, so a manual
# decision is never overwritten by a later run.

suppressPackageStartupMessages({
    library(data.table)
    library(glue)
})

project_directory <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/"
fine_mapping_directory <- file.path(project_directory, "data/fine_mapping")
overrides_file <- file.path(project_directory, "scripts/fine_mapping/locus_overrides.txt")

WIDER_WINDOW <- 2000

# the extended MHC on hg19, taken wide on purpose: a locus anywhere in it inherits
# the region's LD complexity
MHC_CHROMOSOME <- 6
MHC_START <- 25e6
MHC_END <- 35e6

OVERRIDE_COLUMNS <- c("config_key", "locus_id", "window_size", "L", "reason")

MHC_REASON <- paste(
    "MHC: held at 1 Mb. A wider window overlaps the neighbouring MHC loci, so one",
    "signal would be reported under several locus IDs, and the LD matrix becomes too",
    "large to decompose. Reported with the window-edge caveat instead."
)
EDGE_REASON <- "a credible set sat within 100 kb of the 1 Mb window edge"


in_mhc <- function(loci) {
    loci$chr == MHC_CHROMOSOME & loci$pos >= MHC_START & loci$pos <= MHC_END
}


main <- function() {
    loci <- fread(file.path(fine_mapping_directory, "fine_mapping_loci.txt"))
    credible_sets <- fread(file.path(fine_mapping_directory, "credible_set_summary.txt"))

    existing <- if (file.exists(overrides_file)) {
        fread(overrides_file, colClasses = list(character = c("config_key", "locus_id", "reason")))
    } else {
        data.table(config_key = character(), locus_id = character(),
                   window_size = integer(), L = integer(), reason = character())
    }

    near_edge <- unique(credible_sets[near_window_edge == TRUE,
                                      .(config_key = cancer, locus_id)])
    near_edge <- loci[near_edge, on = .(config_key, locus_id), nomatch = NULL]

    mhc <- near_edge[in_mhc(near_edge), .(config_key, locus_id,
                                          window_size = 1000L, L = NA_integer_,
                                          reason = MHC_REASON)]
    wider <- near_edge[!in_mhc(near_edge), .(config_key, locus_id,
                                             window_size = as.integer(WIDER_WINDOW),
                                             L = NA_integer_, reason = EDGE_REASON)]

    new_rows <- rbind(mhc, wider)
    new_rows <- new_rows[!existing, on = .(config_key, locus_id)]

    overrides <- rbind(existing[, ..OVERRIDE_COLUMNS], new_rows)
    setorder(overrides, config_key, locus_id)
    fwrite(overrides, overrides_file, sep = "\t", na = "")

    message(glue(
        "{nrow(near_edge)} loci had a credible set near the window edge: ",
        "{nrow(mhc)} in the MHC (held at 1 Mb), {nrow(wider) } widened to {WIDER_WINDOW} kb."
    ))
    message(glue("{nrow(new_rows)} rows added; {nrow(existing)} kept; ",
                 "{nrow(overrides)} in {overrides_file}"))

    if (nrow(new_rows) > 0) {
        message("\nResults for these loci must be deleted before they are rerun:")
        print(new_rows[, .(config_key, locus_id, window_size)])
    }
}

main()
