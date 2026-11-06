# Extract one GWAS's fine-mapping windows from its summary statistics.
#
# Reading a GWAS once and writing only the variants its loci need keeps the
# per-locus jobs from re-reading tens of millions of rows apiece. Windows are
# extracted at margin_factor times the half-window the loci table asks for, so a
# locus can be rerun with a wider window without returning to the full file.

suppressPackageStartupMessages({
    library(data.table)
    library(glue)
    library(yaml)
})

project_directory <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/"
source(file.path(project_directory, "scripts/fine_mapping/fine_mapping.R"))


parse_args <- function(arguments = commandArgs(trailingOnly = TRUE)) {
    if (length(arguments) %% 2 != 0) stop("Arguments must be --name value pairs.")
    names_and_values <- split(arguments, rep(c("name", "value"), length(arguments) / 2))
    values <- as.list(names_and_values$value)
    names(values) <- sub("^--", "", names_and_values$name)
    values
}


strip_chr <- function(chromosomes) sub("^chr", "", as.character(chromosomes), ignore.case = TRUE)


main <- function() {
    args <- parse_args()
    margin_factor <- as.numeric(args$margin_factor)

    loci <- fread(args$loci)
    loci <- loci[config_key == args$config_key]
    if (nrow(loci) == 0) stop(glue("No loci for {args$config_key} in {args$loci}."))
    message(glue("{args$config_key}: {nrow(loci)} loci"))

    config <- load_gwas_config(args$config_key, metadata_yaml = args$metadata, quiet = TRUE)
    sumstats <- load_sumstats(config)
    message(glue("Loaded {nrow(sumstats)} variants from {config$input}"))

    # one interval per locus, widened by margin_factor
    half_window <- loci$window_size * 1000 / 2 * margin_factor
    windows <- data.table(
        chr = strip_chr(loci$chr),
        start = loci$pos - half_window,
        end = loci$pos + half_window
    )

    sumstats[, chr_key := strip_chr(chr)]
    setkey(windows, chr, start, end)
    sumstats[, `:=`(start = pos, end = pos)]
    in_window <- foverlaps(sumstats, windows, by.x = c("chr_key", "start", "end"), nomatch = NULL)

    extracted <- unique(in_window, by = "rsID")
    extracted <- extracted[, .(rsID, chr, pos, gwas_ea, gwas_oa, beta, se)]
    message(glue("Extracted {nrow(extracted)} variants across {nrow(loci)} windows"))

    dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
    fwrite(extracted, args$output, sep = "\t", na = "NA")
    message(glue("Wrote {args$output}"))
}

main()
