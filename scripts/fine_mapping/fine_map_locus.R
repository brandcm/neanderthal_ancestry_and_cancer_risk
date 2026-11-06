# Fine-map one locus.
#
# Takes the windows its GWAS's variants were extracted into rather than the full
# summary statistics, and writes one result object. The credible set tables are
# assembled afterwards by summarize_fine_mapping_results.R, so nothing here writes
# to a shared file and the loci can run in parallel.

suppressPackageStartupMessages({
    library(data.table)
    library(glue)
    library(susieR)
    library(yaml)
})

project_directory <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/"
source(file.path(project_directory, "scripts/fine_mapping/fine_mapping.R"))

options(warn = 1)


parse_args <- function(arguments = commandArgs(trailingOnly = TRUE)) {
    if (length(arguments) %% 2 != 0) stop("Arguments must be --name value pairs.")
    names_and_values <- split(arguments, rep(c("name", "value"), length(arguments) / 2))
    values <- as.list(names_and_values$value)
    names(values) <- sub("^--", "", names_and_values$name)
    values
}


main <- function() {
    args <- parse_args()

    windows <- fread(args$windows, colClasses = list(character = c("rsID", "gwas_ea", "gwas_oa")))
    introgressed <- load_introgressed_variants()

    result <- fine_mapping(
        introgressed,
        config_key = args$config_key,
        ancestry_group = args$ancestry_group,
        locus_id = args$locus_id,
        n_samples = as.integer(args$n_samples),
        lead_rsid = args$lead_rsid,
        window_size = as.integer(args$window_size),
        L = as.integer(args$L),
        sumstats = windows,
        save_results = FALSE,
        main_file = NULL
    )

    # Neither is needed downstream and both are large: the LD matrix can be
    # rebuilt from the panel, and kriging's ggplot carries its own copy of the
    # locus data plus enclosing environments (>2 GB at the largest loci).
    # Only kriging$conditional_dist is read, by the logLR diagnostic.
    result$ld <- NULL
    result$kriging$plot <- NULL

    dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, args$output)
    message(glue("Wrote {args$output}"))
}

main()
