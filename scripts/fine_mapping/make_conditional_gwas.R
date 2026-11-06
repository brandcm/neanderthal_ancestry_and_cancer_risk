#!/usr/bin/env Rscript
# Turn a conditional analysis result into a GWAS-shaped file the locuszoom
# function can read.
#
#   ./run_rscript.sh make_conditional_gwas.R \
#       --input  ../../data/fine_mapping/conditional/<file>.txt \
#       --chr    12 \
#       --output ../../data/fine_mapping/conditional/pthlh_conditioned.txt.gz
#
# generate_locuszoom_plot_with_annotations() expects a gzipped table with
# chromosome, position, rsID and a p-value, and plots -log10(p). The conditional
# analysis produces z-scores, so the p-value is recovered as a two-sided normal
# tail: p = 2 * pnorm(-|z|).
#
# Variants whose conditional z is NA are dropped. These are variants in near
# perfect LD with the conditioned-on variant, where the conditional statistic is
# undefined rather than zero -- plotting them at p = 1 would wrongly suggest the
# signal was tested and found absent.

suppressPackageStartupMessages({
    library(data.table)
    library(glue)
})

parse_args <- function(arguments = commandArgs(trailingOnly = TRUE)) {
    if (length(arguments) %% 2 != 0) stop("Arguments must be --name value pairs.")
    names_and_values <- split(arguments, rep(c("name", "value"), length(arguments) / 2))
    values <- as.list(names_and_values$value)
    names(values) <- sub("^--", "", names_and_values$name)
    values
}

args <- parse_args()
for (required in c("input", "output")) {
    if (is.null(args[[required]])) stop(glue("--{required} is required."))
}

conditional <- fread(args$input, na.strings = c("NA", ""))
message(glue("Read {nrow(conditional)} rows from {basename(args$input)}."))
message(glue("Columns: {paste(names(conditional), collapse = ', ')}"))

`%||%` <- function(a, b) if (is.null(a)) b else a

pick <- function(candidates, label) {
    found <- intersect(candidates, names(conditional))
    if (length(found) == 0) {
        stop(glue("No {label} column found. Looked for: {paste(candidates, collapse = ', ')}. ",
                  "Pass the right name with --{label}_column."))
    }
    found[1]
}

rsid_column <- args$rsid_column %||% pick(c("rsID", "SNP", "rsid"), "rsid")
pos_column  <- args$pos_column  %||% pick(c("pos", "BP", "position"), "pos")
z_column    <- args$z_column    %||% pick(c("z_conditional", "Z_conditional", "z_cond"), "z")

message(glue("Using rsID = {rsid_column}, position = {pos_column}, z = {z_column}."))

output <- conditional[, .(
    SNP = get(rsid_column),
    BP  = as.integer(get(pos_column)),
    Z   = as.numeric(get(z_column))
)]

n_undefined <- sum(is.na(output$Z))
output <- output[!is.na(Z)]
if (n_undefined > 0) {
    message(glue("Dropped {n_undefined} variants whose conditional z is undefined ",
                 "(near-perfect LD with the conditioned-on variant)."))
}

chromosome <- args$chr
if (is.null(chromosome)) {
    if ("chr" %in% names(conditional)) {
        chromosome <- as.character(unique(conditional$chr))[1]
    } else {
        stop("No chr column found; pass --chr.")
    }
}
output[, CHR := chromosome]

# two-sided normal tail; log.p keeps very small p-values from underflowing to 0
output[, P := 2 * pnorm(-abs(Z))]
underflow <- output$P == 0
if (any(underflow)) {
    output[underflow, P := .Machine$double.xmin]
    message(glue("{sum(underflow)} p-values underflowed to zero and were floored; ",
                 "max |z| is {round(max(abs(output$Z)), 1)}."))
}

setcolorder(output, c("CHR", "BP", "SNP", "P", "Z"))
setorder(output, BP)

# quote = FALSE matters: fwrite quotes character columns by default, and the
# plotting function filters the file with awk and joins rsIDs against a .bim,
# neither of which expects quotes
fwrite(output, args$output, sep = "\t", na = "NA", quote = FALSE, compress = "gzip")
message(glue("Wrote {nrow(output)} variants to {args$output}"))
message(glue("Strongest remaining signal: |z| = {round(max(abs(output$Z)), 2)} ",
             "(p = {signif(min(output$P), 3)})"))
