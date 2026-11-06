#!/usr/bin/env Rscript
# Refit the loci where susie_rss did not converge.
#
# These loci share a signature: L = 10, many credible sets, a large fraction of
# them single-variant, and no severe kriging flags. That is not LD mismatch and
# it is not weak signal -- it is one real signal being split across several
# components while the model oscillates between assignments. The remedies are
# refine = TRUE (restart from several initialisations and keep the best ELBO)
# and a smaller L.
#
# One locus and L per invocation:
#
#   ./run_rscript.sh refit_nonconverged.R --config_key <key> \
#       --ancestry_group <anc> --locus_id <locus> --L <L>
#
# With --list, prints the task table and exits without fitting. With
# --list_output <file>, writes the same table as TSV, which is how the Snakemake
# checkpoint discovers which refits to run. --task_id still selects a task by
# row number, for running one by hand from the printed table.
#
# Results go to data/fine_mapping/refits/<config_key>/<locus>_<anc>_L<L>.rds and
# never overwrite the originals.

suppressPackageStartupMessages({
    library(data.table)
    library(glue)
    library(susieR)
    library(yaml)
})

project_directory <- "/wynton/group/capra/projects/neanderthal_ancestry_and_cancer_risk/"
source(file.path(project_directory, "scripts/fine_mapping/fine_mapping.R"))

options(warn = 1)

fine_mapping_directory <- file.path(project_directory, "data/fine_mapping")
refit_directory <- file.path(fine_mapping_directory, "refits")

# L values to try, in the order tasks should run. The smaller fit is both more
# likely to converge and cheaper, so it comes first.
L_VALUES <- c(5, 10)


parse_args <- function(arguments = commandArgs(trailingOnly = TRUE)) {
    values <- list()
    index <- 1
    while (index <= length(arguments)) {
        name <- sub("^--", "", arguments[index])
        if (name == "list") {
            values$list <- TRUE
            index <- index + 1
        } else {
            values[[name]] <- arguments[index + 1]
            index <- index + 2
        }
    }
    values
}


build_tasks <- function() {
    summary_table <- fread(file.path(fine_mapping_directory, "locus_summary.txt"),
                           na.strings = c("NA", ""))
    loci <- fread(file.path(fine_mapping_directory, "fine_mapping_loci.txt"),
                  na.strings = c("NA", ""))

    unconverged <- summary_table[converged == FALSE,
                                 .(config_key = cancer, ancestry_group, locus_id,
                                   original_L = L, n_variants)]

    tasks <- merge(
        unconverged,
        loci[, .(config_key, ancestry_group, locus_id, lead_rsid,
                 n_samples, window_size)],
        by = c("config_key", "ancestry_group", "locus_id")
    )
    if (nrow(tasks) != nrow(unconverged)) {
        warning(glue("{nrow(unconverged)} unconverged loci but {nrow(tasks)} matched ",
                     "rows in fine_mapping_loci.txt; check for duplicate locus ids."))
    }

    # the window size actually used, which may differ from the loci file when an
    # override applied
    overrides_file <- file.path(project_directory, "scripts/fine_mapping/locus_overrides.txt")
    if (file.exists(overrides_file)) {
        overrides <- fread(overrides_file, na.strings = c("NA", ""))
        if (all(c("config_key", "locus_id", "window_size") %in% names(overrides))) {
            tasks <- merge(tasks, overrides[, .(config_key, locus_id,
                                                override_window = window_size)],
                           by = c("config_key", "locus_id"), all.x = TRUE)
            tasks[!is.na(override_window), window_size := override_window]
            tasks[, override_window := NULL]
        }
    }

    # cheapest first, so a truncated run still finishes the easy ones
    setorder(tasks, n_variants)
    tasks <- tasks[rep(seq_len(.N), each = length(L_VALUES))]
    tasks[, L := rep(L_VALUES, times = .N / length(L_VALUES))]
    tasks[, task_id := .I]
    tasks[]
}


# A task is addressed by the locus it refits rather than by its row number, so
# that the mapping does not shift when build_tasks() is re-evaluated against a
# results directory that has changed in the meantime. --task_id is kept for
# running one task by hand from the printed table.
select_task <- function(tasks, args) {
    if (!is.null(args$config_key)) {
        for (required in c("ancestry_group", "locus_id", "L")) {
            if (is.null(args[[required]])) {
                stop(glue("--config_key given without --{required}."))
            }
        }
        task <- tasks[config_key == args$config_key &
                          ancestry_group == args$ancestry_group &
                          locus_id == args$locus_id &
                          L == as.integer(args$L)]
        if (nrow(task) != 1) {
            stop(glue("Expected one task for {args$config_key} {args$locus_id} ",
                      "{args$ancestry_group} L={args$L}; found {nrow(task)}."))
        }
        return(task)
    }

    task_id <- as.integer(args$task_id)
    if (is.na(task_id) || task_id < 1 || task_id > nrow(tasks)) {
        stop(glue("Pass --config_key, --ancestry_group, --locus_id and --L, or ",
                  "--task_id between 1 and {nrow(tasks)}."))
    }
    tasks[task_id]
}


main <- function() {
    args <- parse_args()
    tasks <- build_tasks()

    if (!is.null(args$list_output)) {
        fwrite(tasks, args$list_output, sep = "\t", na = "NA")
        message(glue("Wrote {nrow(tasks)} tasks to {args$list_output}"))
        return(invisible(NULL))
    }

    if (isTRUE(args$list)) {
        print(tasks[, .(task_id, config_key = substr(config_key, 1, 30),
                        ancestry_group, locus_id, n_variants,
                        window_size, original_L, L)], nrows = 100)
        message(glue("\n{nrow(tasks)} tasks."))
        return(invisible(NULL))
    }

    task <- select_task(tasks, args)

    message(glue("{task$config_key} {task$ancestry_group} ",
                 "{task$locus_id} at L = {task$L}, refine = TRUE ",
                 "({task$n_variants} variants, window {task$window_size} kb)"))

    windows_file <- file.path(fine_mapping_directory, "windows",
                              glue("{task$config_key}_windows.txt.gz"))
    windows <- fread(windows_file,
                     colClasses = list(character = c("rsID", "gwas_ea", "gwas_oa")))
    introgressed <- load_introgressed_variants()

    result <- fine_mapping(
        introgressed,
        config_key      = task$config_key,
        ancestry_group  = task$ancestry_group,
        locus_id        = task$locus_id,
        n_samples       = as.integer(task$n_samples),
        lead_rsid       = task$lead_rsid,
        window_size     = as.integer(task$window_size),
        L               = as.integer(task$L),
        refine          = TRUE,
        sumstats        = windows,
        save_results    = FALSE,
        main_file       = NULL
    )

    result$ld <- NULL
    result$kriging$plot <- NULL

    output_directory <- file.path(refit_directory, task$config_key)
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
    output_file <- file.path(
        output_directory,
        glue("{task$locus_id}_{task$ancestry_group}_L{task$L}.rds")
    )
    saveRDS(result, output_file)
    message(glue("Wrote {output_file}"))
    message(glue("converged: {isTRUE(result$converged)}; credible sets: ",
                 "{if (is.null(result$cs_summary)) 0 else nrow(result$cs_summary)}"))
}

main()
