#!/usr/bin/env Rscript
# =============================================================================
# Selection Stability of High-Dimensional Mediation (HIMA)
# =============================================================================
#
# Regenerates the simulated dataset of 02_high_dimensional_mediation.R
# (n = 300, p = 200, 5 active mediators, seed 6489) and refits HIMA on B = 200
# random subsamples of half the observations, the subsample size of stability
# selection (Meinshausen and Buhlmann, 2010). For each penalty (MCP as in the
# code example, the lasso, and the de-biased lasso, which is the default of
# hima() in HIMA 2.3.4) it records how often each mediator is selected.
#
# Outputs in ../results/stability/:
#   selection-frequencies.csv   one row per penalty and mediator
#   stability.log               summary per penalty
#
# Run from code/ (on Linux: OPENBLAS_NUM_THREADS=1 nice Rscript ...).
# =============================================================================

source("parallel_setup.R")
suppressPackageStartupMessages(loadNamespace("HIMA"))

n_subsamples <- 200L
penalties <- c("MCP", "lasso", "DBlasso")

# Same draws, in the same order, as 02_high_dimensional_mediation.R.
make_data <- function() {
    set.seed(6489)
    n <- 300
    p <- 200
    c1 <- rnorm(n)
    c2 <- rnorm(n)
    a <- rbinom(n, 1, plogis(0.3 * c1 + 0.3 * c2))
    alpha <- `[<-`(rep(0, p), 1:5, c(0.8, 0.7, 0.6, 0.5, 0.5))
    beta <- `[<-`(rep(0, p), 1:5, c(0.6, 0.5, 0.5, 0.4, 0.4))
    m <- a %*% t(alpha) + outer(c1, rep(0.2, p)) + matrix(rnorm(n * p), n, p)
    y <- as.numeric(0.5 * a + m %*% beta + 0.3 * c1 + 0.3 * c2 + rnorm(n))
    list(pheno = data.frame(Y = y, A = a, C1 = c1, C2 = c2),
         M = `colnames<-`(m, paste0("M", seq_len(p))))
}

# IDs of the mediators HIMA selects on one subsample. hima() returns a list
# whose ID element is empty when no mediator is selected.
selected <- function(dat, rows, penalty) {
    fit <- HIMA::hima(Y ~ A + C1 + C2, data.pheno = dat[["pheno"]][rows, ],
                      data.M = dat[["M"]][rows, ], penalty = penalty, scale = TRUE)
    as.character(fit[["ID"]])
}

main <- function() {
    out_dir <- file.path("..", "results", "stability")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    dat <- make_data()
    n <- nrow(dat[["pheno"]])
    ids <- colnames(dat[["M"]])
    seeds <- par_seeds(6489L, n_subsamples)

    # PSOCK workers (Windows) need the selection function exported, and the
    # fork backend ignores the list.
    cfg <- setup_parallel(n_cores = 2L, export = "selected")
    on.exit(cleanup_parallel(cfg), add = TRUE)

    per_penalty <- Map(function(penalty) {
        sel <- par_map(cfg, function(s) {
            set.seed(s)
            selected(dat, sort(sample.int(n, n %/% 2L)), penalty)
        }, seeds) |> par_values()
        counts <- table(factor(unlist(sel), levels = ids))
        data.frame(penalty = penalty, mediator = ids, active = ids %in% paste0("M", 1:5),
                   frequency = as.numeric(counts) / n_subsamples,
                   mean_selected = mean(lengths(sel)))
    }, penalties)
    freq <- Reduce(rbind, per_penalty)
    write.csv(freq, file.path(out_dir, "selection-frequencies.csv"), row.names = FALSE)

    # Percentages come from the integer counts and are rounded half up, since with 200 subsamples
    # many fall on a half and binary floating point would round them inconsistently.
    summary_lines <- Map(function(d) {
        active <- d[d[["active"]], ]
        noise <- d[!d[["active"]], ]
        counts <- round(active[["frequency"]] * n_subsamples)
        pct <- (200 * counts + n_subsamples) %/% (2 * n_subsamples)
        sprintf(paste("%-8s M1-M5 selected (%% of subsamples): %s | noise mediators selected at least",
                      "once: %d, max noise frequency %.3f | mean selected per subsample %.2f"),
                d[["penalty"]][1L], paste(sprintf("%d", pct), collapse = " "),
                sum(noise[["frequency"]] > 0), max(noise[["frequency"]]), d[["mean_selected"]][1L])
    }, per_penalty)
    lines <- c(sprintf("HIMA %s, %d subsamples of n/2 = %d", as.character(packageVersion("HIMA")),
                       n_subsamples, n %/% 2L), unlist(summary_lines))
    writeLines(lines, file.path(out_dir, "stability.log"))
    writeLines(lines)
}

main()
