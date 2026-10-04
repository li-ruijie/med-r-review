#!/usr/bin/env Rscript
# =============================================================================
# Robust Option of bmem and bmemLavaan with an Exposure-Mediator Product Term
# =============================================================================
#
# bmem has no interaction syntax, and bmemLavaan fits lavaan models from sample
# moments, so the colon operator, which lavaan 0.7-2 fits for observed
# variables, fails there. The script therefore supplies the product of the
# exposure and the mediator as one more observed variable, which all three
# packages fit, and the natural indirect effects are then defined from the path
# labels. Maximum likelihood for this model equals least squares for each
# equation. bmem and bmemLavaan also offer a robust option, which downweights
# observations by their distance in all model variables (Yuan and Zhang, 2012),
# the product included. This script checks that option against maximum
# likelihood (lavaan). With a binary exposure and the product term, the robust
# option of bmemLavaan stops with a singular-matrix error, which the output
# records.
#
# Data: n = 1,000, exposure X binary (p = 0.5) or standard normal,
#   M = 0.5 X + e_M,  Y = 0.4 X + 0.3 M + g X M + e_Y,  e_M, e_Y ~ N(0, 1),
# with g = 0.2 (the model includes the product XM) or g = 0 (it does not),
# and with or without outliers (Y + 10 for 5% of the observations). True
# values: a = 0.5, PNIE = a b = 0.15, and TNIE = a (b + g) = 0.15 + 0.5 g for
# X from 0 to 1. 200 replications per design. Only point estimates are
# needed, so the bootstrap of bmem and bmemLavaan runs two replicates.
#
# Outputs in ../results/bmem-check/:
#   bmem-check.csv   mean estimate, Monte Carlo SE, true value, bias, number of
#                    failed fits, and the first error message of those fits per
#                    design, estimator, and parameter
#
# Run from code/ (on Linux: OPENBLAS_NUM_THREADS=1 Rscript ... > log).
# =============================================================================

source("parallel_setup.R")
invisible(suppressPackageStartupMessages(Map(loadNamespace, c("lavaan", "bmemLavaan", "bmem", "sem"))))

n_obs <- 1000L
n_reps <- 200L
params <- c("a", "b", "g", "pnie", "tnie")

designs <- expand.grid(exposure = c("continuous", "binary"), product = c(TRUE, FALSE),
                       outliers = c(0, 0.05), stringsAsFactors = FALSE)

gen <- function(design) {
    x <- if (design[["exposure"]] == "binary") rbinom(n_obs, 1, 0.5) else rnorm(n_obs)
    g <- if (design[["product"]]) 0.2 else 0
    m <- 0.5 * x + rnorm(n_obs)
    y <- 0.4 * x + 0.3 * m + g * x * m + rnorm(n_obs)
    shifted <- sample.int(n_obs, round(design[["outliers"]] * n_obs))
    y <- `[<-`(y, shifted, y[shifted] + 10)
    data.frame(x = x, m = m, y = y, xm = x * m)
}

truth <- function(design) {
    g <- if (design[["product"]]) 0.2 else 0
    c(a = 0.5, b = 0.3, g = g, pnie = 0.15, tnie = 0.15 + 0.5 * g)
}

# lavaan syntax, used by lavaan and bmemLavaan.
lavaan_model <- function(product) {
    if (product) {
        "m ~ a*x\n y ~ c*x + b*m + g*xm\n pnie := a*b\n tnie := a*(b + g)"
    } else {
        "m ~ a*x\n y ~ c*x + b*m\n pnie := a*b\n tnie := a*b"
    }
}

# RAM paths for bmem (sem package), the same model with the exogenous
# variances and covariance free.
ram_model <- function(product) {
    paths <- c("x -> m, a, NA", "x -> y, c, NA", "m -> y, b, NA", "m <-> m, vm, NA", "y <-> y, vy, NA",
               "x <-> x, vx, NA")
    extra <- c("xm -> y, g, NA", "xm <-> xm, vxm, NA", "x <-> xm, cxxm, NA")
    sem::specifyModel(text = paste(c(paths, if (product) extra), collapse = "\n"), quiet = TRUE)
}

# Point estimates of the five parameters (g = 0 without the product term).
lavaan_fit <- function(d, product) {
    pe <- lavaan::parameterEstimates(lavaan::sem(lavaan_model(product), data = d))
    est <- pe[["est"]][match(params, pe[["label"]])]
    `[<-`(est, is.na(est), 0)
}

# The robust weights depend on every column of the data passed, so the
# product column is passed only when the model contains it.
model_vars <- function(product) if (product) c("x", "m", "y", "xm") else c("x", "m", "y")

bmem_lavaan_fit <- function(d, product) {
    invisible(capture.output(fit <- bmemLavaan::bmem(d[, model_vars(product)], lavaan_model(product),
                                                     method = "list", boot = 2L, robust = TRUE)))
    est <- fit[["ci"]][, "estimate"][params]
    `[<-`(unname(est), is.na(est), 0)
}

# bmem applies the robust option only with two-stage maximum likelihood
# (method = "tsml", its default), which equals maximum likelihood on complete
# data.
bmem_fit <- function(d, product) {
    indirect <- if (product) c("a*b", "a*b+a*g") else c("a*b", "a*b")
    invisible(capture.output(fit <- bmem::bmem(d[, model_vars(product)], ram_model(product), indirect,
                                               method = "tsml", boot = 2L, robust = TRUE)))
    ci <- fit[["ci"]][, "estimate"]
    unname(c(ci[["a"]], ci[["b"]], if (product) ci[["g"]] else 0, ci[[indirect[[1L]]]],
             ci[[indirect[[2L]]]]))
}

estimators <- list(lavaan = lavaan_fit, bmemLavaan = bmem_lavaan_fit, bmem = bmem_fit)

# Every estimator on one dataset. A fit that fails is recorded as NA, with
# its error message.
check_rep <- function(design) {
    d <- gen(design)
    Map(function(name, f) {
        fit <- tryCatch(list(est = f(d, design[["product"]]), message = NA_character_),
                        error = function(e) list(est = rep(NA_real_, length(params)),
                                                 message = gsub("[[:space:]]+", " ", conditionMessage(e))))
        data.frame(design[rep(1L, length(params)), ], estimator = name, parameter = params,
                   est = fit[["est"]], message = fit[["message"]], row.names = NULL)
    }, names(estimators), estimators) |>
        Reduce(f = rbind)
}

main <- function() {
    out_dir <- file.path("..", "results", "bmem-check")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    tasks <- rep(seq_len(nrow(designs)), each = n_reps)
    seeds <- par_seeds(6489L, length(tasks))

    cfg <- setup_parallel(export = ls(globalenv()))
    on.exit(cleanup_parallel(cfg), add = TRUE)
    reps <- par_map(cfg, function(i, s) {
        set.seed(s)
        check_rep(designs[i, ])
    }, tasks, seeds) |>
        par_values() |>
        Reduce(f = rbind)

    cells <- split(reps, list(reps[["exposure"]], reps[["product"]], reps[["outliers"]],
                              reps[["estimator"]], reps[["parameter"]]), drop = TRUE)
    out <- Map(function(d) {
        tv <- truth(d[1L, ])[[d[["parameter"]][1L]]]
        est <- d[["est"]][!is.na(d[["est"]])]
        mean_est <- if (length(est) > 0L) mean(est) else NA_real_
        msg <- d[["message"]][!is.na(d[["message"]])]
        data.frame(d[1L, c("exposure", "product", "outliers", "estimator", "parameter")],
                   reps = length(est), failed = nrow(d) - length(est), est = mean_est,
                   mcse = if (length(est) > 1L) sd(est) / sqrt(length(est)) else NA_real_, truth = tv,
                   bias = mean_est - tv, message = if (length(msg) > 0L) msg[[1L]] else NA_character_,
                   row.names = NULL)
    }, cells) |>
        Reduce(f = rbind)
    out <- out[order(out[["exposure"]], -out[["product"]], out[["outliers"]], out[["estimator"]],
                     match(out[["parameter"]], params)), ]
    write.csv(out, file.path(out_dir, "bmem-check.csv"), row.names = FALSE)
    cat(sprintf("lavaan %s, bmemLavaan %s, bmem %s, n = %d, %d replications per design\n",
                as.character(packageVersion("lavaan")), as.character(packageVersion("bmemLavaan")),
                as.character(packageVersion("bmem")), n_obs, n_reps))
    print(format(out[, names(out) != "message"], digits = 3), row.names = FALSE)
    failures <- unique(out[!is.na(out[["message"]]), c("exposure", "product", "outliers", "estimator",
                                                        "message")])
    cat("\nFirst error message of each design and estimator with failed fits:\n")
    print(failures, row.names = FALSE, right = FALSE)
}

timing <- system.time(main())
cat(sprintf("Elapsed: %.1f s\n", timing[["elapsed"]]))
