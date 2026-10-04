#!/usr/bin/env Rscript
# =============================================================================
# Large-Sample Check of the medflex Specification for a Continuous Exposure
# =============================================================================
#
# 07_simulation.R fits medflex with a continuous exposure through a
# conditional natural effect model, because marginal effects from exposure
# weights (xFit) were biased in large samples. This script records that check.
# For the four designs with a continuous exposure and normal errors (outcome
# continuous or binary, interaction absent or present) it draws n = 20,000
# observations from the data-generating model of 07 and fits both approaches
# (weighting and imputation) under two specifications:
#   marginal      Y ~ A0 * A1 with xFit = glm(A ~ C1 + C2), the documented
#                 route to marginal effects
#   conditional   the specification of 07, with exposure-by-covariate terms
#                 and the covariates centred at C = (0, 0.5)
# The marginal estimates are compared with the marginal true values and the
# conditional ones with the true values at C = (0, 0.5), which equal the
# marginal ones for a continuous outcome. Only point estimates are needed, so
# the natural effect models use robust standard errors instead of the
# bootstrap. With a binary outcome the robust standard errors of marginal
# weighting can be singular, which stops neModel(). Such a fit is counted as
# failed and left out of the means, and its estimates, recovered from the same
# model with a two-replicate bootstrap, are recorded with the error message.
#
# Outputs in ../results/medflex-check/:
#   medflex-check.csv          mean estimate, Monte Carlo SE, true value, bias,
#                              and number of failed fits per design, approach,
#                              specification, and effect
#   medflex-check-failed.csv   replication, estimate, and error message of each
#                              failed fit
#
# Run from code/ (on Linux: OPENBLAS_NUM_THREADS=1 Rscript ... > log).
# =============================================================================

source("07_simulation.R")

n_check <- 20000L
n_reps_check <- 20L

# fit_nemodel() of 07 with robust standard errors (se = "robust") or a
# two-replicate bootstrap (se = "bootstrap") instead of 1,000 resamples.
# medflex re-evaluates the stored call, so, as in 07, the formula and the family
# are spliced into it, and the expanded data and the exposure model are held by
# name in the environment in which it runs.
fit_nemodel_se <- function(formula, family, exp_data, xfit, se) {
    env <- list2env(list(expData = exp_data, xFit = xfit), parent = globalenv())
    f <- as.call(as.list(formula))
    call <- if (is.null(xfit)) {
        bquote(medflex::neModel(.(f), family = .(family), expData = expData, se = .(se), nBoot = 2L))
    } else {
        bquote(medflex::neModel(.(f), family = .(family), expData = expData, xFit = xFit,
                                se = .(se), nBoot = 2L))
    }
    eval(call, env)
}

# The robust fit, or, when its standard errors are singular, the bootstrap fit
# of the same expanded data with the error message. The bootstrap runs on a
# copy of the random-number state, so the later fits of the replication are
# those of a run in which every robust fit succeeds.
fit_nemodel_check <- function(formula, family, exp_data, xfit) {
    tryCatch(list(fit = fit_nemodel_se(formula, family, exp_data, xfit, "robust"), message = NA_character_),
             error = function(e) {
                 state <- get(".Random.seed", envir = globalenv())
                 on.exit(assign(".Random.seed", state, envir = globalenv()))
                 list(fit = fit_nemodel_se(formula, family, exp_data, xfit, "bootstrap"),
                      message = gsub("[[:space:]]+", " ", conditionMessage(e)))
             })
}

# TNIE and PNDE of one specification fitted to one dataset.
check_fit <- function(d, scn, approach, spec) {
    fam <- outcome_family(scn)
    cc <- prm[["c_cond"]]
    conditional <- spec == "conditional"
    dd <- if (conditional) {
        shifted <- `[[<-`(d, "C1", value = d[["C1"]] - cc[["C1"]])
        `[[<-`(shifted, "C2", value = d[["C2"]] - cc[["C2"]])
    } else {
        d
    }
    imp_formula <- if (conditional) Y ~ A * M + A:C1 + A:C2 + C1 + C2 else Y ~ A * M + C1 + C2
    nem_formula <- if (conditional) Y ~ A0 * A1 + A0:C1 + A0:C2 + C1 + C2 else Y ~ A0 * A1
    exp_data <- if (approach == "weight") {
        medflex::neWeight(fit_glm(M ~ A + C1 + C2, gaussian(), dd), data = dd)
    } else {
        medflex::neImpute(fit_glm(imp_formula, fam, dd), data = dd)
    }
    xfit <- if (conditional) NULL else fit_glm(A ~ C1 + C2, gaussian(), dd)
    nem <- fit_nemodel_check(nem_formula, fam, exp_data, xfit)
    eff <- medflex::neEffdecomp(nem[["fit"]], xRef = c(0, 1))
    list(est = unname(coef(eff)[c("total indirect effect", "pure direct effect")]),
         message = nem[["message"]])
}

# Every approach and specification on one dataset (replication r) of one
# design.
check_rep <- function(scn, r) {
    d <- gen_data(scn)
    grid <- expand.grid(approach = c("weight", "impute"), spec = c("marginal", "conditional"),
                        stringsAsFactors = FALSE)
    Map(function(approach, spec) {
        fit <- check_fit(d, scn, approach, spec)
        data.frame(scenario = scn[["id"]], rep = r, outcome = scn[["outcome"]],
                   interaction = scn[["interaction"]], approach = approach, spec = spec,
                   effect = c("tnie", "pnde"), est = fit[["est"]], message = fit[["message"]])
    }, grid[["approach"]], grid[["spec"]]) |>
        Reduce(f = rbind)
}

# True value of a cell from the true values `tv` of its design: marginal for
# the marginal specification, at C = (0, 0.5) for the conditional one (the
# same for a continuous outcome).
true_value <- function(tv, outcome, spec, effect) {
    scale <- if (outcome == "continuous") "md" else if (spec == "marginal") "logor" else "clogor"
    tv[tv[["scale"]] == scale & tv[["effect"]] == effect, "truth"]
}

main <- function() {
    out_dir <- file.path("..", "results", "medflex-check")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    grid <- make_scenarios()
    keep <- grid[["exposure"]] == "continuous" & grid[["errors"]] == "normal" & grid[["n"]] == 1000L
    designs <- Map(function(i) `[[<-`(as.list(grid[i, ]), "n", n_check), which(keep))
    tasks <- rep(designs, each = n_reps_check)
    task_reps <- rep(seq_len(n_reps_check), times = length(designs))
    seeds <- par_seeds(master_seed, length(tasks))

    cfg <- setup_parallel(export = ls(globalenv()))
    on.exit(cleanup_parallel(cfg), add = TRUE)
    reps <- par_map(cfg, function(scn, r, s) {
        set.seed(s)
        check_rep(scn, r)
    }, tasks, task_reps, seeds) |>
        par_values() |>
        Reduce(f = rbind)

    truths <- `names<-`(Map(truth_values, designs), mapply(function(scn) scn[["id"]], designs))
    cells <- split(reps, list(reps[["scenario"]], reps[["approach"]], reps[["spec"]], reps[["effect"]]),
                   drop = TRUE)
    out <- Map(function(d) {
        truth <- true_value(truths[[d[["scenario"]][1L]]], d[["outcome"]][1L], d[["spec"]][1L],
                            d[["effect"]][1L])
        est <- d[["est"]][is.na(d[["message"]])]
        data.frame(d[1L, c("scenario", "outcome", "interaction", "approach", "spec", "effect")],
                   reps = length(est), failed = nrow(d) - length(est), est = mean(est),
                   mcse = sd(est) / sqrt(length(est)), truth = truth, bias = mean(est) - truth)
    }, cells) |>
        Reduce(f = rbind)
    out <- out[order(out[["outcome"]], out[["interaction"]], out[["approach"]], out[["spec"]],
                     out[["effect"]]), ]
    write.csv(out, file.path(out_dir, "medflex-check.csv"), row.names = FALSE)
    failed <- reps[!is.na(reps[["message"]]), ]
    failed <- failed[order(failed[["scenario"]], failed[["rep"]], failed[["approach"]], failed[["spec"]],
                           failed[["effect"]]), ]
    write.csv(failed, file.path(out_dir, "medflex-check-failed.csv"), row.names = FALSE)
    cat(sprintf("medflex %s, n = %d, %d replications per design\n",
                as.character(packageVersion("medflex")), n_check, n_reps_check))
    print(format(out, digits = 3), row.names = FALSE)
    cat("\nFailed fits, left out of the means above:\n")
    print(format(failed, digits = 3), row.names = FALSE, right = FALSE)
}

timing <- system.time(main())
cat(sprintf("Elapsed: %.1f s\n", timing[["elapsed"]]))
