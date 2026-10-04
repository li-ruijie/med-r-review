#!/usr/bin/env Rscript
# =============================================================================
# Simulation: Natural Effect Estimators Across Mediation Packages
# =============================================================================
#
# Compares the total natural indirect effect (TNIE) and pure natural direct
# effect (PNDE) estimators of mediation, medflex, CMAverse, regmedint,
# robmed, and lavaan over 48 scenarios: outcome type (continuous, binary) x
# exposure type (binary, continuous) x exposure-mediator interaction (absent,
# present) x error distribution (normal, scaled t with 3 df) x sample size
# (100, 300, 1000). Every resampling method uses 1,000 resamples.
#
# Usage (from code/, launched by run_simulation.sh, which pins BLAS to
# one thread per process):
#   Rscript 07_simulation.R pilot   2 replications per scenario
#   Rscript 07_simulation.R check   32 replications per scenario
#   Rscript 07_simulation.R full    1,000 replications per scenario
#   Rscript 07_simulation.R full medflex-weight
#                                   re-runs the named methods in every
#                                   finished scenario and replaces their rows,
#                                   with the seeds a full run gives them
#
# Outputs in ../results/simulation/<mode>/:
#   truth.rds            true value of every estimand in every scenario
#   prevalence.rds       P(Y = 1) of every binary-outcome scenario
#   scenario-<id>.rds    one row per replication, method, and effect
#   run.log              progress and timing
#
# A scenario whose output file exists is skipped, so an interrupted run
# resumes where it stopped.
# =============================================================================

source("parallel_setup.R")

packages <- c("mediation", "medflex", "CMAverse", "regmedint", "robmed", "lavaan")
invisible(Map(function(p) suppressPackageStartupMessages(loadNamespace(p)), packages))

n_resamples <- 1000L
master_seed <- 6489L

# Data-generating parameters (intercept first, then in the order of the
# variable names).
prm <- list(
    a = c(intercept = -0.2, C1 = 0.4, C2 = 0.4),
    m = c(intercept = 0, A = 0.5, C1 = 0.3, C2 = 0.3),
    y_continuous = c(intercept = 0, A = 0.4, M = 0.4, C1 = 0.3, C2 = 0.3),
    y_binary = c(intercept = -1, A = 0.4, M = 0.4, C1 = 0.3, C2 = 0.3),
    interaction = 0.3,
    c_cond = c(C1 = 0, C2 = 0.5)
)

# ── Scenarios ────────────────────────────────────────────────────────────────

make_scenarios <- function() {
    grid <- expand.grid(
        outcome = c("continuous", "binary"),
        exposure = c("binary", "continuous"),
        interaction = c(FALSE, TRUE),
        errors = c("normal", "t3"),
        n = c(100L, 300L, 1000L),
        stringsAsFactors = FALSE
    )
    `[[<-`(grid, "id", value = sprintf("s%02d", seq_len(nrow(grid))))
}

# ── Data generation ──────────────────────────────────────────────────────────

draw_errors <- function(n, errors) {
    if (errors == "normal") rnorm(n) else rt(n, df = 3) / sqrt(3)
}

error_density <- function(x, errors) {
    if (errors == "normal") dnorm(x) else sqrt(3) * dt(sqrt(3) * x, df = 3)
}

mediator_mean <- function(a, c1, c2) {
    b <- prm[["m"]]
    b[["intercept"]] + b[["A"]] * a + b[["C1"]] * c1 + b[["C2"]] * c2
}

outcome_lp <- function(scn, a, m, c1, c2) {
    b <- prm[[paste0("y_", scn[["outcome"]])]]
    g <- if (scn[["interaction"]]) prm[["interaction"]] else 0
    b[["intercept"]] + b[["A"]] * a + b[["M"]] * m + g * a * m + b[["C1"]] * c1 + b[["C2"]] * c2
}

gen_data <- function(scn) {
    n <- scn[["n"]]
    c1 <- rnorm(n)
    c2 <- rbinom(n, 1L, 0.5)
    lp_a <- prm[["a"]][["intercept"]] + prm[["a"]][["C1"]] * c1 + prm[["a"]][["C2"]] * c2
    a <- if (scn[["exposure"]] == "binary") rbinom(n, 1L, plogis(lp_a)) else lp_a + rnorm(n)
    m <- mediator_mean(a, c1, c2) + draw_errors(n, scn[["errors"]])
    lp_y <- outcome_lp(scn, a, m, c1, c2)
    y <- if (scn[["outcome"]] == "continuous") {
        lp_y + draw_errors(n, scn[["errors"]])
    } else {
        rbinom(n, 1L, plogis(lp_y))
    }
    data.frame(C1 = c1, C2 = c2, A = a, M = m, Y = y)
}

# ── True values ──────────────────────────────────────────────────────────────

# E[Y(a, M(astar)) | C1 = c1, C2 = c2]. The continuous outcome is linear in
# M, so the mediator mean can be plugged in. The binary outcome integrates
# the outcome probability over the mediator error distribution.
nested_mean <- function(scn, a, astar, c1, c2) {
    mu <- mediator_mean(astar, c1, c2)
    if (scn[["outcome"]] == "continuous") return(outcome_lp(scn, a, mu, c1, c2))
    integrand <- function(e) plogis(outcome_lp(scn, a, mu + e, c1, c2)) * error_density(e, scn[["errors"]])
    integrate(integrand, -Inf, Inf, rel.tol = 1e-10)[["value"]]
}

# E[Y(a, M(astar))], averaging over C1 ~ N(0, 1) and C2 ~ Bernoulli(0.5).
marginal_mean <- function(scn, a, astar) {
    over_c1 <- function(c2) {
        f <- function(x) mapply(function(c1) nested_mean(scn, a, astar, c1, c2) * dnorm(c1), x)
        integrate(f, -Inf, Inf, rel.tol = 1e-10)[["value"]]
    }
    mean(c(over_c1(0), over_c1(1)))
}

truth_values <- function(scn) {
    effects <- c("tnie", "pnde")
    p11 <- marginal_mean(scn, 1, 1)
    p10 <- marginal_mean(scn, 1, 0)
    p00 <- marginal_mean(scn, 0, 0)
    rd <- data.frame(effect = effects, truth = c(p11 - p10, p10 - p00))
    if (scn[["outcome"]] == "continuous") return(cbind(rd, scale = "md"))
    cc <- prm[["c_cond"]]
    q <- function(a, astar) nested_mean(scn, a, astar, cc[["C1"]], cc[["C2"]])
    q11 <- q(1, 1)
    q10 <- q(1, 0)
    q00 <- q(0, 0)
    rbind(
        cbind(rd, scale = "rd"),
        data.frame(effect = effects, truth = c(qlogis(p11) - qlogis(p10), qlogis(p10) - qlogis(p00)),
                   scale = "logor"),
        data.frame(effect = effects, truth = c(qlogis(q11) - qlogis(q10), qlogis(q10) - qlogis(q00)),
                   scale = "clogor")
    )
}

# P(Y = 1) of a binary-outcome scenario, averaging over the covariates, the
# exposure given the covariates, and the mediator given both. Under the
# data-generating model E[Y | A = a, C] is the nested mean with astar = a.
outcome_prevalence <- function(scn) {
    given_c <- function(c1, c2) {
        lp_a <- prm[["a"]][["intercept"]] + prm[["a"]][["C1"]] * c1 + prm[["a"]][["C2"]] * c2
        if (scn[["exposure"]] == "binary") {
            return(plogis(lp_a) * nested_mean(scn, 1, 1, c1, c2) +
                       plogis(-lp_a) * nested_mean(scn, 0, 0, c1, c2))
        }
        f <- function(x) mapply(function(a) nested_mean(scn, a, a, c1, c2) * dnorm(a - lp_a), x)
        integrate(f, -Inf, Inf, rel.tol = 1e-8)[["value"]]
    }
    over_c1 <- function(c2) {
        f <- function(x) mapply(function(c1) given_c(c1, c2) * dnorm(c1), x)
        integrate(f, -Inf, Inf, rel.tol = 1e-8)[["value"]]
    }
    mean(c(over_c1(0), over_c1(1)))
}

# ── Estimators ───────────────────────────────────────────────────────────────
#
# Each returns list(est, se, lower, upper), every element ordered (TNIE, PNDE).

outcome_family <- function(scn) {
    if (scn[["outcome"]] == "binary") binomial("logit") else gaussian()
}

# mediation (bootstrap) and medflex re-evaluate the stored model call, so the
# family object is spliced into the call rather than referenced by a local
# name that no longer exists when the call is re-run.
fit_glm <- function(formula, family, data) {
    do.call("glm", list(formula = formula, family = family, data = data))
}

fit_mediation <- function(d, scn, boot) {
    mfit <- lm(M ~ A + C1 + C2, data = d)
    yfit <- fit_glm(Y ~ A * M + C1 + C2, outcome_family(scn), d)
    med <- mediation::mediate(mfit, yfit, treat = "A", mediator = "M", control.value = 0,
                              treat.value = 1, sims = n_resamples, boot = boot)
    list(est = c(med[["d1"]], med[["z0"]]),
         se = c(sd(med[["d1.sims"]]), sd(med[["z0.sims"]])),
         lower = c(med[["d1.ci"]][[1L]], med[["z0.ci"]][[1L]]),
         upper = c(med[["d1.ci"]][[2L]], med[["z0.ci"]][[2L]]))
}

# neModel() is called by name on the expanded data and the exposure model, so
# it runs in an environment that holds them under those names, with the
# family spliced into the call. The formula is spliced as a bare call, so
# medflex creates it in its own frame. Without an exposure model, medflex
# evaluates the weights, attr(expData, "weights"), in the formula's
# environment, and a formula that carried its own expData would make every
# bootstrap refit use the weights of the original sample. Its standard errors
# are the default bootstrap ones for both approaches.
fit_nemodel <- function(formula, family, exp_data, xfit) {
    env <- list2env(list(expData = exp_data, xFit = xfit), parent = globalenv())
    f <- as.call(as.list(formula))
    call <- if (is.null(xfit)) {
        bquote(medflex::neModel(.(f), family = .(family), expData = expData, nBoot = .(n_resamples)))
    } else {
        bquote(medflex::neModel(.(f), family = .(family), expData = expData, xFit = xFit,
                                nBoot = .(n_resamples)))
    }
    eval(call, env)
}

# Binary exposure: coded as a factor, as medflex requires for a categorical
# exposure, and made marginal through exposure weights (xFit) in a saturated
# natural effect model. Continuous exposure: exposure weighting through xFit
# was biased in large samples (12_medflex_check.R), so the conditional natural effect model
# with exposure-by-covariate terms is fitted instead, with the covariates
# centred at c_cond so that its effects are those at c_cond. The imputation
# model then carries the same exposure-by-covariate terms, as medflex asks.
fit_medflex <- function(d, scn, approach) {
    fam <- outcome_family(scn)
    binary_a <- scn[["exposure"]] == "binary"
    cc <- prm[["c_cond"]]
    dd <- if (binary_a) {
        `[[<-`(d, "A", value = factor(d[["A"]], levels = c(0, 1)))
    } else {
        shifted <- `[[<-`(d, "C1", value = d[["C1"]] - cc[["C1"]])
        `[[<-`(shifted, "C2", value = d[["C2"]] - cc[["C2"]])
    }
    imp_formula <- if (binary_a) Y ~ A * M + C1 + C2 else Y ~ A * M + A:C1 + A:C2 + C1 + C2
    nem_formula <- if (binary_a) Y ~ A0 * A1 else Y ~ A0 * A1 + A0:C1 + A0:C2 + C1 + C2
    exp_data <- if (approach == "weight") {
        medflex::neWeight(fit_glm(M ~ A + C1 + C2, gaussian(), dd), data = dd)
    } else {
        medflex::neImpute(fit_glm(imp_formula, fam, dd), data = dd)
    }
    xfit <- if (binary_a) fit_glm(A ~ C1 + C2, binomial("logit"), dd) else NULL
    nem <- fit_nemodel(nem_formula, fam, exp_data, xfit)
    eff <- if (binary_a) medflex::neEffdecomp(nem) else medflex::neEffdecomp(nem, xRef = c(0, 1))
    nm <- c("total indirect effect", "pure direct effect")
    ci <- confint(eff)
    # vcov() of the effects, K cov(B) K', cancels catastrophically when a few
    # bootstrap coefficients are huge, so the standard error comes from the
    # bootstrap replicates of each effect, as the interval of confint() does.
    boot_eff <- eff[["linfct"]][nm, , drop = FALSE] %*% t(nem[["bootRes"]][["t"]])
    list(est = unname(coef(eff)[nm]),
         se = unname(apply(boot_eff, 1L, sd)),
         lower = unname(ci[nm, 1L]),
         upper = unname(ci[nm, 2L]))
}

fit_cmaverse <- function(d, scn) {
    binary <- scn[["outcome"]] == "binary"
    fit <- CMAverse::cmest(data = d, model = "rb", outcome = "Y", exposure = "A", mediator = "M",
                           basec = c("C1", "C2"), EMint = TRUE, mreg = list("linear"),
                           yreg = if (binary) "logistic" else "linear", astar = 0, a = 1,
                           mval = list(0), estimation = "imputation", inference = "bootstrap",
                           nboot = n_resamples)
    # Ratio-scale effects are reported on the log scale, with the standard
    # error carried over by the delta method.
    nm <- if (binary) c("Rtnie", "Rpnde") else c("tnie", "pnde")
    pe <- unname(fit[["effect.pe"]][nm])
    se <- unname(fit[["effect.se"]][nm])
    lo <- unname(fit[["effect.ci.low"]][nm])
    hi <- unname(fit[["effect.ci.high"]][nm])
    if (binary) list(est = log(pe), se = se / pe, lower = log(lo), upper = log(hi))
    else list(est = pe, se = se, lower = lo, upper = hi)
}

fit_regmedint <- function(d, scn) {
    binary <- scn[["outcome"]] == "binary"
    fit <- regmedint::regmedint(data = d, yvar = "Y", avar = "A", mvar = "M", cvar = c("C1", "C2"),
                                a0 = 0, a1 = 1, m_cde = 0, c_cond = unname(prm[["c_cond"]]),
                                mreg = "linear", yreg = if (binary) "logistic" else "linear",
                                interaction = TRUE, casecontrol = FALSE)
    nm <- c("tnie", "pnde")
    ci <- confint(fit)
    list(est = unname(coef(fit)[nm]),
         se = unname(sqrt(diag(vcov(fit)))[nm]),
         lower = unname(ci[nm, 1L]),
         upper = unname(ci[nm, 2L]))
}

fit_robmed <- function(d) {
    tst <- robmed::test_mediation(d, x = "A", y = "Y", m = "M", covariates = c("C1", "C2"),
                                  test = "boot", R = n_resamples, method = "regression",
                                  robust = TRUE)
    # robmed reports no standard error for the bootstrap indirect effect, so
    # only its interval is used, and the bootstrap standard error of the direct
    # effect is not extracted either.
    nm <- c("Indirect", "Direct")
    ci <- confint(tst)
    list(est = unname(coef(tst)[nm]),
         se = c(NA_real_, NA_real_),
         lower = unname(ci[nm, 1L]),
         upper = unname(ci[nm, 2L]))
}

lavaan_model <- "
    M ~ a * A + C1 + C2
    Y ~ b * M + cp * A + C1 + C2
    tnie := a * b
"

fit_lavaan <- function(d, boot) {
    fit <- if (boot) {
        lavaan::sem(lavaan_model, data = d, se = "bootstrap", bootstrap = n_resamples)
    } else {
        lavaan::sem(lavaan_model, data = d)
    }
    pe <- lavaan::parameterEstimates(fit, boot.ci.type = "perc")
    rows <- match(c("tnie", "cp"), pe[["label"]])
    list(est = pe[["est"]][rows], se = pe[["se"]][rows],
         lower = pe[["ci.lower"]][rows], upper = pe[["ci.upper"]][rows])
}

method_specs <- function(scn) {
    binary <- scn[["outcome"]] == "binary"
    medflex_scale <- if (!binary) "md" else if (scn[["exposure"]] == "binary") "logor" else "clogor"
    spec <- function(name, scale, fun) list(name = name, scale = scale, fun = fun)
    common <- list(
        spec("mediation-qb", if (binary) "rd" else "md", function(d) fit_mediation(d, scn, FALSE)),
        spec("mediation-boot", if (binary) "rd" else "md", function(d) fit_mediation(d, scn, TRUE)),
        spec("medflex-weight", medflex_scale, function(d) fit_medflex(d, scn, "weight")),
        spec("medflex-impute", medflex_scale, function(d) fit_medflex(d, scn, "impute")),
        spec("cmaverse-gcomp", if (binary) "logor" else "md", function(d) fit_cmaverse(d, scn)),
        spec("regmedint-delta", if (binary) "clogor" else "md", function(d) fit_regmedint(d, scn))
    )
    continuous_only <- list(
        spec("robmed-mm", "md", fit_robmed),
        spec("lavaan-delta", "md", function(d) fit_lavaan(d, FALSE)),
        spec("lavaan-boot", "md", function(d) fit_lavaan(d, TRUE))
    )
    if (binary) common else c(common, continuous_only)
}

# ── Replications ─────────────────────────────────────────────────────────────

# Runs fun() with console output, messages, and warnings suppressed (CMAverse
# and medflex print progress bars).
quietly <- function(fun) {
    sink(nullfile())
    on.exit(sink(), add = TRUE)
    suppressWarnings(suppressMessages(fun()))
}

run_method <- function(spec, d, seed) {
    set.seed(seed)
    t0 <- proc.time()[["elapsed"]]
    res <- tryCatch(quietly(function() spec[["fun"]](d)), error = function(e) e)
    seconds <- proc.time()[["elapsed"]] - t0
    if (!inherits(res, "error") && !all(is.finite(res[["est"]]))) {
        res <- simpleError("non-finite point estimate")
    }
    failed <- inherits(res, "error")
    pick <- function(field) if (failed) c(NA_real_, NA_real_) else as.numeric(res[[field]])
    data.frame(method = spec[["name"]], effect = c("tnie", "pnde"), scale = spec[["scale"]],
               est = pick("est"), se = pick("se"), lower = pick("lower"), upper = pick("upper"),
               seconds = seconds, error = if (failed) conditionMessage(res) else NA_character_)
}

spec_names <- function(scn) mapply(function(spec) spec[["name"]], method_specs(scn))

# Runs every method, or only those named in `methods`. Each method's seed
# depends on its position among all the methods, so a re-run of one method
# reproduces what a full run gives it.
run_rep <- function(scn, rep, seed, methods = NULL) {
    set.seed(seed)
    d <- gen_data(scn)
    specs <- method_specs(scn)
    method_seeds <- (seed + seq_along(specs)) %% .Machine[["integer.max"]]
    keep <- is.null(methods) | spec_names(scn) %in% methods
    rows <- Map(function(spec, s) run_method(spec, d, s), specs[keep], method_seeds[keep]) |>
        Reduce(f = rbind)
    cbind(scenario = scn[["id"]], rep = rep, rows)
}

# With `methods`, the new rows replace those methods' rows in the scenario's
# existing output, in the row order of a full run.
scenario_rows <- function(rows, scn, path, methods) {
    new <- rows[rows[["scenario"]] == scn[["id"]], ]
    if (is.null(methods)) return(new)
    old <- readRDS(path)
    all <- rbind(old[!(old[["method"]] %in% methods), ], new)
    `rownames<-`(all[order(all[["rep"]], match(all[["method"]], spec_names(scn))), ], NULL)
}

log_line <- function(path, ...) {
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...), "\n", file = path, append = TRUE)
}

# Runs every replication of a batch of scenarios in one parallel map, then
# saves one file per scenario.
run_batch <- function(batch, n_reps, seeds, cfg, out_dir, methods = NULL) {
    tasks <- expand.grid(k = seq_along(batch), rep = seq_len(n_reps))
    task_ids <- mapply(function(k) batch[[k]][["id"]], tasks[["k"]])
    task_seeds <- seeds[["seed"]][match(paste(task_ids, tasks[["rep"]]),
                                        paste(seeds[["id"]], seeds[["rep"]]))]
    t0 <- Sys.time()
    res <- par_map(cfg, function(k, rep, seed) run_rep(batch[[k]], rep, seed, methods),
                   tasks[["k"]], tasks[["rep"]], task_seeds)
    rows <- res |> par_values() |> Reduce(f = rbind)
    invisible(Map(function(scn) {
        path <- file.path(out_dir, sprintf("scenario-%s.rds", scn[["id"]]))
        saveRDS(scenario_rows(rows, scn, path, methods), path)
    }, batch))
    log_line(file.path(out_dir, "run.log"), "%s done in %.1f min",
             paste(Map(function(scn) scn[["id"]], batch), collapse = ","),
             as.numeric(difftime(Sys.time(), t0, units = "mins")))
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
    mode <- if (length(args) > 0L) args[[1L]] else "pilot"
    methods <- if (length(args) > 1L) args[-1L] else NULL
    n_reps <- switch(mode, pilot = 2L, check = 32L, full = 1000L,
                     stop("mode must be 'pilot', 'check', or 'full'"))
    out_dir <- file.path("..", "results", "simulation", mode)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    log_path <- file.path(out_dir, "run.log")

    grid <- make_scenarios()
    scenarios <- Map(function(i) as.list(grid[i, ]), seq_len(nrow(grid)))

    truth_path <- file.path(out_dir, "truth.rds")
    if (!file.exists(truth_path)) {
        truth <- scenarios |>
            Map(f = function(scn) cbind(scenario = scn[["id"]], truth_values(scn))) |>
            Reduce(f = rbind)
        saveRDS(truth, truth_path)
    }
    prevalence_path <- file.path(out_dir, "prevalence.rds")
    if (!file.exists(prevalence_path)) {
        binary <- Filter(function(scn) scn[["outcome"]] == "binary", scenarios)
        prevalence <- data.frame(scenario = mapply(function(scn) scn[["id"]], binary),
                                 prevalence = mapply(outcome_prevalence, binary))
        saveRDS(prevalence, prevalence_path)
    }

    # One seed per (scenario, replication), fixed by the master seed and the
    # number of replications, so a resumed run reproduces the same data.
    seeds <- data.frame(id = rep(grid[["id"]], each = n_reps),
                        rep = rep(seq_len(n_reps), times = nrow(grid)),
                        seed = par_seeds(master_seed, nrow(grid) * n_reps))

    done <- function(scn) file.exists(file.path(out_dir, sprintf("scenario-%s.rds", scn[["id"]])))
    todo <- if (is.null(methods)) {
        Filter(Negate(done), scenarios)
    } else {
        Filter(function(scn) done(scn) && any(spec_names(scn) %in% methods), scenarios)
    }

    # PSOCK workers (Windows) see only what is exported, so every top-level
    # function and constant goes to them. The fork backend ignores the list.
    cfg <- setup_parallel(export = ls(globalenv()))
    on.exit(cleanup_parallel(cfg), add = TRUE)
    log_line(log_path, "start %s: %d of %d scenarios to run, %d replications each, %d cores%s",
             mode, length(todo), length(scenarios), n_reps, cfg[["n_cores"]],
             if (is.null(methods)) "" else paste0(", re-running ", paste(methods, collapse = ", ")))

    per_batch <- max(1L, ceiling(2 * cfg[["n_cores"]] / n_reps))
    batches <- split(todo, ceiling(seq_along(todo) / per_batch))
    invisible(Map(function(batch) run_batch(batch, n_reps, seeds, cfg, out_dir, methods), batches))
    log_line(log_path, "finished %s", mode)
}

# Run only as a script, so that 12_medflex_check.R can source the functions.
if (sys.nframe() == 0L) main()
