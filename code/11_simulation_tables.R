#!/usr/bin/env Rscript
# =============================================================================
# Simulation Tables: Summary and Full Results
# =============================================================================
#
# Writes the LaTeX rows of the simulation tables from the performance summary
# of 08_simulation_summary.R.
#
#   summary-rows.tex        Summary at n = 1,000 over four scenario groups
#   scenario-rows.tex       Every scenario, estimator, and scale
#   runtime-rows.tex        Mean runtime per fit by estimator and sample size
#
# Output: ../results/tables/, and the summary statistics of the simulation, each
# with its definition, in
# ../results/simulation/full/summary.log
#
# Run from code/ after 08_simulation_summary.R full.
# =============================================================================

estimators <- list(
    list(method = "mediation-qb", label = "\\textbf{mediation}, quasi-Bayesian", binary = "rd"),
    list(method = "mediation-boot", label = "\\textbf{mediation}, bootstrap", binary = "rd"),
    list(method = "medflex-weight", label = "\\textbf{medflex}, weighting", binary = c("logor", "clogor")),
    list(method = "medflex-impute", label = "\\textbf{medflex}, imputation", binary = c("logor", "clogor")),
    list(method = "cmaverse-gcomp", label = "\\textbf{CMAverse}, g-computation", binary = "logor"),
    list(method = "regmedint-delta", label = "\\textbf{regmedint}, closed form", binary = "clogor"),
    list(method = "robmed-mm", label = "\\textbf{robmed}, MM-estimation", binary = character(0)),
    list(method = "lavaan-delta", label = "\\textbf{lavaan}, delta method", binary = character(0)),
    list(method = "lavaan-boot", label = "\\textbf{lavaan}, bootstrap", binary = character(0))
)

scale_labels <- c(md = "MD", rd = "RD", logor = "log OR", clogor = "cond.\\ log OR")

# Fixed decimals with a typographic minus sign. Values that round to zero
# print without a sign, and values of 1,000 or more in absolute size (the
# degenerate medflex weighting fits under t errors) print in scientific form.
num <- function(x, digits) {
    if (is.na(x)) return("--")
    if (abs(x) >= 1e3) {
        e <- floor(log10(abs(x)))
        return(sprintf("\\(%s%.1f\\mathrm{e}%d\\)", if (x < 0) "-" else "", abs(x) / 10^e, e))
    }
    rounded <- round(x, digits)
    shown <- if (rounded == 0) 0 else rounded
    sub("^-", "\\\\(-\\\\)", formatC(shown, format = "f", digits = digits))
}

# Rows of `perf` that match every named column value in `by`.
select_rows <- function(perf, by) {
    keep <- Map(function(col, val) perf[[col]] %in% val, names(by), by) |> Reduce(f = `&`)
    perf[keep, ]
}

# Mean of `measure` over the selected rows (the two exposure types).
group_mean <- function(perf, by, measure) {
    d <- select_rows(perf, by)
    if (nrow(d) == 0L) NA_real_ else mean(d[[measure]])
}

# Median interval width relative to mediation (quasi-Bayesian) in the same
# scenario and effect, averaged over the selected scenarios.
relative_width <- function(perf, by) {
    d <- select_rows(perf, by)
    ref <- select_rows(perf, `[[<-`(by, "method", "mediation-qb"))
    m <- merge(d, ref, by = c("scenario", "effect"), suffixes = c("", "_ref"))
    mean(m[["width_median"]] / m[["width_median_ref"]])
}

# Scenario groups of the summary table, all at n = 1,000 and averaged over
# the binary and continuous exposures.
groups <- list(
    correct = list(outcome = "continuous", interaction = FALSE, errors = "normal", width = TRUE),
    interaction = list(outcome = "continuous", interaction = TRUE, errors = "normal", width = FALSE),
    t3 = list(outcome = "continuous", interaction = FALSE, errors = "t3", width = TRUE),
    binary = list(outcome = "binary", interaction = TRUE, errors = "normal", width = FALSE)
)

summary_cells <- function(perf, est, effect, g) {
    scales <- if (g[["outcome"]] == "binary") est[["binary"]] else "md"
    by <- list(method = est[["method"]], effect = effect, scale = scales, outcome = g[["outcome"]],
               interaction = g[["interaction"]], errors = g[["errors"]], n = 1000L)
    cells <- c(num(group_mean(perf, by, "bias"), 3), num(group_mean(perf, by, "coverage"), 2))
    if (g[["width"]]) c(cells, num(relative_width(perf, by), 2)) else cells
}

summary_row <- function(perf, est, effect, first) {
    cells <- groups |> Map(f = function(g) summary_cells(perf, est, effect, g)) |> unlist()
    # The row breaks after the effect so that the source lines stay short.
    sprintf("%s & %s\n    & %s \\\\", if (first) est[["label"]] else "", toupper(effect),
            paste(cells, collapse = " & "))
}

summary_rows <- function(perf) {
    estimators |>
        Map(f = function(est) c(summary_row(perf, est, "tnie", TRUE),
                                summary_row(perf, est, "pnde", FALSE))) |>
        unlist()
}

scenario_header <- function(s) {
    sprintf(paste0("\\multicolumn{8}{@{}l}{\\textit{%s: %s outcome, %s exposure, %s,\n",
                   "    %s errors, \\(n = %s\\)}} \\\\*"),
            s[["scenario"]], s[["outcome"]], s[["exposure"]],
            if (s[["interaction"]]) "interaction" else "no interaction",
            if (s[["errors"]] == "t3") "\\(t_3\\)" else "normal", format(s[["n"]], big.mark = "{,}"))
}

effect_cells <- function(d) {
    if (nrow(d) == 0L) return(rep("--", 3L))
    c(sprintf("%s (%s)", num(d[["bias"]], 3), num(d[["bias_mcse"]], 3)), num(d[["coverage"]], 3),
      num(d[["width_median"]], 3))
}

result_row <- function(perf_s, est, scale) {
    d <- perf_s[perf_s[["method"]] == est[["method"]] & perf_s[["scale"]] == scale, ]
    sprintf("%s & %s\n    & %s\n    & %s \\\\", est[["label"]], scale_labels[[scale]],
            paste(effect_cells(d[d[["effect"]] == "tnie", ]), collapse = " & "),
            paste(effect_cells(d[d[["effect"]] == "pnde", ]), collapse = " & "))
}

scenario_rows <- function(perf, scenario) {
    perf_s <- perf[perf[["scenario"]] == scenario, ]
    pairs <- estimators |>
        Map(f = function(est) {
            scales <- intersect(c("md", "rd", "logor", "clogor"),
                                perf_s[perf_s[["method"]] == est[["method"]], "scale"])
            Map(function(sc) result_row(perf_s, est, sc), scales)
        }) |>
        unlist()
    # Every row but the last ends with \\* so that longtable keeps each scenario on one page.
    rows <- c(paste0(head(pairs, -1L), "*"), tail(pairs, 1L))
    c(scenario_header(perf_s[1L, ]), rows, "\\addlinespace")
}

runtime_rows <- function(perf) {
    estimators |>
        Map(f = function(est) {
            secs <- c(100L, 300L, 1000L) |>
                Map(f = function(n) {
                    mean(perf[perf[["method"]] == est[["method"]] & perf[["n"]] == n, "seconds"])
                }) |>
                unlist()
            sprintf("%s & %s \\\\", est[["label"]], paste(mapply(num, secs, 2L), collapse = " & "))
        }) |>
        unlist()
}

# ── Summary statistics ──────────────────────────────────────────────────

omits_interaction <- c("robmed-mm", "lavaan-delta", "lavaan-boot")
counterfactual <- c("mediation-qb", "mediation-boot", "medflex-weight", "medflex-impute", "cmaverse-gcomp",
                    "regmedint-delta")

# A cell counts as correctly specified when the estimator's working models
# contain the data-generating ones. The exception is medflex with a continuous
# exposure and a binary outcome, whose logit-linear natural effect model only
# approximates the nested means. Under t errors only the continuous outcome is
# counted, since with a binary outcome mediation and CMAverse draw mediators
# from a normal model. medflex imputation would also qualify there, exactly
# with a binary exposure and nearly with a continuous one. Leaving out the
# binary exposure changes no summary statistic, but the continuous exposure
# would raise the largest n = 1,000 bias on the conditional log OR to 0.0081
# and the n = 100 standard-error ratios to 18.2 (root mean square) and 324.05
# (largest), and lower the smallest median ratio to 1.1014, so these
# statistics are scoped to normal errors. medflex
# weighting does not qualify under t errors, since its weights use the normal
# mediator density. robmed and lavaan omit the interaction, and the closed
# form of regmedint relies on the rare-outcome approximation, which the binary
# outcome of the design does not meet.
is_correct <- function(perf) {
    m <- perf[["method"]]
    (perf[["errors"]] == "normal" | perf[["outcome"]] == "continuous") &
        !(m %in% omits_interaction & perf[["interaction"]]) &
        !(m == "regmedint-delta" & perf[["outcome"]] == "binary") &
        !(m == "medflex-weight" & perf[["errors"]] == "t3")
}

cell_id <- function(d) {
    sprintf("%s %s %s %s", d[["scenario"]], d[["method"]], toupper(d[["effect"]]), d[["scale"]])
}

# Smallest and largest value of `measure` over the rows of `d`, with the
# cells that attain them.
extremes <- function(d, measure) {
    if (nrow(d) == 0L) return("no cells")
    lo <- which.min(d[[measure]])
    hi <- which.max(d[[measure]])
    sprintf("%.4f (%s) to %.4f (%s)", d[[measure]][[lo]], cell_id(d[lo, ]), d[[measure]][[hi]],
            cell_id(d[hi, ]))
}

# Adds the ratio of the model standard error (the root mean square of the
# replicates' standard errors, 08) to the empirical standard error.
se_ratio <- function(d) `[[<-`(d, "se_ratio", value = d[["mod_se"]] / d[["emp_se"]])

# The same ratio with the median of the replicates' standard errors (08).
median_se_ratio <- function(d) `[[<-`(d, "se_ratio", value = d[["se_median"]] / d[["emp_se"]])

# The same ratio with the largest of the replicates' standard errors (08).
max_se_ratio <- function(d) `[[<-`(d, "se_ratio", value = d[["se_max"]] / d[["emp_se"]])

# Largest absolute value of `measure`, with its cell.
largest_abs <- function(d, measure) {
    i <- which.max(abs(d[[measure]]))
    sprintf("%.4f (%s)", d[[measure]][[i]], cell_id(d[i, ]))
}

# Median widths of methods `a` and `b` in the same scenario, effect, and
# scale, with ratio = width of b / width of a.
paired_widths <- function(d, a, b) {
    keys <- c("scenario", "effect", "scale", "n")
    da <- d[d[["method"]] == a, c(keys, "width_median")]
    db <- d[d[["method"]] == b, c(keys, "width_median")]
    m <- merge(da, db, by = keys, suffixes = c("_a", "_b"))
    m <- `[[<-`(m, "method", value = rep(sprintf("%s/%s", b, a), nrow(m)))
    `[[<-`(m, "ratio", value = m[["width_median_b"]] / m[["width_median_a"]])
}

pairs_of <- function(d, pairs) {
    pairs |> Map(f = function(p) paired_widths(d, p[[1L]], p[[2L]])) |> Reduce(f = rbind)
}

# Largest ratio of median widths among the estimators of one cell.
cell_spread <- function(d) {
    hi <- which.max(d[["width_median"]])
    lo <- which.min(d[["width_median"]])
    data.frame(scenario = d[["scenario"]][[1L]], effect = d[["effect"]][[1L]], scale = d[["scale"]][[1L]],
               method = sprintf("%s/%s", d[["method"]][[hi]], d[["method"]][[lo]]),
               ratio = d[["width_median"]][[hi]] / d[["width_median"]][[lo]])
}

stat <- function(key, definition, value) {
    c(sprintf("%s: %s", key, paste(value, collapse = "; ")), sprintf("    %s", definition), "")
}

pick <- function(d, ...) select_rows(d, list(...))

# Largest ratio of median widths among the correct estimators of each cell
# with at least two of them.
width_spread <- function(d) {
    split(d, list(d[["scenario"]], d[["effect"]], d[["scale"]]), drop = TRUE) |>
        Filter(f = function(cell) nrow(cell) >= 2L) |>
        Map(f = cell_spread) |>
        Reduce(f = rbind)
}

# Median TNIE width of lavaan over each counterfactual estimator, continuous
# outcome, no interaction.
lavaan_ratios <- function(d, comparators) {
    pairs <- Reduce(c, Map(function(cf) list(c(cf, "lavaan-delta"), c(cf, "lavaan-boot")), comparators))
    pairs_of(pick(d, outcome = "continuous", interaction = FALSE, effect = "tnie"), pairs)
}

by_exposure <- function(d, effect, measure) {
    Map(function(e) sprintf("%s %s", e, extremes(pick(d, exposure = e, effect = effect), measure)),
        c("binary", "continuous")) |>
        unlist()
}

summary_statistics <- function(perf, prevalence) {
    ok <- perf[is_correct(perf), ]
    ok_nw <- ok[ok[["method"]] != "medflex-weight", ]
    ok_ls <- ok_nw[ok_nw[["method"]] != "robmed-mm", ]
    normal_ok <- ok[ok[["errors"]] == "normal", ]
    normal_nw <- normal_ok[normal_ok[["method"]] != "medflex-weight", ]
    tnie_rows <- perf[perf[["effect"]] == "tnie", ]
    failures <- tapply(tnie_rows[["failed"]], tnie_rows[["method"]], sum)
    switches <- pairs_of(ok, list(c("mediation-qb", "mediation-boot"), c("lavaan-delta", "lavaan-boot")))
    rd100 <- switches[["scale"]] == "rd" & switches[["n"]] == 100L
    bias_by_scale <- function(d) {
        Map(function(sc) sprintf("%s %s", sc, largest_abs(pick(d, n = 1000L, scale = sc), "bias")),
            c("md", "rd", "logor", "clogor")) |>
            unlist()
    }
    logor_vs <- pairs_of(pick(normal_ok, scale = "logor"),
                         list(c("cmaverse-gcomp", "medflex-impute"), c("cmaverse-gcomp", "medflex-weight")))
    omitted <- pick(perf, method = omits_interaction, interaction = TRUE, outcome = "continuous", n = 1000L)
    regmed <- pick(perf, method = "regmedint-delta", outcome = "binary", errors = "normal")
    weight_t <- pick(perf, method = "medflex-weight", outcome = "continuous", errors = "t3",
                     effect = "tnie")
    weight_cc <- pick(perf, method = "medflex-weight", exposure = "continuous", outcome = "continuous",
                      errors = "normal")
    weight_cc <- se_ratio(weight_cc)
    weight_cc <- `[[<-`(weight_cc, "bias_se", value = weight_cc[["bias"]] / weight_cc[["emp_se"]])
    # The quasi-Bayesian point estimate of mediation is the mean of the
    # simulated effects, and the bootstrap one the effect on the original data.
    med100 <- pick(perf, method = c("mediation-qb", "mediation-boot"), outcome = "binary",
                   errors = "normal", n = 100L, effect = "tnie", scale = "rd")
    med100 <- `[[<-`(med100, "rel_bias", value = med100[["bias"]] / med100[["truth"]])
    med100 <- `[[<-`(med100, "bias_z", value = med100[["bias"]] / med100[["bias_mcse"]])
    by_mediation <- function(measure) {
        Map(function(m) sprintf("%s %s", m, extremes(pick(med100, method = m), measure)),
            c("mediation-qb", "mediation-boot")) |>
            unlist()
    }
    weight_bin <- max_se_ratio(pick(perf, method = "medflex-weight", outcome = "binary", errors = "normal"))
    weight_bb <- pick(weight_bin, exposure = "binary")
    weight_bb <- weight_bb[weight_bb[["se_ratio"]] > 1000, ]
    weight_bc100 <- pick(weight_bin, exposure = "continuous", n = 100L, effect = "pnde")
    robmed_t <- pairs_of(pick(perf, outcome = "continuous", interaction = FALSE, errors = "t3"),
                         list(c("lavaan-delta", "robmed-mm"), c("lavaan-boot", "robmed-mm")))
    robmed_normal <- pairs_of(pick(perf, outcome = "continuous", interaction = FALSE, errors = "normal"),
                              list(c("lavaan-delta", "robmed-mm"), c("lavaan-boot", "robmed-mm")))
    summary_cells <- pick(perf, n = 1000L, outcome = "continuous", interaction = FALSE) |>
        rbind(pick(perf, n = 1000L, outcome = "continuous", interaction = TRUE, errors = "normal")) |>
        rbind(pick(perf, n = 1000L, outcome = "binary", interaction = TRUE, errors = "normal"))
    summary_mcse <- summary_cells[!(summary_cells[["method"]] == "medflex-weight" &
                                    summary_cells[["errors"]] == "t3"), ]
    runtime <- Map(function(m) {
        secs <- mapply(function(n) mean(perf[perf[["method"]] == m & perf[["n"]] == n, "seconds"]),
                       c(100L, 300L, 1000L))
        sprintf("%s %.2f to %.2f", m, min(secs), max(secs))
    }, unique(perf[["method"]]))
    c(
        "Summary statistics of the simulation. A cell is one scenario, estimator,",
        "effect, and scale. Correct cells: see is_correct() in 11_simulation_tables.R.",
        "",
        stat("outcome prevalence",
             "P(Y = 1) of the binary-outcome scenarios, by numerical integration (07_simulation.R).",
             sprintf("%.4f (%s) to %.4f (%s)", min(prevalence[["prevalence"]]),
                     prevalence[["scenario"]][[which.min(prevalence[["prevalence"]])]],
                     max(prevalence[["prevalence"]]),
                     prevalence[["scenario"]][[which.max(prevalence[["prevalence"]])]])),
        stat("failed fits", "Failed fits per estimator over all scenarios and replications.",
             sprintf("%s %d", names(failures), as.integer(failures))),
        stat("extreme estimates",
             "Cells with estimates of absolute size 1,000 or more, and their number in each cell.",
             Map(function(i) {
                 cell <- perf[i, ]
                 sprintf("%s %d", cell_id(cell), cell[["extreme"]])
             }, which(perf[["extreme"]] > 0)) |> unlist()),
        stat("bias, correct, n = 1,000", "Largest absolute bias of a correct cell at n = 1,000, by scale.",
             bias_by_scale(ok)),
        stat("bias, correct, n = 1,000, no weighting", "As above, without medflex weighting.",
             bias_by_scale(ok_nw)),
        stat("coverage, correct", "Coverage of the correct cells, all sample sizes.",
             extremes(ok, "coverage")),
        stat("coverage, correct, no weighting", "As above, without medflex weighting.",
             extremes(ok[ok[["method"]] != "medflex-weight", ], "coverage")),
        stat("coverage, delta TNIE, n = 100",
             "Coverage of the delta-method TNIE intervals in correct cells.",
             extremes(pick(ok, n = 100L, effect = "tnie", method = c("lavaan-delta", "regmedint-delta")),
                      "coverage")),
        stat("coverage, resampling TNIE, n = 100",
             paste("Coverage of the bootstrap and quasi-Bayesian TNIE intervals in correct cells with a",
                   "continuous outcome (the scenarios of the delta-method statistic), without medflex",
                   "weighting."),
             extremes(pick(ok, n = 100L, effect = "tnie", outcome = "continuous",
                           method = c("mediation-qb", "mediation-boot", "medflex-impute", "cmaverse-gcomp",
                                      "lavaan-boot", "robmed-mm")), "coverage")),
        stat("mediation, binary, n = 100, TNIE RD relative bias",
             paste("Bias over the true value of the quasi-Bayesian and bootstrap point estimates of",
                   "mediation, binary outcome, normal errors, n = 100."),
             by_mediation("rel_bias")),
        stat("mediation, binary, n = 100, TNIE RD bias over MCSE",
             "As above, bias over its Monte Carlo SE.", by_mediation("bias_z")),
        stat("coverage, medflex imputation, cond. log OR, n = 100",
             paste("Coverage of medflex imputation on the conditional log OR scale at n = 100, correct",
                   "cells, normal errors."),
             extremes(pick(ok, method = "medflex-impute", scale = "clogor", n = 100L), "coverage")),
        stat("medflex imputation, cond. log OR, n = 100, SE ratio",
             paste("Root mean square of the bootstrap standard errors over the empirical standard error,",
                   "TNIE then PNDE, correct cells, normal errors."),
             Map(function(e) extremes(se_ratio(pick(ok, method = "medflex-impute", scale = "clogor",
                                                    n = 100L, effect = e)), "se_ratio"),
                 c("tnie", "pnde")) |> unlist()),
        stat("medflex imputation, cond. log OR, n = 100, median SE ratio",
             "As above, with the median of the bootstrap standard errors.",
             Map(function(e) extremes(median_se_ratio(pick(ok, method = "medflex-impute", scale = "clogor",
                                                           n = 100L, effect = e)), "se_ratio"),
                 c("tnie", "pnde")) |> unlist()),
        stat("medflex imputation, cond. log OR, n = 100, largest SE ratio",
             "As above, with the largest bootstrap standard error.",
             Map(function(e) extremes(max_se_ratio(pick(ok, method = "medflex-impute", scale = "clogor",
                                                        n = 100L, effect = e)), "se_ratio"),
                 c("tnie", "pnde")) |> unlist()),
        stat("medflex imputation, cond. log OR, n = 100, large SE share",
             "As above, share of replications with a bootstrap SE above three times the empirical SE.",
             Map(function(e) extremes(pick(ok, method = "medflex-impute", scale = "clogor", n = 100L,
                                           effect = e), "se_over3"),
                 c("tnie", "pnde")) |> unlist()),
        stat("interval method, width change",
             paste("Median width with the bootstrap over the quasi-Bayesian (mediation) or delta-method",
                   "(lavaan) interval, correct cells, without RD at n = 100."),
             extremes(switches[!rd100, ], "ratio")),
        stat("interval method, width change, RD n = 100", "As above, RD cells at n = 100.",
             extremes(switches[rd100, ], "ratio")),
        stat("same estimand, width spread",
             paste("Widest over narrowest median width among the correct estimators of one scenario,",
                   "effect, and scale, normal errors."),
             extremes(width_spread(normal_ok), "ratio")),
        stat("same estimand, width spread, no weighting", "As above, without medflex weighting.",
             extremes(width_spread(normal_nw), "ratio")),
        stat("lavaan TNIE width",
             paste("Median TNIE width of lavaan over that of each counterfactual estimator, continuous",
                   "outcome, no interaction, normal errors."),
             extremes(lavaan_ratios(normal_ok, counterfactual), "ratio")),
        stat("lavaan TNIE width, no weighting", "As above, without medflex weighting.",
             extremes(lavaan_ratios(normal_ok, setdiff(counterfactual, "medflex-weight")), "ratio")),
        stat("lavaan TNIE width, medflex imputation", "As above, against medflex imputation only.",
             extremes(lavaan_ratios(normal_ok, "medflex-impute"), "ratio")),
        stat("medflex and CMAverse, log OR width",
             paste("Median width of medflex over CMAverse on the marginal log OR scale, correct cells,",
                   "normal errors."),
             Map(function(m) extremes(logor_vs[logor_vs[["method"]] == m, ], "ratio"),
                 unique(logor_vs[["method"]])) |> unlist()),
        stat("omitted interaction, TNIE bias",
             paste("robmed and lavaan, continuous outcome with interaction, n = 1,000, normal and t",
                   "errors, by exposure."),
             by_exposure(omitted, "tnie", "bias")),
        stat("omitted interaction, TNIE coverage", "As above, coverage.",
             by_exposure(omitted, "tnie", "coverage")),
        stat("omitted interaction, PNDE bias", "As above, PNDE bias.",
             by_exposure(omitted, "pnde", "bias")),
        stat("omitted interaction, PNDE coverage", "As above, PNDE coverage.",
             by_exposure(omitted, "pnde", "coverage")),
        stat("regmedint, binary, no interaction, bias",
             "Closed form, binary outcome, normal errors, all n.",
             largest_abs(pick(regmed, interaction = FALSE), "bias")),
        stat("regmedint, binary, no interaction, coverage", "As above, coverage.",
             extremes(pick(regmed, interaction = FALSE), "coverage")),
        stat("regmedint, binary, interaction, PNDE bias, n = 1,000",
             "Closed form, PNDE, cond. log OR, normal errors.",
             extremes(pick(regmed, interaction = TRUE, effect = "pnde", n = 1000L), "bias")),
        stat("regmedint, binary, interaction, PNDE truth", "True cond. log OR of the PNDE.",
             extremes(pick(regmed, interaction = TRUE, effect = "pnde", n = 1000L), "truth")),
        stat("regmedint, binary, interaction, PNDE coverage, n = 1,000", "As above, coverage.",
             extremes(pick(regmed, interaction = TRUE, effect = "pnde", n = 1000L), "coverage")),
        stat("medflex weighting, coverage, normal errors",
             "All cells with normal errors, by exposure, TNIE then PNDE.",
             c(by_exposure(pick(perf, method = "medflex-weight", errors = "normal"), "tnie", "coverage"),
               by_exposure(pick(perf, method = "medflex-weight", errors = "normal"), "pnde", "coverage"))),
        stat("medflex weighting, continuous exposure and outcome, coverage",
             "Normal errors, every n, with or without the interaction.", extremes(weight_cc, "coverage")),
        stat("medflex weighting, continuous exposure and outcome, SE ratio",
             paste("Root mean square of the bootstrap standard errors over the empirical standard",
                   "error, as above."),
             extremes(weight_cc, "se_ratio")),
        stat("medflex weighting, continuous exposure and outcome, TNIE bias over empirical SE",
             "As above, TNIE only.", extremes(pick(weight_cc, effect = "tnie"), "bias_se")),
        stat("medflex weighting, binary outcome, largest SE ratio",
             paste("Largest bootstrap standard error over the empirical standard error, normal errors,",
                   "by exposure, every n and both effects."),
             Map(function(e) sprintf("%s %s", e, extremes(pick(weight_bin, exposure = e), "se_ratio")),
                 c("binary", "continuous")) |> unlist()),
        stat("medflex weighting, binary outcome, binary exposure, SE ratio above 1,000",
             "As above, every binary-exposure cell whose largest ratio exceeds 1,000.",
             sprintf("%s %.4g", cell_id(weight_bb), weight_bb[["se_ratio"]])),
        stat("medflex weighting, binary outcome, continuous exposure, n = 100, PNDE large SE share",
             "Share of replications with a bootstrap SE above three times the empirical SE, normal errors.",
             extremes(weight_bc100, "se_over3")),
        stat("medflex weighting, binary outcome, continuous exposure, n = 100, PNDE coverage",
             "As above, coverage.", extremes(weight_bc100, "coverage")),
        stat("medflex weighting, t errors, empirical SE",
             "TNIE, continuous outcome, mean over the four t-error scenarios at each n.",
             mapply(function(n) sprintf("n = %s: %.4f", format(n, big.mark = ","),
                                        mean(pick(weight_t, n = n)[["emp_se"]])),
                    c(100L, 300L, 1000L))),
        stat("medflex weighting, t errors, coverage", "All cells with t errors.",
             extremes(pick(perf, method = "medflex-weight", errors = "t3"), "coverage")),
        stat("robmed over lavaan, t errors, width",
             "Median width of robmed over lavaan, continuous outcome, no interaction, t errors.",
             extremes(robmed_t, "ratio")),
        stat("coverage, t errors, correct, no weighting",
             paste("Coverage of the correct cells under t errors (continuous outcome), without medflex",
                   "weighting and robmed."),
             extremes(pick(ok_ls, errors = "t3", outcome = "continuous"), "coverage")),
        stat("bias, t errors, correct, no weighting", "As above, largest absolute bias.",
             largest_abs(pick(ok_ls, errors = "t3", outcome = "continuous"), "bias")),
        stat("robmed, t errors, coverage",
             "Coverage of robmed under t errors (continuous outcome, no interaction).",
             extremes(pick(perf, method = "robmed-mm", interaction = FALSE, errors = "t3"), "coverage")),
        stat("robmed over lavaan, normal errors, width",
             "Median width of robmed over lavaan, continuous outcome, no interaction, normal errors.",
             extremes(robmed_normal, "ratio")),
        stat("runtime", "Mean seconds per fit, smallest and largest over the three sample sizes.",
             unlist(runtime)),
        stat("summary table, bias MCSE",
             "Cells of the summary table, without medflex weighting under t errors.",
             extremes(summary_mcse, "bias_mcse")),
        stat("summary table, coverage MCSE", "As above.", extremes(summary_mcse, "coverage_mcse"))
    )
}

main <- function() {
    sim_dir <- file.path("..", "results", "simulation", "full")
    perf <- read.csv(file.path(sim_dir, "performance.csv"), stringsAsFactors = FALSE)
    out_dir <- file.path("..", "results", "tables")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    writeLines(summary_rows(perf), file.path(out_dir, "summary-rows.tex"))
    scenarios <- sort(unique(perf[["scenario"]]))
    writeLines(scenarios |> Map(f = function(s) scenario_rows(perf, s)) |> unlist(),
               file.path(out_dir, "scenario-rows.tex"))
    writeLines(runtime_rows(perf), file.path(out_dir, "runtime-rows.tex"))
    prevalence <- readRDS(file.path(sim_dir, "prevalence.rds"))
    writeLines(summary_statistics(perf, prevalence), file.path(sim_dir, "summary.log"))
    message(sprintf("Wrote the summary, %d scenarios, and the runtime table to %s, and summary.log to %s",
                    length(scenarios), out_dir, sim_dir))
}

main()
