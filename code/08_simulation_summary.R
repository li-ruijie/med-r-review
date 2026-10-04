#!/usr/bin/env Rscript
# =============================================================================
# Simulation Summary: Performance Measures with Monte Carlo Standard Errors
# =============================================================================
#
# Reads the per-scenario output of 07_simulation.R and computes, for every
# scenario, method, effect (TNIE, PNDE), and scale, the performance measures
# of Morris, White, and Crowther (2019). Bias, empirical SE, RMSE, coverage of
# the nominal 95% interval, and mean interval width come with their Monte
# Carlo standard errors. The model SE (the root mean square of the replicates'
# standard errors), the median and the largest of those standard errors, the
# share of them above three times the empirical SE, the median interval width,
# the number of failed fits, the number of estimates of absolute size 1,000 or
# more, and the mean runtime come without them.
# The median width is reported because single replicates can produce extreme
# intervals (a bootstrap resample with near-separation, or extreme weights)
# that make the mean width uninformative.
#
# Usage (from code/):
#   Rscript 08_simulation_summary.R check   or   full
#
# Output: ../results/simulation/<mode>/performance.csv
# =============================================================================

scenario_grid <- function() {
    grid <- expand.grid(
        outcome = c("continuous", "binary"),
        exposure = c("binary", "continuous"),
        interaction = c(FALSE, TRUE),
        errors = c("normal", "t3"),
        n = c(100L, 300L, 1000L),
        stringsAsFactors = FALSE
    )
    `[[<-`(grid, "scenario", value = sprintf("s%02d", seq_len(nrow(grid))))
}

# Performance of one (scenario, method, effect, scale) cell. Failed fits are
# excluded from every measure except the number of failed fits and the runtime.
performance <- function(d) {
    ok <- d[is.na(d[["error"]]), ]
    k <- nrow(ok)
    err <- ok[["est"]] - ok[["truth"]]
    covered <- ok[["lower"]] <= ok[["truth"]] & ok[["truth"]] <= ok[["upper"]]
    width <- ok[["upper"]] - ok[["lower"]]
    emp_se <- sd(ok[["est"]])
    rmse <- sqrt(mean(err^2))
    no_se <- all(is.na(ok[["se"]]))
    mod_se <- if (no_se) NA_real_ else sqrt(mean(ok[["se"]]^2, na.rm = TRUE))
    se_median <- if (no_se) NA_real_ else median(ok[["se"]], na.rm = TRUE)
    se_over3 <- if (no_se) NA_real_ else mean(ok[["se"]] > 3 * emp_se, na.rm = TRUE)
    se_max <- if (no_se) NA_real_ else max(ok[["se"]], na.rm = TRUE)
    cover <- mean(covered)
    data.frame(
        scenario = d[["scenario"]][1L], method = d[["method"]][1L], effect = d[["effect"]][1L],
        scale = d[["scale"]][1L], truth = d[["truth"]][1L],
        reps = nrow(d), failed = nrow(d) - k, extreme = sum(abs(ok[["est"]]) >= 1e3),
        bias = mean(err), bias_mcse = emp_se / sqrt(k),
        emp_se = emp_se, emp_se_mcse = emp_se / sqrt(2 * (k - 1)),
        mod_se = mod_se, se_median = se_median, se_over3 = se_over3, se_max = se_max,
        rmse = rmse, rmse_mcse = sqrt(var(err^2) / k) / (2 * rmse),
        coverage = cover, coverage_mcse = sqrt(cover * (1 - cover) / k),
        width = mean(width), width_mcse = sd(width) / sqrt(k), width_median = median(width),
        seconds = mean(d[["seconds"]])
    )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
    mode <- if (length(args) > 0L) args[[1L]] else "check"
    dir <- file.path("..", "results", "simulation", mode)
    files <- list.files(dir, pattern = "^scenario-s[0-9]+\\.rds$", full.names = TRUE)
    res <- files |> Map(f = readRDS) |> Reduce(f = rbind)
    truth <- readRDS(file.path(dir, "truth.rds"))
    res <- merge(res, truth, by = c("scenario", "effect", "scale"))

    cells <- split(res, list(res[["scenario"]], res[["method"]], res[["effect"]], res[["scale"]]),
                   drop = TRUE)
    perf <- cells |> Map(f = performance) |> Reduce(f = rbind) |> merge(scenario_grid(), by = "scenario")
    perf <- perf[order(perf[["scenario"]], perf[["effect"]], perf[["method"]]), ]
    write.csv(perf, file.path(dir, "performance.csv"), row.names = FALSE)
    message(sprintf("%d scenarios, %d cells written to %s", length(files), nrow(perf),
                    file.path(dir, "performance.csv")))
}

main()
