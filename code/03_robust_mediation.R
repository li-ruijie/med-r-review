# =============================================================================
# Robust Mediation Analysis with MM-Estimators (robmed)
# =============================================================================
#
# Package: robmed
# Companion: 04_sem_mediation.R (manymome/lavaan) analyses the same path
#            model, enabling direct comparison of robust vs ML estimation.
#
# PURPOSE
# -------
# Demonstrate robust mediation analysis using MM-estimators, which provide
# resistance to outliers and influential observations. This script analyses
# the same path model and dataset as 04_sem_mediation.R (manymome/lavaan),
# enabling direct comparison of estimation strategies on identical data.
#
# DATASET
# -------
# Brader et al.'s "framing" experiment (mediation package).
#   Path:       treat -> emo -> p_harm (continuous outcome)
#   Covariates: age, education, gender (female), income
#
# This script uses the continuous outcome p_harm (perceived harm of
# immigration) rather than the binary cong_mesg used in script 01, because
# MM-estimators require continuous dependent variables. The treatment and
# mediator (treat, emo) are the same, making this a natural complement to
# the foundational comparison.
#
# PACKAGE OVERVIEW
# ----------------
# robmed (Alfons et al., 2022)
#   Implements mediation analysis using MM-estimators (robust regression).
#   MM-estimation is a two-stage procedure:
#     1. S-estimator: finds an initial estimate with high breakdown point
#        (~50%), meaning it remains reliable even if up to half the data
#        are outliers.
#     2. M-estimator: refines the S-estimate for high asymptotic efficiency
#        under normality, so that little precision is lost when the data
#        are well behaved.
#
#   The test_mediation() function fits the full mediation model:
#     - a path:  treat -> emo      (mediator model)
#     - b path:  emo -> p_harm     (controlling for treat)
#     - c' path: treat -> p_harm   (direct effect, controlling for emo)
#     - ab:      a * b             (indirect effect, bootstrap CI)
#
#   USE CASE: When data may contain outliers or influential observations
#   that could bias OLS/ML estimates. Common scenarios:
#     - Survey data with extreme responses or careless respondents
#     - Behavioural experiments with reaction-time outliers
#     - Observational studies with measurement error or data entry errors
#     - Any setting that calls for estimation robust to a few anomalous
#       observations without ad-hoc outlier removal
#   Example: A psychology experiment measuring response times as mediators
#   of a priming effect on decision accuracy, where a few participants
#   have extreme RTs due to inattention. Rather than trimming at an
#   arbitrary threshold (e.g., 2 SD), robmed downweights those cases
#   automatically while preserving the full sample.
#
# COMPARISON WITH 04_sem_mediation.R (manymome/lavaan)
# ----------------------------------------------------
# Both scripts analyse: treat -> emo -> p_harm with the same covariates.
#
# SIMILARITIES:
#   - Same estimand: indirect effect ab = a * b
#   - Same kind of interval: percentile bootstrap for ab (the fast-and-robust
#     bootstrap in robmed, the nonparametric bootstrap in manymome)
#   - Same dataset and path specification
#   - Both produce path coefficients (a, b, c') and indirect effect
#
# DIFFERENCES:
#   - Estimation: MM-estimators (robmed) vs ML (lavaan)
#   - Outlier handling: robmed downweights extreme observations. lavaan
#     (ML) assumes multivariate normality and treats all observations
#     equally, making it sensitive to outliers
#   - Model framework: regression-based (robmed) vs SEM (lavaan)
#   - Flexibility: lavaan supports latent variables, multiple mediators,
#     and complex path models. robmed handles simple, parallel, and serial
#     mediators with observed variables only
#   - When data are clean: both give similar results, and the MM-estimator
#     gives up a little efficiency under normality in exchange for its
#     resistance to outliers
#   - When data contain outliers: robmed estimates will be more stable
#     while lavaan estimates may be pulled toward extreme values
# =============================================================================

# --- Log output --------------------------------------------------------------
.log_file <- {
  f <- grep("--file=", commandArgs(), value = TRUE)
  if (length(f)) sub("\\.R$", ".log", sub("--file=", "", f)) else NULL
}
if (!is.null(.log_file)) sink(.log_file, split = TRUE)

cat(sprintf("Script: %s\n", "03_robust_mediation.R"))
cat(sprintf("Date:   %s\n", Sys.time()))
cat(strrep("=", 72), "\n\n")

suppressPackageStartupMessages({
  library("mediation")  # for the framing dataset
  library("robmed")
})

# --- Data preparation --------------------------------------------------------
data("framing", package = "mediation")
framing$female <- as.numeric(framing$gender == "female")
framing$educ_num <- as.numeric(framing$educ)

# --- robmed: MM-estimator with bootstrap CI -----------------------------------
# test_mediation() fits the mediation model using MM-estimators (robust = TRUE)
# and constructs bootstrap CIs for the indirect effect.
# R = 10000 bootstrap replicates for stable inference.
#
# With robust = FALSE, robmed uses OLS instead of MM-estimation. Comparing
# robust = TRUE vs FALSE on the same data reveals the influence of outliers.
set.seed(6489)
t_robmed <- system.time({
  rob_mm <- test_mediation(framing, x = "treat", y = "p_harm", m = "emo",
    covariates = c("age", "educ_num", "female", "income"),
    test = "boot", robust = TRUE, R = 10000)
})
summary(rob_mm)

# --- Extract path coefficients and indirect effect ----------------------------
# coef() returns the means of the bootstrap distribution by default. The
# estimates on the original data (type = "data") are extracted instead, so
# that they compare like with like with the maximum likelihood estimates of
# 04_sem_mediation.R. The CI of the indirect effect is the percentile
# bootstrap interval.
ab_est <- coef(rob_mm, parm = "indirect", type = "data")
ab_ci  <- confint(rob_mm, parm = "indirect")
a_est  <- coef(rob_mm, parm = "a", type = "data")
b_est  <- coef(rob_mm, parm = "b", type = "data")
cp_est <- coef(rob_mm, parm = "direct", type = "data")

# Robustness weights of the two regressions. The boot_test_mediation object
# stores the fit to the original data in $fit, with the mediator regression
# (lmrob) in $fit$fit_mx and the outcome regression in $fit$fit_ymx.
rob_wt_m <- weights(rob_mm$fit$fit_mx, type = "robustness")
rob_wt <- weights(rob_mm$fit$fit_ymx, type = "robustness")
n_outlier <- sum(rob_wt == 0)

# --- Pretty-printed results --------------------------------------------------
cat("ROBUST MEDIATION: MM-ESTIMATOR BOOTSTRAP (robmed)\n")
cat(strrep("=", 72), "\n")
cat("Data:     framing (mediation package), n =", nrow(framing), "\n")
cat("Path:     treat -> emo -> p_harm\n")
cat("Method:   MM-estimator (robust = TRUE), 10000 bootstrap replicates\n")
cat("Seed:     6489\n\n")

cat(sprintf("TIMING:  robmed %.1fs\n\n", t_robmed["elapsed"]))
cat("PATH COEFFICIENTS\n")
cat(strrep("-", 72), "\n")
cat(sprintf("  a  path (treat -> emo):      %7.3f\n", a_est))
cat(sprintf("  b  path (emo -> p_harm):     %7.3f\n", b_est))
cat(sprintf("  c' path (direct effect):     %7.3f\n", cp_est))
cat(sprintf("  ab      (indirect effect):   %7.3f  95%% CI [%.3f, %.3f]\n",
  ab_est, ab_ci[1], ab_ci[2]))
cat(sprintf("  Minimum weight, mediator model:            %.2f\n",
  min(rob_wt_m)))
cat(sprintf("  Observations with weight 0, outcome model: %d of %d\n",
  n_outlier, length(rob_wt)))
cat("\n")

cat("INTERPRETATION\n")
cat(strrep("-", 72), "\n")
cat("  The indirect effect ab is the product of the a and b paths.\n")
cat("  MM-estimators downweight outliers, providing estimates that are\n")
cat("  resistant to influential observations while maintaining high\n")
cat("  efficiency when data are well-behaved.\n")
cat("  \n")
cat("  Compare these estimates with 04_sem_mediation.R (ML via lavaan).\n")
cat("  The difference between the two shows how much the observations\n")
cat("  that MM-estimation downweights move the ML estimates.\n")
cat(strrep("=", 72), "\n")

if (!is.null(.log_file)) sink()
