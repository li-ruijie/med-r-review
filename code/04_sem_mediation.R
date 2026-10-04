# =============================================================================
# SEM-Based Mediation Analysis (manymome + lavaan)
# =============================================================================
#
# Packages: manymome, lavaan
# Companion: 03_robust_mediation.R (robmed) analyses the same path model,
#            enabling direct comparison of ML vs robust estimation.
#
# PURPOSE
# -------
# Demonstrate SEM-based mediation analysis using lavaan for model specification
# and manymome for computing indirect effects with bootstrap CIs. This script
# analyses the same path model and dataset as 03_robust_mediation.R (robmed),
# enabling direct comparison of maximum-likelihood vs robust estimation.
#
# DATASET
# -------
# Brader et al.'s "framing" experiment (mediation package).
#   Path:       treat -> emo -> p_harm (continuous outcome)
#   Covariates: age, education, gender (female), income
#
# Same path and data as 03_robust_mediation.R for paired comparison.
#
# PACKAGE OVERVIEW
# ----------------
# lavaan (Rosseel, 2012)
#   The standard R package for structural equation modelling (SEM). Specifies
#   path models using a concise formula syntax with labelled parameters.
#   Estimation via maximum likelihood (ML). Supports latent variables,
#   multiple groups, missing data (FIML), and complex models.
#
# manymome (Cheung & Cheung, 2024)
#   A post-estimation package that computes indirect effects, conditional
#   indirect effects, and standardised effects from fitted lavaan (or lm)
#   models. Uses the nonparametric bootstrap, or Monte Carlo simulation for
#   lavaan models, for inference on indirect effects. Designed for
#   models with "many" mediators and/or moderators.
#
#   USE CASE: When the mediation model is part of a larger SEM or involves
#   multiple mediators, serial mediation, or moderated mediation. The lavaan
#   + manymome combination offers:
#     - Flexible path model specification (lavaan syntax with labels)
#     - Latent variable support (measurement models + structural paths)
#     - Bootstrap inference for indirect effects (manymome)
#     - Conditional indirect effects at different moderator levels
#   Example: A health psychology study testing whether a behavioural
#   intervention reduces BMI through two serial mediators (self-efficacy
#   then exercise frequency), with treatment fidelity as a moderator.
#   lavaan specifies the full path model, and manymome computes each
#   specific indirect path and tests conditional effects at high and low
#   fidelity.
#
#   Another example: An organisational behaviour study with latent constructs
#   (job satisfaction measured by 5 items, burnout measured by 3 items)
#   where lavaan handles the measurement model and manymome computes the
#   indirect effect through the latent mediator.
#
# COMPARISON WITH 03_robust_mediation.R (robmed)
# -----------------------------------------------
# Both scripts analyse: treat -> emo -> p_harm with the same covariates.
#
# SIMILARITIES:
#   - Same estimand: indirect effect ab = a * b
#   - Same kind of interval: percentile bootstrap for ab (the nonparametric
#     bootstrap in manymome, the fast-and-robust bootstrap in robmed)
#   - Same dataset and path specification
#
# DIFFERENCES:
#   - Estimation: ML (lavaan) vs MM-estimators (robmed)
#   - Assumptions: lavaan (ML) assumes multivariate normality, and robmed
#     does not
#   - Model framework: SEM (lavaan) vs regression (robmed)
#   - Flexibility: lavaan handles latent variables, multiple mediators,
#     serial mediation, and moderated mediation. robmed handles simple,
#     parallel, and serial mediators with observed variables only
#   - Architecture: lavaan + manymome separates model fitting from effect
#     computation (one fitted model, many queries), and robmed's
#     test_mediation() combines both steps
#   - When data are clean: both give similar ab estimates, since the robust
#     and ML estimators target the same coefficients under normality
#   - When data contain outliers: robmed is more resistant
#
# BOOTSTRAP INTERVALS FOR BOTH
#   The indirect effect ab = a * b is a product of two coefficients. Even
#   if a and b are normally distributed, their product follows a skewed
#   distribution. The bootstrap approximates the sampling distribution of ab
#   directly, whereas the Sobel test and the delta method treat ab as
#   normal.
# =============================================================================

# --- Log output --------------------------------------------------------------
.log_file <- {
  f <- grep("--file=", commandArgs(), value = TRUE)
  if (length(f)) sub("\\.R$", ".log", sub("--file=", "", f)) else NULL
}
if (!is.null(.log_file)) sink(.log_file, split = TRUE)

cat(sprintf("Script: %s\n", "04_sem_mediation.R"))
cat(sprintf("Date:   %s\n", Sys.time()))
cat(strrep("=", 72), "\n\n")

suppressPackageStartupMessages({
  library("mediation")  # for the framing dataset
  library("manymome")
  library("lavaan")
})

# --- Data preparation --------------------------------------------------------
data("framing", package = "mediation")
framing$female <- as.numeric(framing$gender == "female")
framing$educ_num <- as.numeric(framing$educ)

# --- lavaan: specify and fit the SEM path model -------------------------------
# The model specifies two equations:
#   1. Mediator equation: emo = a * treat + covariates
#   2. Outcome equation:  p_harm = b * emo + cp * treat + covariates
# indirect_effect() of manymome finds the path from the variable names, so the
# labels (a, b, cp) only name the coefficients.
model <- "
  emo ~ a * treat + age + educ_num + female + income
  p_harm ~ b * emo + cp * treat + age + educ_num + female + income
"
fit <- sem(model, data = framing)

# --- manymome: bootstrap indirect effect -------------------------------------
# do_boot() generates bootstrap replicates of the model parameters. This is
# done once, and the stored replicates can be reused for any number of
# indirect_effect() queries, an architectural advantage over robmed, where
# bootstrapping is coupled to a single mediation test.
#
# indirect_effect() traces the path treat -> emo -> p_harm and computes
# the product ab with a percentile bootstrap CI.
set.seed(6489)
t_manymome <- system.time({
  boot_out <- do_boot(fit, R = 10000, seed = 6489)
})
ind <- indirect_effect(x = "treat", y = "p_harm", m = "emo",
  fit = fit, boot_ci = TRUE, boot_out = boot_out)
print(ind)

# --- Extract results ---------------------------------------------------------
pe <- parameterEstimates(fit)
a_est  <- pe$est[pe$lhs == "emo"    & pe$rhs == "treat"]
b_est  <- pe$est[pe$lhs == "p_harm" & pe$rhs == "emo"]
ab_est <- coef(ind)
ab_ci  <- confint(ind)

# --- Pretty-printed results --------------------------------------------------
cat("SEM MEDIATION: MANYMOME + LAVAAN (ML + BOOTSTRAP)\n")
cat(strrep("=", 72), "\n")
cat("Data:     framing (mediation package), n =", nrow(framing), "\n")
cat("Path:     treat -> emo -> p_harm\n")
cat("Method:   ML estimation (lavaan), 10000 bootstrap replicates (manymome)\n")
cat("Seed:     6489\n\n")

cat(sprintf("TIMING:  manymome bootstrap %.1fs\n\n", t_manymome["elapsed"]))
cat("PATH COEFFICIENTS\n")
cat(strrep("-", 72), "\n")
cat(sprintf("  a  path (treat -> emo):      %7.3f\n", a_est))
cat(sprintf("  b  path (emo -> p_harm):     %7.3f\n", b_est))
cat(sprintf("  ab      (indirect effect):   %7.3f  95%% CI [%.3f, %.3f]\n",
  ab_est, ab_ci[1], ab_ci[2]))
cat("\n")

cat("INTERPRETATION\n")
cat(strrep("-", 72), "\n")
cat("  Compare the ML-based indirect effect from lavaan and manymome with\n")
cat("  the MM-estimates of 03_robust_mediation.R. The difference between\n")
cat("  the two shows how much the observations that MM-estimation\n")
cat("  downweights move the ML estimates.\n")
cat("  \n")
cat("  The lavaan + manymome architecture separates model fitting from\n")
cat("  effect computation: one set of bootstrap replicates (do_boot) can\n")
cat("  be queried for multiple indirect paths, conditional effects, or\n")
cat("  standardised effects without re-bootstrapping.\n")
cat(strrep("=", 72), "\n")

if (!is.null(.log_file)) sink()
