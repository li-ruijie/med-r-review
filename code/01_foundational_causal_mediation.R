# =============================================================================
# Cross-Package Comparison of Foundational Causal Mediation Packages
# =============================================================================
#
# Compares: mediation, medflex, CMAverse, regmedint
#
# PURPOSE
# -------
# These four packages all implement causal mediation analysis under the
# counterfactual/potential-outcomes framework, estimating natural indirect
# effects (NIE) and natural direct effects (NDE). This script applies all four
# to the same data and covariates to show how they differ in estimation
# strategy, effect scale, and uncertainty quantification.
#
# DATASET
# -------
# Brader et al.'s "framing" experiment (bundled with the mediation package).
# A randomised experiment crossing two factors of a news story on immigration,
# its tone (positive or negative) and the ethnicity of the featured immigrant
# (European or Latino). treat is the product of the two, so treat = 1 marks
# the one cell in which both equal 1 (the negative story about a Latino
# immigrant, 68 subjects), and treat = 0 the other three cells (197 subjects).
# The mediator is emotional response (emo, continuous), and the outcome is
# whether the subject sent a congressional message supporting stricter
# immigration policy (cong_mesg, binary).
#
#   Path:       treat -> emo -> cong_mesg
#   Covariates: age, education, gender (female), income
#
# PACKAGE OVERVIEW
# ----------------
# mediation (Tingley et al., 2014)
#   Quasi-Bayesian Monte Carlo simulation for CIs. Requires separate mediator
#   and outcome model objects. Reports on the probability-difference scale.
#   USE CASE: Default choice for simple single-mediator analysis. Well-
#   documented, widely cited, and approachable for applied researchers.
#   Example: A political scientist testing whether framing effects on policy
#   preferences operate through emotional arousal.
#
# medflex (Steen et al., 2017)
#   Natural effect models via weighting (neWeight) or imputation (neImpute),
#   with bootstrap standard errors by default and sandwich standard errors as
#   an option. Conditional estimates come from a natural effect model with
#   covariates, and marginal ones from a model without them, weighted by an
#   exposure model (xFit) unless the exposure is randomised. A binary exposure
#   must be coded as a factor, since medflex treats a numeric exposure as
#   continuous.
#   USE CASE: Natural effect models stated explicitly, or a comparison of the
#   weighting and imputation approaches, with marginal and conditional
#   estimates from a single framework.
#   Example: An epidemiologist comparing population-average and conditional
#   mediation effects in an observational study, checking robustness across
#   both estimation strategies.
#
# CMAverse (Shi et al., 2021)
#   VanderWeele's regression-based approach with full 4-way decomposition
#   (CDE, intref, intmed, PIE). Supports closed-form (delta method) and
#   g-computation (bootstrap). Built-in sensitivity analysis via E-values.
#   USE CASE: Full decomposition of exposure-mediator interaction, combined
#   with integrated sensitivity analysis, when both closed-form and
#   g-computation estimates and an assessment of unmeasured confounding are
#   needed.
#   Example: A clinical researcher decomposing a treatment effect into
#   mediated, interaction, and direct components while quantifying how
#   strong unmeasured confounding would need to be to explain the result.
#
# regmedint (Yoshida et al., 2026)
#   Valeri & VanderWeele's closed-form regression-based formulas with
#   exposure-mediator interaction. Conditional effects at user-specified
#   covariate values. Delta-method CIs only.
#   USE CASE: Precise regression-based mediation at specific covariate
#   profiles. Clean, minimal implementation of the parametric approach.
#   Example: A pharmacoepidemiologist estimating mediation effects for a
#   patient with specific demographic characteristics (age, sex, BMI).
#
# KEY COMPARISONS
# ---------------
# 1. Marginal vs conditional effects:
#    Non-collapsibility of the odds ratio means marginal (population-average)
#    and conditional (at covariate means) estimates differ even absent
#    confounding. This script estimates both to demonstrate the distinction.
#    The closed-form conditional OR of CMAverse and regmedint relies on the
#    rare-outcome approximation, so the conditional OR is also computed
#    exactly at the covariate means by numerical integration over the mediator.
#
# 2. CI methods: quasi-Bayesian (mediation), bootstrap (medflex, CMAverse
#    g-computation, hand-coded g-computation), delta method (CMAverse closed
#    form, regmedint).
#
# 3. Sensitivity analysis: medsens (mediation) varies the correlation rho
#    between the error terms of the two models and needs, for a binary
#    outcome, a probit outcome model without the exposure-mediator
#    interaction. cmsens (CMAverse) reports E-values.
#
# 4. Agreement: estimators that share an estimand and a formula agree
#    closely, and the comparison separates differences of estimand (marginal
#    or conditional) from differences of estimator.
#
# SIMILARITIES
# ------------
# - All four estimate natural indirect and direct effects under the same
#   identification assumptions (sequential ignorability).
# - All use the framing data with the same covariates, and every outcome
#   model includes the treat:emo interaction (medflex weighting through its
#   saturated natural effect model).
# - All handle binary outcomes via logistic regression.
# - When using closed-form parametric methods, CMAverse and regmedint
#   produce identical point estimates (same underlying formulas).
#
# DIFFERENCES
# -----------
# - Scale: mediation reports probability differences. The others report
#   odds ratios.
# - Marginal vs conditional: medflex and CMAverse can do both, and regmedint
#   is conditional only (the marginal effect needs g-computation by hand).
# - CI method: quasi-Bayes or bootstrap (mediation), bootstrap or sandwich (medflex),
#   delta or bootstrap (CMAverse), delta (regmedint).
# - Decomposition: CMAverse offers 4-way, and the others offer 2-way.
# - Sensitivity: medsens needs a probit outcome model without the
#   interaction for a binary outcome, and the CMAverse E-values need no
#   refit, although cmsens() does not cover survival outcomes.
# =============================================================================

# --- Log output --------------------------------------------------------------
# When run via Rscript, output is written to a .log file alongside this
# script. When the script is sourced interactively, output goes to the
# console only.
.log_file <- {
  f <- grep("--file=", commandArgs(), value = TRUE)
  if (length(f)) sub("\\.R$", ".log", sub("--file=", "", f)) else NULL
}
if (!is.null(.log_file)) sink(.log_file, split = TRUE)

cat(sprintf("Script: %s\n", "01_foundational_causal_mediation.R"))
cat(sprintf("Date:   %s\n", Sys.time()))
cat(strrep("=", 72), "\n\n")

suppressPackageStartupMessages({
  library("mediation")
  library("medflex")
  library("CMAverse")
  library("regmedint")
})

# --- Data preparation --------------------------------------------------------
# The framing dataset is a randomised experiment. We recode gender to a binary
# indicator and education to numeric for use as covariates in the models.
data("framing", package = "mediation")
framing$female <- as.numeric(framing$gender == "female")
framing$educ_num <- as.numeric(framing$educ)
# medflex treats a numeric exposure as continuous, so it gets a factor copy.
framing$treatf <- factor(framing$treat)

# Fit shared mediator and outcome models. The mediator model is linear
# (emo is continuous). The outcome model is logistic (cong_mesg is binary)
# and includes the treat:emo interaction, as the other packages below do.
med.fit <- lm(emo ~ treat + age + educ_num + female + income,
  data = framing)
out.fit <- glm(cong_mesg ~ emo * treat + age + educ_num + female + income,
  data = framing, family = binomial("logit"))

# --- mediation: quasi-Bayesian Monte Carlo -----------------------------------
# The mediation package draws from the sampling distribution of the model
# parameters via simulation and computes the ACME for each draw. This yields
# CIs on the probability-difference scale, not odds ratios.
#
# SIMILARITY: Like CMAverse and regmedint, it uses a two-model approach
# (mediator model + outcome model). Unlike medflex, it does not expand the
# dataset through counterfactual imputation/weighting.
#
# DIFFERENCE: Only package here that reports on the probability-difference
# scale by default, while the others report on the odds-ratio scale.
set.seed(6489)
t_mediation <- system.time({
  med.out <- mediate(med.fit, out.fit, treat = "treat",
    mediator = "emo", sims = 10000)
})
summary(med.out)

# --- medflex: natural effect models (weighting) ------------------------------
# medflex expands the data with two hypothetical exposure values, treatf0 for
# the direct path and treatf1 for the indirect path. The weighting approach
# (neWeight) weights each row by the ratio of the mediator densities under
# the two values, taken from the mediator model, and the natural effect model
# (neModel) is a weighted GLM on the expanded data. With the binary exposure
# coded as a factor, the model in treatf0 * treatf1 is saturated, and its
# treatf11 coefficient is the log OR of the pure natural indirect effect. The
# exposure is randomised, so the marginal model needs no exposure model.
#
# SIMILARITY with imputation below: Both target the same natural effects on
# the log-odds scale, and both use the default bootstrap standard errors,
# which refit the working model on each resample.
#
# DIFFERENCE from mediation: medflex explicitly constructs counterfactual
# data, making the natural effect model a standard GLM on expanded data.
med.glm <- glm(emo ~ treatf + age + educ_num + female + income,
  data = framing, family = gaussian())
expData.wt <- neWeight(med.glm)
set.seed(6489)
nem.wt <- neModel(cong_mesg ~ treatf0 * treatf1,
  family = binomial("logit"), expData = expData.wt, nBoot = 10000)
summary(neEffdecomp(nem.wt))

# --- medflex: natural effect models (imputation) -----------------------------
# The imputation approach (neImpute) imputes the nested counterfactual
# outcomes Y(a, M(a*)) from the outcome model, whose formula must list the
# exposure before the mediator. The outcome model includes the treat:emo
# interaction, as CMAverse (EMint = TRUE) and regmedint (interaction = TRUE)
# do.
out.reord <- glm(cong_mesg ~ treatf * emo + age + educ_num + female + income,
  data = framing, family = binomial("logit"))
expData.imp <- neImpute(out.reord)
set.seed(6489)
nem.imp <- neModel(cong_mesg ~ treatf0 * treatf1,
  family = binomial("logit"), expData = expData.imp, nBoot = 10000)
summary(neEffdecomp(nem.imp))

# --- CMAverse: conditional (delta-method CIs) --------------------------------
# With estimation = "paramfunc" and inference = "delta", CMAverse computes
# closed-form conditional effects using the delta method for CIs.
#
# SIMILARITY with regmedint: Both implement the same Valeri & VanderWeele
# closed-form formulas, so their conditional estimates should be identical.
#
# DIFFERENCE from regmedint: CMAverse also supports g-computation and has
# built-in sensitivity analysis via E-values (cmsens).
set.seed(6489)
cma.out <- cmest(data = framing, model = "rb",
  outcome = "cong_mesg", exposure = "treat", mediator = "emo",
  basec = c("age", "educ_num", "female", "income"),
  EMint = TRUE, mreg = list("linear"), yreg = "logistic",
  astar = 0, a = 1, mval = list(0),
  estimation = "paramfunc", inference = "delta")
summary(cma.out)

# --- regmedint: conditional (delta-method CIs) --------------------------------
# regmedint is a clean implementation of Valeri & VanderWeele's formulas.
# Conditional effects are evaluated at user-specified covariate values
# (here, the sample means). Delta-method CIs.
#
# SIMILARITY with CMAverse (paramfunc): Identical formulas, so point
# estimates match exactly when covariates are set to the same values.
reg.out <- regmedint(data = framing,
  yvar = "cong_mesg", avar = "treat", mvar = "emo",
  cvar = c("age", "educ_num", "female", "income"),
  a0 = 0, a1 = 1, m_cde = 0,
  c_cond = colMeans(framing[, c("age", "educ_num", "female", "income")]),
  mreg = "linear", yreg = "logistic", interaction = TRUE)
reg.s <- summary(reg.out)
print(reg.s)

# --- medflex: conditional (covariates in neModel) ----------------------------
# Including covariates in the neModel makes the estimates conditional on
# those covariates, so the same data expansion supports both marginal and
# conditional estimation. The covariates enter as main effects only, so the
# model assumes the same conditional odds ratio at every covariate value.
set.seed(6489)
nem.wt.cond <- neModel(cong_mesg ~ treatf0 * treatf1 +
  age + educ_num + female + income,
  family = binomial("logit"), expData = expData.wt, nBoot = 10000)

set.seed(6489)
nem.imp.cond <- neModel(cong_mesg ~ treatf0 * treatf1 +
  age + educ_num + female + income,
  family = binomial("logit"), expData = expData.imp, nBoot = 10000)

# --- CMAverse: marginal (g-computation, bootstrap) ---------------------------
# G-computation (estimation = "imputation") produces marginal effects by
# averaging over the empirical covariate distribution. Bootstrap CIs account
# for the additional uncertainty from the averaging step.
#
# DIFFERENCE from conditional above: g-computation marginalises over
# covariates, giving population-average effects. On the OR scale, these
# differ from conditional effects due to non-collapsibility.
set.seed(6489)
cma.marg <- cmest(data = framing, model = "rb",
  outcome = "cong_mesg", exposure = "treat", mediator = "emo",
  basec = c("age", "educ_num", "female", "income"),
  EMint = TRUE, mreg = list("linear"), yreg = "logistic",
  astar = 0, a = 1, mval = list(0),
  estimation = "imputation", inference = "bootstrap", nboot = 10000)
cma.marg.s <- summary(cma.marg)
print(cma.marg.s)

# --- Hand-coded g-computation from the regmedint models (bootstrap CI) -------
# regmedint reports conditional effects only, so the marginal effect is
# computed by hand from the mediator and outcome models that regmedint fits.
# For each bootstrap replicate:
#   1. Resample the data with replacement
#   2. Fit mediator and outcome models
#   3. Simulate counterfactual mediator values under treat=1 and treat=0
#   4. Predict outcomes under treat=0 with each mediator distribution
#   5. Compute the marginal OR from the averaged predictions
#   6. Compute the conditional OR at the covariate means exactly
#
# This gives a marginal PNIE comparable to CMAverse g-computation and medflex.
#
# The conditional OR at the covariate means c is the OR of
# P(Y(0, M(1)) = 1 | c) and P(Y(0, M(0)) = 1 | c), each the logistic outcome
# model integrated over the normal mediator distribution. The closed form of
# CMAverse and regmedint approximates this OR when the outcome is rare, and
# the exact value shows the error of the approximation at the prevalence of
# cong_mesg. With treat = 0 the treat:emo term drops out of the outcome model.
covs <- c("age", "educ_num", "female", "income")
c_bar <- colMeans(framing[, covs])
cond_or <- function(m_fit, y_fit) {
  bm <- coef(m_fit); by <- coef(y_fit)
  lp0 <- by[["(Intercept)"]] + sum(by[covs] * c_bar)
  mu <- function(a) bm[["(Intercept)"]] + bm[["treat"]] * a + sum(bm[covs] * c_bar)
  q <- function(a) integrate(function(m) plogis(lp0 + by[["emo"]] * m) *
    dnorm(m, mu(a), sigma(m_fit)), -Inf, Inf)$value
  (q(1) / (1 - q(1))) / (q(0) / (1 - q(0)))
}
set.seed(6489)
B <- 10000; n <- nrow(framing)
boot_p1 <- numeric(B); boot_p0 <- numeric(B); boot_ors <- numeric(B)
boot_cond_ors <- numeric(B)
t_gcomp <- system.time({
for (b in seq_len(B)) {
  idx <- sample(n, replace = TRUE)
  d_boot <- framing[idx, ]
  m_boot <- lm(emo ~ treat + age + educ_num + female + income, data = d_boot)
  y_boot <- glm(cong_mesg ~ emo * treat + age + educ_num + female + income,
    data = d_boot, family = binomial("logit"))
  M_a1 <- predict(m_boot, newdata = transform(d_boot, treat = 1)) +
    rnorm(n, 0, sigma(m_boot))
  M_a0 <- predict(m_boot, newdata = transform(d_boot, treat = 0)) +
    rnorm(n, 0, sigma(m_boot))
  boot_p1[b] <- mean(predict(y_boot,
    newdata = transform(d_boot, treat = 0, emo = M_a1), type = "response"))
  boot_p0[b] <- mean(predict(y_boot,
    newdata = transform(d_boot, treat = 0, emo = M_a0), type = "response"))
  boot_ors[b] <- (boot_p1[b] / (1 - boot_p1[b])) /
    (boot_p0[b] / (1 - boot_p0[b]))
  boot_cond_ors[b] <- cond_or(m_boot, y_boot)
}
})
# The exact conditional OR needs no Monte Carlo draws, so its point estimate
# comes from the models fitted to the original sample.
exact_cf_or <- cond_or(
  lm(emo ~ treat + age + educ_num + female + income, data = framing),
  glm(cong_mesg ~ emo * treat + age + educ_num + female + income,
    data = framing, family = binomial("logit")))

# --- Sensitivity analysis: medsens (requires probit) -------------------------
# medsens varies the correlation (rho) between the error terms of the
# mediator and outcome models, and reports the rho at which the ACME crosses
# zero.
#
# NOTE: For a binary outcome, medsens requires a probit link (not logit) and
# an outcome model without the treat:emo interaction, so the outcome model
# is refitted without it. This is a practical limitation compared to
# CMAverse's E-value approach, which needs no refit for a logistic outcome model.
out.probit <- glm(cong_mesg ~ emo + treat + age + educ_num +
  female + income, data = framing, family = binomial("probit"))
set.seed(6489)
med.probit <- mediate(med.fit, out.probit, treat = "treat",
  mediator = "emo", sims = 10000)
sens.out <- medsens(med.probit, rho.by = 0.05, sims = 10000)

# --- Sensitivity analysis: CMAverse E-values ---------------------------------
# E-values quantify the minimum strength of unmeasured confounding (on the
# risk-ratio scale) needed to reduce the observed effect to the null. Unlike
# medsens, this does not require refitting the model with a different link.
# For a logistic outcome model cmsens() applies the risk-ratio formula to the
# odds ratio itself, which presumes a rare outcome, so at the prevalence of
# cong_mesg (0.33) its E-values overstate the confounding strength needed.
set.seed(6489)
cma.boot <- cmest(data = framing, model = "rb",
  outcome = "cong_mesg", exposure = "treat", mediator = "emo",
  basec = c("age", "educ_num", "female", "income"),
  EMint = TRUE, mreg = list("linear"), yreg = "logistic",
  astar = 0, a = 1, mval = list(0),
  estimation = "imputation", inference = "bootstrap", nboot = 10000)
evals <- cmsens(cma.boot, sens = "uc")

# --- Extract results (before printing, to avoid progress bar noise) ----------
# neEffdecomp() reports the pure natural indirect effect on the log OR scale,
# with the default bootstrap interval of neModel().
or_fn <- function(m) exp(coef(neEffdecomp(m))[["pure indirect effect"]])
ci_fn <- function(m) {
  invisible(capture.output(v <- confint(neEffdecomp(m))["pure indirect effect", ]))
  exp(v)
}
reg.pnie <- coef(reg.s)["pnie", ]
s <- summary(med.out)
ss <- summary(sens.out)

# Marginal estimates
marg_wt_or  <- or_fn(nem.wt);       marg_wt_ci  <- ci_fn(nem.wt)
marg_imp_or <- or_fn(nem.imp);      marg_imp_ci <- ci_fn(nem.imp)
cma.rpnie <- cma.marg.s$summarydf["Rpnie", ]
cma_or  <- cma.rpnie[["Estimate"]]
cma_cil <- cma.rpnie[["95% CIL"]]; cma_ciu <- cma.rpnie[["95% CIU"]]
p1_bar <- mean(boot_p1); p0_bar <- mean(boot_p0)
reg_gcomp_or  <- (p1_bar / (1 - p1_bar)) / (p0_bar / (1 - p0_bar))
reg_gcomp_cil <- quantile(boot_ors, 0.025)
reg_gcomp_ciu <- quantile(boot_ors, 0.975)

# Conditional estimates
cond_wt_or  <- or_fn(nem.wt.cond);  cond_wt_ci  <- ci_fn(nem.wt.cond)
cond_imp_or <- or_fn(nem.imp.cond); cond_imp_ci <- ci_fn(nem.imp.cond)
cma_cond.s <- summary(cma.out)
cma_cond_rpnie <- cma_cond.s$summarydf["Rpnie", ]
cma_cf_or  <- cma_cond_rpnie[["Estimate"]]
cma_cf_cil <- cma_cond_rpnie[["95% CIL"]]; cma_cf_ciu <- cma_cond_rpnie[["95% CIU"]]
reg_cf_or  <- exp(reg.pnie["est"])
reg_cf_cil <- exp(reg.pnie["lower"]); reg_cf_ciu <- exp(reg.pnie["upper"])
exact_cf_cil <- quantile(boot_cond_ors, 0.025)
exact_cf_ciu <- quantile(boot_cond_ors, 0.975)
prevalence <- mean(framing$cong_mesg)

# --- Pretty-printed results --------------------------------------------------
cat("CROSS-PACKAGE COMPARISON: FOUNDATIONAL CAUSAL MEDIATION\n")
cat(strrep("=", 72), "\n")
cat("Packages: mediation, medflex, CMAverse, regmedint\n")
cat("Data:     framing (mediation package), n =", nrow(framing), "\n")
cat("Path:     treat -> emo -> cong_mesg\n")
cat("Seed:     6489, B = 10000 bootstrap/simulation replicates\n\n")

cat("MARGINAL (population-average) natural indirect effect (OR scale)\n")
cat(strrep("-", 72), "\n")
cat(sprintf("  %-28s  OR = %.3f  95%% CI [%.3f, %.3f]  %s\n",
  "medflex weighting",     marg_wt_or,  marg_wt_ci[1],  marg_wt_ci[2],
  "Bootstrap"))
cat(sprintf("  %-28s  OR = %.3f  95%% CI [%.3f, %.3f]  %s\n",
  "medflex imputation",    marg_imp_or, marg_imp_ci[1], marg_imp_ci[2],
  "Bootstrap"))
cat(sprintf("  %-28s  OR = %.3f  95%% CI [%.3f, %.3f]  %s\n",
  "CMAverse g-computation", cma_or,     cma_cil,        cma_ciu,
  "Bootstrap"))
cat(sprintf("  %-28s  OR = %.3f  95%% CI [%.3f, %.3f]  %s\n",
  "Hand-coded g-computation", reg_gcomp_or, reg_gcomp_cil, reg_gcomp_ciu,
  "Bootstrap"))
cat("\n")

cat("CONDITIONAL (at covariate means) natural indirect effect (OR scale)\n")
cat(strrep("-", 72), "\n")
cat(sprintf("  %-28s  OR = %.3f  95%% CI [%.3f, %.3f]  %s\n",
  "medflex weighting",     cond_wt_or,  cond_wt_ci[1],  cond_wt_ci[2],
  "Bootstrap"))
cat(sprintf("  %-28s  OR = %.3f  95%% CI [%.3f, %.3f]  %s\n",
  "medflex imputation",    cond_imp_or, cond_imp_ci[1], cond_imp_ci[2],
  "Bootstrap"))
cat(sprintf("  %-28s  OR = %.3f  95%% CI [%.3f, %.3f]  %s\n",
  "CMAverse closed-form",  cma_cf_or,   cma_cf_cil,     cma_cf_ciu,
  "Delta method"))
cat(sprintf("  %-28s  OR = %.3f  95%% CI [%.3f, %.3f]  %s\n",
  "regmedint closed-form",  reg_cf_or,  reg_cf_cil,     reg_cf_ciu,
  "Delta method"))
cat(sprintf("  %-28s  OR = %.3f  95%% CI [%.3f, %.3f]  %s\n",
  "Hand-coded exact", exact_cf_or, exact_cf_cil, exact_cf_ciu,
  "Bootstrap"))
cat(sprintf("  Prevalence of cong_mesg: %.3f\n", prevalence))
cat("\n")

cat("PROBABILITY-DIFFERENCE SCALE\n")
cat(strrep("-", 72), "\n")
cat(sprintf("  %-28s  ACME = %.3f  95%% CI [%.3f, %.3f]  %s\n",
  "mediation quasi-Bayes", s$d0, s$d0.ci[1], s$d0.ci[2], "Simulation"))
cat("\n")

cat("TIMING\n")
cat(strrep("-", 72), "\n")
cat(sprintf("  mediation (quasi-Bayes, 10000 sims): %.1fs\n", t_mediation["elapsed"]))
cat(sprintf("  Hand-coded g-computation (10000 boot): %.1fs\n", t_gcomp["elapsed"]))
cat("\n")

cat("SENSITIVITY ANALYSIS\n")
cat(strrep("-", 72), "\n")
print(ss)
cat(sprintf("  medsens rho at ACME = 0:  %.2f\n", ss$err.cr.d[1]))
cat("  CMAverse E-values:\n")
print(evals)
cat("\n")

cat("INTERPRETATION\n")
cat(strrep("-", 72), "\n")
cat(sprintf(paste0("  CMAverse and regmedint share the closed-form formulas and give the\n",
  "  same conditional OR (%.3f and %.3f).\n"), cma_cf_or, reg_cf_or))
cat(sprintf(paste0("  The closed form relies on the rare-outcome approximation. At an\n",
  "  outcome prevalence of %.2f the exact conditional OR is %.3f.\n"), prevalence, exact_cf_or))
cat(sprintf(paste0("  CMAverse g-computation and the hand-coded g-computation differ by\n",
  "  %.1f%% on the marginal scale.\n"), 100 * abs(cma_or / reg_gcomp_or - 1)))
cat("  The mediation package reports on the probability-difference scale, so\n")
cat("  its estimate is not comparable with the odds ratios without\n")
cat("  transformation.\n")
cat(strrep("=", 72), "\n")

if (!is.null(.log_file)) sink()
