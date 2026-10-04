# =============================================================================
# High-Dimensional Mediation: Frequentist (HIMA) vs Bayesian (bama)
# =============================================================================
#
# Compares: HIMA, bama
#
# PURPOSE
# -------
# When the number of potential mediators is large (e.g., genomic, epigenomic,
# or proteomic data), standard single-mediator methods are inadequate. This
# script compares two approaches to high-dimensional mediation: a frequentist
# penalized regression method (HIMA) and a Bayesian shrinkage method (bama).
#
# SIMULATED DATA
# --------------
# n = 300 observations, p = 200 candidate mediators, of which 5 are truly
# active (nonzero alpha and beta paths). A binary exposure A is generated
# from a logistic model with two confounders C1 and C2. The outcome Y is
# continuous and depends on A, the five active mediators, and the confounders.
#
# True indirect effects (alpha * beta) for the 5 active mediators:
#   M1: 0.48, M2: 0.35, M3: 0.30, M4: 0.20, M5: 0.20
#
# The remaining 195 mediators have zero indirect effects. This setup tests
# each method's ability to distinguish true signals from noise.
#
# PACKAGE OVERVIEW
# ----------------
# HIMA (Zhang et al., 2016)
#   Frequentist penalized regression for variable selection among high-
#   dimensional mediators. Uses MiniMax Concave Penalty (MCP) to select
#   mediators, then applies joint significance testing with Bonferroni
#   correction. The MCP penalty has an oracle property, giving nearly
#   unbiased estimates for strong signals and exact zeros for noise.
#   USE CASE: Fast screening in genomic/epigenomic studies where thousands
#   of CpG sites or genes are candidate mediators. Ideal for discovery-phase
#   analyses where speed matters more than full posterior inference.
#   Example: An EWAS study testing which DNA methylation sites mediate the
#   effect of smoking on lung cancer risk, screening 450K CpG sites.
#
# bama (Song et al., 2020)
#   Bayesian Sparse Linear Mixed Model (BSLMM) for high-dimensional
#   mediation. Uses a continuous shrinkage prior that produces posterior
#   inclusion probabilities (PIPs) and credible intervals. Requires MCMC.
#   USE CASE: When full Bayesian inference is needed, including uncertainty
#   quantification via credible intervals and evidence ranking via PIPs.
#   PIPs provide a continuous measure of evidence (0 to 1), avoiding the
#   binary significant/not-significant decision of HIMA.
#   Example: A proteomics study ranking 500 plasma proteins by their
#   posterior probability of mediating the effect of a drug on patient
#   outcomes, where effect sizes are small and a ranking is more useful
#   than a hard significance cutoff.
#
# KEY COMPARISONS
# ---------------
# 1. Variable selection:
#    HIMA selects mediators via penalised regression + Bonferroni p-values.
#    bama ranks mediators via posterior inclusion probabilities (PIPs).
#    With 5 true mediators among 200, both should identify the strongest
#    signals but may differ on borderline cases (M4, M5 with ab = 0.20).
#
# 2. Speed:
#    HIMA is much faster (under a second here) because it solves a penalised
#    optimisation problem. bama requires MCMC sampling (about a minute here)
#    but provides richer inference. The timing comparison quantifies this trade-off.
#
# 3. Preprocessing:
#    bama's BSLMM prior is scale-dependent, so outcome, mediators, and
#    covariates are standardised and the binary exposure is left unscaled.
#    HIMA (scale = TRUE) standardises the exposure, mediators, and
#    covariates internally but not the outcome. The two sets of estimates
#    are therefore on different scales, and both are converted back to the
#    original scale (per unit of A, in units of Y) before they are compared
#    with the true values.
#
# 4. Output:
#    HIMA returns point estimates and p-values for selected mediators ONLY.
#    bama returns posterior summaries for ALL mediators, allowing ranking.
#
# SIMILARITIES
# ------------
# - Both handle the high-dimensional setting (p can exceed n).
# - Both estimate the indirect effect as the product alpha * beta.
# - Both use the same data-generating process for fair comparison.
# - Both aim to separate true signals from noise mediators.
#
# DIFFERENCES
# -----------
# - Framework: frequentist (HIMA) vs Bayesian (bama).
# - Inference: p-values + Bonferroni (HIMA) vs PIPs + credible intervals.
# - Speed: under a second (HIMA) vs about a minute (bama) in this example.
# - Sparsity: hard thresholding via penalty (HIMA) vs continuous shrinkage.
# - Coverage: HIMA returns only selected mediators, and bama returns all p.
# =============================================================================

# --- Log output --------------------------------------------------------------
.log_file <- {
  f <- grep("--file=", commandArgs(), value = TRUE)
  if (length(f)) sub("\\.R$", ".log", sub("--file=", "", f)) else NULL
}
if (!is.null(.log_file)) sink(.log_file, split = TRUE)

cat(sprintf("Script: %s\n", "02_high_dimensional_mediation.R"))
cat(sprintf("Date:   %s\n", Sys.time()))
cat(strrep("=", 72), "\n\n")

suppressPackageStartupMessages({
  library("HIMA")
  library("bama")
})

# --- Simulated data ----------------------------------------------------------
# Construct a realistic high-dimensional mediation scenario:
# - Binary exposure A depends on confounders C1, C2 via a logistic model
# - 200 mediators: only M1-M5 have nonzero paths from A (alpha) and to Y (beta)
# - Continuous outcome Y depends on A (direct effect = 0.5), the five active
#   mediators (indirect via beta), and the confounders
# - True alpha*beta ranges from 0.48 (M1, strongest) to 0.20 (M4-M5, weakest)
set.seed(6489)
n <- 300; p <- 200
C1 <- rnorm(n); C2 <- rnorm(n)
A <- rbinom(n, 1, plogis(0.3 * C1 + 0.3 * C2))
alpha <- rep(0, p); alpha[1:5] <- c(0.8, 0.7, 0.6, 0.5, 0.5)
beta  <- rep(0, p); beta[1:5]  <- c(0.6, 0.5, 0.5, 0.4, 0.4)
M <- A %*% t(alpha) + outer(C1, rep(0.2, p)) + matrix(rnorm(n * p), n, p)
Y <- as.numeric(0.5 * A + M %*% beta + 0.3 * C1 + 0.3 * C2 + rnorm(n))
pheno <- data.frame(Y = Y, A = A, C1 = C1, C2 = C2)

# --- HIMA: penalized regression (MCP) ----------------------------------------
# For a continuous outcome HIMA screens the mediators by their association
# with the outcome, selects among them with MCP in the outcome model,
# estimates the exposure-mediator path of each selected mediator, and keeps
# those whose joint (max-p) test passes a Bonferroni correction over the
# selected mediators. The MCP penalty shrinks weak signals to exactly zero.
#
# scale = TRUE tells HIMA to standardise the exposure, mediators, and
# covariates internally (not the outcome), so that the penalty treats all
# mediators alike.
t_hima <- system.time({
  hima.out <- hima(Y ~ A + C1 + C2, data.pheno = pheno,
    data.M = M, penalty = "MCP", scale = TRUE)
})
head(hima.out)

# --- bama: Bayesian shrinkage (BSLMM) ----------------------------------------
# bama uses MCMC to sample from the posterior of a Bayesian sparse linear
# mixed model. The BSLMM prior uses a two-component normal mixture: one
# narrow component concentrated near zero and one wide component for true
# effects, providing continuous shrinkage rather than hard thresholding.
#
# IMPORTANT: Standardise outcome, mediators, and covariates because the
# BSLMM prior is scale-dependent. The binary exposure A is left unscaled
# because standardising a 0/1 variable distorts the treatment contrast.
#
# bama's C1 argument holds the covariates of the outcome model and C2 those
# of the mediator model. bama adds no intercept, so each matrix starts with a
# column of ones, as in the package examples. The confounders C1 and C2 enter
# both models.
Y.std <- as.numeric(scale(Y))
M.std <- scale(M)
C.mat.std <- cbind(1, scale(cbind(C1, C2)))
t_bama <- system.time({
  bama.out <- bama(Y = Y.std, A = A, M = M.std,
    C1 = C.mat.std, C2 = C.mat.std,
    method = "BSLMM", burnin = 3000, ndraws = 10000, seed = 6489)
})
s.bama <- summary(bama.out)
print(s.bama)

# summary() reports the mediator-outcome coefficient beta.m with the
# posterior inclusion probability (the share of draws in which both paths
# are in the large-effect component). The indirect effect alpha.a * beta.m is
# summarised from the posterior draws instead.
ab_draws <- bama.out[["alpha.a"]] * bama.out[["beta.m"]]
ab_quantile <- function(prob) {
  mapply(function(j) quantile(ab_draws[, j], prob), seq_len(ncol(ab_draws)))
}

# --- Pretty-printed results --------------------------------------------------
true_ab <- (alpha * beta)[1:5]

# Reorder HIMA results by mediator index for aligned comparison with bama.
# HIMA only returns selected (significant) mediators, so IDs may be a subset.
ids <- as.integer(hima.out$ID)
h_ab <- hima.out[["alpha*beta"]]
h_pv <- hima.out[["p-value"]]
hima_ord <- order(ids)
ids <- ids[hima_ord]; h_ab <- h_ab[hima_ord]; h_pv <- h_pv[hima_ord]

# Back to the original scale (per unit of A, in units of Y). HIMA scaled A
# but not Y, so its alpha*beta is per standard deviation of A. bama received
# the raw A and a standardised Y, so its estimates are in standard deviations
# of Y. The scaling of the mediators cancels in the product alpha*beta.
h_ab <- h_ab / sd(A)
b_est <- colMeans(ab_draws) * sd(Y)
b_lo <- ab_quantile(0.025) * sd(Y)
b_hi <- ab_quantile(0.975) * sd(Y)

# Least squares with the five active mediators known (the oracle model), the
# benchmark for both packages: alpha from each mediator model, and beta from
# one outcome model with the five mediators.
ls_alpha <- mapply(function(j) coef(lm(M[, j] ~ A + C1 + C2))[["A"]], 1:5)
ls_beta <- unname(coef(lm(Y ~ A + M[, 1:5] + C1 + C2))[3:7])
ls_ab <- ls_alpha * ls_beta

cat("HIGH-DIMENSIONAL MEDIATION: HIMA vs BAMA\n")
cat(strrep("=", 82), "\n")
cat("Simulation: n = 300, p = 200 candidate mediators, 5 truly active\n")
cat("Seed: 6489\n\n")

cat(sprintf("TIMING:  HIMA %.1fs  |  bama %.1fs\n",
  t_hima["elapsed"], t_bama["elapsed"]))
cat("  (bama is slower due to MCMC sampling but provides richer inference)\n\n")

cat("COMPARISON TABLE: TRUE vs ESTIMATED INDIRECT EFFECTS (active mediators)\n")
cat("Original scale: per unit of A, in units of Y\n")
cat(strrep("-", 82), "\n")
cat(sprintf("%-8s  %8s  %10s  %12s  %10s  %18s  %6s\n",
  "Mediator", "True ab", "HIMA ab", "HIMA p", "bama est",
  "bama 95% CrI", "PIP"))
cat(strrep("-", 82), "\n")
for (i in 1:5) {
  # For HIMA, check whether this mediator was selected, and mark it if not
  hima_idx <- which(ids == i)
  if (length(hima_idx)) {
    cat(sprintf("M%-7d  %8.2f  %10.3f  %12.2e  %10.3f  [%6.3f, %6.3f]  %6.3f\n",
      i, true_ab[i], h_ab[hima_idx], h_pv[hima_idx],
      b_est[i], b_lo[i], b_hi[i], s.bama$pip[i]))
  } else {
    cat(sprintf("M%-7d  %8.2f  %10s  %12s  %10.3f  [%6.3f, %6.3f]  %6.3f\n",
      i, true_ab[i], "not sel.", "---",
      b_est[i], b_lo[i], b_hi[i], s.bama$pip[i]))
  }
}
cat(strrep("-", 82), "\n\n")

cat("LEAST SQUARES ON THE FIVE ACTIVE MEDIATORS (oracle benchmark)\n")
cat(strrep("-", 72), "\n")
cat(sprintf("%-8s  %8s  %12s  %12s\n", "Mediator", "LS ab", "HIMA - LS", "bama - LS"))
for (i in 1:5) {
  hima_idx <- which(ids == i)
  h_diff <- if (length(hima_idx)) sprintf("%12.3f", h_ab[hima_idx] - ls_ab[i]) else
    sprintf("%12s", "not sel.")
  cat(sprintf("M%-7d  %8.3f  %s  %12.3f\n", i, ls_ab[i], h_diff, b_est[i] - ls_ab[i]))
}
cat(strrep("-", 72), "\n\n")

cat("SUMMARY\n")
cat(strrep("-", 72), "\n")
cat(sprintf("  HIMA detected: %d mediators (Bonferroni p < 0.05)\n", length(ids)))
cat(sprintf("  bama PIP range (true mediators M1-M5): %.3f -- %.3f\n",
  min(s.bama$pip[1:5]), max(s.bama$pip[1:5])))
cat(sprintf("  bama PIP range (noise mediators M6-M200): %.3f -- %.3f\n",
  min(s.bama$pip[6:p]), max(s.bama$pip[6:p])))
cat("\n")

cat("INTERPRETATION\n")
cat(strrep("-", 72), "\n")
cat("  HIMA selects all five active mediators (Bonferroni p < 0.05) and no\n")
cat("  noise mediator. bama gives a PIP above 0.5 to M1, M2, M3, and M5 and\n")
cat("  a PIP of zero to three decimals to every noise mediator. M4, whose\n")
cat("  exposure-mediator path bama shrinks most, has a PIP below 0.1. HIMA\n")
cat("  is faster, and bama reports posterior uncertainty through credible\n")
cat("  intervals and PIPs.\n")
cat(strrep("=", 82), "\n")

if (!is.null(.log_file)) sink()
