# Mediation Analysis in R: Methods, Packages, and Gaps

Companion code for the review paper "Mediation analysis in R: methods,
packages, and gaps" (submitted to *Statistical Methods in Medical
Research*).

The paper reviews 25 R packages for mediation analysis across six
methodological families and identifies four ecosystem gaps. The scripts
in this repository reproduce the cross-package comparisons discussed in
the paper and its appendix.

## Scripts

All scripts use `set.seed(6489)` and 10,000 bootstrap/simulation
replicates for reproducibility. Each script logs output to a `.log` file
when run via `Rscript`.

| Script | Packages | Description |
|--------|----------|-------------|
| `code/00_install_packages.R` | -- | Installs all CRAN and GitHub dependencies. Run once before the analysis scripts. |
| `code/01_foundational_causal_mediation.R` | mediation, medflex, CMAverse, regmedint | Cross-package comparison on the `framing` dataset (binary outcome). Marginal vs conditional NIE, quasi-Bayesian vs bootstrap vs delta-method CIs, sensitivity analysis (medsens rho, CMAverse E-values). |
| `code/02_high_dimensional_mediation.R` | HIMA, bama | Frequentist penalised regression vs Bayesian shrinkage on simulated data (n = 300, p = 200, 5 active mediators). Variable selection, PIPs, and speed comparison. |
| `code/03_robust_mediation.R` | robmed | MM-estimator bootstrap mediation on the `framing` dataset (continuous outcome `p_harm`). Companion to script 04 for robust vs ML comparison. |
| `code/04_sem_mediation.R` | manymome, lavaan | SEM-based mediation with bootstrap indirect effects on the same path model as script 03. Companion to script 03 for ML vs robust comparison. |

## Requirements

R >= 4.5.0 with the following packages:

- **mediation** (>= 4.5.1)
- **medflex** (>= 0.6-11)
- **CMAverse** (>= 0.1.0) -- installed from GitHub: `remotes::install_github("BS1125/CMAverse")`
- **regmedint** (>= 1.0.2)
- **HIMA** (>= 2.3.3)
- **bama** (>= 1.3.1)
- **robmed** (>= 1.3.0)
- **manymome** (>= 0.3.3)
- **lavaan** (>= 0.6-19)

Run the install script to set up all dependencies at once:

```bash
cd code
Rscript 00_install_packages.R
```

This installs all CRAN packages (including **devtools**) and
**CMAverse** from GitHub. Already-installed packages are skipped.

## Running the scripts

```bash
cd code
Rscript 01_foundational_causal_mediation.R
Rscript 02_high_dimensional_mediation.R
Rscript 03_robust_mediation.R
Rscript 04_sem_mediation.R
```

Each script writes a `.log` file alongside itself (e.g.,
`code/01_foundational_causal_mediation.log`).

## Platform

Results in the paper were produced on a 64-bit Linux platform (AMD Ryzen
9 9950X3D, 32 threads) under R 4.5.0. Results may differ under
alternative package versions or random-number generator states.

## Licence

AGPL-3.0-or-later. See [LICENCE](LICENCE).
