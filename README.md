# Mediation Analysis in R: A Review of Packages and Gaps

Code and results for a review of R packages for causal mediation analysis. The repository holds a
reproducible search of CRAN for mediation packages, cross-package code examples on shared datasets, a
simulation study of natural effect estimators, a selection-stability analysis of high-dimensional
mediation, and two checks of package specifications.

## Contents

| Path | Purpose |
|------|---------|
| `code/00_install_packages.R` | CRAN, Bioconductor, and GitHub dependencies (CMAverse pinned) |
| `code/01_foundational_causal_mediation.R` | mediation, medflex, CMAverse, and regmedint on `framing` |
| `code/02_high_dimensional_mediation.R` | HIMA and bama on simulated data (n = 300, p = 200) |
| `code/03_robust_mediation.R` | robmed MM-estimation on the `framing` data |
| `code/04_sem_mediation.R` | manymome and lavaan on the same path model as script 03 |
| `code/05_cran_search.R` | CRAN snapshot of 1 October 2026 and keyword and task-view candidates |
| `code/06_cran_inventory.R` | Joins the screening decisions with 12-month download counts |
| `code/07_simulation.R` | Simulation of nine estimator and interval combinations over 48 scenarios |
| `code/08_simulation_summary.R` | Performance measures with Monte Carlo standard errors |
| `code/09_selection_stability.R` | HIMA selection frequencies over half-samples |
| `code/10_inventory_table.R` | LaTeX rows of the inventory table |
| `code/11_simulation_tables.R` | LaTeX rows of the simulation tables and the summary statistics |
| `code/12_medflex_check.R` | Large-sample check of the medflex specification for a continuous exposure |
| `code/13_bmem_check.R` | Robust options of bmem and bmemLavaan with an exposure-mediator product term |
| `code/parallel_setup.R` | Parallel helpers (fork on Linux, PSOCK on Windows) and seeds |
| `code/run_simulation.sh` | Launches script 07 detached, with BLAS pinned to one thread |
| `code/cran-screening-2026-10-01.csv` | Screening of 133 search candidates and four added packages |
| `results/cran-search/` | CRAN database snapshot, candidates, inventory, and search log |
| `results/simulation/full/` | Per-scenario results, true values, measures, and summary |
| `results/stability/` | Selection frequencies and log |
| `results/medflex-check/`, `results/bmem-check/` | Outputs of scripts 12 and 13 |
| `results/tables/` | Table rows written by scripts 10 and 11 |

Scripts 01 to 04 write a `.log` beside themselves when run with `Rscript`, the logs of scripts 12 and 13
come from redirecting their output, and the repository includes all six logs.

## Requirements

R 4.5.0 or later. Install the dependencies once, from `code/`:

```bash
cd code
Rscript 00_install_packages.R
```

## Running the scripts

Run every script from `code/`, since each resolves its paths from there. Scripts 07, 09, 12, and 13
run their parallel work through `parallel_setup.R`, with one seed per task drawn from the master seed
6489, so the results do not depend on the number of workers.

```bash
cd code
Rscript 01_foundational_causal_mediation.R
Rscript 05_cran_search.R
Rscript 06_cran_inventory.R
Rscript 10_inventory_table.R
OPENBLAS_NUM_THREADS=1 Rscript 09_selection_stability.R
cd ..
bash code/run_simulation.sh full
cd code
Rscript 08_simulation_summary.R full
Rscript 11_simulation_tables.R
```

Scripts 02 to 04 run like script 01, and scripts 12 and 13 as
`OPENBLAS_NUM_THREADS=1 Rscript 12_medflex_check.R > 12_medflex_check.log 2>&1`. The full simulation
takes several hours on 31 cores. `run_simulation.sh` returns at once and the simulation runs detached,
so run script 08 after `results/simulation/full/run.log` ends with "finished full".
`05_cran_search.R` downloads the current CRAN database, so a rerun at a later date finds a different
snapshot than the one stored in `results/cran-search/`.

## Platform

The results were produced under R 4.5.0 on 64-bit Linux with a 16-core processor (32 hardware
threads). Other package versions or random-number generator states may give different results.

## Licence

AGPL-3.0-or-later. See [LICENCE](LICENCE).
