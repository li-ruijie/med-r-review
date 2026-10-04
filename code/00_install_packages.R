# =============================================================================
# Install Required Packages
# =============================================================================
#
# Installs all R packages needed by scripts 01-13. Run once before running
# the analysis scripts. Dependencies are installed automatically.
#
# CRAN packages:
#   mediation, medflex, regmedint, HIMA, bama, robmed, manymome, lavaan
#   (code examples, simulation, and selection stability), jsonlite (download
#   counts for the CRAN inventory), bmem, bmemLavaan, and sem (check of the
#   robust options, 13)
#
# Bioconductor package (installed via BiocManager):
#   qvalue (required by HDMT, which HIMA imports, and not on CRAN)
#
# GitHub packages (installed via devtools):
#   CMAverse (BS1125/CMAverse, pinned to the commit used for the results)
# =============================================================================

# Rscript has no CRAN mirror set by default.
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  "devtools",
  "mediation",
  "medflex",
  "regmedint",
  "HIMA",
  "bama",
  "robmed",
  "manymome",
  "lavaan",
  "jsonlite",
  "bmem",
  "bmemLavaan",
  "sem"
)

if (!requireNamespace("qvalue", quietly = TRUE)) {
  cat("Installing qvalue from Bioconductor...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("qvalue")
}

to_install <- cran_pkgs[!vapply(cran_pkgs, requireNamespace,
  logical(1), quietly = TRUE)]

if (length(to_install)) {
  cat("Installing CRAN packages:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install)
} else {
  cat("All CRAN packages already installed.\n")
}

if (!requireNamespace("CMAverse", quietly = TRUE)) {
  cat("Installing CMAverse from GitHub...\n")
  devtools::install_github("BS1125/CMAverse@b2ce0598ed362ddfa0db8fff5486ee43d3e73f54")
} else {
  cat("CMAverse already installed.\n")
}
