# =============================================================================
# Install Required Packages
# =============================================================================
#
# Installs all R packages needed by scripts 01-04. Run once before running
# the analysis scripts. Dependencies are installed automatically.
#
# CRAN packages:
#   mediation, medflex, regmedint, HIMA, bama, robmed, manymome, lavaan
#
# GitHub packages (installed via devtools):
#   CMAverse (BS1125/CMAverse)
# =============================================================================

cran_pkgs <- c(
  "devtools",
  "mediation",
  "medflex",
  "regmedint",
  "HIMA",
  "bama",
  "robmed",
  "manymome",
  "lavaan"
)

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
  devtools::install_github("BS1125/CMAverse")
} else {
  cat("CMAverse already installed.\n")
}
