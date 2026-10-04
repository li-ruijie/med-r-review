#!/usr/bin/env Rscript
# =============================================================================
# Inventory Table of the Screened Mediation Packages
# =============================================================================
#
# Writes the LaTeX rows of the inventory table (longtable body) from the
# inventory produced by 06_cran_inventory.R: every screened package with its
# scope, role, twelve-month download count, and screening note. Included
# packages come first, alphabetically, then the excluded candidates.
#
# Output: ../results/tables/inventory-rows.tex
#
# Run from code/.
# =============================================================================

latex_escape <- function(x) {
    specials <- c("\\\\" = "\\\\textbackslash{}", "&" = "\\\\&", "%" = "\\\\%", "\\$" = "\\\\$",
                  "#" = "\\\\#", "_" = "\\\\_", "\\{" = "\\\\{", "\\}" = "\\\\}",
                  "~" = "\\\\textasciitilde{}", "\\^" = "\\\\textasciicircum{}")
    Reduce(function(s, k) gsub(k, specials[[k]], s), names(specials), x)
}

# Straight double quotes in the notes become LaTeX opening and closing quotes.
latex_quotes <- function(x) gsub('"([^"]*)"', "``\\1''", x)

# The note starts a new source line and is wrapped, so that no line of the
# generated file exceeds 108 characters.
format_row <- function(r) {
    downloads <- if (is.na(r[["downloads"]])) "--" else format(r[["downloads"]], big.mark = ",")
    note <- latex_quotes(latex_escape(r[["note"]])) |> strwrap(width = 100) |> paste(collapse = "\n    ")
    sprintf("\\textbf{%s} & %s & %s & %s\n    & %s \\\\", latex_escape(r[["package"]]),
            if (r[["screen"]] == "include") r[["scope"]] else "excluded",
            if (r[["screen"]] == "include") r[["role"]] else "--",
            downloads, note)
}

main <- function() {
    inv <- read.csv(file.path("..", "results", "cran-search", "inventory-2026-10-01.csv"),
                    stringsAsFactors = FALSE)
    inv <- inv[order(inv[["screen"]] != "include", tolower(inv[["package"]])), ]
    rows <- Map(function(i) format_row(inv[i, ]), seq_len(nrow(inv))) |> unlist()
    out_dir <- file.path("..", "results", "tables")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    writeLines(rows, file.path(out_dir, "inventory-rows.tex"))
    message(sprintf("%d rows (%d included, %d excluded)", length(rows),
                    sum(inv[["screen"]] == "include"), sum(inv[["screen"]] == "exclude")))
}

main()
