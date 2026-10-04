#!/usr/bin/env Rscript
# =============================================================================
# CRAN Search for Mediation Packages
# =============================================================================
#
# Snapshots the CRAN package database and keeps every package whose title or
# description mentions mediation ("mediat") or indirect effects. Packages
# named in the mediation paragraph of the CRAN Task View on Causal Inference
# are added as a cross-check, read from a pinned commit of the view's source.
#
# Outputs (written to ../results/cran-search/, stamped with the search date):
#   cran-db-<date>.rds      full CRAN package database snapshot
#   candidates-<date>.csv   keyword and task-view candidates for screening
#   search-<date>.log       search metadata and counts
#
# Run from code/.
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

ctv_commit <- "9630cf2127f347ad3a2b05d9cfce85c149fa1fc5"
ctv_url <- paste0(
    "https://raw.githubusercontent.com/cran-task-views/CausalInference/",
    ctv_commit, "/CausalInference.md"
)
# Word-initial so that "intermediate" and "immediate" do not match.
keyword_regex <- "\\bmediat|\\bindirect[ -]effect"

squish <- function(x) gsub("\\s+", " ", trimws(x))

# Lines of the task-view bullet that starts with "*Mediation analysis*", up to
# the next top-level bullet.
ctv_mediation_packages <- function(lines) {
    start <- grep("^-\\s+\\*Mediation analysis\\*", lines)
    bullets <- grep("^-\\s+", lines)
    end <- min(bullets[bullets > start]) - 1L
    paragraph <- paste(lines[start:end], collapse = " ")
    hits <- regmatches(paragraph, gregexpr('pkg\\("[^"]+"', paragraph))[[1]]
    unique(sub('pkg\\("', "", sub('"$', "", hits)))
}

main <- function() {
    out_dir <- file.path("..", "results", "cran-search")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    stamp <- format(Sys.time(), "%Y-%m-%d", tz = "UTC")

    db <- tools::CRAN_package_db()
    db <- db[!duplicated(db[["Package"]]), ]
    saveRDS(db, file.path(out_dir, paste0("cran-db-", stamp, ".rds")))

    title <- squish(db[["Title"]])
    description <- squish(db[["Description"]])
    in_title <- grepl(keyword_regex, title, ignore.case = TRUE, perl = TRUE)
    in_description <- grepl(keyword_regex, description, ignore.case = TRUE, perl = TRUE)

    ctv_pkgs <- ctv_mediation_packages(readLines(ctv_url, warn = FALSE))
    in_ctv <- db[["Package"]] %in% ctv_pkgs

    keep <- in_title | in_description | in_ctv
    candidates <- data.frame(
        package = db[["Package"]][keep],
        version = db[["Version"]][keep],
        published = db[["Published"]][keep],
        title = title[keep],
        description = description[keep],
        keyword_title = in_title[keep],
        keyword_description = in_description[keep],
        task_view = in_ctv[keep]
    )
    candidates <- candidates[order(tolower(candidates[["package"]])), ]
    write.csv(candidates, file.path(out_dir, paste0("candidates-", stamp, ".csv")),
              row.names = FALSE)

    meta <- c(
        sprintf("Search date (UTC):            %s", stamp),
        sprintf("R version:                    %s", R.version.string),
        sprintf("CRAN mirror:                  %s", getOption("repos")[["CRAN"]]),
        sprintf("Packages in CRAN database:    %d", nrow(db)),
        sprintf("Keyword regex (PCRE, ignore case): %s", keyword_regex),
        sprintf("Keyword hits, title:          %d", sum(in_title)),
        sprintf("Keyword hits, description:    %d", sum(in_description)),
        sprintf("Keyword hits, either field:   %d", sum(in_title | in_description)),
        sprintf("Task view commit:             %s", ctv_commit),
        sprintf("Task view mediation packages: %d (%s)", length(ctv_pkgs),
                paste(sort(ctv_pkgs), collapse = ", ")),
        sprintf("Task view packages not on CRAN: %s",
                paste(setdiff(ctv_pkgs, db[["Package"]]), collapse = ", ")),
        sprintf("Task view packages without keyword hit: %s",
                paste(sort(setdiff(ctv_pkgs, db[["Package"]][in_title | in_description])),
                      collapse = ", ")),
        sprintf("Candidates (union):           %d", nrow(candidates))
    )
    writeLines(meta, file.path(out_dir, paste0("search-", stamp, ".log")))
    writeLines(meta)
}

main()
