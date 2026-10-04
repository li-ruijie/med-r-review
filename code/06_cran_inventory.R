#!/usr/bin/env Rscript
# =============================================================================
# Inventory of Screened Mediation Packages
# =============================================================================
#
# Joins the CRAN search candidates (05_cran_search.R) with the screening
# decisions in cran-screening-<date>.csv and adds total CRAN downloads from
# the RStudio mirror logs (cranlogs) for the twelve months that end the day
# before the search date.
#
# Output: ../results/cran-search/inventory-<date>.csv
#
# Run from code/.
# =============================================================================

search_date <- as.Date("2026-10-01")
window <- c(seq(search_date, by = "-1 year", length.out = 2L)[2L], search_date - 1L)
cranlogs_url <- "https://cranlogs.r-pkg.org/downloads/total/%s:%s/%s"

fetch_downloads <- function(pkgs) {
    url <- sprintf(cranlogs_url, window[1L], window[2L], paste(pkgs, collapse = ","))
    res <- jsonlite::fromJSON(url)
    data.frame(package = res[["package"]], downloads = res[["downloads"]])
}

main <- function() {
    stamp <- format(search_date)
    out_dir <- file.path("..", "results", "cran-search")
    candidates <- read.csv(file.path(out_dir, paste0("candidates-", stamp, ".csv")))
    screening <- read.csv(paste0("cran-screening-", stamp, ".csv"))

    unscreened <- setdiff(candidates[["package"]], screening[["package"]])
    if (length(unscreened) > 0L) {
        stop("Candidates without a screening decision: ", paste(unscreened, collapse = ", "))
    }

    on_cran <- setdiff(screening[["package"]], c("CMAverse", "medltmle"))
    chunks <- split(on_cran, ceiling(seq_along(on_cran) / 40L))
    downloads <- chunks |> Map(f = fetch_downloads) |> Reduce(f = rbind)

    # Version, date, and title come from the snapshot so that packages added
    # outside the search (lavaan, brms) carry them too.
    db <- readRDS(file.path(out_dir, paste0("cran-db-", stamp, ".rds")))
    cran_fields <- data.frame(
        package = db[["Package"]],
        version = db[["Version"]],
        published = db[["Published"]],
        title = gsub("\\s+", " ", trimws(db[["Title"]]))
    )

    inventory <- screening |>
        merge(cran_fields, by = "package", all.x = TRUE) |>
        merge(candidates[, c("package", "keyword_title", "keyword_description", "task_view")],
              by = "package", all.x = TRUE) |>
        merge(downloads, by = "package", all.x = TRUE)
    inventory <- inventory[order(inventory[["screen"]] != "include", -inventory[["downloads"]],
                                 tolower(inventory[["package"]])), ]
    write.csv(inventory, file.path(out_dir, paste0("inventory-", stamp, ".csv")),
              row.names = FALSE)

    message(sprintf("Download window: %s to %s", window[1L], window[2L]))
    message(sprintf("Screened: %d, included: %d, excluded: %d", nrow(inventory),
                    sum(inventory[["screen"]] == "include"),
                    sum(inventory[["screen"]] == "exclude")))
}

main()
