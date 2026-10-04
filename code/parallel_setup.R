# ══════════════════════════════════════════════════════════════
#  parallel_setup.R, a portable parallel backend for R
#
#  Works on Linux/macOS (fork via mcmapply) and Windows (PSOCK
#  via clusterMap). Source this file at the top of any script
#  that needs parallelism.
#
#  Primary API: par_map (parallel Map over several iterables)
#               par_values (unwrap par_map records, stop on error)
#               par_mapply (par_map + simplify)
#               par_map_chunked (one data-frame chunk per worker)
#               par_seeds (generate one seed per element)
#               par_benchmark (sequential versus parallel timings)
#               print.par_cfg (one-line summary of a config)
#
#  lapply/sapply and their parallel wrappers are intentionally
#  excluded. Use Map/mapply (and par_map/par_mapply) exclusively.
#
#  Style: functional (no $, no side-effect assignment, no for).
#  All extraction via [["name"]]. All modification via
#  functional replacement forms (`[[<-`, `class<-`, etc.).
#
#  Usage:
#    source("parallel_setup.R")
#    cfg <- setup_parallel(n_cores = 4)
#    seeds <- par_seeds(6489L, length(items))
#    results <- par_map(cfg, function(s, x) {
#        set.seed(s)
#        # ... work ...
#    }, seeds, items)
#    values <- par_values(results)   # stops if any element failed
#    cleanup_parallel(cfg)
#
#  Worker initialisation (PSOCK + fork):
#    cfg <- setup_parallel(
#        n_cores = 4,
#        packages = c("lavaan", "semTools"),
#        init_worker = function() {
#            Sys.setenv(OPENBLAS_NUM_THREADS = "1")
#        }
#    )
#    # On PSOCK: init_worker runs via clusterCall after packages
#    # On fork:  init_worker runs once in the parent before forking
#    #           (forked children inherit the parent's environment)
# ══════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(parallel))

# ── Backend detection ─────────────────────────────────────────

.is_fork_safe <- function() {

    # fork works on Unix-like systems but NOT inside RStudio on
    # macOS (risk of GUI deadlocks). WSL counts as Linux here.
    os <- .Platform[["OS.type"]]
    if (os != "unix") return(FALSE)

    # Respect the user/system override
    env_val <- Sys.getenv("R_PARALLEL_BACKEND", unset = "")
    if (nzchar(env_val)) return(tolower(env_val) == "fork")

    # Disable fork inside RStudio on macOS to avoid crashes
    if (Sys.info()[["sysname"]] == "Darwin" &&
        nzchar(Sys.getenv("RSTUDIO", unset = ""))) {
        return(FALSE)
    }

    TRUE
}

# ── Seed generation ──────────────────────────────────────────

#' Generate n unique seeds, one per element, from a master seed.
#'
#' Uses the default RNG (Mersenne-Twister). Saves and restores the
#' global RNG state so calling par_seeds has no side effect on the
#' caller's random stream.
#'
#' @param seed  Integer master seed.
#' @param n     Number of unique seeds to generate.
#' @return Integer vector of length n.
par_seeds <- function(seed, n) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        old_seed <- get(".Random.seed", envir = .GlobalEnv)
        on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
    } else {
        on.exit(rm(".Random.seed", envir = .GlobalEnv))
    }
    set.seed(seed)
    sample.int(.Machine[["integer.max"]], size = n)
}

# ── Setup ─────────────────────────────────────────────────────

#' Initialise the parallel backend.
#'
#' @param n_cores     Number of worker cores. NULL = detectCores() - 1.
#' @param packages    Character vector of packages to load on workers.
#'                    PSOCK: loaded via library() on each worker.
#'                    Fork: loaded in the parent (children inherit).
#' @param export      Character vector of object names to export to PSOCK
#'                    workers. Looked up in the calling environment.
#'                    Ignored for fork (children share parent memory).
#' @param init_worker A function (no arguments) to run on each worker
#'                    for arbitrary initialisation (e.g., Sys.setenv,
#'                    options, .libPaths). PSOCK: runs via clusterCall.
#'                    Fork: runs once in the parent before forking.
#'                    NULL (default) = no extra initialisation.
#' @return A par_cfg list with [["type"]], [["n_cores"]], [["cl"]].
#'
#' @details Seeds are NOT managed by setup_parallel. Use par_seeds()
#'   to generate one seed per element, then pass them as an iterable to
#'   par_map. Each element calls set.seed() with its own seed, so the
#'   results do not depend on the number of workers. This uses the
#'   default RNG (Mersenne-Twister) throughout, with no L'Ecuyer-CMRG.
setup_parallel <- function(n_cores     = NULL,
                           packages    = character(0),
                           export      = character(0),
                           init_worker = NULL) {

    max_cores <- parallel::detectCores(logical = TRUE)
    n_cores <- if (is.null(n_cores)) max(1L, max_cores - 1L) else n_cores
    n_cores <- n_cores |> as.integer() |> min(max_cores) |> max(1L)

    use_fork <- .is_fork_safe() && n_cores > 1L

    backend <- if (n_cores == 1L) "sequential"
               else if (use_fork) "fork"
               else "psock"

    cfg <- `class<-`(list(
        type    = backend,
        n_cores = n_cores,
        cl      = NULL
    ), "par_cfg")

    if (backend == "sequential") {
        message("[parallel] Running sequentially (1 core).")
        if (length(packages) > 0L) {
            invisible(Map(
                function(p) library(p, character.only = TRUE),
                packages
            ))
        }
        if (!is.null(init_worker)) init_worker()
        return(cfg)
    }

    if (backend == "fork") {
        message(sprintf("[parallel] Using fork backend with %d cores.",
                        n_cores))
        if (length(packages) > 0L) {
            invisible(Map(
                function(p) library(p, character.only = TRUE),
                packages
            ))
        }
        if (!is.null(init_worker)) init_worker()
        return(cfg)
    }

    # PSOCK path
    message(sprintf("[parallel] Using PSOCK backend with %d cores.",
                    n_cores))
    # outfile = "" keeps worker stdout and stderr on the parent's console.
    # The default sends them to the null device, which silences the
    # .progress counter and every message() a worker emits.
    cl <- parallel::makeCluster(n_cores, outfile = "")

    # Load packages on workers
    if (length(packages) > 0L) {
        parallel::clusterCall(cl, function(pkgs) {
            invisible(Map(
                function(p) library(p, character.only = TRUE),
                pkgs
            ))
        }, pkgs = packages)
    }

    # Arbitrary worker initialisation (env vars, options, etc.)
    if (!is.null(init_worker)) {
        parallel::clusterCall(cl, init_worker)
    }

    # Export objects from the caller's environment
    if (length(export) > 0L) {
        envir <- parent.frame()
        parallel::clusterExport(cl, varlist = export, envir = envir)
    }

    `class<-`(`[[<-`(cfg, "cl", cl), "par_cfg")
}

# ── Cleanup ───────────────────────────────────────────────────

#' Stop the PSOCK cluster if one exists. Safe to call multiple times.
cleanup_parallel <- function(cfg) {
    if (!is.null(cfg[["cl"]])) {
        tryCatch(
            parallel::stopCluster(cfg[["cl"]]),
            error = function(e) NULL
        )
    }
    invisible(NULL)
}

# ── Internal: wrap FUN for safe execution + progress ──────────

# par_map passes every wrapper the element index as .par_idx. The progress
# wrapper names its marker file after it, the plain wrapper drops it.
#
# FUN is forced here, in the parent. PSOCK serialises the wrapper with every
# element, and an unforced FUN would be sent as a promise together with the
# frame of par_map, which holds every element.
.make_safe_fun <- function(FUN) {
    force(FUN)
    function(..., .par_idx) {
        tryCatch(
            list(error = FALSE, value = FUN(...)),
            error = function(e) list(error = TRUE, message = conditionMessage(e))
        )
    }
}

.make_progress_fun <- function(FUN, prog_dir, n) {
    safe_fun <- .make_safe_fun(FUN)
    step <- max(1L, n %/% 20L)
    function(..., .par_idx) {
        res <- safe_fun(...)
        # Each completion creates one marker file, named by the element
        # index, in a directory the parent made. Creating a uniquely named
        # file is atomic across processes, whereas concurrent appends to one
        # shared file overwrite each other on Windows and lose a count, so
        # the file count is exact and the last finisher sees all n markers.
        tryCatch({
            file.create(file.path(prog_dir, .par_idx))
            done <- length(list.files(prog_dir))
            if (done %% step == 0L || done == n) {
                message(sprintf("[parallel] %d / %d (%.0f%%)",
                                done, n, 100 * done / n))
            }
        }, error = function(e) NULL)
        res
    }
}

.prepare_fun <- function(FUN, n, .progress) {
    if (.progress) {
        prog_dir <- tempfile("par_progress_")
        dir.create(prog_dir)
        list(
            run_fun  = .make_progress_fun(FUN, prog_dir, n),
            prog_dir = prog_dir
        )
    } else {
        list(
            run_fun  = .make_safe_fun(FUN),
            prog_dir = NULL
        )
    }
}

# ── Primary dispatch: par_map (parallel Map) ──────────────────

#' Parallel Map, the primary parallel dispatch.
#'
#' Applies FUN element-wise over one or more iterables, using
#' mcmapply (fork), clusterMap (PSOCK), or Map (sequential).
#'
#' For reproducibility, generate one seed per element with par_seeds()
#' and pass them as the first iterable:
#'
#'   seeds <- par_seeds(6489L, length(items))
#'   par_map(cfg, function(s, x) { set.seed(s); ... }, seeds, items)
#'
#' @param cfg       A par_cfg object from setup_parallel().
#' @param FUN       Function to apply. Receives one element from each
#'                  iterable as positional arguments.
#' @param ...       One or more vectors/lists to iterate over in
#'                  parallel (passed as the dots of Map/mcmapply/clusterMap).
#' @param .progress Logical. If TRUE, print a progress counter.
#' @return A list the same length as the first iterable. Each element
#'   is a record: list(error = FALSE, value = <FUN result>) on success,
#'   list(error = TRUE, message = <condition message>) on failure.
#'   Pass the list to par_values() to unwrap it.
par_map <- function(cfg, FUN, ..., .progress = FALSE) {

    dots <- list(...)
    n <- length(dots[[1L]])
    prep <- .prepare_fun(FUN, n, .progress)
    if (!is.null(prep[["prog_dir"]])) {
        on.exit(unlink(prep[["prog_dir"]], recursive = TRUE), add = TRUE)
    }

    results <- switch(cfg[["type"]],
        fork = {
            parallel::mcmapply(
                prep[["run_fun"]], ..., .par_idx = seq_len(n),
                mc.cores = cfg[["n_cores"]],
                mc.set.seed = FALSE,
                SIMPLIFY = FALSE
            )
        },
        psock = {
            # run_fun is a closure over FUN, prog_dir, and n, and
            # clusterMap serialises it with its environment, so nothing
            # needs exporting here.
            parallel::clusterMap(
                cfg[["cl"]],
                prep[["run_fun"]], ..., .par_idx = seq_len(n),
                SIMPLIFY = FALSE
            )
        },
        sequential = Map(prep[["run_fun"]], ..., .par_idx = seq_len(n))
    )

    results
}

#' Unwrap par_map records into their values.
#'
#' Stops with the first failure's message, naming the element, when any
#' record has error = TRUE. Names are kept.
#'
#' @param results The list returned by par_map().
#' @return A list of FUN results in the original order.
par_values <- function(results) {
    failed <- results |>
        Map(f = function(r) isTRUE(r[["error"]])) |>
        unlist(use.names = FALSE)
    if (any(failed)) {
        first <- which(failed)[[1L]]
        stop(sprintf("par_map: element %d of %d failed: %s",
                     first, length(results), results[[first]][["message"]]),
             call. = FALSE)
    }
    Map(function(r) r[["value"]], results)
}

#' Parallel mapply, like par_map but simplifying the result.
#' A failed element becomes NA rather than stopping.
par_mapply <- function(cfg, FUN, ..., .progress = FALSE) {
    res <- par_map(cfg, FUN, ..., .progress = .progress)
    vals <- Map(function(r) {
        if (isTRUE(r[["error"]])) NA else r[["value"]]
    }, res)
    tryCatch(simplify2array(vals), error = function(e) vals)
}

# ── Convenience wrappers ──────────────────────────────────────

#' Chunked parallel Map. Splits a data frame into n_cores
#' chunks, sends one chunk per worker, and Reduce(rbind)s the
#' results. It suits PSOCK, where it minimises serialisation.
#'
#' FUN receives a data.frame chunk, followed by any extra arguments
#' given in ..., and must return a data.frame.
par_map_chunked <- function(cfg, df, FUN, ..., .progress = FALSE) {
    n <- nrow(df)
    idx <- ceiling(seq_len(n) * cfg[["n_cores"]] / n)
    chunks <- split(df, idx)

    results <- par_map(cfg, .chunk_fun(FUN, list(...)), chunks,
                       .progress = .progress)

    # Unwrap safe_fun structure and Reduce(rbind)
    results |>
        Map(f = function(r) {
            if (isTRUE(r[["error"]])) {
                warning("Chunk failed: ", r[["message"]])
                NULL
            } else {
                r[["value"]]
            }
        }) |>
        Filter(f = Negate(is.null)) |>
        Reduce(f = rbind)
}

# Built outside par_map_chunked, with its arguments forced, so that its
# environment holds only FUN and the extra arguments. PSOCK serialises the
# function with every chunk, and a closure made inside par_map_chunked, or
# an unforced argument, would carry df and all of its chunks.
.chunk_fun <- function(FUN, extra) {
    force(FUN)
    force(extra)
    function(chunk) do.call(FUN, c(list(chunk), extra))
}

# ── Utilities ─────────────────────────────────────────────────

#' Print method for par_cfg
print.par_cfg <- function(x, ...) {
    cat(sprintf("Parallel config: %s backend, %d core(s)\n",
                x[["type"]], x[["n_cores"]]))
    invisible(x)
}

#' Quick benchmark: compare sequential vs parallel on a test function.
#' Returns a data.frame of timings.
par_benchmark <- function(cfg, X, FUN, ..., times = 3L) {
    seq_cfg <- `class<-`(list(
        type = "sequential", n_cores = 1L, cl = NULL
    ), "par_cfg")

    timings <- Map(function(run) {
        t_seq <- system.time(par_map(seq_cfg, FUN, X, ...))[["elapsed"]]
        t_par <- system.time(par_map(cfg, FUN, X, ...))[["elapsed"]]
        data.frame(
            backend = c("sequential", cfg[["type"]]),
            run     = c(run, run),
            elapsed = c(t_seq, t_par)
        )
    }, seq_len(times)) |> Reduce(f = rbind)

    agg <- aggregate(elapsed ~ backend, data = timings, FUN = function(x) {
        c(mean = mean(x), sd = sd(x))
    })

    message("[benchmark] Mean elapsed times:")
    invisible(Map(function(i) {
        message(sprintf("  %-12s %.3fs (sd %.3fs)",
                        agg[["backend"]][i], agg[["elapsed"]][i, "mean"],
                        agg[["elapsed"]][i, "sd"]))
    }, seq_len(nrow(agg))))

    timings
}
