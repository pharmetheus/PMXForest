#' Start the parallel workers for a `%dopar%` loop
#'
#' On macOS and Linux, doParallel forks this session: each worker is a copy of
#' it and shares its loaded packages. On Windows it cannot fork, so the workers
#' are fresh R sessions, which load packages by searching `.libPaths()`. A
#' package this session loaded from a folder that is not on `.libPaths()` -
#' PMXRenv's versioned library, or `library(lib.loc = )` - then resolves to a
#' different installed copy on the workers, which runs a different version's
#' code: an internal function it lacks fails with "could not find function",
#' and one it has silently computes with the old code.
#'
#' Fresh workers therefore get the folders this session loaded its packages
#' from first on their library path, load `pkgs` before any task reaches them,
#' and are checked against the session: any package whose version differs
#' stops the run, naming both versions and where they came from.
#'
#' @param ncores Number of workers.
#' @param pkgs Packages the loop needs on the workers.
#' @param fresh Start fresh worker sessions instead of forking. Defaults to
#'   `TRUE` on Windows, where doParallel cannot fork. The
#'   `PMXForest.freshWorkers` option sets it elsewhere, for testing.
#'
#' @return A function that stops the workers, to be called on exit. Its
#'   `cluster` attribute holds the workers, for tests.
#' @noRd
.forestStartWorkers <- function(ncores, pkgs = "PMXForest",
                                fresh = getOption("PMXForest.freshWorkers", .Platform$OS.type == "windows")) {
  if (!fresh) {
    doParallel::registerDoParallel(cores = ncores)
    return(function() doParallel::stopImplicitCluster())
  }

  ## Where and at which version this session loaded each package. The version
  ## is the one in memory, not the one now on disk: the two differ once a
  ## library has been updated under a running session.
  ## base has no namespace path and always comes with R itself.
  ns <- setdiff(loadedNamespaces(), "base")
  nsPath <- vapply(ns, function(p) getNamespaceInfo(p, "path"), character(1))
  nsVersion <- vapply(ns, function(p) as.character(getNamespaceVersion(p)), character(1))

  ## The folders to search first: those of `pkgs`, then of everything else
  ## loaded. R's own library is on every worker's path already.
  first <- c(nsPath[intersect(pkgs, ns)], nsPath)
  dirs <- setdiff(unique(dirname(first)), c(.Library, .Library.site))

  cl <- parallel::makePSOCKcluster(ncores)
  ready <- FALSE
  on.exit(if (!ready) parallel::stopCluster(cl), add = TRUE)

  ## Run on every worker. Its environment must not be this package's
  ## namespace: serialising the function would carry a reference to it, and
  ## the worker would load PMXForest - from its own library path - before the
  ## body could put this session's folders first.
  setup <- function(dirs, pkgs) {
    .libPaths(c(dirs, .libPaths()))
    failed <- character(0)
    for (p in pkgs) {
      msg <- tryCatch(
        {
          loadNamespace(p)
          NULL
        },
        error = function(e) conditionMessage(e)
      )
      if (!is.null(msg)) failed[p] <- msg
    }
    loaded <- setdiff(loadedNamespaces(), "base")
    list(
      failed = failed,
      version = vapply(loaded, function(p) as.character(getNamespaceVersion(p)), character(1)),
      path = vapply(loaded, function(p) getNamespaceInfo(p, "path"), character(1))
    )
  }
  environment(setup) <- baseenv()
  workers <- parallel::clusterCall(cl, setup, dirs, pkgs)

  failed <- unlist(lapply(workers, `[[`, "failed"))
  if (length(failed) > 0) {
    failed <- failed[!duplicated(names(failed))]
    stop("The parallel workers could not load ",
      paste0(names(failed), " (", failed, ")", collapse = "; "),
      ". Use ncores = 1, or install the package where a new R session finds it.",
      call. = FALSE
    )
  }

  clash <- character(0)
  for (w in workers) {
    common <- intersect(names(w$version), ns)
    for (p in common[w$version[common] != nsVersion[common]]) {
      clash[p] <- paste0(
        p, ": this session uses ", nsVersion[[p]], " from ", nsPath[[p]],
        ", the workers loaded ", w$version[[p]], " from ", w$path[[p]]
      )
    }
  }
  if (length(clash) > 0) {
    stop("The parallel workers could not load the package versions this session uses.\n",
      paste0("  ", clash, collapse = "\n"),
      "\nA new R session finds a different copy, for example when a package was ",
      "loaded from a folder that is not on .libPaths(), or was updated after this ",
      "session loaded it. Restart R so the session and the installed packages ",
      "agree, or use ncores = 1.",
      call. = FALSE
    )
  }

  ready <- TRUE
  doParallel::registerDoParallel(cl)
  structure(function() {
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }, cluster = cl)
}
