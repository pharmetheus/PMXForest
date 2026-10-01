## Fresh worker sessions (what doParallel starts on Windows) load packages by
## searching .libPaths(). A package the session loaded from a folder that is
## not on .libPaths() - PMXRenv's versioned library, library(lib.loc = ) -
## then resolves to a different copy on the workers. These tests reproduce that
## with a one-function package installed twice: an old version in a
## "qualified" library the workers see, and a newer one in a "versioned"
## folder only the session loads from.

## Install forestToyPkg `version` into library folder `lib`.
installToyPkg <- function(lib, version) {
  src <- file.path(withr::local_tempdir(), "forestToyPkg")
  dir.create(file.path(src, "R"), recursive = TRUE)
  writeLines(c(
    "Package: forestToyPkg", paste0("Version: ", version),
    "Title: Toy Package", "Description: A toy package for tests.",
    "Author: PMXForest tests", "Maintainer: PMXForest tests <noreply@example.com>",
    "License: GPL-3"
  ), file.path(src, "DESCRIPTION"))
  writeLines("export(toyVersion)", file.path(src, "NAMESPACE"))
  writeLines(sprintf('toyVersion <- function() "%s"', version), file.path(src, "R", "toy.R"))
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  out <- system2(file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--no-docs", "--no-test-load", paste0("--library=", shQuote(lib)), shQuote(src)),
    stdout = TRUE, stderr = TRUE
  )
  if (!file.exists(file.path(lib, "forestToyPkg", "DESCRIPTION"))) {
    stop("could not install forestToyPkg ", version, ":\n", paste(out, collapse = "\n"))
  }
  invisible(lib)
}

## The qualified/versioned setup: workers see `qualified` (1.0.0) through
## R_LIBS; the session loads 2.0.0 from `versioned`, which is not on
## .libPaths().
versionedSetup <- function(env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = env)
  qualified <- installToyPkg(file.path(root, "qualified"), "1.0.0")
  versioned <- installToyPkg(file.path(root, "versioned", "foresttoypkg", "2.0.0"), "2.0.0")
  withr::local_envvar(R_LIBS = qualified, .local_envir = env)
  if (isNamespaceLoaded("forestToyPkg")) unloadNamespace("forestToyPkg")
  loadNamespace("forestToyPkg", lib.loc = versioned)
  withr::defer(if (isNamespaceLoaded("forestToyPkg")) unloadNamespace("forestToyPkg"), envir = env)
  list(qualified = qualified, versioned = versioned)
}

test_that("fresh workers load a package from the folder the session loaded it from", {
  skip_on_cran()
  paths <- versionedSetup()
  expect_equal(forestToyPkg::toyVersion(), "2.0.0")
  expect_false(any(startsWith(paths$versioned, .libPaths())))

  ## The premise: an ordinary fresh worker finds the qualified copy.
  cl <- parallel::makePSOCKcluster(1)
  plain <- parallel::clusterEvalQ(cl, forestToyPkg::toyVersion())[[1]]
  parallel::stopCluster(cl)
  expect_equal(plain, "1.0.0")

  ## With .forestStartWorkers() they get the session's copy.
  stopWorkers <- .forestStartWorkers(2, pkgs = "forestToyPkg", fresh = TRUE)
  withr::defer(stopWorkers())
  got <- foreach::foreach(i = 1:2, .packages = "forestToyPkg") %dopar% forestToyPkg::toyVersion()
  expect_equal(unlist(got), c("2.0.0", "2.0.0"))
})

test_that("setting up the workers does not itself load PMXForest there", {
  skip_on_cran()
  ## Were the setup function to carry PMXForest's namespace with it, each worker
  ## would load PMXForest from its own library path before the setup could put
  ## the session's folders first - the very mismatch this exists to prevent.
  paths <- versionedSetup()
  stopWorkers <- .forestStartWorkers(1, pkgs = "forestToyPkg", fresh = TRUE)
  withr::defer(stopWorkers())
  loaded <- parallel::clusterEvalQ(attr(stopWorkers, "cluster"), loadedNamespaces())[[1]]
  expect_true("forestToyPkg" %in% loaded)
  expect_false("PMXForest" %in% loaded)
})

test_that("a worker that cannot load the session's version stops with both versions named", {
  skip_on_cran()
  paths <- versionedSetup()

  ## The session keeps 2.0.0 in memory, but the folder it came from now holds
  ## 1.5.0 - what happens when a library is updated under a running session.
  ## R reads a function body from disk on first use, so use it before the
  ## files change underneath it.
  expect_equal(forestToyPkg::toyVersion(), "2.0.0")
  unlink(file.path(paths$versioned, "forestToyPkg"), recursive = TRUE)
  installToyPkg(paths$versioned, "1.5.0")
  expect_equal(forestToyPkg::toyVersion(), "2.0.0")

  err <- expect_error(
    .forestStartWorkers(2, pkgs = "forestToyPkg", fresh = TRUE),
    "forestToyPkg"
  )
  expect_match(conditionMessage(err), "2.0.0", fixed = TRUE)
  expect_match(conditionMessage(err), "1.5.0", fixed = TRUE)
  expect_match(conditionMessage(err), "ncores = 1", fixed = TRUE)
})

test_that("a package the workers cannot load at all is reported as such", {
  skip_on_cran()
  expect_error(
    .forestStartWorkers(1, pkgs = "forestNoSuchPkg", fresh = TRUE),
    "forestNoSuchPkg"
  )
})

test_that("forked workers are registered as before", {
  skip_on_os("windows")
  stopWorkers <- .forestStartWorkers(2, fresh = FALSE)
  withr::defer(stopWorkers())
  expect_equal(foreach::getDoParName(), "doParallelMC")
  expect_equal(foreach::getDoParWorkers(), 2)
})

test_that("getForestDFSCM gives the same result on fresh workers as on one core", {
  skip_on_cran()
  ## Fresh workers load the installed PMXForest, so the comparison is only
  ## meaningful when this session runs an installed copy - not the source tree
  ## through load_all(), as the unit-test job does.
  skip_if_not(
    file.exists(file.path(getNamespaceInfo("PMXForest", "path"), "Meta", "package.rds")),
    "PMXForest is not loaded from an installed library"
  )
  dfData <- read.csv(system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest"))
  extFile <- system.file("extdata", "SimVal/run7.ext", package = "PMXForest")
  covFile <- system.file("extdata", "SimVal/run7.cov", package = "PMXForest")
  dfCovs <- setupDfCovs(dfData, covariates = c("WT", "FOOD"), idVar = "ID")
  set.seed(1)
  dfSamples <- getSamples(covFile, extFile, n = 20)
  paramFunction <- function(thetas, df, ...) {
    TVCL <- thetas[4]
    if (any(names(df) == "WT") && df$WT != -99) TVCL <- thetas[4] * (df$WT / 75)^thetas[2]
    list(CL = TVCL)
  }
  forest <- function(ncores) {
    getForestDFSCM(
      dfCovs = dfCovs, functionList = list(paramFunction), functionListName = "CL",
      noBaseThetas = 14, dfParameters = dfSamples, ncores = ncores
    )
  }
  withr::local_options(PMXForest.freshWorkers = TRUE)
  expect_equal(forest(2), forest(1))
})
