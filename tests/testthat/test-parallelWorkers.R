## Fresh worker sessions (what doParallel starts on Windows) load packages by
## searching .libPaths(). A package the session loaded from a folder that is
## not on .libPaths() - PMXRenv's versioned library, library(lib.loc = ) -
## then resolves to a different copy on the workers. These tests reproduce that
## with a one-function package installed twice: an old version in a
## "qualified" library the workers see, and a newer one in a "versioned"
## folder only the session loads from.

## Install a toy package `name` at `version` into library folder `lib`, with
## the extra DESCRIPTION lines `fields`, the NAMESPACE lines `namespace` and
## the R code `code`.
installToySource <- function(lib, name, version, fields, namespace, code) {
  src <- file.path(withr::local_tempdir(), name)
  dir.create(file.path(src, "R"), recursive = TRUE)
  writeLines(c(
    paste0("Package: ", name), paste0("Version: ", version),
    "Title: Toy Package", "Description: A toy package for tests.",
    "Author: PMXForest tests", "Maintainer: PMXForest tests <noreply@example.com>",
    "License: GPL-3", fields
  ), file.path(src, "DESCRIPTION"))
  writeLines(namespace, file.path(src, "NAMESPACE"))
  writeLines(code, file.path(src, "R", "toy.R"))
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  out <- system2(file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--no-docs", "--no-test-load", paste0("--library=", shQuote(lib)), shQuote(src)),
    stdout = TRUE, stderr = TRUE
  )
  if (!file.exists(file.path(lib, name, "DESCRIPTION"))) {
    stop("could not install ", name, " ", version, ":\n", paste(out, collapse = "\n"))
  }
  invisible(lib)
}

## Install forestToyPkg `version` into library folder `lib`.
installToyPkg <- function(lib, version) {
  installToySource(
    lib, "forestToyPkg", version, NULL,
    "export(toyVersion)", sprintf('toyVersion <- function() "%s"', version)
  )
}

## Install forestToyTop into library folder `lib`. forestToyTop imports
## forestToyPkg, at version `needs` or later unless `needs` is NULL, as PMXFrem
## imports PMXForest, and its depVersion() returns the version of forestToyPkg
## it runs with. R CMD INSTALL loads forestToyPkg from library folder `depLib`.
installToyTop <- function(lib, needs, depLib) {
  imports <- if (is.null(needs)) "forestToyPkg" else paste0("forestToyPkg (>= ", needs, ")")
  withr::with_envvar(
    c(R_LIBS = depLib),
    installToySource(
      lib, "forestToyTop", "1.0.0", paste0("Imports: ", imports),
      c("import(forestToyPkg)", "export(depVersion)"), "depVersion <- function() toyVersion()"
    )
  )
}

## The layout of an old PMXForest next to PMXFrem, with a newer PMXForest in
## another library: the library `top` holds forestToyTop and forestToyPkg
## 1.0.0, and the session loads forestToyPkg 2.0.0 from the library `newer`.
## forestToyTop requires forestToyPkg `needs`, or any version if NULL.
dependencySetup <- function(needs, env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = env)
  newer <- installToyPkg(file.path(root, "newer"), "2.0.0")
  top <- installToyTop(file.path(root, "top"), needs = needs, depLib = newer)
  installToyPkg(top, "1.0.0")
  unloadToys <- function() {
    for (p in c("forestToyTop", "forestToyPkg")) if (isNamespaceLoaded(p)) unloadNamespace(p)
  }
  unloadToys()
  withr::defer(unloadToys(), envir = env)
  loadNamespace("forestToyPkg", lib.loc = newer)
  loadNamespace("forestToyTop", lib.loc = top)
  list(top = top, newer = newer)
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

  ## With startParallelWorkers() they get the session's copy.
  stopWorkers <- startParallelWorkers(2, pkgs = "forestToyPkg", fresh = TRUE)
  withr::defer(stopWorkers())
  got <- foreach::foreach(i = 1:2, .packages = "forestToyPkg") %dopar% forestToyPkg::toyVersion()
  expect_equal(unlist(got), c("2.0.0", "2.0.0"))
})

test_that("an older copy of a dependency next to a package in pkgs is not loaded on the workers", {
  skip_on_cran()
  paths <- dependencySetup(needs = "2.0.0")
  expect_equal(forestToyTop::depVersion(), "2.0.0")

  ## The premise: a fresh session that searches the library of forestToyTop
  ## first finds forestToyPkg 1.0.0 there and cannot load forestToyTop.
  tryLoad <- function(dirs) {
    .libPaths(c(dirs, .libPaths()))
    tryCatch(
      {
        loadNamespace("forestToyTop")
        "loaded"
      },
      error = function(e) conditionMessage(e)
    )
  }
  environment(tryLoad) <- baseenv()
  cl <- parallel::makePSOCKcluster(1)
  plain <- parallel::clusterCall(cl, tryLoad, c(paths$top, paths$newer))[[1]]
  parallel::stopCluster(cl)
  expect_match(plain, "'forestToyPkg' 1.0.0 is being loaded, but >= 2.0.0 is required", fixed = TRUE)

  ## With startParallelWorkers() the workers load forestToyPkg 2.0.0, which
  ## this session uses, and then forestToyTop.
  stopWorkers <- startParallelWorkers(2, pkgs = "forestToyTop", fresh = TRUE)
  withr::defer(stopWorkers())
  got <- parallel::clusterEvalQ(attr(stopWorkers, "cluster"), forestToyTop::depVersion())
  expect_equal(unlist(got), c("2.0.0", "2.0.0"))
})

test_that("a dependency is loaded before the package that imports it, so the session's copy is used", {
  skip_on_cran()
  ## Without a version requirement, loading forestToyTop first would load the
  ## forestToyPkg 1.0.0 beside it without an error, and the workers would run
  ## an older forestToyPkg than the session.
  paths <- dependencySetup(needs = NULL)
  expect_equal(forestToyTop::depVersion(), "2.0.0")

  stopWorkers <- startParallelWorkers(2, pkgs = "forestToyTop", fresh = TRUE)
  withr::defer(stopWorkers())
  got <- parallel::clusterEvalQ(attr(stopWorkers, "cluster"), forestToyTop::depVersion())
  expect_equal(unlist(got), c("2.0.0", "2.0.0"))
})

test_that("startParallelWorkers() is exported for packages that build on PMXForest", {
  skip_on_cran()
  ## PMXFrem starts its workers with this function rather than a copy of it,
  ## so it must be in NAMESPACE - load_all() alone would reach it either way.
  expect_true("export(startParallelWorkers)" %in% readLines(system.file("NAMESPACE", package = "PMXForest")))

  ## Called the way another package calls it, for that package's own copy.
  paths <- versionedSetup()
  stopWorkers <- PMXForest::startParallelWorkers(1, pkgs = "forestToyPkg", fresh = TRUE)
  withr::defer(stopWorkers())
  got <- parallel::clusterEvalQ(attr(stopWorkers, "cluster"), forestToyPkg::toyVersion())[[1]]
  expect_equal(got, "2.0.0")
})

test_that("setting up the workers does not itself load PMXForest there", {
  skip_on_cran()
  ## Were the setup function to carry PMXForest's namespace with it, each worker
  ## would load PMXForest from its own library path before the setup could put
  ## the session's folders first - the very mismatch this exists to prevent.
  paths <- versionedSetup()
  stopWorkers <- startParallelWorkers(1, pkgs = "forestToyPkg", fresh = TRUE)
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
    startParallelWorkers(2, pkgs = "forestToyPkg", fresh = TRUE),
    "forestToyPkg"
  )
  expect_match(conditionMessage(err), "2.0.0", fixed = TRUE)
  expect_match(conditionMessage(err), "1.5.0", fixed = TRUE)
  expect_match(conditionMessage(err), "ncores = 1", fixed = TRUE)
})

test_that("a package the workers cannot load at all is reported as such", {
  skip_on_cran()
  expect_error(
    startParallelWorkers(1, pkgs = "forestNoSuchPkg", fresh = TRUE),
    "forestNoSuchPkg"
  )
})

test_that("forked workers are registered as before", {
  skip_on_os("windows")
  stopWorkers <- startParallelWorkers(2, fresh = FALSE)
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
