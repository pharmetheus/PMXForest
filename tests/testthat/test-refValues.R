## Reference-value resolution shared by setupDfCovs() and setupDfRefRow().

modFile <- system.file("extdata", "SimVal/run7.mod", package = "PMXForest")

simData <- function() {
  read.csv(system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv",
                       package = "PMXForest"))
}

## 24 distinct weights, so WT is comfortably above the minLevels = 10 threshold
## and is treated as continuous.
mockData <- data.frame(
  ID   = 1:24,
  WT   = seq(60, 106, by = 2),
  SEX  = rep(c(1, 1, 1, 2), length.out = 24),          # mode 1
  RACE = rep(c(1, 2, 2, 3, 2, 2), length.out = 24)     # lowest 1, mode 2
)

test_that("a bare scalar applies to every covariate", {
  r <- refResolve(mockData, c("WT", "SEX"), contRef = "mean", idVar = "ID")
  expect_equal(r$WT$value, signif(mean(mockData$WT), 3))
  expect_equal(r$SEX$value, 1)                       # categorical falls back to mode
})

test_that("a list gives per-covariate settings with a default", {
  d <- cbind(mockData, AGE = seq(20, 65, length.out = 24))
  r <- refResolve(d, c("WT", "AGE"), idVar = "ID",
                  contRef = list(WT = 75, default = "mean"))
  expect_equal(r$WT$value, 75)
  expect_match(r$WT$source, "supplied directly")
  expect_equal(r$AGE$value, signif(mean(d$AGE), 3))
})

test_that("an explicitly supplied value is not rounded", {
  r <- refResolve(mockData, "WT", contRef = list(WT = 123.4567), idVar = "ID",
                  nsig = 3)
  expect_equal(r$WT$value, 123.4567)
})

test_that("categorical settings cover mode, lowest and an explicit level", {
  expect_equal(refResolve(mockData, "RACE", catRef = "mode", idVar = "ID")$RACE$value, 2)
  expect_equal(refResolve(mockData, "RACE", catRef = "lowest", idVar = "ID")$RACE$value, 1)
  expect_equal(refResolve(mockData, "RACE", catRef = list(RACE = 3),
                          idVar = "ID")$RACE$value, 3)
})

test_that("\"model\" reads the reference out of the control stream", {
  d <- simData()
  r <- suppressWarnings(refResolve(
    d, c("WT", "FORM", "SEX"), contRef = "model", catRef = "model",
    model = modFile, idVar = "ID"
  ))
  expect_equal(r$WT$value, 75)      # normalisation constant, not the data median
  expect_equal(r$FORM$value, 1)     # "; Most common" branch, not the data mode
  expect_equal(r$SEX$value, 1)
  expect_match(r$WT$source, "normalisation constant")
})

test_that("\"model\" also accepts a createParamFunction() result", {
  d   <- simData()
  gen <- suppressWarnings(createParamFunction(modFile, parameters = "CL",
                                              quiet = TRUE))
  fromPath <- suppressWarnings(refResolve(d, "WT", contRef = "model",
                                          model = modFile, idVar = "ID"))
  fromObj  <- refResolve(d, "WT", contRef = "model", model = gen, idVar = "ID")
  expect_equal(fromObj$WT$value, fromPath$WT$value)
})

test_that("\"model\" without a model, or for an underivable covariate, errors", {
  d <- simData()
  expect_error(refResolve(d, "WT", contRef = "model", idVar = "ID"),
               "`model` is required")
  # CRCL is not referenced in run7.mod's $PK at all
  expect_error(refResolve(d, "CRCL", contRef = "model", model = modFile,
                          idVar = "ID"),
               "No reference value for covariate CRCL")
  expect_error(refResolve(d, "WT", contRef = "model", model = 42, idVar = "ID"),
               "control stream path")
})

test_that("a setting inappropriate to the covariate type is rejected", {
  expect_error(refResolve(mockData, "SEX", catRef = "mean", idVar = "ID"),
               "not a reference for the categorical")
  expect_error(refResolve(mockData, "WT", contRef = "mode", idVar = "ID"),
               "not a reference for the continuous")
  expect_error(refResolve(mockData, "SEX", catRef = "median", idVar = "ID"),
               "not a reference for the categorical")
  expect_error(refResolve(mockData, "WT", contRef = "nonsense", idVar = "ID"),
               "Unknown reference setting")
  expect_error(refResolve(mockData, "WT", contRef = c(1, 2), idVar = "ID"),
               "single value")
})

test_that("references are computed on deduplicated data", {
  dup <- rbind(mockData, mockData[rep(1, 40), ])   # subject 1 forty times over
  expect_equal(refResolve(dup, "WT", contRef = "median", idVar = "ID")$WT$value,
               refResolve(mockData, "WT", contRef = "median", idVar = "ID")$WT$value)
})

test_that("refLevels is deprecated but still honoured", {
  expect_warning(
    dfr <- setupDfCovs(mockData, covariates = "WT", additionalCovs = "RACE",
                       refLevels = list(RACE = 2)),
    "`refLevels` is deprecated"
  )
  expect_true("RACE_1" %in% names(dfr))
  expect_false("RACE_2" %in% names(dfr))

  expect_warning(getCovStats(mockData, "RACE", refLevels = list(RACE = 2)),
                 "deprecated")
  expect_warning(setupCovExpressionsList(mockData, "RACE", minSubjects = 1,
                                         includeReference = FALSE,
                                         refLevels = list(RACE = 2)),
                 "deprecated")
})

test_that("supplying both catRef and refLevels is an error", {
  expect_error(
    suppressWarnings(getCovStats(mockData, "RACE", refLevels = list(RACE = 2),
                                 catRef = list(RACE = 1))),
    "not both"
  )
})

test_that("setupDfRefRow with model references removes the reference mismatch", {
  # The bug this was written for: setupDfRefRow put WT at the data median (85.4)
  # while $PK normalises at 75, so every row where WT was inactive sat at
  # (75/85.4)^theta instead of 1. Same for FORM, at the data mode rather than
  # the "; Most common" level.
  d      <- simData()
  covs   <- c("WT", "SEX", "FOOD", "FORM")
  dfCovs <- setupDfCovs(d, covariates = covs, idVar = "ID")

  dataRef  <- setupDfRefRow(dfCovs, d, covs, idVar = "ID")
  modelRef <- setupDfRefRow(dfCovs, d, covs, idVar = "ID", contRef = "model",
                            catRef = "model", model = modFile)

  expect_equal(dataRef$WT, 85.4)      # unchanged default behaviour
  expect_equal(dataRef$FORM, 0)
  expect_equal(modelRef$WT, 75)       # the model's normalisation weight
  expect_equal(modelRef$FORM, 1)      # the model's reference formulation
})

test_that("a model-derived reference row makes the inactive rows sit at 1", {
  d      <- simData()
  covs   <- c("WT", "SEX", "FOOD", "FORM")
  dfCovs <- setupDfCovs(d, covariates = covs, idVar = "ID")
  dfRef  <- setupDfRefRow(dfCovs, d, covs, idVar = "ID", contRef = "model",
                          catRef = "model", model = modFile)

  gen <- createParamFunction(modFile, parameters = c("CL", "V"),
                             covRef = list(GENO1 = 0, GENO3 = 0), quiet = TRUE)
  fun <- eval(parse(text = gen$code))

  ext    <- getExt(system.file("extdata", "SimVal/run7.ext", package = "PMXForest"))
  thetas <- as.numeric(ext[ext$ITERATION == -1000000000, 2:15])

  ref <- unlist(fun(thetas, dfRef[, setdiff(names(dfRef), "COVARIATEGROUPS"),
                                  drop = FALSE]))
  # A row varying only SEX leaves WT and FORM inactive; with the reference taken
  # from the model, CL and V there must equal the reference exactly.
  sexRow <- dfCovs[dfCovs$COVARIATEGROUPS == "SEX", , drop = FALSE][1, , drop = FALSE]
  got    <- unlist(fun(thetas, sexRow[, setdiff(names(sexRow), "COVARIATEGROUPS"),
                                      drop = FALSE]))
  expect_equal(got[["CL"]] / ref[["CL"]], 1, tolerance = 1e-12)
  expect_equal(got[["V"]]  / ref[["V"]],  1, tolerance = 1e-12)
})

test_that("catRef = \"model\" and \"mode\" resolve the encoding level too", {
  d <- simData()
  # "model": GENO's reference level comes from the control stream. run7.mod uses
  # the GENO1..GENO4 dummies rather than raw GENO, so ask about a covariate the
  # $PK block does name.
  st <- suppressWarnings(getCovStats(d, "RACEL", idVar = "ID", catRef = "mode"))
  # mode of RACEL is the dropped level, so its dummy is absent
  modeLev <- refMode(d$RACEL[!duplicated(d$ID)])
  expect_false(paste0("RACEL_", modeLev) %in% names(st$RACEL))

  st2 <- getCovStats(d, "RACEL", idVar = "ID", catRef = "lowest")
  expect_false("RACEL_1" %in% names(st2$RACEL))
})

test_that("catRef = \"model\" for the encoding level errors when underivable", {
  d <- simData()
  expect_error(
    getCovStats(d, "RACEL", idVar = "ID", catRef = "model", model = modFile),
    "No reference level for covariate RACEL"
  )
  expect_error(getCovStats(d, "RACEL", idVar = "ID", catRef = "nonsense"),
               "Unknown reference setting")
})

test_that("an inferred model reference warns from setupDfRefRow too", {
  d      <- simData()
  covs   <- c("WT", "GENO1")     # GENO1's reference is inferred, not stated
  dfCovs <- setupDfCovs(d, covariates = covs, idVar = "ID")
  expect_warning(
    setupDfRefRow(dfCovs, d, covs, idVar = "ID", contRef = "model",
                  catRef = "model", model = modFile),
    "inferred from the model"
  )
})

test_that("a covariate absent from the data is reported by name", {
  expect_error(refResolve(mockData, "NOPE", idVar = "ID"),
               "Covariate NOPE is not present in the data")
})

test_that("a single-level covariate is handled", {
  d <- cbind(mockData, ONE = 1)
  expect_equal(refResolve(d, "ONE", catRef = "mode", idVar = "ID")$ONE$value, 1)
})

test_that("model reference errors surface clearly", {
  expect_error(refModelValues(list(a = 1), -99), "createParamFunction")
  f <- withr::local_tempfile(fileext = ".mod")
  writeLines(c("$PROBLEM x", "$INPUT ID DV"), f)
  expect_error(refModelValues(f, -99), "No \\$PK record")
})

test_that("\"lowest\" is rejected for a continuous covariate", {
  expect_error(refResolve(mockData, "WT", contRef = "lowest", idVar = "ID"),
               "not a reference for the continuous")
})

test_that("a malformed per-covariate setting is rejected", {
  # A length-1 logical cannot reach refSpec()'s length check - it is caught
  # further on by refResolve()'s type test, whose message happens to contain
  # the same words. Use a length-2 value so the length check is what fires.
  expect_error(
    refResolve(mockData, "WT", contRef = list(WT = c(70, 80)), idVar = "ID"),
    "given for WT must be a single value"
  )
  expect_error(
    refResolve(mockData, "WT", contRef = list(default = c(70, 80)), idVar = "ID"),
    "given as `default` must be a single value"
  )
  # the bare-vector spelling of the same mistake is still rejected
  expect_error(
    refResolve(mockData, "WT", contRef = c(70, 80), idVar = "ID"),
    "given as a vector must be a single value"
  )
  # and the original type-error path still works
  expect_error(
    refResolve(mockData, "WT", contRef = list(WT = TRUE), idVar = "ID"),
    "must be a single value"
  )
})
