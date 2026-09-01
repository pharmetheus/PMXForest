## Mock data: 240 subjects, one row each plus a duplicated row for subject 1.
mock_data <- local({
  n  <- 240
  df <- data.frame(
    ID   = seq_len(n),
    WT   = seq(50, 120, length.out = n),
    AGE  = rep(c(20, 30, 40, 50, 60, 70), length.out = n),
    SEX  = rep(c(1, 2), length.out = n),
    GENO = rep(c(1, 2, 3, 4), length.out = n),
    # period 4 so SEX & FOOD combinations are all populated
    FOOD = rep(c(0, 1, 1, 0), length.out = n),
    # cycling ramps so CRCL / CRP are not collinear with WT / ID order
    CRCL = rep(seq(40, 158, length.out = 24), length.out = n),
    CRP  = rep(seq(0.137, 88.42, length.out = 24), length.out = n)
  )
  rbind(df, df[1, ]) # duplicate record for subject 1
})

expr_strings <- function(out) {
  vapply(out$covExpressionsList, as.character, character(1), USE.NAMES = FALSE)
}

test_that("return structure is the two getForestDFemp arguments", {
  out <- setupCovExpressionsList(mock_data, c("WT", "SEX", "GENO"))

  expect_named(out, c("covExpressionsList", "cdfCovsNames"))
  expect_type(out$covExpressionsList, "list")
  expect_true(all(vapply(out$covExpressionsList, is.expression, logical(1))))
  expect_equal(length(out$covExpressionsList), length(out$cdfCovsNames))
  expect_equal(names(out$covExpressionsList),
               c("WT", "WT", "SEX", "SEX", "GENO", "GENO", "GENO", "GENO"))
})

test_that("binary covariate emits both sorted levels", {
  out <- setupCovExpressionsList(mock_data, "SEX")
  expect_equal(expr_strings(out), c("SEX == 1", "SEX == 2"))
  expect_equal(names(out$covExpressionsList), c("SEX", "SEX"))
})

test_that("multi-level categorical includes every level by default", {
  out <- setupCovExpressionsList(mock_data, "GENO")
  expect_equal(expr_strings(out),
               c("GENO == 1", "GENO == 2", "GENO == 3", "GENO == 4"))
})

test_that("includeReference = FALSE drops the reference level", {
  lowest <- setupCovExpressionsList(mock_data, "GENO", includeReference = FALSE)
  expect_equal(expr_strings(lowest), c("GENO == 2", "GENO == 3", "GENO == 4"))

  named <- setupCovExpressionsList(mock_data, "GENO", includeReference = FALSE,
                                   catRef = list(GENO = 2))
  expect_equal(expr_strings(named), c("GENO == 1", "GENO == 3", "GENO == 4"))

  expect_error(
    setupCovExpressionsList(mock_data, "GENO", includeReference = FALSE,
                            catRef = list(GENO = 9)),
    "not present in the data"
  )
})

test_that("continuous quantile split uses the probs quantiles", {
  v  <- mock_data$WT[!duplicated(mock_data$ID)]
  qs <- signif(stats::quantile(v, probs = c(0.1, 0.9), names = FALSE), 3)

  out <- setupCovExpressionsList(mock_data, "WT", probs = c(0.1, 0.9))
  expect_equal(expr_strings(out),
               c(paste0("WT < ", qs[1]), paste0("WT >= ", qs[2])))
  expect_equal(out$cdfCovsNames, c(paste0("WT <", qs[1]), paste0("WT >=", qs[2])))

  wide <- setupCovExpressionsList(mock_data, "WT", probs = c(0.25, 0.75))
  expect_false(identical(expr_strings(out), expr_strings(wide)))
})

test_that("continuous median split partitions the subjects", {
  v <- mock_data$WT[!duplicated(mock_data$ID)]
  m <- signif(stats::median(v), 3)

  out <- setupCovExpressionsList(mock_data, "WT", contSplit = "median")
  expect_equal(expr_strings(out), c(paste0("WT < ", m), paste0("WT >= ", m)))

  lo <- sum(v <  m)
  hi <- sum(v >= m)
  expect_equal(lo + hi, length(v))
})

test_that("nsig rounds the threshold in expression and label", {
  v  <- mock_data$CRP[!duplicated(mock_data$ID)]
  m2 <- signif(stats::median(v), 2)
  m3 <- signif(stats::median(v), 3)
  expect_false(identical(m2, m3)) # rounding is actually visible here

  out <- setupCovExpressionsList(mock_data, "CRP", contSplit = "median", nsig = 2)
  expect_equal(expr_strings(out), c(paste0("CRP < ", m2), paste0("CRP >= ", m2)))
  expect_equal(out$cdfCovsNames, c(paste0("CRP <", m2), paste0("CRP >=", m2)))
})

test_that("categorical additionalCov: own rows plus a condition on the others", {
  out <- setupCovExpressionsList(mock_data, c("WT", "SEX"), contSplit = "median",
                                 additionalCovs = list(FOOD = 1))

  expect_equal(names(out$covExpressionsList),
               c("WT", "WT", "SEX", "SEX", "FOOD", "FOOD"))
  es <- expr_strings(out)
  # every primary expression carries the FOOD condition
  expect_true(all(grepl("& FOOD == 1$", es[1:4])))
  # FOOD's own rows are not self-conditioned
  expect_equal(es[5:6], c("FOOD == 0", "FOOD == 1"))
})

test_that("continuous additionalCov with prob places the condition at the quantile", {
  v <- mock_data$CRCL[!duplicated(mock_data$ID)]
  p <- signif(stats::quantile(v, probs = 0.5, names = FALSE), 3)

  out <- setupCovExpressionsList(mock_data, "WT", contSplit = "median",
                                 additionalCovs = list(CRCL = list(prob = 0.5, dir = "gt")))
  es <- expr_strings(out)
  expect_true(all(grepl(paste0("& CRCL > ", p, "$"), es[grepl("^WT", es)])))
  # CRCL own rows present, not self-conditioned
  expect_true(any(grepl("^CRCL < ", es)))
})

test_that("continuous additionalCov with value and dir = lt", {
  out <- setupCovExpressionsList(mock_data, "WT", contSplit = "median",
                                 additionalCovs = list(CRCL = list(value = 118, dir = "lt")))
  es <- expr_strings(out)
  expect_true(all(grepl("& CRCL < 118$", es[grepl("^WT", es)])))
})

test_that("multiple additionalCovs cross-condition each other but not themselves", {
  out <- setupCovExpressionsList(mock_data, "WT", contSplit = "median",
                                 additionalCovs = list(FOOD = 1,
                                                       CRCL = list(prob = 0.25, dir = "gt")))
  es <- setNames(expr_strings(out), names(out$covExpressionsList))

  wt <- es[names(es) == "WT"]
  expect_true(all(grepl("& FOOD == 1 & CRCL > ", wt)))

  food <- es[names(es) == "FOOD"]
  expect_true(all(grepl("& CRCL > ", food)))
  expect_false(any(grepl("FOOD == 1 & FOOD", food)))

  crcl <- es[names(es) == "CRCL"]
  expect_true(all(grepl("& FOOD == 1$", crcl)))
  expect_false(any(grepl("CRCL > .* & CRCL", crcl)))

  expect_equal(unique(names(out$covExpressionsList)), c("WT", "FOOD", "CRCL"))
})

test_that("probs and minSubjects are validated", {
  expect_error(setupCovExpressionsList(mock_data, "WT", probs = 0.5),
               "length 2")
  expect_error(setupCovExpressionsList(mock_data, "WT", minSubjects = 0),
               "single number")
  expect_error(setupCovExpressionsList(mock_data, "WT", minSubjects = c(1, 2)),
               "single number")
})

test_that("additionalCovs argument is validated", {
  expect_error(setupCovExpressionsList(mock_data, "WT", additionalCovs = list(1)),
               "named list")
  expect_error(
    setupCovExpressionsList(mock_data, "WT",
                            additionalCovs = list(CRCL = 100)),
    "continuous covariate"
  )
  expect_error(
    setupCovExpressionsList(mock_data, "WT",
                            additionalCovs = list(CRCL = list(prob = 0.5))),
    "dir"
  )
  expect_error(
    setupCovExpressionsList(mock_data, "WT",
                            additionalCovs = list(CRCL = list(prob = 1.5, dir = "gt"))),
    "in \\(0, 1\\)"
  )
  expect_error(
    setupCovExpressionsList(mock_data, "WT",
                            additionalCovs = list(CRCL = list(prob = 0.5, value = 1, dir = "gt"))),
    "exactly one of"
  )
  expect_error(
    setupCovExpressionsList(mock_data, "WT",
                            additionalCovs = list(SEX = list(prob = 0.5, dir = "gt"))),
    "categorical covariate"
  )
  expect_error(
    setupCovExpressionsList(mock_data, "WT", additionalCovs = list(FOOD = 9)),
    "not present in the data"
  )
  expect_error(
    setupCovExpressionsList(mock_data, c("WT", "FOOD"),
                            additionalCovs = list(FOOD = 1)),
    "both primary and additional"
  )
})

test_that("minSubjects stops when an expression selects too few subjects", {
  expect_error(
    setupCovExpressionsList(mock_data, "WT", probs = c(0.02, 0.98)),
    "minSubjects"
  )
  # a restrictive additionalCov condition that empties a primary subset
  expect_error(
    setupCovExpressionsList(mock_data, "WT", probs = c(0.05, 0.95),
                            additionalCovs = list(CRCL = list(value = 157, dir = "gt"))),
    "minSubjects"
  )
  # lowering the floor lets the same call through
  out <- setupCovExpressionsList(mock_data, "WT", probs = c(0.02, 0.98),
                                 minSubjects = 1)
  expect_type(out$covExpressionsList, "list")

  # a continuous additionalCov condition that empties the primary subsets
  expect_error(
    setupCovExpressionsList(mock_data, "SEX",
                            additionalCovs = list(WT = list(value = 118, dir = "gt"))),
    "minSubjects"
  )
})

test_that("labels follow the documented terse format", {
  out <- setupCovExpressionsList(mock_data, c("SEX", "GENO"))
  expect_equal(out$cdfCovsNames, c("SEX 1", "SEX 2", "GENO 1", "GENO 2",
                                   "GENO 3", "GENO 4"))
})

test_that("duplicate subject records do not skew the quantiles or the level set", {
  no_dup <- mock_data[!duplicated(mock_data$ID), ]
  expect_equal(setupCovExpressionsList(mock_data, c("WT", "GENO")),
               setupCovExpressionsList(no_dup, c("WT", "GENO")))
})

test_that("missing values are excluded, including a custom missVal", {
  d <- mock_data
  d$AGE[d$ID %in% 1:6] <- -99
  out <- setupCovExpressionsList(d, "AGE")
  expect_false(any(grepl("-99", expr_strings(out))))

  d2 <- mock_data
  d2$AGE[d2$ID %in% 1:6] <- -999
  out2 <- setupCovExpressionsList(d2, "AGE", missVal = -999)
  expect_false(any(grepl("-999", expr_strings(out2))))
})

test_that("NA in a covariate column is dropped, not turned into a level", {
  d <- mock_data
  d$GENO[d$ID %in% 1:3] <- NA
  out <- setupCovExpressionsList(d, "GENO")
  expect_false(any(grepl("NA", expr_strings(out))))
  expect_equal(expr_strings(out),
               c("GENO == 1", "GENO == 2", "GENO == 3", "GENO == 4"))
})

test_that("covariate errors surface clearly", {
  expect_error(setupCovExpressionsList(mock_data, "NOPE"),
               "present in the data")

  d <- mock_data
  d$AGE <- -99
  expect_error(setupCovExpressionsList(d, "AGE"), "only missing values")
})

test_that("a single-level covariate warns and emits one row", {
  d <- mock_data
  d$FORM <- 1
  expect_warning(out <- setupCovExpressionsList(d, "FORM", minSubjects = 1),
                 "only one non-missing level")
  expect_equal(expr_strings(out), "FORM == 1")
})

test_that("tibble and data.frame inputs behave identically", {
  skip_if_not_installed("tibble")
  expect_equal(
    setupCovExpressionsList(tibble::as_tibble(mock_data), c("WT", "SEX", "GENO")),
    setupCovExpressionsList(mock_data, c("WT", "SEX", "GENO"))
  )
})

test_that("output feeds straight into getForestDFemp()", {
  dfData <- read.csv(
    system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")
  )
  dfDataEmp <- dfData[!duplicated(dfData$ID), ]
  dfSamples <- getSamples(
    system.file("extdata", "SimVal/run7.cov", package = "PMXForest"),
    system.file("extdata", "SimVal/run7.ext", package = "PMXForest"),
    n = 15
  )
  pf  <- function(thetas, df, ...) {
    if (df$WT != -99) return(list(CL = thetas[4] * (df$WT / 75)^thetas[2]))
    list(CL = thetas[4])
  }
  out <- setupCovExpressionsList(dfData, c("WT", "SEX", "GENO"), idVar = "ID")

  res <- getForestDFemp(
    dfData             = dfDataEmp,
    covExpressionsList = out$covExpressionsList,
    cdfCovsNames       = out$cdfCovsNames,
    functionList       = list(pf),
    functionListName   = "CL",
    noBaseThetas       = 14,
    dfParameters       = dfSamples,
    ncores             = 1
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), length(out$covExpressionsList))
  expect_equal(as.character(unique(res$GROUPNAME)), c("WT", "SEX", "GENO"))
  # the two WT tail rows summarise different subjects -> different POINT
  wt <- res[res$GROUPNAME == "WT", ]
  expect_false(isTRUE(all.equal(wt$POINT[1], wt$POINT[2])))
})

test_that("median split keeps every empirical subset non-empty in the round trip", {
  dfData <- read.csv(
    system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")
  )
  dfDataEmp <- dfData[!duplicated(dfData$ID), ]
  dfSamples <- getSamples(
    system.file("extdata", "SimVal/run7.cov", package = "PMXForest"),
    system.file("extdata", "SimVal/run7.ext", package = "PMXForest"),
    n = 10
  )
  pf  <- function(thetas, df, ...) list(CL = thetas[4])
  out <- setupCovExpressionsList(dfData, c("WT", "AGE"), contSplit = "median",
                                 idVar = "ID")

  expect_error(
    getForestDFemp(
      dfData             = dfDataEmp,
      covExpressionsList = out$covExpressionsList,
      cdfCovsNames       = out$cdfCovsNames,
      functionList       = list(pf),
      functionListName   = "CL",
      noBaseThetas       = 14,
      dfParameters       = dfSamples,
      ncores             = 1
    ),
    NA
  )
})

test_that("deterministic reference row from setupDfRefRow() drives getForestDFemp()", {
  dfData <- read.csv(
    system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")
  )
  dfDataEmp <- dfData[!duplicated(dfData$ID), ]
  dfSamples <- getSamples(
    system.file("extdata", "SimVal/run7.cov", package = "PMXForest"),
    system.file("extdata", "SimVal/run7.ext", package = "PMXForest"),
    n = 10
  )
  pf     <- function(thetas, df, ...) list(CL = thetas[4])
  covs   <- c("WT", "SEX")
  dfCovs <- setupDfCovs(dfData, covariates = covs, idVar = "ID")
  dfRef  <- setupDfRefRow(dfCovs, dfData, covariates = covs, singleRef = TRUE,
                          idVar = "ID")
  out    <- setupCovExpressionsList(dfData, covariates = covs, idVar = "ID")

  res <- getForestDFemp(
    dfData             = dfDataEmp,
    covExpressionsList = out$covExpressionsList,
    cdfCovsNames       = out$cdfCovsNames,
    functionList       = list(pf),
    functionListName   = "CL",
    noBaseThetas       = 14,
    dfParameters       = dfSamples,
    dfRefRow           = dfRef,
    ncores             = 1
  )
  expect_true(all(res$REFROW == "YES"))
})
