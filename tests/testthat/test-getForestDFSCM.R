test_that("getForestDFSCM handles parallel setup and auto-naming", {
  # 1. Setup mock parameters and covariates
  df_params <- data.frame(THETA1 = c(10, 11), THETA2 = c(2, 2.1))

  # Ensure multiple columns to prevent indexing issues
  df_covs <- data.frame(
    WT = c(70, 100),
    AGE = c(25, 50),
    COVARIATEGROUPS = c("Weight", "Age"),
    stringsAsFactors = FALSE
  )

  # Function must use 'thetas' to match internal call
  p_func <- function(thetas, df, ...) thetas[1] * (df$WT / 70)

  # 2. Test ncores = 1 and auto-naming (Lines 123-139)
  res_auto <- getForestDFSCM(
    dfCovs = df_covs,
    cdfCovsNames = NULL,
    functionList = list(p_func),
    functionListName = "CL",
    noBaseThetas = 2,
    dfParameters = df_params,
    ncores = 1
  )

  expect_s3_class(res_auto, "data.frame")
  expect_s3_class(res_auto$GROUPNAME, "factor") # Line 210
  expect_true(any(grepl("WT=70", res_auto$COVNAME))) # Line 183

  # 3. Test Custom Reference Row (Lines 102-104)
  df_ref <- data.frame(WT = 70, AGE = 25)
  res_ref <- getForestDFSCM(
    dfCovs = df_covs,
    dfRefRow = df_ref,
    functionList = list(p_func),
    functionListName = "CL",
    noBaseThetas = 2,
    dfParameters = df_params
  )

  expect_equal(unique(res_ref$REFROW), "YES") # Line 206
})

test_that("getForestDFSCM keeps relative CI endpoints ordered when paramFunction is negative", {
  # Uncertainty samples with a clear spread (row 1 = final estimates)
  df_params <- data.frame(
    THETA1 = c(10, 5, 7, 9, 11, 13, 15),
    THETA2 = 1
  )

  df_covs <- data.frame(WT = c(70, 100), AGE = c(25, 25), stringsAsFactors = FALSE)
  df_ref <- data.frame(WT = 70, AGE = 25)

  # Function returns a NEGATIVE value; reference value is therefore negative too
  p_func_neg <- function(thetas, df, ...) -thetas[1] * (df$WT / 70)
  p_func_pos <- function(thetas, df, ...) thetas[1] * (df$WT / 70)

  res_neg <- getForestDFSCM(
    dfCovs = df_covs,
    dfRefRow = df_ref,
    functionList = list(p_func_neg),
    functionListName = "CL",
    noBaseThetas = 2,
    dfParameters = df_params,
    probs = c(0.05, 0.95)
  )

  res_pos <- getForestDFSCM(
    dfCovs = df_covs,
    dfRefRow = df_ref,
    functionList = list(p_func_pos),
    functionListName = "CL",
    noBaseThetas = 2,
    dfParameters = df_params,
    probs = c(0.05, 0.95)
  )

  # Absolute quantile columns are always ascending (Q1 = lower prob, Q2 = upper prob)
  expect_true(all(res_neg$Q1 <= res_neg$Q2))

  # The relative CI columns must also stay ascending: Q1_REL_* is the lower limit,
  # Q2_REL_* the upper limit, since positions 1 and 2 of `probs` are used as the
  # plotted uncertainty. This holds for the positive function ...
  expect_true(all(res_pos$Q1_REL_REFFUNC <= res_pos$Q2_REL_REFFUNC))
  expect_true(all(res_pos$Q1_REL_REFFINAL <= res_pos$Q2_REL_REFFINAL))

  # ... and must equally hold for the negative function. Dividing the ascending
  # absolute quantiles by a negative reference reverses their order, so without
  # the fix Q1_REL_* ends up above Q2_REL_* and the forest CI is drawn reversed.
  expect_true(all(res_neg$Q1_REL_REFFUNC <= res_neg$Q2_REL_REFFUNC))
  expect_true(all(res_neg$Q1_REL_REFFINAL <= res_neg$Q2_REL_REFFINAL))
})

test_that("getForestDFSCM parallel logic and input conversion", {
  df_params <- data.frame(THETA1 = 10, THETA2 = 2)
  p_func <- function(thetas, df, ...) thetas[1]

  # 4. Test ncores=2 (Lines 114-120)
  # Use manual names to skip the auto-naming helper if desired
  df_covs_par <- data.frame(WT = 70, AGE = 25)
  res_par <- getForestDFSCM(
    dfCovs = df_covs_par,
    cdfCovsNames = "ManualName",
    functionList = list(p_func),
    noBaseThetas = 2,
    dfParameters = df_params,
    ncores = 2
  )

  expect_equal(as.character(res_par$COVNAME[1]), "ManualName")

  # 5. Test input conversion from list (Line 79)
  # FIX: Provide at least two covariates in the list to ensure ncol > 0
  cov_list <- list("WT" = c(70, 80), "AGE" = c(30, 40))
  res_list <- getForestDFSCM(
    dfCovs = cov_list,
    functionList = list(p_func),
    noBaseThetas = 2,
    dfParameters = df_params
  )

  expect_s3_class(res_list, "data.frame")
  expect_gt(nrow(res_list), 0)
})

# --- oneHot argument ---

test_that("getForestDFSCM oneHot matches manually pre-encoded dfCovs", {
  df_params <- data.frame(THETA1 = c(10, 11, 9, 12), THETA2 = c(0.3, 0.35, 0.28, 0.31))

  # dfCovs with a raw multi-level GENO column, reference level 2
  df_covs_raw <- data.frame(
    GENO = c(1, 3, 4),
    COVARIATEGROUPS = "GENO",
    stringsAsFactors = FALSE
  )

  p_func <- function(thetas, df, ...) {
    x <- thetas[1]
    if (isTRUE(df$GENO_1 == 1)) x <- x * (1 + thetas[2])
    if (isTRUE(df$GENO_3 == 1)) x <- x * (1 - thetas[2])
    if (isTRUE(df$GENO_4 == 1)) x <- x * (1 + 2 * thetas[2])
    x
  }

  res_onehot <- getForestDFSCM(
    dfCovs = df_covs_raw, functionList = list(p_func), functionListName = "CL",
    noBaseThetas = 2, dfParameters = df_params, oneHot = list(GENO = list(ref = 2))
  )

  df_covs_manual <- oneHotEncode(df_covs_raw,
    spec = list(GENO = list(ref = 2)),
    dropOriginal = TRUE
  )
  res_manual <- getForestDFSCM(
    dfCovs = df_covs_manual, functionList = list(p_func), functionListName = "CL",
    noBaseThetas = 2, dfParameters = df_params
  )

  expect_equal(res_onehot, res_manual)
})

test_that("getForestDFSCM oneHot = NULL leaves the result unchanged", {
  df_params <- data.frame(THETA1 = c(10, 11), THETA2 = c(2, 2.1))
  df_covs <- data.frame(
    WT = c(70, 100), AGE = c(25, 50),
    COVARIATEGROUPS = c("Weight", "Age"), stringsAsFactors = FALSE
  )
  p_func <- function(thetas, df, ...) thetas[1] * (df$WT / 70)

  a <- getForestDFSCM(
    dfCovs = df_covs, functionList = list(p_func),
    functionListName = "CL", noBaseThetas = 2, dfParameters = df_params
  )
  b <- getForestDFSCM(
    dfCovs = df_covs, functionList = list(p_func),
    functionListName = "CL", noBaseThetas = 2, dfParameters = df_params,
    oneHot = NULL
  )
  expect_equal(a, b)
})

test_that("getForestDFSCM oneHotSep controls the dummy column separator", {
  df_params <- data.frame(THETA1 = c(10, 11), THETA2 = c(0.3, 0.32))
  df_covs <- data.frame(
    GENO = c(1, 3, 4), COVARIATEGROUPS = "GENO",
    stringsAsFactors = FALSE
  )
  p_func <- function(thetas, df, ...) {
    x <- thetas[1]
    if (isTRUE(df$GENO1 == 1)) x <- x * (1 + thetas[2])
    x
  }
  res <- getForestDFSCM(
    dfCovs = df_covs, functionList = list(p_func),
    functionListName = "CL", noBaseThetas = 2,
    dfParameters = df_params,
    oneHot = list(GENO = list(ref = 2)), oneHotSep = ""
  )
  ref <- getForestDFSCM(
    dfCovs = oneHotEncode(df_covs,
      spec = list(GENO = list(ref = 2)),
      sep = "", dropOriginal = TRUE
    ),
    functionList = list(p_func), functionListName = "CL", noBaseThetas = 2,
    dfParameters = df_params
  )
  expect_equal(res, ref)
})

# --- coverage: dfRefRow validation and oneHot + data.frame dfRefRow ---

test_that("getForestDFSCM errors on a dfRefRow with an invalid number of rows", {
  df_params <- data.frame(THETA1 = c(10, 11), THETA2 = c(2, 2.1))
  df_covs <- data.frame(
    WT = c(70, 100), COVARIATEGROUPS = "WT",
    stringsAsFactors = FALSE
  )
  bad_ref <- data.frame(WT = c(70, 80, 90)) # 3 rows, dfCovs has 2

  expect_error(
    getForestDFSCM(
      dfCovs = df_covs, dfRefRow = bad_ref,
      functionList = list(function(thetas, df, ...) thetas[1]),
      functionListName = "CL", noBaseThetas = 2, dfParameters = df_params
    ),
    "number of reference rows"
  )
})

test_that("functionList may hold several functions, named end to end", {
  # functionList is a list because it can carry more than one function. Each is
  # evaluated for every dfCovs row and every parameter vector, and the returns
  # are concatenated in order, so functionListName names all of them end to end.
  df_params <- data.frame(THETA1 = c(10, 11), THETA2 = c(2, 2.1))
  df_covs <- data.frame(
    WT = c(70, 100), COVARIATEGROUPS = "WT",
    stringsAsFactors = FALSE
  )

  pk <- function(thetas, df, ...) list(thetas[1], thetas[2])
  exposure <- function(thetas, df, ...) list(80 / thetas[1])

  both <- getForestDFSCM(
    dfCovs = df_covs,
    functionList = list(pk, exposure),
    functionListName = c("CL", "V", "AUC"),
    noBaseThetas = 2, dfParameters = df_params, ncores = 1
  )

  # PARAMETER is a factor, which is what preserves the functionListName order
  # downstream; compare the levels rather than the codes.
  expect_equal(as.character(unique(both$PARAMETER)), c("CL", "V", "AUC"))
  expect_equal(levels(both$PARAMETER), c("CL", "V", "AUC"))

  # each function sees the same thetas and row as it would alone
  alone1 <- getForestDFSCM(
    dfCovs = df_covs, functionList = list(pk),
    functionListName = c("CL", "V"),
    noBaseThetas = 2, dfParameters = df_params, ncores = 1
  )
  alone2 <- getForestDFSCM(
    dfCovs = df_covs, functionList = list(exposure),
    functionListName = "AUC",
    noBaseThetas = 2, dfParameters = df_params, ncores = 1
  )
  expect_equal(
    both$POINT[both$PARAMETER == "CL"], alone1$POINT[alone1$PARAMETER == "CL"]
  )
  expect_equal(
    both$POINT[both$PARAMETER == "AUC"], alone2$POINT[alone2$PARAMETER == "AUC"]
  )
})

test_that("functionListName must cover every return of every function", {
  df_params <- data.frame(THETA1 = c(10, 11), THETA2 = c(2, 2.1))
  df_covs <- data.frame(WT = 70, COVARIATEGROUPS = "WT", stringsAsFactors = FALSE)
  pk <- function(thetas, df, ...) list(thetas[1], thetas[2])
  exposure <- function(thetas, df, ...) list(80 / thetas[1])

  # three returns across the two functions, only two names
  expect_error(
    getForestDFSCM(
      dfCovs = df_covs, functionList = list(pk, exposure),
      functionListName = c("CL", "V"),
      noBaseThetas = 2, dfParameters = df_params, ncores = 1
    )
  )
})

test_that("getForestDFSCM oneHot also encodes a data.frame dfRefRow", {
  df_params <- data.frame(THETA1 = c(10, 11, 9, 12), THETA2 = c(0.3, 0.35, 0.28, 0.31))
  df_covs <- data.frame(
    GENO = c(1, 3, 4), COVARIATEGROUPS = "GENO",
    stringsAsFactors = FALSE
  )
  # A single-row reference carries one level, so its dummy columns can only come
  # from the encoding derived for dfCovs. GENO 1 is a level the parameter
  # function reacts to, so an unencoded reference row would be visibly wrong.
  df_ref <- data.frame(GENO = 1) # raw column, needs encoding too

  p_func <- function(thetas, df, ...) {
    x <- thetas[1]
    if (isTRUE(df$GENO_1 == 1)) x <- x * (1 + thetas[2])
    if (isTRUE(df$GENO_3 == 1)) x <- x * (1 - thetas[2])
    list(CL = x)
  }
  spec <- list(GENO = list(ref = 2))

  res_onehot <- getForestDFSCM(
    dfCovs = df_covs, dfRefRow = df_ref, functionList = list(p_func),
    functionListName = "CL", noBaseThetas = 2, dfParameters = df_params, oneHot = spec
  )
  res_manual <- getForestDFSCM(
    dfCovs = oneHotEncode(df_covs, spec = spec, dropOriginal = TRUE),
    # The dummy columns dfCovs gets, stated outright.
    dfRefRow = data.frame(GENO_1 = 1, GENO_3 = 0, GENO_4 = 0),
    functionList = list(p_func), functionListName = "CL",
    noBaseThetas = 2, dfParameters = df_params
  )
  # The reference row really was encoded: GENO 1 carries the THETA2 effect, so
  # the reference is not the bare thetas[1] an unencoded row would give.
  expect_equal(
    unique(res_onehot$REFFUNC),
    median(df_params$THETA1 * (1 + df_params$THETA2))
  )
  expect_equal(res_onehot, res_manual)
  expect_equal(unique(res_onehot$REFROW), "YES")
})

test_that("getForestDFSCM accepts a tibble dfCovs / dfRefRow (issue: tibble `[` does not drop)", {
  skip_if_not_installed("tibble")

  df_params <- data.frame(THETA1 = c(10, 11, 9), THETA2 = c(2, 2.1, 1.9))
  df_covs <- createInputForestData(list(WT = c(60, 115), AGE = c(30, 60)))
  df_ref <- data.frame(WT = 75, AGE = 45)
  p_func <- function(thetas, df, ...) list(CL = thetas[1] * (df$WT / 70))

  ref <- getForestDFSCM(df_covs,
    dfRefRow = df_ref, functionList = list(p_func),
    functionListName = "CL", noBaseThetas = 2, dfParameters = df_params
  )

  tbl <- getForestDFSCM(tibble::as_tibble(df_covs),
    dfRefRow = tibble::as_tibble(df_ref),
    functionList = list(p_func), functionListName = "CL",
    noBaseThetas = 2, dfParameters = df_params
  )

  expect_equal(tbl, ref)
})

test_that("getForestDFSCM works with a single covariate (old ticket: column named 'dfCovs[i, ]')", {
  df_params <- data.frame(THETA1 = c(10, 11, 9), THETA2 = c(2, 2.1, 1.9))
  df_covs <- createInputForestData(list(LIVSTBL = c(0, 1))) # one covariate only
  p_func <- function(thetas, df, ...) list(CL = thetas[1])

  res <- getForestDFSCM(df_covs,
    functionList = list(p_func), functionListName = "CL",
    noBaseThetas = 2, dfParameters = df_params
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  # The covariate column keeps its name, not "dfCovs[i, ]"
  expect_true("LIVSTBL" %in% names(res))
  expect_false(any(grepl("dfCovs", names(res))))
  expect_equal(as.character(unique(res$GROUPNAME)), "LIVSTBL")
  expect_setequal(res$COVNAME, c("LIVSTBL=0", "LIVSTBL=1"))

  # supplied cdfCovsNames are used
  res2 <- getForestDFSCM(df_covs,
    cdfCovsNames = c("Normal", "Impaired"),
    functionList = list(p_func), functionListName = "CL",
    noBaseThetas = 2, dfParameters = df_params
  )
  expect_setequal(res2$COVNAME, c("Normal", "Impaired"))
})

# --- dfCovs against the covariates a generated parameter function reads ----

## run7, reduced to CL, FREL and V. $PK reads SEX, GENO4, FORM, FOOD, WT, GENO1
## and GENO3; GENO2 is the reference genotype.
run7Gen <- function() {
  createParamFunction(
    system.file("extdata", "SimVal/run7.mod", package = "PMXForest"),
    parameters = c("CL", "FREL", "V"),
    covRef = list(GENO1 = 0, GENO3 = 0), quiet = TRUE
  )
}
run7Fun <- function() eval(parse(text = paste(run7Gen()$code, collapse = "\n")))

## The final estimates of run7, twice: getForestDFSCM() needs a sample.
run7Final <- function() {
  ext <- getExt(system.file("extdata", "SimVal/run7.ext", package = "PMXForest"))
  fin <- ext[ext$ITERATION == -1000000000, grep("^THETA", names(ext))]
  fin[c(1, 1), ]
}

run7Covs <- function(...) {
  dfData <- read.csv(
    system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")
  )
  setupDfCovs(dfData,
    covariates = c("WT", "SEX", "FOOD", "FORM", "GENO"),
    idVar = "ID", ...
  )
}

forestRun7 <- function(dfCovs, functionList = list(run7Fun()),
                       functionListName = c("CL", "FREL", "V"), ...) {
  getForestDFSCM(
    dfCovs = dfCovs, functionList = functionList,
    functionListName = functionListName, noBaseThetas = 14,
    dfParameters = run7Final(), ...
  )
}

test_that("indicator columns named apart from $PK are reported, not silently flat", {
  ## setupDfCovs() names the GENO indicators GENO_2, GENO_3, GENO_4; run7 reads
  ## GENO1, GENO3, GENO4. Every GENO row then comes out as "no effect".
  w <- expect_warning(
    res <- forestRun7(run7Covs()),
    class = "pmxCovariateMismatch"
  )
  expect_equal(w$unread, c("GENO_2", "GENO_3", "GENO_4"))
  expect_equal(w$absent, c("GENO4", "GENO1", "GENO3"))
  expect_match(conditionMessage(w), "GENO_2, GENO_3, GENO_4", fixed = TRUE)
  expect_match(conditionMessage(w), "GENO4, GENO1, GENO3", fixed = TRUE)
  expect_false(any(res$COVEFF[grepl("GENO", res$COVNAME)]))
})

test_that("matching covariate names give no warning", {
  ## sep = "" and the reference genotype of run7 give exactly the columns $PK
  ## reads; COVARIATEGROUPS is bookkeeping and is not reported.
  dfCovs <- run7Covs(sep = "", catRef = list(GENO = 2))
  expect_true("COVARIATEGROUPS" %in% names(dfCovs))
  expect_no_warning(res <- forestRun7(dfCovs))
  expect_true(any(res$COVEFF[grepl("GENO1=1", res$COVNAME)]))
})

test_that("a wrong reference genotype is reported", {
  ## Without catRef, setupDfCovs() leaves out GENO1 instead of GENO2: the GENO2
  ## column is unread and GENO1 is absent, so the UM row comes out flat.
  w <- expect_warning(
    res <- forestRun7(run7Covs(sep = "")),
    class = "pmxCovariateMismatch"
  )
  expect_equal(w$unread, "GENO2")
  expect_equal(w$absent, "GENO1")
  expect_match(conditionMessage(w), "catRef in setupDfCovs()", fixed = TRUE)
  expect_false(any(res$COVEFF[grepl("GENO2=1", res$COVNAME)]))
})

test_that("a column no function reads gives no warning on its own", {
  dfCovs <- run7Covs(sep = "", catRef = list(GENO = 2))
  dfCovs$AGE <- -99
  dfCovs$AGE[1] <- 40
  expect_no_warning(forestRun7(dfCovs))
})

test_that("a covariate dfCovs lacks gives no warning on its own", {
  dfCovs <- run7Covs(sep = "", catRef = list(GENO = 2))
  dfCovs$SEX <- NULL
  expect_no_warning(forestRun7(dfCovs))
})

test_that("an unread column and an absent covariate are reported together", {
  ## Both may be intended, as in the quick start; the pair is reported anyway,
  ## because a misspelled name gives the same pair.
  dfCovs <- run7Covs(sep = "", catRef = list(GENO = 2))
  dfCovs$AGE <- -99
  dfCovs$AGE[1] <- 40
  dfCovs$SEX <- NULL
  w <- expect_warning(forestRun7(dfCovs), class = "pmxCovariateMismatch")
  expect_equal(w$unread, "AGE")
  expect_equal(w$absent, "SEX")
})

## run7Covs() with one raw GENO column instead of the indicator columns, for
## getForestDFSCM() to encode.
run7RawGeno <- function() {
  dfCovs <- run7Covs(sep = "", catRef = list(GENO = 2))
  geno <- dfCovs$COVARIATEGROUPS == "GENO"
  lev <- ifelse(dfCovs$GENO1 == 1, 1,
    ifelse(dfCovs$GENO3 == 1, 3, ifelse(dfCovs$GENO4 == 1, 4, 2))
  )
  dfCovs$GENO <- ifelse(geno, lev, -99)
  dfCovs[c("GENO1", "GENO3", "GENO4")] <- NULL
  dfCovs
}

test_that("the comparison is made after the one-hot encoding of getForestDFSCM", {
  dfCovs <- run7RawGeno()

  expect_no_warning(
    forestRun7(dfCovs, oneHot = list(GENO = list(ref = 2)), oneHotSep = "")
  )
  w <- expect_warning(
    forestRun7(dfCovs, oneHot = list(GENO = list(ref = 2)), oneHotSep = "_"),
    class = "pmxCovariateMismatch"
  )
  expect_equal(w$unread, c("GENO_1", "GENO_3", "GENO_4"))
  expect_equal(w$absent, c("GENO4", "GENO1", "GENO3"))
})

test_that("a hand-written parameter function is not checked", {
  handWritten <- function(thetas, df, ...) list(CL = thetas[4])
  expect_no_warning(
    forestRun7(run7Covs(), functionList = list(handWritten), functionListName = "CL")
  )
})

test_that("next to a hand-written function, no check is made", {
  ## The hand-written function may read any column, so no column can be called
  ## unread, and an absent covariate alone is not reported.
  handWritten <- function(thetas, df, ...) list(KA = thetas[6])
  expect_no_warning(
    forestRun7(run7Covs(),
      functionList = list(run7Fun(), handWritten),
      functionListName = c("CL", "FREL", "V", "KA")
    )
  )
})

# --- dfRefRow against dfCovs -------------------------------------------------

test_that("a column of dfRefRow that dfCovs lacks is reported, and shifts every row", {
  ## wt for WT: the reference subject weighs 75 kg, the reference of $PK, and
  ## not 60 kg, so every row is divided by the wrong reference.
  dfCovs <- run7Covs(sep = "", catRef = list(GENO = 2))
  expect_no_warning(ok <- forestRun7(dfCovs, dfRefRow = data.frame(WT = 60)))
  w <- expect_warning(
    typo <- forestRun7(dfCovs, dfRefRow = data.frame(wt = 60)),
    class = "pmxCovariateMismatch"
  )
  expect_equal(w$refOnly, "wt")
  expect_match(conditionMessage(w), "dfRefRow has columns that dfCovs lacks: wt.", fixed = TRUE)
  expect_match(conditionMessage(w), "setupDfRefRow()", fixed = TRUE)
  expect_equal(typo$REFFUNC, forestRun7(dfCovs)$REFFUNC)
  expect_false(isTRUE(all.equal(typo$REFFUNC, ok$REFFUNC)))
})

test_that("dfRefRow is checked next to a hand-written parameter function too", {
  handWritten <- function(thetas, df, ...) list(CL = thetas[4])
  w <- expect_warning(
    forestRun7(run7Covs(sep = "", catRef = list(GENO = 2)),
      functionList = list(handWritten), functionListName = "CL",
      dfRefRow = data.frame(wt = 60)
    ),
    class = "pmxCovariateMismatch"
  )
  expect_equal(w$refOnly, "wt")
})

test_that("a dfRefRow from setupDfRefRow() gives no warning", {
  ## The reference row has the columns of dfCovs, COVARIATEGROUPS included,
  ## which getForestDFSCM() takes out of dfCovs before the comparison.
  dfCovs <- run7Covs(sep = "", catRef = list(GENO = 2))
  dfData <- read.csv(
    system.file("extdata", "SimVal/DAT-1-MI-PMX-2.csv", package = "PMXForest")
  )
  dfRefRow <- setupDfRefRow(dfCovs, dfData,
    covariates = c("WT", "SEX", "FOOD", "FORM", "GENO"), idVar = "ID"
  )
  expect_true("COVARIATEGROUPS" %in% names(dfRefRow))
  expect_no_warning(forestRun7(dfCovs, dfRefRow = dfRefRow))
})

test_that("dfRefRow is compared with dfCovs after the one-hot encoding", {
  ## A raw GENO column in dfRefRow becomes the indicator columns of dfCovs,
  ## and indicator columns already in dfRefRow match those of the encoded
  ## dfCovs, though not the raw GENO column dfCovs was given.
  for (ref in list(data.frame(GENO = 1, WT = 60), data.frame(GENO1 = 1, GENO3 = 0, GENO4 = 0, WT = 60))) {
    expect_no_warning(
      forestRun7(run7RawGeno(),
        dfRefRow = ref,
        oneHot = list(GENO = list(ref = 2)), oneHotSep = ""
      )
    )
  }
})
