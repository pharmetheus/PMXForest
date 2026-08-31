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
  p_func <- function(thetas, df, ...) thetas[1] * (df$WT/70)

  # 2. Test ncores = 1 and auto-naming (Lines 123-139)
  res_auto <- getForestDFSCM(dfCovs = df_covs,
                             cdfCovsNames = NULL,
                             functionList = list(p_func),
                             functionListName = "CL",
                             noBaseThetas = 2,
                             dfParameters = df_params,
                             ncores = 1)

  expect_s3_class(res_auto, "data.frame")
  expect_s3_class(res_auto$GROUPNAME, "factor") # Line 210
  expect_true(any(grepl("WT=70", res_auto$COVNAME))) # Line 183

  # 3. Test Custom Reference Row (Lines 102-104)
  df_ref <- data.frame(WT = 70, AGE = 25)
  res_ref <- getForestDFSCM(dfCovs = df_covs,
                            dfRefRow = df_ref,
                            functionList = list(p_func),
                            functionListName = "CL",
                            noBaseThetas = 2,
                            dfParameters = df_params)

  expect_equal(unique(res_ref$REFROW), "YES") # Line 206
})

test_that("getForestDFSCM keeps relative CI endpoints ordered when paramFunction is negative", {
  # Uncertainty samples with a clear spread (row 1 = final estimates)
  df_params <- data.frame(
    THETA1 = c(10, 5, 7, 9, 11, 13, 15),
    THETA2 = 1
  )

  df_covs <- data.frame(WT = c(70, 100), AGE = c(25, 25), stringsAsFactors = FALSE)
  df_ref  <- data.frame(WT = 70, AGE = 25)

  # Function returns a NEGATIVE value; reference value is therefore negative too
  p_func_neg <- function(thetas, df, ...) -thetas[1] * (df$WT / 70)
  p_func_pos <- function(thetas, df, ...)  thetas[1] * (df$WT / 70)

  res_neg <- getForestDFSCM(dfCovs           = df_covs,
                            dfRefRow         = df_ref,
                            functionList     = list(p_func_neg),
                            functionListName = "CL",
                            noBaseThetas     = 2,
                            dfParameters     = df_params,
                            probs            = c(0.05, 0.95))

  res_pos <- getForestDFSCM(dfCovs           = df_covs,
                            dfRefRow         = df_ref,
                            functionList     = list(p_func_pos),
                            functionListName = "CL",
                            noBaseThetas     = 2,
                            dfParameters     = df_params,
                            probs            = c(0.05, 0.95))

  # Absolute quantile columns are always ascending (Q1 = lower prob, Q2 = upper prob)
  expect_true(all(res_neg$Q1 <= res_neg$Q2))

  # The relative CI columns must also stay ascending: Q1_REL_* is the lower limit,
  # Q2_REL_* the upper limit, since positions 1 and 2 of `probs` are used as the
  # plotted uncertainty. This holds for the positive function ...
  expect_true(all(res_pos$Q1_REL_REFFUNC  <= res_pos$Q2_REL_REFFUNC))
  expect_true(all(res_pos$Q1_REL_REFFINAL <= res_pos$Q2_REL_REFFINAL))

  # ... and must equally hold for the negative function. Dividing the ascending
  # absolute quantiles by a negative reference reverses their order, so without
  # the fix Q1_REL_* ends up above Q2_REL_* and the forest CI is drawn reversed.
  expect_true(all(res_neg$Q1_REL_REFFUNC  <= res_neg$Q2_REL_REFFUNC))
  expect_true(all(res_neg$Q1_REL_REFFINAL <= res_neg$Q2_REL_REFFINAL))
})

test_that("getForestDFSCM parallel logic and input conversion", {
  df_params <- data.frame(THETA1 = 10, THETA2 = 2)
  p_func <- function(thetas, df, ...) thetas[1]

  # 4. Test ncores=2 (Lines 114-120)
  # Use manual names to skip the auto-naming helper if desired
  df_covs_par <- data.frame(WT = 70, AGE = 25)
  res_par <- getForestDFSCM(dfCovs = df_covs_par,
                            cdfCovsNames = "ManualName",
                            functionList = list(p_func),
                            noBaseThetas = 2,
                            dfParameters = df_params,
                            ncores = 2)

  expect_equal(as.character(res_par$COVNAME[1]), "ManualName")

  # 5. Test input conversion from list (Line 79)
  # FIX: Provide at least two covariates in the list to ensure ncol > 0
  cov_list <- list("WT" = c(70, 80), "AGE" = c(30, 40))
  res_list <- getForestDFSCM(dfCovs = cov_list,
                             functionList = list(p_func),
                             noBaseThetas = 2,
                             dfParameters = df_params)

  expect_s3_class(res_list, "data.frame")
  expect_gt(nrow(res_list), 0)
})

# --- oneHot argument ---

test_that("getForestDFSCM oneHot matches manually pre-encoded dfCovs", {
  df_params <- data.frame(THETA1 = c(10, 11, 9, 12), THETA2 = c(0.3, 0.35, 0.28, 0.31))

  # dfCovs with a raw multi-level GENO column, reference level 2
  df_covs_raw <- data.frame(
    GENO            = c(1, 3, 4),
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

  df_covs_manual <- oneHotEncode(df_covs_raw, spec = list(GENO = list(ref = 2)),
                                 dropOriginal = TRUE)
  res_manual <- getForestDFSCM(
    dfCovs = df_covs_manual, functionList = list(p_func), functionListName = "CL",
    noBaseThetas = 2, dfParameters = df_params
  )

  expect_equal(res_onehot, res_manual)
})

test_that("getForestDFSCM oneHot = NULL leaves the result unchanged", {
  df_params <- data.frame(THETA1 = c(10, 11), THETA2 = c(2, 2.1))
  df_covs <- data.frame(WT = c(70, 100), AGE = c(25, 50),
                        COVARIATEGROUPS = c("Weight", "Age"), stringsAsFactors = FALSE)
  p_func <- function(thetas, df, ...) thetas[1] * (df$WT / 70)

  a <- getForestDFSCM(dfCovs = df_covs, functionList = list(p_func),
                      functionListName = "CL", noBaseThetas = 2, dfParameters = df_params)
  b <- getForestDFSCM(dfCovs = df_covs, functionList = list(p_func),
                      functionListName = "CL", noBaseThetas = 2, dfParameters = df_params,
                      oneHot = NULL)
  expect_equal(a, b)
})

test_that("getForestDFSCM oneHotSep controls the dummy column separator", {
  df_params <- data.frame(THETA1 = c(10, 11), THETA2 = c(0.3, 0.32))
  df_covs   <- data.frame(GENO = c(1, 3, 4), COVARIATEGROUPS = "GENO",
                          stringsAsFactors = FALSE)
  p_func <- function(thetas, df, ...) {
    x <- thetas[1]
    if (isTRUE(df$GENO1 == 1)) x <- x * (1 + thetas[2])
    x
  }
  res <- getForestDFSCM(dfCovs = df_covs, functionList = list(p_func),
                        functionListName = "CL", noBaseThetas = 2,
                        dfParameters = df_params,
                        oneHot = list(GENO = list(ref = 2)), oneHotSep = "")
  ref <- getForestDFSCM(
    dfCovs = oneHotEncode(df_covs, spec = list(GENO = list(ref = 2)),
                          sep = "", dropOriginal = TRUE),
    functionList = list(p_func), functionListName = "CL", noBaseThetas = 2,
    dfParameters = df_params
  )
  expect_equal(res, ref)
})

# --- coverage: dfRefRow validation and oneHot + data.frame dfRefRow ---

test_that("getForestDFSCM errors on a dfRefRow with an invalid number of rows", {
  df_params <- data.frame(THETA1 = c(10, 11), THETA2 = c(2, 2.1))
  df_covs   <- data.frame(WT = c(70, 100), COVARIATEGROUPS = "WT",
                          stringsAsFactors = FALSE)
  bad_ref   <- data.frame(WT = c(70, 80, 90))   # 3 rows, dfCovs has 2

  expect_error(
    getForestDFSCM(dfCovs = df_covs, dfRefRow = bad_ref,
                   functionList = list(function(thetas, df, ...) thetas[1]),
                   functionListName = "CL", noBaseThetas = 2, dfParameters = df_params),
    "number of reference rows"
  )
})

test_that("getForestDFSCM oneHot also encodes a data.frame dfRefRow", {
  df_params <- data.frame(THETA1 = c(10, 11, 9, 12), THETA2 = c(0.3, 0.35, 0.28, 0.31))
  df_covs   <- data.frame(GENO = c(1, 3, 4), COVARIATEGROUPS = "GENO",
                          stringsAsFactors = FALSE)
  df_ref    <- data.frame(GENO = 2)             # raw column, needs encoding too

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
    dfCovs   = oneHotEncode(df_covs, spec = spec, dropOriginal = TRUE),
    dfRefRow = oneHotEncode(df_ref, spec = spec, dropOriginal = TRUE),
    functionList = list(p_func), functionListName = "CL",
    noBaseThetas = 2, dfParameters = df_params
  )
  expect_equal(res_onehot, res_manual)
  # The reference row was encoded: reference genotype 2 -> all GENO_* dummies 0
  expect_equal(unique(res_onehot$REFROW), "YES")
})

test_that("getForestDFSCM accepts a tibble dfCovs / dfRefRow (issue: tibble `[` does not drop)", {
  skip_if_not_installed("tibble")

  df_params <- data.frame(THETA1 = c(10, 11, 9), THETA2 = c(2, 2.1, 1.9))
  df_covs   <- createInputForestData(list(WT = c(60, 115), AGE = c(30, 60)))
  df_ref    <- data.frame(WT = 75, AGE = 45)
  p_func    <- function(thetas, df, ...) list(CL = thetas[1] * (df$WT / 70))

  ref <- getForestDFSCM(df_covs, dfRefRow = df_ref, functionList = list(p_func),
                        functionListName = "CL", noBaseThetas = 2, dfParameters = df_params)

  tbl <- getForestDFSCM(tibble::as_tibble(df_covs),
                        dfRefRow = tibble::as_tibble(df_ref),
                        functionList = list(p_func), functionListName = "CL",
                        noBaseThetas = 2, dfParameters = df_params)

  expect_equal(tbl, ref)
})

test_that("getForestDFSCM works with a single covariate (old ticket: column named 'dfCovs[i, ]')", {
  df_params <- data.frame(THETA1 = c(10, 11, 9), THETA2 = c(2, 2.1, 1.9))
  df_covs   <- createInputForestData(list(LIVSTBL = c(0, 1)))   # one covariate only
  p_func    <- function(thetas, df, ...) list(CL = thetas[1])

  res <- getForestDFSCM(df_covs, functionList = list(p_func), functionListName = "CL",
                        noBaseThetas = 2, dfParameters = df_params)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  # The covariate column keeps its name, not "dfCovs[i, ]"
  expect_true("LIVSTBL" %in% names(res))
  expect_false(any(grepl("dfCovs", names(res))))
  expect_equal(as.character(unique(res$GROUPNAME)), "LIVSTBL")
  expect_setequal(res$COVNAME, c("LIVSTBL=0", "LIVSTBL=1"))

  # supplied cdfCovsNames are used
  res2 <- getForestDFSCM(df_covs, cdfCovsNames = c("Normal", "Impaired"),
                         functionList = list(p_func), functionListName = "CL",
                         noBaseThetas = 2, dfParameters = df_params)
  expect_setequal(res2$COVNAME, c("Normal", "Impaired"))
})
