test_that("getForestDFemp handles observed references and expressions", {
  # 1. Setup Mock Data
  df_params <- data.frame(THETA1 = c(10, 11), THETA2 = c(2, 2.1))
  df_data <- data.frame(ID = 1:4, WT = c(60, 70, 80, 90), SEX = c(1, 1, 2, 2))

  # Expression list: grouped by Sex
  ls_expr <- list("Sex" = expression(SEX == 1), "Sex" = expression(SEX == 2))

  # FIX: Function must use 'thetas' to match internal call at Line 120
  p_func <- function(thetas, df, ...) thetas[1] * (df$WT / 70)

  # 2. Test Auto-Reference from Data (Lines 155-159)
  # NULL dfRefRow triggers metricFunction summary of typical predictions
  res_auto <- getForestDFemp(dfData = df_data,
                             covExpressionsList = ls_expr,
                             noBaseThetas = 2,
                             dfParameters = df_params,
                             functionList = list(p_func),
                             dfRefRow = NULL)

  expect_s3_class(res_auto, "data.frame")
  expect_equal(unique(as.character(res_auto$REFROW)), "NO") # Line 248

  # 3. Test Expression-based Reference (Lines 163-173)
  # Exercises the is.expression(dfRefRow[[1]]) branch
  ref_expr <- list(expression(WT < 80))
  res_exp_ref <- getForestDFemp(dfData = df_data,
                                covExpressionsList = ls_expr,
                                noBaseThetas = 2,
                                dfParameters = df_params,
                                functionList = list(p_func),
                                dfRefRow = ref_expr)

  expect_equal(unique(as.character(res_exp_ref$REFROW)), "YES")

  # 4. Test parallelization setup (Lines 114-116)
  res_par <- getForestDFemp(dfData = df_data,
                            covExpressionsList = ls_expr,
                            noBaseThetas = 2,
                            dfParameters = df_params,
                            functionList = list(p_func),
                            ncores = 2)
  expect_s3_class(res_par, "data.frame")
})

test_that("getForestDFemp keeps relative CI endpoints ordered when paramFunction is negative", {
  # Uncertainty samples with a clear spread (row 1 = final estimates)
  df_params <- data.frame(
    THETA1 = c(10, 5, 7, 9, 11, 13, 15),
    THETA2 = 1
  )

  df_data <- data.frame(
    ID  = 1:6,
    WT  = c(55, 62, 68, 78, 85, 95),
    SEX = c(1, 1, 1, 2, 2, 2)
  )

  ls_expr <- list("WT" = expression(WT < 70), "WT" = expression(WT >= 70))
  df_ref  <- data.frame(WT = 70, SEX = 1)

  # Function returns a NEGATIVE value; the reference value is therefore negative too
  p_func_neg <- function(thetas, df, ...) -thetas[1] * (df$WT / 70)
  p_func_pos <- function(thetas, df, ...)  thetas[1] * (df$WT / 70)

  res_neg <- getForestDFemp(dfData             = df_data,
                            covExpressionsList = ls_expr,
                            noBaseThetas       = 2,
                            dfParameters       = df_params,
                            functionList       = list(p_func_neg),
                            functionListName   = "CL",
                            dfRefRow           = df_ref,
                            probs              = c(0.05, 0.95))

  res_pos <- getForestDFemp(dfData             = df_data,
                            covExpressionsList = ls_expr,
                            noBaseThetas       = 2,
                            dfParameters       = df_params,
                            functionList       = list(p_func_pos),
                            functionListName   = "CL",
                            dfRefRow           = df_ref,
                            probs              = c(0.05, 0.95))

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

test_that("getForestDFemp error handling", {
  # FIX: Ensure 'thetas' is used even in error-test functions
  p_func <- function(thetas, df, ...) thetas[1]
  df_params <- data.frame(THETA1 = 10, THETA2 = 2)
  df_data <- data.frame(ID = 1:2, WT = c(70, 80), SEX = c(1, 2))
  ls_expr <- list(expression(WT == 70))

  # Input length mismatch (Line 92)
  expect_error(
    getForestDFemp(df_data, ls_expr, cdfCovsNames = c("L1", "L2"),
                   noBaseThetas = 2, dfParameters = df_params, functionList = list(p_func)),
    "cdfCovsNames should have the same length"
  )

  # No data for subset (Line 167)
  expect_error(
    getForestDFemp(df_data, list(expression(WT == 999)),
                   noBaseThetas = 2, dfParameters = df_params, functionList = list(p_func)),
    "no available data for subset"
  )
})


test_that("getForestDFemp covers data frame reference logic (Lines 142-148)", {
  # 1. Setup Mock Data
  df_params <- data.frame(THETA1 = 10, THETA2 = 2)
  df_data <- data.frame(ID = 1:2, WT = c(70, 80), SEX = c(1, 2))
  ls_expr <- list("WT" = expression(WT == 70), "WT" = expression(WT == 80))

  # FIX: p_func must return a list of numeric values to avoid "length zero" errors
  p_func <- function(thetas, df, ...) {
    # internalCalc can pass rows as vectors; coerce to list for $ subsetting
    df_l <- as.list(df)
    # Ensure a numeric value is returned even if WT is missing in the slice
    val <- if(!is.null(df_l$WT)) thetas[1] * (as.numeric(df_l$WT) / 70) else thetas[1]
    return(list(val))
  }

  # 2. Case A: Single-row data frame reference (Hits Lines 142-143)
  # Satisfies (m==1) and initializes valbase
  df_ref_single <- data.frame(WT = 70)
  res_single <- getForestDFemp(dfData = df_data,
                               covExpressionsList = ls_expr,
                               noBaseThetas = 2,
                               dfParameters = df_params,
                               functionList = list(p_func),
                               functionListName = "CL", # Length 1 matches list(val)
                               dfRefRow = df_ref_single)

  expect_s3_class(res_single, "data.frame")

  # 3. Case B: Multi-row data frame reference (Hits Lines 142-148)
  # Satisfies nrow(dfRefRow) > 1 and traverses the n/l assignment loops
  df_ref_multi <- data.frame(WT = c(70, 75))
  res_multi <- getForestDFemp(dfData = df_data,
                              covExpressionsList = ls_expr,
                              noBaseThetas = 2,
                              dfParameters = df_params,
                              functionList = list(p_func),
                              functionListName = "CL",
                              dfRefRow = df_ref_multi)

  expect_s3_class(res_multi, "data.frame")
  # Verifies that both expressions were processed
  expect_equal(nrow(res_multi), 2)
})

# --- oneHot argument ---

test_that("getForestDFemp oneHot matches manually pre-encoded dfData", {
  df_params <- data.frame(THETA1 = c(10, 11, 9, 12), THETA2 = c(0.3, 0.35, 0.28, 0.31))

  df_data <- data.frame(
    ID   = 1:8,
    WT   = c(60, 70, 80, 90, 65, 75, 72, 68),
    GENO = c(1, 2, 3, 4, 2, 3, 1, 4)
  )

  ls_expr <- list(
    "GENO" = expression(GENO == 1),
    "GENO" = expression(GENO == 3),
    "GENO" = expression(GENO == 4)
  )

  p_func <- function(thetas, df, ...) {
    x <- thetas[1]
    if (isTRUE(df$GENO_1 == 1)) x <- x * (1 + thetas[2])
    if (isTRUE(df$GENO_3 == 1)) x <- x * (1 - thetas[2])
    if (isTRUE(df$GENO_4 == 1)) x <- x * (1 + 2 * thetas[2])
    x
  }

  spec <- list(GENO = list(ref = 2))

  res_onehot <- getForestDFemp(
    dfData = df_data, covExpressionsList = ls_expr, noBaseThetas = 2,
    dfParameters = df_params, functionList = list(p_func), oneHot = spec
  )

  res_manual <- getForestDFemp(
    dfData = oneHotEncode(df_data, spec = spec), covExpressionsList = ls_expr,
    noBaseThetas = 2, dfParameters = df_params, functionList = list(p_func)
  )

  expect_equal(res_onehot, res_manual)
})

test_that("getForestDFemp oneHot = NULL leaves the result unchanged", {
  df_params <- data.frame(THETA1 = c(10, 11), THETA2 = c(2, 2.1))
  df_data <- data.frame(ID = 1:4, WT = c(60, 70, 80, 90), SEX = c(1, 1, 2, 2))
  ls_expr <- list("Sex" = expression(SEX == 1), "Sex" = expression(SEX == 2))
  p_func <- function(thetas, df, ...) thetas[1] * (df$WT / 70)

  a <- getForestDFemp(dfData = df_data, covExpressionsList = ls_expr, noBaseThetas = 2,
                      dfParameters = df_params, functionList = list(p_func))
  b <- getForestDFemp(dfData = df_data, covExpressionsList = ls_expr, noBaseThetas = 2,
                      dfParameters = df_params, functionList = list(p_func), oneHot = NULL)
  expect_equal(a, b)
})

test_that("getForestDFemp oneHotSep controls the dummy column separator", {
  df_params <- data.frame(THETA1 = c(10, 11), THETA2 = c(0.3, 0.32))
  df_data <- data.frame(ID = 1:8, WT = c(60,70,80,90,65,75,72,68),
                        GENO = c(1, 2, 3, 4, 2, 3, 1, 4))
  ls_expr <- list("GENO" = expression(GENO == 1), "GENO" = expression(GENO == 3))
  p_func <- function(thetas, df, ...) {
    x <- thetas[1]
    if (isTRUE(df$GENO1 == 1)) x <- x * (1 + thetas[2])
    if (isTRUE(df$GENO3 == 1)) x <- x * (1 - thetas[2])
    x
  }
  spec <- list(GENO = list(ref = 2))
  res <- getForestDFemp(dfData = df_data, covExpressionsList = ls_expr,
                        functionList = list(p_func), noBaseThetas = 2,
                        dfParameters = df_params, oneHot = spec, oneHotSep = "")
  ref <- getForestDFemp(dfData = oneHotEncode(df_data, spec = spec, sep = ""),
                        covExpressionsList = ls_expr, functionList = list(p_func),
                        noBaseThetas = 2, dfParameters = df_params)
  expect_equal(res, ref)
})

# --- coverage: dfRefRow validation, expression-ref errors, oneHot + data.frame ref ---

test_that("getForestDFemp errors on a data.frame dfRefRow with an invalid number of rows", {
  df_params <- data.frame(THETA1 = c(10, 11), THETA2 = c(2, 2.1))
  df_data   <- data.frame(ID = 1:4, WT = c(60, 70, 80, 90), SEX = c(1, 1, 2, 2))
  ls_expr   <- list("Sex" = expression(SEX == 1), "Sex" = expression(SEX == 2))
  bad_ref   <- data.frame(WT = c(70, 80, 90))   # 3 rows, expr list has 2

  expect_error(
    getForestDFemp(dfData = df_data, covExpressionsList = ls_expr, dfRefRow = bad_ref,
                   functionList = list(function(thetas, df, ...) thetas[1]),
                   noBaseThetas = 2, dfParameters = df_params),
    "number of reference rows/expressions"
  )
})

test_that("getForestDFemp errors when an expression reference selects no subjects", {
  df_params <- data.frame(THETA1 = c(10, 11), THETA2 = c(2, 2.1))
  df_data   <- data.frame(ID = 1:4, WT = c(60, 70, 80, 90), SEX = c(1, 1, 2, 2))
  ls_expr   <- list("Sex" = expression(SEX == 1), "Sex" = expression(SEX == 2))

  expect_error(
    getForestDFemp(dfData = df_data, covExpressionsList = ls_expr,
                   dfRefRow = list(expression(WT > 1000)),
                   functionList = list(function(thetas, df, ...) thetas[1] * (df$WT / 70)),
                   noBaseThetas = 2, dfParameters = df_params),
    "no available data for reference subset"
  )
})

test_that("getForestDFemp oneHot also encodes a data.frame dfRefRow", {
  df_params <- data.frame(THETA1 = c(10, 11, 9, 12), THETA2 = c(0.3, 0.35, 0.28, 0.31))
  df_data   <- data.frame(ID = 1:8, WT = c(60, 70, 80, 90, 65, 75, 72, 68),
                          GENO = c(1, 2, 3, 4, 2, 3, 1, 4))
  ls_expr   <- list("GENO" = expression(GENO == 1), "GENO" = expression(GENO == 3))
  df_ref    <- data.frame(GENO = 2, WT = 70)   # >1 column so row-subsetting stays a data.frame

  p_func <- function(thetas, df, ...) {
    x <- thetas[1]
    if (isTRUE(df$GENO_1 == 1)) x <- x * (1 + thetas[2])
    if (isTRUE(df$GENO_3 == 1)) x <- x * (1 - thetas[2])
    list(CL = x)
  }
  spec <- list(GENO = list(ref = 2))

  res_onehot <- getForestDFemp(
    dfData = df_data, covExpressionsList = ls_expr, dfRefRow = df_ref,
    functionList = list(p_func), functionListName = "CL",
    noBaseThetas = 2, dfParameters = df_params, oneHot = spec
  )
  res_manual <- getForestDFemp(
    dfData   = oneHotEncode(df_data, spec = spec),
    covExpressionsList = ls_expr,
    dfRefRow = oneHotEncode(df_ref, spec = spec),
    functionList = list(p_func), functionListName = "CL",
    noBaseThetas = 2, dfParameters = df_params
  )
  expect_equal(res_onehot, res_manual)
})

test_that("getForestDFemp uses supplied cdfCovsNames for the row labels", {
  df_params <- data.frame(THETA1 = c(10, 11), THETA2 = c(2, 2.1))
  df_data   <- data.frame(ID = 1:4, WT = c(60, 70, 80, 90), SEX = c(1, 1, 2, 2))
  ls_expr   <- list("Sex" = expression(SEX == 1), "Sex" = expression(SEX == 2))

  res <- getForestDFemp(
    dfData = df_data, covExpressionsList = ls_expr,
    cdfCovsNames = c("Men", "Women"),
    functionList = list(function(thetas, df, ...) thetas[1] * (df$WT / 70)),
    functionListName = "CL", noBaseThetas = 2, dfParameters = df_params
  )
  expect_setequal(unique(res$COVNAME), c("Men", "Women"))
})
