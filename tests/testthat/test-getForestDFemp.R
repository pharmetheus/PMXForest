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
