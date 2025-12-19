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
