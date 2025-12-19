test_that("getSamples works correctly for .cov input", {
  # Lock RNG for cross-version stability
  suppressWarnings(RNGversion("3.5.0"))

  covFile <- system.file("extdata", "SimVal/run7.cov", package = "PMXForest")
  extFile <- system.file("extdata", "SimVal/run7.ext", package = "PMXForest")

  # Error handling coverage
  expect_error(getSamples(10), "input needs to be a character string")
  expect_error(getSamples("test.txt"), "either a .cov or .csv extension")

  set.seed(123)
  n_samples <- 20
  tmp <- getSamples(covFile, extFile = extFile, n = n_samples)

  # Robustness: Check structure and values
  expect_s3_class(tmp, "data.frame")
  # .cov path adds the final estimate as the first row, so total is n + 1
  expect_equal(nrow(tmp), n_samples + 1)

  # Check that the first row matches the final estimates from the .ext file
  # In run7.ext, THETA1 is roughly 1.0
  expect_equal(tmp[1, 1], 1.0, tolerance = 0.1)
})

test_that("getSamples handles CSV with and without 'n' (Bootstrap/SIR)", {
  suppressWarnings(RNGversion("3.5.0"))
  bootFile <- system.file("extdata", "SimVal/bs7.dir/raw_results_run7bs.csv", package = "PMXForest")
  extFile  <- system.file("extdata", "SimVal/run7.ext", package = "PMXForest")

  # Standard Bootstrap (n=NULL)
  tmp0 <- getSamples(bootFile, extFile = extFile)
  expect_s3_class(tmp0, "data.frame")
  expect_true("OBJ" %in% names(tmp0))

  # Small Bootstrap (n=20) - triggers covariance logic in Line 213
  set.seed(123)
  n_samples <- 20
  tmp1 <- getSamples(bootFile, extFile = extFile, n = n_samples)

  # Note: The CSV path with 'n' currently returns exactly n rows
  expect_equal(nrow(tmp1), n_samples)
})

test_that("getSamples handles SIR and Missing Columns", {
  sirFile  <- system.file("extdata", "SimVal/sir7.dir/raw_results_run7.csv", package = "PMXForest")
  extFile  <- system.file("extdata", "SimVal/run7.ext", package = "PMXForest")

  # SIR logic (Line 203)
  tmp_sir <- getSamples(sirFile, extFile = extFile)
  expect_s3_class(tmp_sir, "data.frame")
  expect_true(nrow(tmp_sir) > 0)
})

test_that("getSamples handles TTE models (Missing SIGMA/OMEGA logic)", {
  # Trigger Lines 188-198: Logic for models missing SIGMA or OMEGA
  runno <- "tte_weibull"
  bootFile <- system.file("extdata", "tte", "bootstrap_tte_weibull_n500", paste0("raw_results_", runno, ".csv"), package = "PMXForest")
  extFile <- system.file("extdata", "tte", paste0(runno, ".ext"), package = "PMXForest")

  tmp_tte <- getSamples(bootFile, extFile)
  # Verify dummy SIGMA.1.1. was added to satisfy requirements
  expect_true("SIGMA.1.1." %in% names(tmp_tte))
})

test_that("getSamples handles data frame input without extFile", {
  # Trigger Lines 142-159: Data frame logic
  df_input <- data.frame(THETA1 = c(1, 1.1), THETA2 = c(2, 2.2))

  # Case 1: n is NULL returns input exactly
  res1 <- getSamples(df_input)
  expect_equal(res1, df_input)

  # Case 2: n is provided (MVRNORM sampling from DF)
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(123)
  res2 <- getSamples(df_input, n = 5)
  expect_equal(nrow(res2), 5)
  expect_true("OBJ" %in% names(res2))
})


test_that("getSamples comprehensive coverage", {
  # Robustness Pattern: Lock RNG and define tolerance
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(123)

  extFile  <- system.file("extdata", "SimVal/run7.ext", package = "PMXForest")
  bootFile <- system.file("extdata", "SimVal/bs7.dir/raw_results_run7bs.csv", package = "PMXForest")

  # 1. Trigger Validation Gaps (Lines 74, 80-82)
  expect_error(getSamples("test.txt"), "either a .cov or .csv extension")
  expect_error(getSamples(bootFile, extFile = "wrong.txt"), "does not have a .ext extension")
  expect_error(getSamples(bootFile, extFile = NULL), "Need to provide an .ext file")

  # 2. Trigger Data Frame Reconstruction with sampling (Lines 107-111)
  # We use the correct column indices to match run7.ext (20 parameters)
  df_raw <- read.csv(bootFile, stringsAsFactors = FALSE)
  param_indices <- c(21:34, 40, 35:39)
  df_input <- df_raw[df_raw$ofv != 0, param_indices]

  # This hits the 'else' branch of is.null(n) for data frames
  res_df_samp <- getSamples(df_input, extFile = extFile, n = 5)
  expect_equal(nrow(res_df_samp), 6) # n + 1 (final estimates)
  expect_true("OBJ" %in% names(res_df_samp))

  # 3. Trigger SIR Path (Lines 205-208)
  sirFile <- system.file("extdata", "SimVal/sir7.dir/raw_results_run7.csv", package = "PMXForest")
  res_sir <- getSamples(sirFile, extFile = extFile)
  expect_s3_class(res_sir, "data.frame")
  expect_equal(res_sir$OBJ[1], 0)

  # 4. Trigger Missing OMEGA logic (Lines 194-195)
  # Using the TTE model which lacks OMEGAs in raw results
  tteBoot <- system.file("extdata", "tte", "bootstrap_tte_weibull_n500", "raw_results_tte_weibull.csv", package = "PMXForest")
  tteExt  <- system.file("extdata", "tte", "tte_weibull.ext", package = "PMXForest")
  res_tte <- getSamples(tteBoot, tteExt)
  expect_true("OMEGA.1.1." %in% names(res_tte))

  # 5. Trigger Line 176: Missing raw_results_structure
  # Move bootFile to temp location without its companion files
  tmp_csv <- tempfile(fileext = ".csv")
  file.copy(bootFile, tmp_csv)
  expect_error(getSamples(tmp_csv, extFile = extFile), "does not exist and indexvec is not provided")
  unlink(tmp_csv)
})

test_that("getSamples input validation and data frame checks", {
  # Setup paths and data
  bootFile <- system.file("extdata", "SimVal/bs7.dir/raw_results_run7bs.csv", package = "PMXForest")
  extFile  <- system.file("extdata", "SimVal/run7.ext", package = "PMXForest")

  # 1. Trigger Initial Validation Errors (Lines 74, 77, 80-82)
  # Missing extension
  expect_error(getSamples("testfile"), "needs to have an extension")

  # File does not exist
  expect_error(getSamples("nonexistent.cov"), "Can not find nonexistent.cov")

  # Incorrect .ext extension (Line 80)
  expect_error(getSamples(bootFile, extFile = "test.txt"), "does not have a .ext extension")

  # Nonexistent .ext file (Line 81)
  expect_error(getSamples(bootFile, extFile = "missing.ext"), "Can not find missing.ext")

  # Missing .ext when input is .csv (Line 82)
  expect_error(getSamples(bootFile, extFile = NULL), "Need to provide an .ext file")

  # 2. Trigger Data Frame Path with n = NULL (Lines 87, 92-95)
  # Prepare a numeric data frame to avoid cov() errors later
  df_raw <- read.csv(bootFile, stringsAsFactors = FALSE)
  # Use correct indices to match .ext structure
  param_indices <- c(21:34, 40, 35:39)
  df_input <- df_raw[df_raw$ofv != 0, param_indices]

  # Exercise indexvec subsetting (Line 87) and non-sampling DF branch (Line 92-95)
  # This covers the 'red' lines in the 'if(is.data.frame(input))' block
  res_null <- getSamples(df_input, extFile = extFile, n = NULL)

  expect_s3_class(res_null, "data.frame")
  expect_equal(nrow(res_null), nrow(df_input) + 1) # Final estimate + samples
  expect_true("OBJ" %in% names(res_null))
})
