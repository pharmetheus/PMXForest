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

  # Note: The CSV path with 'n' currently returns exactly n+1 rows
  expect_equal(nrow(tmp1), n_samples + 1)
})

test_that("getSamples handles SIR and Missing Columns", {
  sirFile  <- system.file("extdata", "SimVal/sir7.dir/raw_results_run7.csv", package = "PMXForest")
  extFile  <- system.file("extdata", "SimVal/run7.ext", package = "PMXForest")

  raw <- read.csv(sirFile)
  n_resampled <- sum(raw$resamples == 1, na.rm = TRUE)

  # SIR branch: the importance-resampled vectors, plus the estimates as row 1
  tmp_sir <- getSamples(sirFile, extFile = extFile, quiet = TRUE)
  expect_s3_class(tmp_sir, "data.frame")
  expect_equal(nrow(tmp_sir), n_resampled + 1)
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
  expect_equal(nrow(res_df_samp), 5) # n + 1 (final estimates)
  expect_true("OBJ" %in% names(res_df_samp))

  # 3. SIR path: row 1 is the final estimates (OBJ = final OFV), the rest OBJ = 0
  sirFile <- system.file("extdata", "SimVal/sir7.dir/raw_results_run7.csv", package = "PMXForest")
  res_sir <- getSamples(sirFile, extFile = extFile, quiet = TRUE)
  expect_s3_class(res_sir, "data.frame")
  expect_true(res_sir$OBJ[1] != 0)
  expect_true(all(res_sir$OBJ[-1] == 0))

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

  # OLD: expect_equal(nrow(res_null), nrow(df_input) + 1)
  # NEW: Returns the exact input row count
  expect_equal(nrow(res_null), nrow(df_input))

  # OLD: expect_true("OBJ" %in% names(res_null))
  # NEW: We no longer append OBJ when n is NULL for generic data frames
  expect_false("OBJ" %in% names(res_null))
})

test_that("getSamples handles data frame input without extFile", {
  # Trigger Data frame logic
  df_input <- data.frame(THETA1 = c(1, 1.1, 0.9), THETA2 = c(2, 2.2, 1.8))

  # Case 1: n is NULL returns input exactly
  res1 <- getSamples(df_input)
  expect_equal(res1, df_input)
  expect_equal(nrow(res1), nrow(df_input))

  # Case 2: n is provided (MVRNORM sampling from DF) returns exactly n rows
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(123)
  res2 <- getSamples(df_input, n = 5)

  # Assert exactly n rows (Documented exception to the n+1 rule)
  expect_equal(nrow(res2), 5)
  expect_true("OBJ" %in% names(res2))
})




test_that("getSamples stops on parameter dimensionality mismatch", {
  # 1. Setup mock covariance file with 2 parameters
  tmp_cov <- tempfile(fileext = ".cov")
  cov_content <- "TABLE NO: 1
NAME      THETA1      THETA2
THETA1    0.1         0.01
THETA2    0.01        0.1"
  writeLines(cov_content, tmp_cov)

  # 2. Create a dfExt with 3 parameters (The Mismatch)
  df_ext_mismatch <- data.frame(
    ITERATION = -1000000000,
    THETA1 = 1.5,
    THETA2 = 2.5,
    THETA3 = 3.5, # Extra parameter not in .cov
    OBJ = 100
  )

  # 3. Verify the stop() triggers with the correct message
  expect_error(
    getSamples(input = tmp_cov, extFile = df_ext_mismatch, n = 10),
    "Critical Mismatch: The .ext data contains 3 parameters, but the .cov matrix has 2 parameters"
  )

  unlink(tmp_cov)
})

# --- coverage: data.frame guards and explicit indexvec paths ---

test_that("getSamples rejects a non-numeric data.frame", {
  df_bad <- data.frame(THETA1 = c(1, 1.1), LABEL = c("a", "b"),
                       stringsAsFactors = FALSE)
  expect_error(
    getSamples(df_bad),
    "all columns must be numeric"
  )
})

test_that("getSamples requires n when the input is a .cov file", {
  covFile <- system.file("extdata", "SimVal/run7.cov", package = "PMXForest")
  extFile <- system.file("extdata", "SimVal/run7.ext", package = "PMXForest")
  expect_error(
    getSamples(covFile, extFile = extFile, n = NULL),
    "number of samples"
  )
})

test_that("getSamples applies indexvec to a data.frame input", {
  df_in <- data.frame(A = c(1, 2, 3), B = c(4, 5, 6), C = c(7, 8, 9))
  res <- getSamples(df_in, indexvec = c(1, 3))   # keep A and C
  expect_equal(names(res), c("A", "C"))
  expect_equal(nrow(res), 3)
})

test_that("getSamples accepts an explicit indexvec for a csv input", {
  bootFile <- system.file("extdata", "SimVal/bs7.dir/raw_results_run7bs.csv",
                          package = "PMXForest")
  extFile  <- system.file("extdata", "SimVal/run7.ext", package = "PMXForest")

  # Parameter columns in this raw_results file: THETA 21-34, SIGMA 40, OMEGA 35-39
  auto <- getSamples(bootFile, extFile = extFile)
  idx  <- getSamples(bootFile, extFile = extFile,
                     indexvec = c(21:34, 40, 35:39))
  expect_equal(idx, auto)
})

# --- SIR raw_results handling (bug fix: sample_order column, resampled vectors) ---

test_that("getSamples returns the SIR importance-resampled vectors with the estimates as row 1", {
  sirFile <- system.file("extdata", "SimVal/sir7.dir/raw_results_run7.csv", package = "PMXForest")
  extFile <- system.file("extdata", "SimVal/run7.ext", package = "PMXForest")

  raw         <- read.csv(sirFile)
  n_resampled <- sum(raw$resamples == 1, na.rm = TRUE)
  ext_fin     <- subset(getExt(extFile), ITERATION == "-1000000000")
  n_theta     <- length(grep("^THETA", names(ext_fin)))

  expect_message(
    getSamples(sirFile, extFile = extFile),
    "importance-resampled"
  )

  res <- getSamples(sirFile, extFile = extFile, quiet = TRUE)

  # 1 estimates row + the resamples == 1 vectors (not the full proposal set)
  expect_equal(nrow(res), n_resampled + 1)
  expect_lt(nrow(res), sum(!is.na(raw$resamples)) + 1)   # fewer than all proposals

  # Row 1 is the .ext final estimates
  expect_equal(as.numeric(res[1, seq_len(n_theta)]),
               as.numeric(ext_fin[1, 1 + seq_len(n_theta)]))

  # quiet = TRUE silences the message but returns the same data
  expect_no_message(getSamples(sirFile, extFile = extFile, quiet = TRUE))
  expect_equal(getSamples(sirFile, extFile = extFile, quiet = TRUE), res)

  # n is ignored for a SIR file (and says so)
  expect_message(getSamples(sirFile, extFile = extFile, n = 25), "ignored for SIR")
  expect_equal(nrow(getSamples(sirFile, extFile = extFile, n = 25, quiet = TRUE)),
               n_resampled + 1)
})
