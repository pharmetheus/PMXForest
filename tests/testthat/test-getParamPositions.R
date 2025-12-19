test_that("getParamPositions works correctly with standard files", {
  bootFile  <- system.file("extdata", "SimVal/bs7.dir/raw_results_run7bs.csv", package = "PMXForest")
  rr_struct <- file.path(dirname(bootFile), "raw_results_structure")

  skip_if_not(file.exists(rr_struct))

  res <- getParamPositions(rr_struct)

  # Robustness: Explicit comparison handles R-version differences in list printing
  expect_type(res, "list")
  expect_named(res, c("THETA", "SIGMA", "OMEGA"))

  # Values updated based on the actual SimVal/bs7.dir structure:
  expect_equal(res$THETA, c(20, 14))
  expect_equal(res$SIGMA, c(39, 1))
  expect_equal(res$OMEGA, c(34, 5))
})

test_that("getParamPositions handles incomplete or duplicate files", {
  # 1. Trigger Lines 48-50: File with only THETAs and SIGMAs (Missing OMEGAs)
  tmp_file <- tempfile()
  writeLines(c("theta=1,5", "sigma=10,2"), tmp_file)

  res_incomplete <- getParamPositions(tmp_file)
  expect_equal(res_incomplete$THETA, c(1, 5))
  expect_equal(res_incomplete$SIGMA, c(10, 2))
  expect_null(res_incomplete$OMEGA)

  # 2. Trigger Safety Breaks (Lines 25, 31, 37): Duplicate entries
  tmp_dup <- tempfile()
  writeLines(c("theta=1,5", "theta=6,10"), tmp_dup)

  res_dup <- getParamPositions(tmp_dup)
  # Should break at the second 'theta=' and return the first valid set found
  expect_equal(res_dup$THETA, c(1, 5))

  # 3. Trigger Line 21: Empty File handling
  tmp_empty <- tempfile()
  writeLines(character(0), tmp_empty)
  res_empty <- getParamPositions(tmp_empty)
  expect_null(res_empty)

  unlink(c(tmp_file, tmp_dup, tmp_empty))
})
