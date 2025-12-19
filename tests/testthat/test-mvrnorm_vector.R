test_that("mvrnorm_vector samples correctly", {
  # Robustness Pattern 1: Lock RNG version for cross-version reproducibility
  # This ensures set.seed(dSeed) inside the function behaves identically
  # regardless of the R version running the test.
  suppressWarnings(RNGversion("3.5.0"))

  # Basic checks for fixed values (deterministic)
  expect_equal(mvrnorm_vector(1, sigma = 1, iSampleIndex = 0), 1)
  expect_equal(mvrnorm_vector(c(1, 2), sigma = matrix(c(0, 0, 0, 1), ncol = 2), fixed_mu = c(0, 1)), c(1, 2))

  ## Setup data for sampling tests
  extFile <- getExt(system.file("extdata", "SimVal/run7.ext", package = "PMXForest"))
  dfcov   <- read.table(system.file("extdata", "SimVal/run7.cov", package = "PMXForest"),
                        fill = TRUE, header = TRUE, sep = "",
                        skip = 1, stringsAsFactors = FALSE)
  sigma   <- data.matrix(dfcov[, 2:ncol(dfcov)])

  finPar  <- subset(extFile, ITERATION == "-1000000000")
  mu      <- as.numeric(finPar[, -(c(1, ncol(finPar)))])

  fixedmu <- rep(FALSE, ncol(sigma))
  for (j in 1:ncol(sigma)) {
    fixedmu[j] <- all(sigma[, j] == 0)
  }

  ## Robustness Pattern 2: Explicit seed and tolerance for sampling
  # Use a fixed seed to make the random draw deterministic for the test
  dSeed_val <- 12345

  # Test iSampleIndex = 10
  dfParameters <- as.data.frame(mvrnorm_vector(mu = mu, sigma = sigma,
                                               iSampleIndex = 10,
                                               fixed_mu = fixedmu,
                                               dSeed = dSeed_val))

  expect_equal(nrow(dfParameters), 10)

  # Robustness Pattern 3: Use tolerance for numeric comparisons
  # (Avoids failure due to 1e-15 differences in different R builds)
  # If you later use snapshots, use: expect_snapshot_value(dfParameters, style = "json2", tolerance = 1e-6)

  # Test iSampleIndex = 1
  dfParamSingle <- mvrnorm_vector(mu = mu, sigma = sigma,
                                  iSampleIndex = 1,
                                  fixed_mu = fixedmu,
                                  dSeed = dSeed_val)

  expect_type(dfParamSingle, "double")

  # Optional: Verify the first value matches a known 'golden' value with tolerance
  # expect_equal(dfParamSingle[1], -0.0456, tolerance = 1e-6)

  # Reset RNG to default after test to avoid affecting other tests
  RNGkind(sample.kind = "default")
})

test_that("mvrnorm_vector handles edge cases and default logic", {
  suppressWarnings(RNGversion("3.5.0"))

  # 1. Trigger Line 34: Test with fixed_mu = NULL (Default behavior)
  # Uses the internal rep(0, length(mu)) logic
  res_null_fixed <- mvrnorm_vector(mu = c(1, 1), sigma = diag(2),
                                   fixed_mu = NULL, iSampleIndex = 1)
  expect_length(res_null_fixed, 2)

  # 2. Trigger Line 42: Test with no parameters fixed (Full Sigma)
  # This enters the 'if(!any(fixed_mu!=0))' block where mySigma <- sigma
  res_no_fixed <- mvrnorm_vector(mu = c(0, 0), sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2),
                                 fixed_mu = c(0, 0), iSampleIndex = 1)
  expect_equal(length(res_no_fixed), 2)

  # 3. Trigger Lines 54-57: Matrix reconstruction loop
  # We use iSampleIndex > 1 and a mix of fixed/non-fixed parameters.
  # This ensures we hit the 'else' block and the 'for' loop.
  n_samples <- 5
  mu_val    <- c(10, 20, 30)
  sig_val   <- diag(3)
  fixed_val <- c(0, 1, 0) # Middle is fixed, 1 and 3 are sampled

  res_matrix <- mvrnorm_vector(mu = mu_val, sigma = sig_val,
                               fixed_mu = fixed_val, iSampleIndex = n_samples)

  expect_equal(dim(res_matrix), c(n_samples, 3))
  # Verify the fixed parameter (index 2) remains exactly 20 across all samples
  expect_true(all(res_matrix[, 2] == 20))
})
