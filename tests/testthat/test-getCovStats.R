library(testthat)
library(dplyr)

# 1. SETUP: Create a controlled test data frame
test_data <- tibble::tribble(
  ~ID, ~WT, ~SEX, ~RACE, ~BMI, ~DOSAGE,
  1, 60.5, 0, 1, 22.1, 50,
  1, 61.0, 0, 1, 22.4, 50,
  2, 70.2, 1, 2, 25.3, 100,
  3, 80.8, 0, 3, -99,  150,
  4, 65.1, 1, 1, 23.5, 100,
  5, 90.3, 0, 2, 28.9, 200
)

test_that("Function handles basic continuous and binary covariates correctly", {
  # FIX: Calculate expected value dynamically to be robust to quantile algorithm differences.

  # 1. Manually prepare the data exactly as the function does internally.
  data_for_calc <- test_data %>% distinct(ID, .keep_all = TRUE)

  # 2. Calculate the expected result using this prepared data.
  expected_wt <- signif(quantile(data_for_calc$WT, p = c(0.05, 0.95)), digits = 3)
  expected_sex <- c(0, 1)

  # 3. Run the function.
  stats <- getCovStats(test_data, covariates = c("WT", "SEX"), missVal = -99, minLevels = 4)

  # 4. Compare the function's output to the dynamically calculated expected value.
  expect_equal(stats$WT, expected_wt)
  expect_identical(sort(stats$SEX), expected_sex)
})

test_that("Function correctly handles missing data indicated by `missVal`", {
  # Apply the same robust calculation pattern for BMI.
  data_for_calc <- test_data %>%
    distinct(ID, .keep_all = TRUE) %>%
    filter(BMI != -99)

  expected_bmi <- signif(quantile(data_for_calc$BMI, p = c(0.05, 0.95)), digits = 3)

  stats <- getCovStats(test_data, covariates = "BMI", missVal = -99, minLevels = 3)

  expect_equal(stats$BMI, expected_bmi)
})

test_that("a genuine NA is dropped, like missVal, and not treated as a level", {
  # `x != missVal` is NA where x is NA, and logical-NA row indexing keeps an
  # all-NA row instead of dropping it. The leaked NA used to be counted as a
  # third level, silently switching a binary covariate to the one-hot branch.
  n  <- 40
  df <- data.frame(ID = seq_len(n),
                   SEX = rep(0:1, length.out = n),
                   WT  = seq(50, 120, length.out = n))
  withNA <- df
  withNA$SEX[c(3, 8)] <- NA
  withNA$WT[c(5, 9)]  <- NA

  # binary: still the documented sorted vector, not a nested one-hot list
  expect_equal(getCovStats(withNA, "SEX", idVar = "ID")$SEX, c(0L, 1L))
  expect_type(getCovStats(withNA, "SEX", idVar = "ID")$SEX, "integer")

  # continuous: quantiles are computed, and match dropping the NAs by hand
  expect_equal(
    getCovStats(withNA, "WT", idVar = "ID")$WT,
    signif(quantile(df$WT[-c(5, 9)], p = c(0.05, 0.95)), digits = 3)
  )

  # NA and missVal are treated the same way
  asMissVal <- df
  asMissVal$WT[c(5, 9)] <- -99
  expect_equal(getCovStats(withNA, "WT", idVar = "ID"),
               getCovStats(asMissVal, "WT", idVar = "ID"))
})

test_that("a covariate with no non-missing value is refused, not dropped", {
  # setupDfCovs() used to emit zero rows for such a covariate, so it vanished
  # from the forest plot with no error. refValues() already stopped here.
  df <- data.frame(ID = 1:5, WT = c(60, 70, 80, 90, 100), AGE = rep(-99, 5))
  expect_error(getCovStats(df, "AGE", idVar = "ID"),
               "contains only missing values")
  expect_error(getCovStats(df, c("WT", "AGE"), idVar = "ID"),
               "contains only missing values")
  expect_error(setupDfCovs(df, covariates = c("WT", "AGE"), idVar = "ID"),
               "contains only missing values")
  # all-NA is refused the same way as all-missVal
  df$AGE <- NA_real_
  expect_error(getCovStats(df, "AGE", idVar = "ID"),
               "contains only missing values")
})

test_that("Function correctly handles multi-level categorical covariates", {
  # This test is unchanged as it does not involve quantiles.
  stats <- getCovStats(test_data, covariates = "RACE")

  expected_race <- list(
    RACE_2 = c(0, 1, 0),
    RACE_3 = c(0, 0, 1)
  )

  expect_identical(stats$RACE, expected_race)
})

test_that("Function arguments `probs` and `nsig` are respected", {
  # Apply the same robust calculation pattern here.
  data_for_calc <- test_data %>% distinct(ID, .keep_all = TRUE)

  # Test with non-default quantiles
  stats_probs <- getCovStats(test_data, "WT", probs = c(0.1, 0.9), minLevels = 4)
  expected_wt_probs <- signif(quantile(data_for_calc$WT, p = c(0.1, 0.9)), digits = 3)
  expect_equal(stats_probs$WT, expected_wt_probs)

  # Test with non-default significant digits
  stats_nsig <- getCovStats(test_data, "WT", nsig = 5, minLevels = 4)
  expected_wt_nsig <- signif(quantile(data_for_calc$WT, p = c(0.05, 0.95)), digits = 5)
  expect_equal(stats_nsig$WT, expected_wt_nsig)
})

test_that("Function logic for `minLevels` works as expected", {
  # This test is unchanged.
  stats_cat <- getCovStats(test_data, "DOSAGE")
  expect_true(is.list(stats_cat$DOSAGE))
  expect_length(stats_cat$DOSAGE, 3)

  stats_cont <- getCovStats(test_data, "DOSAGE", minLevels = 3)
  expect_true(is.numeric(stats_cont$DOSAGE))
  expect_length(stats_cont$DOSAGE, 2)
})

test_that("Function handles non-default `idVar` correctly", {
  # Apply the same robust calculation pattern here.
  test_data_subj <- test_data %>% rename(SUBJID = ID)
  data_for_calc <- test_data_subj %>% distinct(SUBJID, .keep_all = TRUE)

  expected_wt <- signif(quantile(data_for_calc$WT, p = c(0.05, 0.95)), digits = 3)

  stats <- getCovStats(test_data_subj, "WT", idVar = "SUBJID", minLevels = 4)

  expect_equal(stats$WT, expected_wt)
})

test_that("Function throws an error for non-existent covariates", {
  # This test is unchanged.
  expect_error(
    getCovStats(test_data, "NonExistentCovariate"),
    "Not all covariates are present in the data."
  )
})

test_that("getCovStats consistently sorts binary covariates regardless of appearance order", {
  # Here, 1 appears before 0 for the binary variable TRT
  mock_data <- data.frame(
    ID = 1:2,
    TRT = c(1, 0)
  )

  stats <- getCovStats(mock_data, covariates = "TRT")

  # The output vector should be strictly sorted: c(0, 1)
  expect_equal(as.numeric(stats$TRT), c(0, 1),
               info = "Binary covariates are not being properly sorted.")
})

# --- refLevels and sep (added for the one-hot encoding refinement) ---

test_that("getCovStats default multi-level output is unchanged by the new arguments", {
  # Regression lock: omitting refLevels/sep must reproduce the historical output.
  stats <- getCovStats(test_data, covariates = "RACE")
  expect_identical(stats$RACE, list(RACE_2 = c(0, 1, 0), RACE_3 = c(0, 0, 1)))
})

test_that("getCovStats refLevels selects a non-lowest reference level", {
  # RACE has levels 1, 2, 3. Make level 2 the reference.
  stats <- getCovStats(test_data, covariates = "RACE", catRef = list(RACE = 2))

  expect_identical(names(stats$RACE), c("RACE_1", "RACE_3"))
  # Level order in the vectors is still sorted 1, 2, 3; the reference (2) is the
  # all-zero row.
  expect_equal(stats$RACE$RACE_1, c(1, 0, 0))
  expect_equal(stats$RACE$RACE_3, c(0, 0, 1))
})

test_that("getCovStats sep controls the one-hot name separator", {
  stats <- getCovStats(test_data, covariates = "RACE", sep = "")
  expect_identical(names(stats$RACE), c("RACE2", "RACE3"))
})

test_that("getCovStats errors on a reference level that is not in the data", {
  expect_error(
    getCovStats(test_data, covariates = "RACE", catRef = list(RACE = 9)),
    "not present in the data"
  )
})
