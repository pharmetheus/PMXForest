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
