# test-createInputForestData.R

library(testthat)

# 1. Test with a simple, non-nested list ("univariate" cases)
test_that("createInputForestData works with a simple list", {

  simple_list <- list(
    SEX = c(1, 2),
    WT = c(70, 80)
  )

  result_df <- createInputForestData(simple_list, iMiss = -99)
  # FIX: Reset row names to NULL to allow for comparison.
  row.names(result_df) <- NULL

  expected_df <- data.frame(
    SEX = c(1, 2, -99, -99),
    WT = c(-99, -99, 70, 80),
    COVARIATEGROUPS = c("SEX", "SEX", "WT", "WT")
  )

  expect_equal(result_df, expected_df)
  expect_equal(nrow(result_df), 4)
  expect_equal(ncol(result_df), 3)
})


# 2. Test with a nested list ("grouped" or "multivariate" case)
test_that("createInputForestData works with a grouped (nested) list", {

  grouped_list <- list(
    GENO = list(
      GENO_2 = c(0, 1, 0),
      GENO_3 = c(0, 0, 1)
    )
  )

  result_df <- createInputForestData(grouped_list, iMiss = -99)
  # FIX: Reset row names to NULL.
  row.names(result_df) <- NULL

  expected_df <- data.frame(
    GENO_2 = c(0, 1, 0),
    GENO_3 = c(0, 0, 1),
    COVARIATEGROUPS = c("GENO", "GENO", "GENO")
  )

  expect_equal(result_df, expected_df)
  expect_true(all(result_df$COVARIATEGROUPS == "GENO"))
})


# 3. Test with a mixed list of simple and grouped covariates
test_that("createInputForestData works with a mixed list", {

  mixed_list <- list(
    FORM = c(0, 1),
    GENO = list(
      GENO_2 = c(0, 1, 0),
      GENO_3 = c(0, 0, 1)
    ),
    WT = c(65, 115)
  )

  result_df <- createInputForestData(mixed_list, iMiss = -99)

  expect_equal(nrow(result_df), 2 + 3 + 2)
  expect_equal(ncol(result_df), 5)

  expect_setequal(
    names(result_df),
    c("FORM", "GENO_2", "GENO_3", "WT", "COVARIATEGROUPS")
  )

  form_row1 <- subset(result_df, COVARIATEGROUPS == "FORM" & FORM == 0)
  expect_equal(form_row1$GENO_2, -99)
  expect_equal(form_row1$WT, -99)

  geno_row2 <- subset(result_df, COVARIATEGROUPS == "GENO")[2, ]
  expect_equal(geno_row2$GENO_2, 1)
  expect_equal(geno_row2$GENO_3, 0)
  expect_equal(geno_row2$FORM, -99)
})


# 4. Test that the iMiss argument is correctly handled
test_that("createInputForestData respects the iMiss argument", {

  list_data <- list(
    SEX = c(1, 2),
    WT = c(70, 80)
  )

  result_df <- createInputForestData(list_data, iMiss = 0)
  # FIX: Reset row names to NULL.
  row.names(result_df) <- NULL

  expected_df <- data.frame(
    SEX = c(1, 2, 0, 0),
    WT = c(0, 0, 70, 80),
    COVARIATEGROUPS = c("SEX", "SEX", "WT", "WT")
  )

  expect_equal(result_df, expected_df)
  expect_false(any(result_df == -99))
})


# 5. Test edge case: empty list input
test_that("createInputForestData handles an empty list gracefully", {
  # With the function now fixed, this test should pass without error.
  result_df <- createInputForestData(list())

  expected_df <- data.frame(COVARIATEGROUPS = character(0))

  result_df$COVARIATEGROUPS <- as.character(result_df$COVARIATEGROUPS)

  expect_equal(result_df, expected_df)
  expect_equal(nrow(result_df), 0)
  expect_equal(ncol(result_df), 1)
})


# 6. Test integration with getCovStats output
test_that("createInputForestData works with output from getCovStats", {
  integration_test_data <- data.frame(
    ID = c(1, 1, 2, 3, 4, 5),
    WT = c(60.5, 61.0, 70.2, 80.8, 65.1, 90.3),
    SEX = c(0, 0, 1, 0, 1, 0),
    RACE = c(1, 1, 2, 3, 1, 2)
  )

  stats_list_from_getCovStats <- getCovStats(
    data = integration_test_data,
    covariates = c("WT", "SEX", "RACE"),
    minLevels = 4
  )

  forest_df <- createInputForestData(stats_list_from_getCovStats, iMiss = -99)

  expect_setequal(names(forest_df), c("WT", "SEX", "RACE_2", "RACE_3", "COVARIATEGROUPS"))
  expect_equal(nrow(forest_df), 2 + 2 + 3)

  group_counts <- table(forest_df$COVARIATEGROUPS)
  expect_equal(as.integer(group_counts["WT"]), 2)
  expect_equal(as.integer(group_counts["SEX"]), 2)
  expect_equal(as.integer(group_counts["RACE"]), 3)

  race_row1 <- subset(forest_df, COVARIATEGROUPS == "RACE")[1, ]
  expect_equal(race_row1$WT, -99)
  expect_equal(race_row1$SEX, -99)
  expect_equal(race_row1$RACE_2, 0)
  expect_equal(race_row1$RACE_3, 0)
})
