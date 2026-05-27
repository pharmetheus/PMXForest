test_that("setupDfRefRow generates correct singleRef = TRUE geometry", {
  mock_data <- data.frame(
    ID = 1:5,
    WT = c(60, 70, 70, 80, 90), # Median 70
    SEX = c(1, 1, 1, 0, 0)      # Mode 1
  )

  df_covs <- setupDfCovs(mock_data, covariates = c("WT", "SEX"), minLevels = 3)

  ref_single <- setupDfRefRow(
    dfCovs = df_covs,
    data = mock_data,
    covariates = c("WT", "SEX"),
    minLevels = 3
  )

  expect_equal(nrow(ref_single), 1, info = "singleRef = TRUE did not return a single row.")
  expect_equal(ref_single$WT, 70, info = "Continuous reference incorrect.")
  expect_equal(ref_single$SEX, 1, info = "Categorical reference incorrect.")
  expect_equal(ref_single$COVARIATEGROUPS, "Reference")
})

test_that("setupDfRefRow generates correct singleRef = FALSE geometry", {
  mock_data <- data.frame(
    ID = 1:5,
    WT = c(60, 70, 70, 80, 90), # Median 70
    SEX = c(1, 1, 1, 0, 0)      # Mode 1
  )

  # WT and SEX have -99s in inactive cells
  df_covs <- setupDfCovs(mock_data, covariates = c("WT", "SEX"), minLevels = 3)

  ref_matrix <- setupDfRefRow(
    dfCovs = df_covs,
    data = mock_data,
    covariates = c("WT", "SEX"),
    singleRef = FALSE,
    minLevels = 3
  )

  expect_equal(nrow(ref_matrix), nrow(df_covs), info = "singleRef = FALSE did not match dfCovs rows.")

  # Extract the WT varying block.
  # WT should be 70 everywhere in this block. SEX should remain -99.
  wt_block <- ref_matrix[ref_matrix$COVARIATEGROUPS == "WT", ]
  expect_true(all(wt_block$WT == 70), info = "Active continuous cells were not overwritten with reference.")
  expect_true(all(wt_block$SEX == -99), info = "Inactive missVal cells were not preserved.")
})
