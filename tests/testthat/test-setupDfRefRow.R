test_that("setupDfRefRow generates correct singleRef = TRUE geometry", {
  mock_data <- data.frame(
    ID = 1:5,
    WT = c(60, 70, 70, 80, 90), # Median 70
    SEX = c(1, 1, 1, 0, 0) # Mode 1
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
    SEX = c(1, 1, 1, 0, 0) # Mode 1
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

test_that("setupDfRefRow honours refLevels so its columns match setupDfCovs", {
  mock_data <- data.frame(
    ID   = 1:8,
    WT   = c(60, 70, 70, 80, 90, 65, 75, 72),
    GENO = c(1, 2, 2, 3, 4, 2, 3, 2) # mode is 2
  )

  df_covs <- setupDfCovs(
    mock_data,
    covariates = c("WT", "GENO"), catRef = list(GENO = 2)
  )

  # refLevels is deprecated; it must still work and forward to catRef.
  expect_warning(
    ref <- setupDfRefRow(
      dfCovs     = df_covs,
      data       = mock_data,
      covariates = c("WT", "GENO"),
      refLevels  = list(GENO = 2)
    ),
    "`refLevels` is deprecated"
  )

  # Column names must line up with df_covs (GENO_1, GENO_3, GENO_4)
  expect_true(all(c("GENO_1", "GENO_3", "GENO_4") %in% names(ref)))
  # Mode genotype is 2 (the reference) -> all GENO dummies at 0
  expect_equal(ref$GENO_1, 0)
  expect_equal(ref$GENO_3, 0)
  expect_equal(ref$GENO_4, 0)
})

test_that("setupDfRefRow errors when a covariate is entirely missing", {
  mock_data <- data.frame(
    ID = 1:4,
    WT = c(60, 70, 80, 90),
    GONE = c(-99, -99, -99, -99)
  )
  df_covs <- setupDfCovs(mock_data, covariates = "WT", idVar = "ID")

  expect_error(
    setupDfRefRow(df_covs,
      data = mock_data, covariates = c("WT", "GONE"),
      idVar = "ID"
    ),
    "contains only missing values"
  )
})
