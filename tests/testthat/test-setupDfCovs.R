test_that("setupDfCovs produces identical output to manual two-step process", {
  # Setup mock data
  mock_data <- data.frame(
    ID = c(1, 2, 3, 4, 5),
    WT = c(60, 70, 80, 90, 100),
    SEX = c(0, 1, 0, 1, 0),
    AGE = c(25, -99, 45, 50, 60)
  )

  covs_to_test <- c("WT", "SEX", "AGE")

  # Manual process
  manual_list <- getCovStats(mock_data, covs_to_test, missVal = -99)
  manual_df <- createInputForestData(manual_list, iMiss = -99)

  # Wrapper process
  wrapper_df <- setupDfCovs(mock_data, covs_to_test, missVal = -99)

  # Assert equivalence
  expect_equal(wrapper_df, manual_df, info = "Wrapper output diverges from manual two-step process.")
})

test_that("setupDfCovs correctly propagates custom missing values", {
  mock_data <- data.frame(
    ID = 1:3,
    WT = c(70, 80, -999) # Using -999 as custom missing
  )

  # If we set missVal to -999, the resulting data frame's non-applicable columns
  # should be filled with -999.
  wrapper_df <- setupDfCovs(mock_data, covariates = "WT", missVal = -999)

  # Since there's only one covariate, there are no non-applicable columns in this specific
  # simple output, but we can check if the underlying stats ignored the -999 properly
  # by ensuring the max weight evaluated isn't -999.
  expect_true(all(wrapper_df$WT != -999), info = "Custom missing value was not properly excluded from statistics.")
})


# --- Tests for conditionalCovs feature in setupDfCovs ---

test_that("setupDfCovs correctly handles continuous conditionalCovs with deduplication", {
  # Mock data: Subject 3 has 9 duplicate rows.
  # Deduplicated AGE vector: c(20, 60, 100) -> 3 unique levels.
  # Median = 60, Mean = 60.
  mock_data <- data.frame(
    ID = c(1, 2, rep(3, 9)),
    WT = c(20, 60, rep(100, 9)),
    AGE = c(20, 60, rep(100, 9))
  )

  # Test contRef = "median"
  # Setting minLevels = 2 forces AGE (3 levels) to be treated as continuous.
  df_median <- setupDfCovs(
    data = mock_data,
    covariates = "WT",
    conditionalCovs = "AGE",
    contRef = "median",
    minLevels = 2,
    idVar = "ID"
  )

  wt_rows <- df_median[df_median$COVARIATEGROUPS == "WT", ]
  expect_true(all(wt_rows$AGE == 60),
    info = "AGE did not correctly backfill with deduplicated median."
  )

  # Test contRef = "mean"
  df_mean <- setupDfCovs(
    data = mock_data,
    covariates = "WT",
    conditionalCovs = "AGE",
    contRef = "mean",
    minLevels = 2,
    idVar = "ID"
  )

  wt_rows_mean <- df_mean[df_mean$COVARIATEGROUPS == "WT", ]
  expect_true(all(wt_rows_mean$AGE == 60),
    info = "AGE did not correctly backfill with deduplicated mean."
  )
})

test_that("setupDfCovs accurately backfills binary categorical conditionalCovs (mode)", {
  mock_data <- data.frame(
    ID = 1:5,
    WT = rep(70, 5),
    FOOD = c(1, 1, 1, 0, 0) # Mode is 1
  )

  df_binary <- setupDfCovs(
    data = mock_data,
    covariates = "WT",
    conditionalCovs = "FOOD",
    minLevels = 3
  )

  wt_rows <- df_binary[df_binary$COVARIATEGROUPS == "WT", ]
  expect_true(all(wt_rows$FOOD == 1),
    info = "Binary conditional covariate did not inherit the baseline mode."
  )
})

test_that("setupDfCovs accurately handles one-hot mapping for multi-level conditionalCovs", {
  mock_data <- data.frame(
    ID = 1:6,
    WT = rep(70, 6),
    # RACE has 3 levels: 1 (ref), 2, 3. Mode is 2.
    RACE = c(1, 2, 2, 2, 3, 3)
  )

  # Min levels set to 4 ensures RACE is treated as multi-level categorical
  df_multi <- setupDfCovs(
    data = mock_data,
    covariates = "WT",
    conditionalCovs = "RACE",
    minLevels = 4
  )

  wt_rows <- df_multi[df_multi$COVARIATEGROUPS == "WT", ]

  # Because Mode is 2, the reference state should have RACE_2 = 1 and RACE_3 = 0.
  expect_true(all(wt_rows$RACE_2 == 1),
    info = "RACE_2 dummy did not correctly inherit the mode reference state."
  )
  expect_true(all(wt_rows$RACE_3 == 0),
    info = "RACE_3 dummy did not correctly inherit the mode reference state."
  )
})

test_that("setupDfCovs throws error if conditionalCovs consists purely of missing values", {
  mock_data <- data.frame(
    ID = 1:3,
    WT = c(70, 80, 90),
    BROKEN_COV = c(-99, -99, -99)
  )

  expect_error(
    setupDfCovs(mock_data, covariates = "WT", conditionalCovs = "BROKEN_COV", missVal = -99),
    regexp = "contains only missing values",
    info = "Function failed to stop when an conditional covariate was entirely missing."
  )
})

test_that("setupDfCovs handles alternative explicit references (useMissVal = FALSE)", {
  mock_data <- data.frame(
    ID = 1:5,
    WT = c(60, 70, 70, 80, 90), # Median = 70
    SEX = c(1, 1, 1, 0, 0) # Mode = 1
  )

  # When useMissVal = FALSE, the primary covariates should have their inactive
  # background states (-99) replaced by their reference values.
  df_alt <- setupDfCovs(
    data = mock_data,
    covariates = c("WT", "SEX"),
    useMissVal = FALSE,
    minLevels = 3
  )

  # Check WT varying rows: background SEX should be the mode (1)
  wt_rows <- df_alt[df_alt$COVARIATEGROUPS == "WT", ]
  expect_true(all(wt_rows$SEX == 1),
    info = "Alternative primary categorical covariate did not backfill with mode."
  )

  # Check SEX varying rows: background WT should be the median (70)
  sex_rows <- df_alt[df_alt$COVARIATEGROUPS == "SEX", ]
  expect_true(all(sex_rows$WT == 70),
    info = "Alternative primary continuous covariate did not backfill with median."
  )
})

# --- refLevels / sep passthrough ---

test_that("setupDfCovs passes refLevels through to the one-hot column names", {
  mock_data <- data.frame(
    ID   = 1:6,
    WT   = c(60, 70, 80, 90, 65, 75),
    GENO = c(1, 2, 3, 4, 2, 3)
  )

  # Default: lowest level (1) is the reference
  df_default <- setupDfCovs(mock_data, covariates = c("WT", "GENO"))
  expect_true(all(c("GENO_2", "GENO_3", "GENO_4") %in% names(df_default)))
  expect_false("GENO_1" %in% names(df_default))

  # refLevels: level 2 is the reference
  df_ref2 <- setupDfCovs(
    mock_data,
    covariates = c("WT", "GENO"), catRef = list(GENO = 2)
  )
  expect_true(all(c("GENO_1", "GENO_3", "GENO_4") %in% names(df_ref2)))
  expect_false("GENO_2" %in% names(df_ref2))
})

test_that("setupDfCovs sep passthrough and conditionalCovs backfill honour catRef", {
  mock_data <- data.frame(
    ID   = 1:8,
    WT   = c(60, 70, 80, 90, 65, 75, 72, 68),
    RACE = c(1, 2, 2, 3, 2, 2, 1, 2) # mode is 2
  )

  df <- setupDfCovs(
    mock_data,
    covariates = "WT",
    conditionalCovs = "RACE",
    catRef = list(RACE = 1),
    sep = "."
  )

  # catRef says RACE = 1 is the reference, so it is both the level dropped when
  # encoding (leaving RACE.2 and RACE.3) and the background state, which is
  # therefore all dummies at 0.
  expect_true(all(c("RACE.2", "RACE.3") %in% names(df)))
  wt_rows <- df$COVARIATEGROUPS == "WT"
  expect_true(all(df$RACE.2[wt_rows] == 0))
  expect_true(all(df$RACE.3[wt_rows] == 0))
})

test_that("without catRef the background still comes from the mode", {
  mock_data <- data.frame(
    ID   = 1:8,
    WT   = c(60, 70, 80, 90, 65, 75, 72, 68),
    RACE = c(1, 2, 2, 3, 2, 2, 1, 2) # lowest is 1, mode is 2
  )

  df <- setupDfCovs(mock_data,
    covariates = "WT", conditionalCovs = "RACE",
    sep = "."
  )

  # Unchanged from earlier versions: encoding drops the lowest level, the
  # background takes the most common one.
  expect_true(all(c("RACE.2", "RACE.3") %in% names(df)))
  wt_rows <- df$COVARIATEGROUPS == "WT"
  expect_true(all(df$RACE.2[wt_rows] == 1))
  expect_true(all(df$RACE.3[wt_rows] == 0))
})
