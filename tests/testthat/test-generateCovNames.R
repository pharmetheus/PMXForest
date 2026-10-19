test_that("generateCovNames handles core continuous and categorical logic", {
  df <- data.frame(
    COVARIATEGROUPS = c("WT", "SEX", "RACEL", "RACEL", "RACEL"),
    WT = c(70, -99, -99, -99, -99),
    SEX = c(-99, 1, -99, -99, -99),
    RACEL_2 = c(-99, -99, 0, 1, 0),
    RACEL_3 = c(-99, -99, 0, 0, 1)
  )

  # 1. Base Logic and Polychotomous Reference Detection
  res <- generateCovNames(df)
  expect_equal(res, c("WT=70", "SEX=1", "RACEL=1", "RACEL=2", "RACEL=3"))

  # 2. Map Overrides (Units)
  res_units <- generateCovNames(df, unit_map = c(WT = "kg"))
  expect_equal(res_units[1], "WT=70 kg")

  # 3. Map Overrides (Labels)
  res_labels <- generateCovNames(df, label_map = c("SEX=1" = "Male", "RACEL=1" = "White"))
  expect_equal(res_labels[2], "Male")
  expect_equal(res_labels[3], "White")

  # 4. Blank overrides are ignored
  res_blank <- generateCovNames(df, label_map = c("SEX=1" = ""))
  expect_equal(res_blank[2], "SEX=1")
})

test_that("generateCovNames throws appropriate errors", {
  df_bad <- data.frame(WT = c(70, 80))
  expect_error(generateCovNames(df_bad), "must contain a 'COVARIATEGROUPS' column")
})

test_that("generateCovNames handles missing or unknown dummy columns gracefully", {
  df_unknown <- data.frame(
    COVARIATEGROUPS = c("MISSINGCOV"),
    WT = 70
  )
  # Should warn and return UNKNOWN
  expect_warning(res <- generateCovNames(df_unknown), "No matching columns found")
  expect_equal(res, "MISSINGCOV=UNKNOWN")
})

test_that("generateCovNames maps labels relationally to dfres", {
  df <- data.frame(
    COVARIATEGROUPS = c("WT", "WT", "SEX"),
    WT = c(70, 90, -99),
    SEX = c(-99, -99, 1)
  )

  # Mock dfres mimicking getForestDFSCM() output (2 parameters per COVNUM scenario)
  dfres_mock <- data.frame(
    COVNUM = c(1, 1, 2, 2, 3, 3),
    PARAMETER = c("CL", "V", "CL", "V", "CL", "V")
  )

  # 1. Successful relational mapping
  mapped_df <- generateCovNames(
    df,
    dfres = dfres_mock,
    unit_map = c(WT = "kg"),
    label_map = c("SEX=1" = "Male")
  )

  expect_s3_class(mapped_df, "data.frame")
  expect_true("COVNAME" %in% names(mapped_df))
  expect_s3_class(mapped_df$COVNAME, "factor")

  # The levels must be unique and strictly preserve the top-to-bottom order of dfCovs
  expected_levels <- c("WT=70 kg", "WT=90 kg", "Male")
  expect_equal(levels(mapped_df$COVNAME), expected_levels)

  # Ensure the labels mapped correctly to their respective COVNUMs
  expect_equal(as.character(mapped_df$COVNAME[mapped_df$COVNUM == 1]), c("WT=70 kg", "WT=70 kg"))
  expect_equal(as.character(mapped_df$COVNAME[mapped_df$COVNUM == 3]), c("Male", "Male"))
})

test_that("generateCovNames print_template = TRUE prints the three maps and returns NULL", {
  df <- data.frame(
    COVARIATEGROUPS = c("WT", "WT", "SEX"),
    WT = c(70, 90, -99),
    SEX = c(-99, -99, 1)
  )
  expect_output(
    out <- generateCovNames(df, print_template = TRUE),
    "group_map <- c\\("
  )
  expect_output(generateCovNames(df, print_template = TRUE), "unit_map <- c\\(")
  expect_output(generateCovNames(df, print_template = TRUE), "label_map <- c\\(")
  # the template lists the distinct groups and base strings
  expect_output(generateCovNames(df, print_template = TRUE), "\"WT=70\"")
  invisible(utils::capture.output(res <- generateCovNames(df, print_template = TRUE)))
  expect_null(res)
})

test_that("generateCovNames applies group_map to GROUPNAME in relational mode", {
  df <- data.frame(
    COVARIATEGROUPS = c("WT", "WT", "SEX"),
    WT = c(70, 90, -99),
    SEX = c(-99, -99, 1)
  )
  dfres_mock <- data.frame(
    COVNUM    = c(1, 1, 2, 2, 3, 3),
    PARAMETER = c("CL", "V", "CL", "V", "CL", "V")
  )

  mapped <- generateCovNames(
    df,
    dfres     = dfres_mock,
    group_map = c(WT = "Body Weight", SEX = "") # "" is ignored -> stays "SEX"
  )

  expect_s3_class(mapped$GROUPNAME, "factor")
  expect_equal(levels(mapped$GROUPNAME), c("Body Weight", "SEX"))
  expect_equal(
    as.character(mapped$GROUPNAME[mapped$COVNUM == 1]),
    c("Body Weight", "Body Weight")
  )
  expect_equal(as.character(mapped$GROUPNAME[mapped$COVNUM == 3]), c("SEX", "SEX"))
})

test_that("generateCovNames strictly enforces COVNUM requirement when dfres is passed", {
  df <- data.frame(
    COVARIATEGROUPS = c("WT"),
    WT = c(70)
  )

  # Missing COVNUM
  dfres_bad <- data.frame(PARAMETER = c("CL", "V"))

  expect_error(
    generateCovNames(df, dfres = dfres_bad),
    "must contain a 'COVNUM' column to securely map"
  )
})
