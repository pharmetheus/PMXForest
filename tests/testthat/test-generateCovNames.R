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

test_that("generateCovNames name_map uses the display name in every row label", {
  df <- data.frame(
    COVARIATEGROUPS = c("WTKGBL", "WTKGBL", "SEXN", "SEXN", "AGE"),
    WTKGBL = c(58.3, 126, -99, -99, -99),
    SEXN = c(-99, -99, 1, 2, -99),
    AGE = c(-99, -99, -99, -99, 65)
  )

  res <- generateCovNames(df,
    name_map = c(WTKGBL = "Weight", SEXN = "Sex"),
    unit_map = c(WTKGBL = "kg", AGE = "years")
  )
  # Continuous rows need no label_map; a covariate missing from name_map keeps
  # its column name, and a row without a unit gets no trailing space
  expect_equal(res, c("Weight = 58.3 kg", "Weight = 126 kg", "Sex = 1", "Sex = 2", "AGE = 65 years"))

  # label_map still takes precedence, and a blank name_map entry is ignored
  res_lab <- generateCovNames(df,
    name_map = c(WTKGBL = "", SEXN = "Sex"),
    label_map = c("SEXN=1" = "Male")
  )
  expect_equal(res_lab, c("WTKGBL = 58.3", "WTKGBL = 126", "Male", "Sex = 2", "AGE = 65"))
})

test_that("generateCovNames label_format controls the name_map labels", {
  df <- data.frame(
    COVARIATEGROUPS = c("WTKGBL", "SEXN"),
    WTKGBL = c(58.3, -99),
    SEXN = c(-99, 1)
  )
  res <- generateCovNames(df,
    name_map = c(WTKGBL = "Weight", SEXN = "Sex"),
    unit_map = c(WTKGBL = "kg"),
    label_format = "{name}: {value} {unit} (P10)"
  )
  expect_equal(res, c("Weight: 58.3 kg (P10)", "Sex: 1 (P10)"))
})

test_that("generateCovNames names the groups from name_map unless group_map says otherwise", {
  df <- data.frame(
    COVARIATEGROUPS = c("WTKGBL", "WTKGBL", "SEXN", "AGE"),
    WTKGBL = c(58.3, 126, -99, -99),
    SEXN = c(-99, -99, 1, -99),
    AGE = c(-99, -99, -99, 65)
  )
  dfres_mock <- data.frame(COVNUM = 1:4, PARAMETER = "CL")

  mapped <- generateCovNames(df,
    dfres = dfres_mock,
    name_map = c(WTKGBL = "Weight", SEXN = "Sex"),
    group_map = c(SEXN = "Sex at birth"),
    unit_map = c(WTKGBL = "kg")
  )
  expect_equal(as.character(mapped$GROUPNAME), c("Weight", "Weight", "Sex at birth", "AGE"))
  expect_equal(levels(mapped$COVNAME), c("Weight = 58.3 kg", "Weight = 126 kg", "Sex = 1", "AGE = 65"))
})

test_that("generateCovNames level_map labels levels per covariate", {
  df <- data.frame(
    COVARIATEGROUPS = c("SEXN", "SEXN", "FOODDN", "FOODDN", "RACE", "RACE", "RACE"),
    SEXN = c(1, 2, -99, -99, -99, -99, -99),
    FOODDN = c(-99, -99, 0, 1, -99, -99, -99),
    RACE_2 = c(-99, -99, -99, -99, 0, 1, 0),
    RACE_3 = c(-99, -99, -99, -99, 0, 0, 1)
  )

  res <- generateCovNames(df,
    name_map = c(SEXN = "Sex", FOODDN = "Food status"),
    level_map = list(
      SEXN = c("1" = "Male", "2" = "Female"),
      FOODDN = c("0" = "Fasted", "1" = ""),
      RACE = c("1" = "White", "3" = "Asian")
    )
  )
  # A blank or missing level falls through to name_map, then to the column name
  expect_equal(res, c("Male", "Female", "Fasted", "Food status = 1", "White", "RACE = 2", "Asian"))

  # Without name_map, an unmapped level keeps the old format
  res_old <- generateCovNames(df, level_map = list(SEXN = c("1" = "Male")))
  expect_equal(res_old[1:3], c("Male", "SEXN=2", "FOODDN=0"))

  # label_map takes precedence over level_map
  res_lab <- generateCovNames(df,
    label_map = c("SEXN=1" = "Men"),
    level_map = list(SEXN = c("1" = "Male", "2" = "Female"))
  )
  expect_equal(res_lab[1:2], c("Men", "Female"))
})

test_that("generateCovNames drops rows whose level_map or label_map entry is NA", {
  df <- data.frame(
    COVARIATEGROUPS = c("WT", "WT", "RACEL1N", "RACEL1N", "RACEL2N", "RACEL2N"),
    WT = c(58, 126, -99, -99, -99, -99),
    RACEL1N = c(-99, -99, 0, 1, -99, -99),
    RACEL2N = c(-99, -99, -99, -99, 0, 1)
  )
  dfres_mock <- data.frame(COVNUM = rep(1:6, each = 2), PARAMETER = c("CL", "V"))

  # The character vector keeps one entry per dfCovs row, NA for a dropped one
  res <- generateCovNames(df, level_map = list(RACEL1N = c("0" = NA, "1" = "White")))
  expect_equal(res, c("WT=58", "WT=126", NA, "White", "RACEL2N=0", "RACEL2N=1"))

  mapped <- generateCovNames(df,
    dfres = dfres_mock,
    name_map = c(WT = "Weight", RACEL1N = "Race", RACEL2N = "Race"),
    level_map = list(
      RACEL1N = c("0" = NA, "1" = "White"),
      RACEL2N = c("0" = NA, "1" = "Black")
    ),
    label_map = c("WT=126" = NA)
  )
  expect_equal(unique(mapped$COVNUM), c(1, 4, 6))
  expect_equal(levels(mapped$COVNAME), c("Weight = 58", "White", "Black"))
  expect_equal(levels(mapped$GROUPNAME), c("Weight", "Race"))
  expect_false(anyNA(mapped$COVNAME))

  # A group with every row dropped leaves no empty factor level behind
  gone <- generateCovNames(df,
    dfres = dfres_mock,
    level_map = list(RACEL2N = c("0" = NA, "1" = NA))
  )
  expect_equal(levels(gone$GROUPNAME), c("WT", "RACEL1N"))
})

test_that("generateCovNames rejects a level_map that is not a named list of named vectors", {
  df <- data.frame(COVARIATEGROUPS = "SEXN", SEXN = 1)
  expect_error(generateCovNames(df, level_map = c("SEXN=1" = "Male")), "level_map must be a named list")
  expect_error(generateCovNames(df, level_map = list(c("1" = "Male"))), "level_map must be a named list")
  expect_error(generateCovNames(df, level_map = list(SEXN = c("Male", "Female"))), "level_map must be a named list")
})

test_that("generateCovNames dummy_groups merges 0/1 dummies into one group with one reference row", {
  df <- data.frame(
    COVARIATEGROUPS = c("WT", rep(c("GENO1N", "GENO3N", "GENO4N"), each = 2), "RACEL1N", "RACEL1N"),
    WT = c(70, rep(-99, 8)),
    GENO1N = c(-99, 0, 1, -99, -99, -99, -99, -99, -99),
    GENO3N = c(-99, -99, -99, 0, 1, -99, -99, -99, -99),
    GENO4N = c(-99, -99, -99, -99, -99, 0, 1, -99, -99),
    RACEL1N = c(rep(-99, 7), 0, 1)
  )
  geno <- list("2D6 genotype" = list(
    ref = "2D6 EM",
    levels = c(GENO1N = "2D6 UM", GENO3N = "2D6 IM", GENO4N = "2D6 PM")
  ))

  # Vector mode: the first 0 row carries the reference, the other 0 rows are NA
  res <- generateCovNames(df, dummy_groups = geno)
  expect_equal(res, c("WT=70", "2D6 EM", "2D6 UM", NA, "2D6 IM", NA, "2D6 PM", "RACEL1N=0", "RACEL1N=1"))

  dfres_mock <- data.frame(
    COVNUM = rep(1:9, each = 2), PARAMETER = c("CL", "V"),
    COVEFF = rep(c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), each = 2)
  )
  mapped <- generateCovNames(df,
    dfres = dfres_mock,
    name_map = c(WT = "Weight", GENO1N = "ignored"),
    group_map = c(GENO3N = "ignored too"),
    dummy_groups = c(geno, list(Race = list(levels = c(RACEL1N = "White"))))
  )
  expect_equal(levels(mapped$GROUPNAME), c("Weight", "2D6 genotype", "Race"))
  expect_equal(levels(mapped$COVNAME), c("Weight = 70", "2D6 EM", "2D6 UM", "2D6 IM", "2D6 PM", "White"))
  # No ref: the reference row is dropped as well
  expect_equal(as.character(mapped$COVNAME[mapped$GROUPNAME == "Race"]), c("White", "White"))
  expect_equal(unique(mapped$COVNUM), c(1, 2, 3, 5, 7, 9))

  # label_map still takes precedence for a single row
  res_lab <- generateCovNames(df, dummy_groups = geno, label_map = c("GENO4N=1" = "Poor"))
  expect_equal(res_lab[7], "Poor")
})

test_that("forestPlot onlySignificant keeps or drops a dummy group as a whole", {
  df <- data.frame(
    COVARIATEGROUPS = c("GENO1N", "GENO1N", "GENO3N", "GENO3N"),
    GENO1N = c(0, 1, -99, -99),
    GENO3N = c(-99, -99, 0, 1)
  )
  # Only GENO1N=1 has an effect, but the merged group is kept whole
  dfres_mock <- data.frame(
    COVNUM = 1:4, PARAMETER = factor("CL"), GROUPNAME = factor(c("GENO1N", "GENO1N", "GENO3N", "GENO3N")),
    COVEFF = c(FALSE, TRUE, FALSE, FALSE), POINT = 1, Q1 = 0.9, Q2 = 1.1,
    POINT_REL_REFFUNC = c(1, 2, 1, 1), Q1_REL_REFFUNC = 0.9, Q2_REL_REFFUNC = 1.1, REFFUNC = 1
  )
  mapped <- generateCovNames(df,
    dfres = dfres_mock,
    dummy_groups = list(G = list(ref = "EM", levels = c(GENO1N = "UM", GENO3N = "IM")))
  )
  plotData <- setupForestPlotData(mapped, onlySignificant = TRUE)
  expect_equal(as.character(unique(plotData$COVNAME)), c("EM", "UM", "IM"))
})

test_that("generateCovNames checks dummy_groups", {
  df <- data.frame(COVARIATEGROUPS = c("GENO1N", "GENO1N"), GENO1N = c(0, 1))
  expect_error(generateCovNames(df, dummy_groups = list(G = c(GENO1N = "UM"))), "dummy_groups must be")
  expect_error(generateCovNames(df, dummy_groups = list(list(levels = c(GENO1N = "UM")))), "dummy_groups must be")
  expect_error(
    generateCovNames(df, dummy_groups = list(G = list(ref = c("a", "b"), levels = c(GENO1N = "UM")))),
    "dummy_groups must be"
  )
  expect_error(
    generateCovNames(df, dummy_groups = list(
      A = list(levels = c(GENO1N = "UM")),
      B = list(levels = c(GENO1N = "UM"))
    )),
    "only one of dummy_groups: GENO1N"
  )
  expect_warning(
    generateCovNames(df, dummy_groups = list(G = list(ref = "EM", levels = c(GENO1N = "UM", GENO9N = "XX")))),
    "not in dfCovs: GENO9N"
  )
})
